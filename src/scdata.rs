/// Cellids is a class that should implement the Rhapsody_cell_keys.r.txt functionality
/// Mainly to translate between R1 read and cell ID.
/// It should also fix some read errors in the cell ids. But that will come a little later

use std::collections::BTreeMap;
use std::collections::HashSet;
use std::io::BufReader;
use std::io::BufRead;

//use crate::geneids::GeneIds;
//use crate::fast_mapper::FastMapper;
use mapping_info::MappingInfo;
use crate::cell_data::GeneUmiHash;
//use crate::ambient_rna_detect::AmbientRnaDetect;
use crate::CellData;
//use crate::cellids::CellIds;
use crate::indexed_genes::IndexedGenes;
use int_to_str::int_to_str::IntToStr;

use std::io::BufWriter;
use std::fs::File;
use std::io::Write;
use std::fs;

use flate2::Compression;
use flate2::write::GzEncoder;
use flate2::read::GzDecoder;

use std::path::PathBuf;
use std::path::Path;

use rayon::prelude::*;
use core::fmt;

#[derive(Debug, PartialEq, Clone)]
pub enum MatrixValueType {
    Integer,
    Real,
    Complex,
    Pattern,
    Unknown(String),
}

// This Scdata needs to copy some of the logics from split2samples - no it actually is totally different
// Here we look for new sample ids and each sample id needs to be a total match to the previousely identified sample id
// Of cause I could also implement something with a whitelist. But that is for the future.
pub struct Scdata{    
    //kmer_size: usize,
    //kmers: BTreeMap<u64, u32>,
    data: [BTreeMap< u64, CellData>; u8::MAX as usize],
    genes_to_print: Vec::<String>,
    //ambient_cell_content: BTreeMap< u64, CellData>,
    checked: bool,
    passing: usize,
    pub genes_with_data: HashSet<usize>,
    pub num_threads:usize,
    //ambient_store:AmbientRnaDetect,
    value_type: MatrixValueType,
}

impl Default for Scdata {
    fn default() -> Self {
        Self::new(1, MatrixValueType::Integer)
    }
}


/*impl fmt::Display for Scdata {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let checked = match self.checked{
            true => "checked".to_string(),
            false => "unchecked".to_string(),
        };
        write!(f, "Scdata ({}) with {} cells and {} genes ({})", self.value_type, self.keys().len(), self.genes_with_data.len(), checked );

    }
}*/

impl fmt::Display for Scdata {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "Scdata Summary")?;
        writeln!(f, "===============")?;
        writeln!(f, "Checked: {}", self.checked)?;
        writeln!(f, "Passing Cells: {}", self.passing)?;
        writeln!(f, "Matrix Value Type: {:?}", self.value_type)?;
        writeln!(f, "Genes to Print: [{} genes]", self.genes_to_print.len())?;
        writeln!(f, "Genes with Data: [{} entries]", self.genes_with_data.len())?;

        let mut total_cells = 0;
        let mut non_empty_buckets = 0;
        for (i, map) in self.data.iter().enumerate() {
            if !map.is_empty() {
                total_cells += map.len();
                non_empty_buckets += 1;
                writeln!(f, "  Bucket {}: {} cells", i, map.len())?;
                if f.alternate() {
                    for (cell_id, cell_data) in map.iter().take(2) { // preview only
                        writeln!(f, "    Cell ID {}:\n{}", cell_id, cell_data)?;
                    }
                    if map.len() > 2 {
                        writeln!(f, "    ... {} more cells", map.len() - 2)?;
                    }
                }
            }
        }
        writeln!(f, "Total Cells: {}", total_cells)?;
        writeln!(f, "Non-empty Buckets: {}", non_empty_buckets)?;

        Ok(())
    }
}

// here the functions
impl Scdata{

    //pub fn new(kmer_size:usize )-> Self {
    pub fn new(num_threads:usize, value_type: MatrixValueType )-> Self {

        const EMPTY_MAP: BTreeMap<u64, CellData> = BTreeMap::new();
        let data = [EMPTY_MAP ;u8::MAX as usize];
        let checked:bool = false;
        let passing = 0;
        let genes_with_data = HashSet::new();
        Self {
            //kmer_size,
            data,
            genes_to_print: Vec::<String>::with_capacity(4),
            //ambient_cell_content: BTreeMap::new(),
            checked,
            passing,
            genes_with_data,
            num_threads,
            //ambient_store:AmbientRnaDetect::new(),
            value_type,
        }
    }

    pub fn value_type( &mut self, value_type: Option<MatrixValueType> ) -> &MatrixValueType {
        if let Some(val) = value_type{
            self.value_type = val;
        }
        &self.value_type
    }

    fn to_key(&self, name: &u64 ) -> usize{
        // 48 = 64 -16
        (*name >> 54) as usize
    }

    fn values(&self) -> impl Iterator<Item = &CellData> {
        self.data.iter().flat_map(|map| map.values())
    }

    pub fn keys(&self) -> Vec<u64> {
        let mut all_keys = Vec::new();
        for map in &self.data {
            all_keys.extend(map.keys().copied());
        }
        all_keys
    }

    pub fn is_empty(&self) -> bool{
        self.data.is_empty()
    }

    pub fn len(&self) -> usize {
        let mut size = 0;
        for map in &self.data {
            size += map.len();
        }
        size
    } 

    fn get_mut(&mut self, key: &u64) -> Option<&mut CellData> {
        let index = self.to_key(key); // Extracting the first u8 of the u64 key
        self.data[index].get_mut(key)
    }

    pub fn get(&self, key: &u64) -> Option<&CellData> {
        let index = self.to_key(key); // Extracting the first u8 of the u64 key
        self.data[index].get(key)
    }

    pub fn keep_only_passing_cells(&mut self) {
        for map in &mut self.data {
            map.retain(|_, cell_data| cell_data.passing);
        }
    }

    pub fn merge(&mut self, other: &Scdata) {
        if other.is_empty() {
            return;
        }

        self.checked = false;
        self.passing = 0;
        self.genes_with_data.clear();

        // Parallelize over the indices of `self.data`
        self.data
        .par_iter_mut()
        .enumerate()
        .for_each(|(index, self_bucket)| {
            // Access the corresponding bucket in `other.data`
            if let Some(other_bucket) = other.data.get(index) {
                for (cell_name, other_cell) in other_bucket {

                    match self_bucket.entry(*cell_name) {
                        std::collections::btree_map::Entry::Occupied(mut entry) => {
                            // If cell exists, merge with existing cell
                            let cell = entry.get_mut();
                            cell.merge(other_cell);
                        }
                        std::collections::btree_map::Entry::Vacant(entry) => {
                            // If cell doesn't exist, insert new cell by copying data from other_cell
                            entry.insert(other_cell.deep_clone()); // Assumes `CellData` implements `Clone`
                        }
                    }
                }
            }
        });
    }

    /// merge two Scdata objects - keep track of the umis!
    pub fn merge_single_thread(&mut self, mut other: Scdata) {
        if ! other.is_empty() {
            // Reset all internal measurements
            self.checked = false;
            self.passing = 0;
            self.genes_with_data.clear();

            for other_cell in other.data.iter_mut().flat_map(|map| map.values_mut()) {
                let index = self.to_key(&other_cell.name); // Extracting the first u8 of the u64 key
                match self.data[index].entry(other_cell.name) {
                    std::collections::btree_map::Entry::Occupied(mut entry) => {
                        // If cell exists, merge with existing cell
                        //println!( "merge cell {} {}", other_cell.name, other_cell);
                        let cell = entry.get_mut();
                        cell.merge(other_cell);
                    }
                    std::collections::btree_map::Entry::Vacant(entry) => {
                        // If cell doesn't exist, insert new cell
                        //println!( "steel cell {} {}", other_cell.name, other_cell);
                        entry.insert(std::mem::take(other_cell)); // Move ownership
                    }
                }
            }
        }
    }

    /// merge two Scdata objects - with a different gene list!!
    pub fn merge_re_id_genes(&mut self, other: Scdata, other_ids: &Vec::<usize> ) {
        if ! other.is_empty() {
            // Reset all internal measurements
            self.checked = false;
            self.passing = 0;
            self.genes_with_data.clear();

            for other_cell in other.data.iter().flat_map(|map| map.values()) {
                let index = self.to_key(&other_cell.name); // Extracting the first u8 of the u64 key
                match self.data[index].entry(other_cell.name) {
                    std::collections::btree_map::Entry::Occupied(mut entry) => {
                        // If cell exists, merge with existing cell
                        //println!( "merge cell {} {}", other_cell.name, other_cell);
                        let cell = entry.get_mut();
                        cell.merge_re_id_genes(other_cell, &other_ids);
                    }
                    std::collections::btree_map::Entry::Vacant(entry) => {
                        // If cell doesn't exist, insert new cell
                        //println!( "steel cell {} {}", other_cell.name, other_cell);
                        let mut cell = CellData::new( other_cell.name.clone() );
                        cell.merge_re_id_genes( other_cell, &other_ids );
                        entry.insert( cell ); // Move ownership
                    }
                }
            }

        }
    }



    /// try_insert now is more slick and insterts the data in 256 different sub-areas.
    pub fn try_insert(&mut self, name: &u64, data: GeneUmiHash, value:f32, report: &mut MappingInfo) -> bool {
        let index = self.to_key(name); // Extracting the first u8 of the u64 key
        self.genes_with_data.insert(data.0);
        self.checked = false;

        let cell_info = self.data[index]
            .entry(*name)
            .or_insert_with(|| CellData::new( *name )); // Insert new cell if not exists

        report.ok_reads += 1;
        if !cell_info.add(data, value) {
            report.pcr_duplicates += 1;
            report.local_dup += 1;
            false
        } else {
            true
        }
    }


    pub fn write (&mut self, file_path: PathBuf, genes:&IndexedGenes, min_count:usize) -> Result< (), &str>{

        let names = genes.get_all_gene_names();
        return self.write_sub( file_path, genes, &names, min_count);
    }

    pub fn write_sub (&mut self, file_path: PathBuf, genes:&IndexedGenes, names: &Vec<String>, min_count:usize) -> Result< (), &str>{
    
        let rs:bool = Path::new( &file_path ).exists();
        if rs && fs::remove_file(  &file_path ).is_ok(){};
        
        let file = match File::create( file_path ){
            Ok(file) => file,
            Err(err) => {
                panic!("Error: {err:#?}");
            }
        };
        let mut writer = BufWriter::new(&file);

        match writeln!( writer, "{}", genes.to_header_n( names ) ){
            Ok(_) => (),
            Err(err) => {
                eprintln!("write error: {err}" );
                return Err::<(), &str>("Header could not be written")
            }
        };

        let mut passed = 0;

        if ! self.checked{
            self.update_genes_to_print( genes, names);
            self.mtx_counts( genes, min_count, self.num_threads );
        }

        //println!("We are here exporting a samples table and want these samples to be included: {:?}", names );
        //println!("And we have these ids for them: {:?} using the offset", genes.ids_for_gene_names( names) );

        for cell_obj in self.values() {
            if ! cell_obj.passing{
                continue;
            }
            let text = cell_obj.to_str( genes, names, 32 );
            //println!("this should contain some info {}", text);
            match writeln!( writer, "{text}" ){
                Ok(_) => passed +=1,
                Err(err) => {
                    eprintln!("write error: {err}");
                    return Err::<(), &str>("cell data could not be written")   
                }
            };
        }

        println!( "dense matrix: {passed} cell written");
        Ok( () )
    }


    /// this will create a path and populate that with 10x kind of files.
    pub fn write_sparse (&mut self, file_path: PathBuf, genes: &IndexedGenes, min_count:usize) -> Result< String, String>{
        let names= genes.get_all_gene_names();
        match self.value_type{
            MatrixValueType::Integer => {
                self.write_sparse_sub( file_path, genes, &names, min_count)
            },
            MatrixValueType::Real => {
                self.write_sparse_sub_real( file_path, genes, &names, min_count)
            },
            _ => panic!("Not suported type")
        }
    }


    /// this utilizes the new f32 value (e.g. mean read quality) as data and writes a real MatrixMarket table
    pub fn write_sparse_sub_real(
        &mut self,
        file_path: PathBuf,
        genes: &IndexedGenes,
        names: &Vec<String>,
        min_count: usize
    ) -> Result<String, String> {
        let rs = Path::new(&file_path).exists();

        self.update_genes_to_print(genes, names);

        if self.genes_to_print.is_empty() {
            let err = format!("No genes to report on - no data written to path {:?}", file_path.to_str());
            eprintln!("{}", err);
            return Ok(err);
        }

        if !rs {
            if let Err(err) = fs::create_dir_all(file_path.clone()) {
                eprintln!("Error?: {err:#?}");
            }
        }

        let _ = fs::remove_file(file_path.join("matrix.mtx.gz"));
        let file = File::create(file_path.join("matrix.mtx.gz"))
            .map_err(|err| format!("Error creating the path?: {err:#?}"))?;
        let file1 = GzEncoder::new(file, Compression::default());
        let mut writer = BufWriter::with_capacity(4096, file1);

        let _ = fs::remove_file(file_path.join("barcodes.tsv.gz"));
        let file_b = File::create(file_path.join("barcodes.tsv.gz"))
            .map_err(|err| format!("Error creating the path?: {err:#?}"))?;
        let file2 = GzEncoder::new(file_b, Compression::default());
        let mut writer_b = BufWriter::with_capacity(4096, file2);

        let _ = fs::remove_file(file_path.join("features.tsv.gz"));
        let file_f = File::create(file_path.join("features.tsv.gz"))
            .map_err(|err| format!("Error creating the path?: {err:?}"))?;
        let file3 = GzEncoder::new(file_f, Compression::default());
        let mut writer_f = BufWriter::with_capacity(4096, file3);

        for name in &self.genes_to_print {
            writeln!(writer_f, "{name}\t{name}\tGene Expression")
                .map_err(|err| {
                    eprintln!("write error: {err}");
                    "feature could not be written".to_string()
                })?;
        }

        self.keep_only_passing_cells();
        let gene_ids = genes.ids_for_gene_names(&self.genes_to_print);
        let i2s = IntToStr::new("AAAA");

        let mut entries = 0;
        let mut cell_id = 0;

        // MatrixMarket header for float values
        writeln!(
            writer,
            "%%MatrixMarket matrix coordinate real general\n{}",
            self.mtx_counts(genes, min_count, self.num_threads)
        )
        .map_err(|err| {
            eprintln!("write error: {err}");
            "Header could not be written".to_string()
        })?;

        for cell_obj in self.values() {
            if !cell_obj.passing {
                continue;
            }

            cell_id += 1;
            //cell_obj.name is a u64

            let cell_name = IntToStr::from_u64( cell_obj.name ).to_string(16);

            writeln!(writer_b, "{}", &cell_name)
                .map_err(|err| {
                    eprintln!("write error: {err}");
                    "cell barcode could not be written".to_string()
                })?;

            // Aggregate mean value for each gene in the cell
            let mut gene_sums: BTreeMap<usize, (f32, usize)> = BTreeMap::new();

            for (gh, value) in &cell_obj.genes {
                if gene_ids.contains(&gh.0) {
                    let entry = gene_sums.entry(gh.0).or_insert((0.0, 0));
                    entry.0 += *value;
                    entry.1 += 1;
                }
            }

            for (i, gene_id) in gene_ids.iter().enumerate() {
                if let Some((sum, count)) = gene_sums.get(gene_id) {
                    if *count > 0 {
                        let mean = sum / *count as f32;
                        writeln!(writer, "{} {cell_id} {:.4}", i + 1, mean)
                            .map_err(|err| {
                                eprintln!("write error: {err}");
                                "cell data could not be written".to_string()
                            })?;
                        entries += 1;
                    }
                }
            }
        }

        let report = format!(
            "sparse Matrix (mean values): {} cell(s), {} gene(s) and {} entries written to path {:?}; ",
            cell_id,
            gene_ids.len(),
            entries,
            file_path.into_os_string().into_string()
        );
        println!("{}", report);
        Ok(report)
    }

    /// this utilizes the UMI counts data and writes an integer MatrixMarket table
    pub fn write_sparse_sub (&mut self, file_path: PathBuf, genes:&IndexedGenes, names: &Vec<String>, min_count:usize) -> Result<String, String>{
            
        let rs = Path::new( &file_path ).exists();

        self.update_genes_to_print( genes, names);

        if self.genes_to_print.is_empty(){
            let err = format!("No genes to report on - no data written to path {:?}", file_path.to_str());
            eprintln!( "{}", err );
            return Ok( err )
        }

        if ! rs {
            match fs::create_dir_all ( file_path.clone() ){
                Ok(_file) => (),
                Err(err) => {
                     eprintln!("Error?: {err:#?}");
                 }
            };
        }

        if fs::remove_file(file_path.join("matrix.mtx.gz") ).is_ok(){};

        let file = match File::create( file_path.join("matrix.mtx.gz") ){
            Ok(file) => file,
            Err(err) => {
                return Err(format!("Error creating the path?: {err:#?}"));
            }
        };
        let file1 = GzEncoder::new(file, Compression::default());
        let mut writer = BufWriter::with_capacity(4096,file1);


        //rs = Path::new( &file_path.clone().join("barcodes.tsv.gz") );
        //if rs {
        //    fs::remove_file( rs );
        //}
        if  fs::remove_file(file_path.join("barcodes.tsv.gz") ).is_ok(){};

        let file_b = match File::create( file_path.join("barcodes.tsv.gz") ){
            Ok(file) => file,
            Err(err) => {
                return Err(format!("Error creating the path?: {err:#?}"));
            }
        };
        let file2 = GzEncoder::new(file_b, Compression::default());
        let mut writer_b = BufWriter::with_capacity(4096,file2);
        match writeln!( writer, "%%MatrixMarket matrix coordinate integer general\n{}", 
            self.mtx_counts( genes, min_count, self.num_threads ) ){
            Ok(_) => (),
            Err(err) => {
                eprintln!("write error: {err}");
                return Err("Header could not be written".to_string());
            }
        };

        if fs::remove_file(file_path.join("features.tsv.gz") ).is_ok(){};
 
        let file_f = match File::create( file_path.join("features.tsv.gz") ){
            Ok(file) => file,
            Err(err) => {
                return Err(format!("Error creating the path?: {err:?}"));
            }
        };
        let file3 = GzEncoder::new(file_f, Compression::default());
        let mut writer_f = BufWriter::with_capacity(4096,file3);

        for name in &self.genes_to_print {
            match writeln!( writer_f, "{name}\t{name}\tGene Expression"  ){
                Ok(_) => (),
                Err(err) => {
                    eprintln!("write error: {err}" );
                    return Err("feature could not be written".to_string())   
                }
            }
        }

        self.keep_only_passing_cells();

        let mut entries = 0;
        let mut cell_id = 0;

        let gene_ids = genes.ids_for_gene_names( &self.genes_to_print );

        let i2s = IntToStr::new("AAAA");

        for cell_obj in self.values() {
            if ! cell_obj.passing {
                continue;
            }
            cell_id += 1;
            // this in fact reduces the cell id to 16 nucleotides (u32) the normal length of a cell id for 10x
            let cell_name =  IntToStr::u8_array_to_str( &cell_obj.name.to_le_bytes() );// IntToStr::from_u64(cell_obj.name).to_string(16); 
            //println!("Cell ID {} became seq {}", &cell_obj.name, &cell_name );
            match writeln!( writer_b,"{}", &cell_name){
                Ok(_) => (),
                Err(err) => {
                    eprintln!("write error: {err}");
                    return Err( "cell barcode could not be written".to_string())   
                }
            };

            for (gene_id, id) in gene_ids.iter().enumerate() {
                let n = cell_obj.n_umi_4_gene_id( id);
                if n > 0 {
                    match writeln!(writer, "{} {cell_id} {n}", gene_id+1) {
                        Ok(_) => { entries += 1; },
                        Err(err) => {
                            eprintln!("write error: {err}");
                            return Err( "cell data could not be written".to_string())
                        }
                    }
                }
            }
        } 
        let report = format!("sparse Matrix: {} cell(s), {} gene(s) and {} entries written to path {:?}; ", cell_id, gene_ids.len(), entries, file_path.into_os_string().into_string());
        println!( "{}",report );
        Ok( report )
    }
    /// Update the gene names for export to sparse
    /// returns the count of cells and the count of total gene values
    pub fn update_genes_to_print( &mut self, genes:&IndexedGenes, names:&Vec<String>) -> [usize; 2] {
        
        let mut entries = 0;

        self.genes_to_print.clear();

        let cell_keys:Vec<u64> = self.keys();
        let chunk_size = cell_keys.len() / self.num_threads +1;

        let gene_data:Vec<(BTreeMap<std::string::String, usize>, usize)> = cell_keys
        .par_chunks(chunk_size)
        .map(|chunk| {
            // Your parallel processing logic here...
            let mut names4sparse:  BTreeMap::<String, usize> = BTreeMap::new();
            let mut n:usize;
            let mut entry = 0;
            let gene_ids = genes.ids_for_gene_names( names );
            for key in chunk {
                if let Some(cell_obj) = self.get(key){
                    if self.checked & ! cell_obj.passing{
                        println!("ignoring cell");
                        continue;
                    }
                    for (id, int_id) in gene_ids.iter().enumerate() {
                    //for name in names {
                        n = *cell_obj.total_reads.get( int_id  ).unwrap_or(&0);
                        if n > 0{
                            entry +=1;
                            if ! names4sparse.contains_key ( &names[id] ){
                                names4sparse.insert( names[id].to_string() , names4sparse.len() + 1 );
                            } 
                        }
                    }
                }
            }
            (names4sparse, entry)
        }).collect();

        // Merge gene data from different chunks
        let mut names4sparse = HashSet::<String>::new();

        for (genelist, n) in &gene_data{
            entries += n;
            for name in genelist.keys() {
                names4sparse.insert( name.to_string() );
            }
        }

        self.genes_to_print = names4sparse.iter().cloned().collect();


        /*if genes.max_id  ==0 && ! names.is_empty() {
            let mut to = 10;
            if names.len() < 10{
                to = names.len() -1;
            }
            eprintln!( "None of the genes have data:\n{} ...", names[0..to].join( ", " ) );
        }
        //else { println!("{} genes requested and {} with data found", names.len(), genes.max_id); }
        if names.len() != genes.max_id{
            // better to run this once more - this does somehow not create the same count if more genes are checked for
            let mut used:Vec<String> = Vec::with_capacity( genes.max_id );
            for name in genes.names4sparse.keys() {
                used.push(name.to_string());
            }
            return self.update_names_4_sparse(genes, &used );
        }*/
        [ self.len(), entries ]
    }


    pub fn mtx_counts(&mut self, genes:&IndexedGenes,  min_count:usize, num_threads: usize ) -> String{
        

        if ! self.checked{

            self.passing= 0;

            //println!("Checking cell for min umi count!");

            // here we should firt check if there is some strange similarity over the cells.
            // do we have overlapping gene_id umi connections?
            // should be VERY rare - let's check that!

            // Split keys into chunks and process them in parallel
            let keys = self.keys();
            let chunk_size = keys.len() / num_threads +1; // You need to implement num_threads() based on your requirement

            let results: Vec<(u64, bool)>  = keys
            .par_chunks(chunk_size) 
            .flat_map(|chunk| {
                let min_count = min_count;

                let mut ret= Vec::<(u64, bool)>::with_capacity(chunk_size);
                for key in chunk {
                    if let Some(cell_obj) = &self.get(key) {
                        let n = cell_obj.n_umi( );
                        ret.push( (*key, n >= min_count) );
                    }
                }
                ret
            })
            .collect();

            let mut bad_cells = 0;
            for ( key, passing ) in results{
                if passing{
                    if let Some(cell_obj) = self.get_mut(&key) {
                        cell_obj.passing = passing;
                    }
                }else {
                    bad_cells +=1;
                }
            }

            println!("Dropping cell with too little counts (n={bad_cells})");
            self.keep_only_passing_cells();

            println!("{} cells have passed the cutoff of {} umi counts per cell.\n\n",self.len(), min_count ); 
            self.checked = true;

        }
        
        let ncell_and_entries = self.update_genes_to_print( genes, &self.genes_to_print.clone() );
        self.passing = ncell_and_entries[0];

        let ret = format!("{} {} {}", &self.genes_to_print.len(), ncell_and_entries[0], ncell_and_entries[1] );
        //println!("mtx_counts -> final return: mtx_counts: {}", ret );
        ret
    }

    pub fn passing_cells( &self) -> usize {
        self.passing
    }

    pub fn n_reads( &mut self, genes:&IndexedGenes, names: &Vec<String> ) -> usize {
        let mut count = 0;

        for cell_obj in self.values() {
            count += cell_obj.n_reads( genes, names )
        }
        count
    }

    pub fn read_matrix_market<P: AsRef<Path>>(path: P) -> Result<(Self, IndexedGenes), String> {

        let base = path.as_ref();

        let mut report = MappingInfo::new( None, 0.0, 0 );

        let barcodes_path = base.join("barcodes.tsv.gz");
        let features_path = base.join("features.tsv.gz");
        let matrix_path = base.join("matrix.mtx.gz");

        // Read barcodes
        let barcodes: Vec<u64> = {
            let f = File::open(barcodes_path).map_err(|e| format!("Cannot open barcodes: {e}"))?;
            let reader = BufReader::new(GzDecoder::new(BufReader::new(f)));

            reader
                .lines()
                .collect::<Result<Vec<_>, _>>()
                .map_err(|e| format!("Failed to read barcodes: {e}"))?
                .into_iter()
                .map(|barcode| {
                    println!("The read barcode: {}", barcode);
                    let encoded = IntToStr::new(&barcode).into_u64();
                    encoded
                })
                .collect()
        };


        let f = File::open(features_path).map_err(|e| format!("Cannot open features: {e}"))?;
        let reader_f = BufReader::new(GzDecoder::new(BufReader::new(f)));

        let mut indexed = IndexedGenes::empty(None);
        for line in reader_f.lines() {
            let line = line.unwrap();
            let fields: Vec<&str> = line.split('\t').collect();
            if let Some(gene_id) = fields.get(0) {
                indexed.add(gene_id);
            } else {
                return Err(format!("Malformed gene line: {line:?}"));
            }
        }

        let file = File::open(matrix_path).map_err(|e| format!("Failed to open file: {e}"))?;
        let reader_mtx = BufReader::new(GzDecoder::new(BufReader::new(file)));


        let mut lines = reader_mtx.lines();
        let mut dims = false;
        let mut header = false;

        let mut scdata = Scdata::new(1, MatrixValueType::Integer );

        // Read entries
        for line in lines {
            let line = line.unwrap();
            if ! header {
                header = true;
                let header_parts: Vec<String> = line
                    .trim()
                    .split_whitespace()
                    .map(|s| s.to_string())
                    .collect();
                if header_parts[3].to_lowercase().as_str() == "real" {
                    scdata.value_type=  MatrixValueType::Real
                }
                continue;
            }

            if ! dims {
                dims = true;
                continue;
            }
            println!("This should be a data line? {}", line);
            let mut parts = line.split_whitespace();
            let row = parts.next().ok_or("Missing row")?.parse::<usize>().map_err(|e| e.to_string())?;
            let col = parts.next().ok_or("Missing col")?.parse::<usize>().map_err(|e| e.to_string())?;
            let val = parts.next().ok_or("Missing value")?.parse::<f32>().map_err(|e| e.to_string())?;
            //pub fn try_insert(&mut self, name: &u64, data: GeneUmiHash, value:f32, report: &mut MappingInfo) -> bool {
            match &scdata.value_type {
                &MatrixValueType::Integer => {
                    for umi in 0..val as u64 {
                        scdata.try_insert( &barcodes[col-1], GeneUmiHash( row-1, umi ), 0.0, &mut report );
                    }
                },
                &MatrixValueType::Real => {
                    scdata.try_insert( &barcodes[col-1], GeneUmiHash( row-1, 1 ), val , &mut report );
                },
                _ => unreachable!()
            };
        }

        Ok(( scdata, indexed ))
    }
}






