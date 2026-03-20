// sparse_matrix/io.rs
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use flate2::Compression;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;

use int_to_str::int_to_str::IntToStr;
use mapping_info::MappingInfo;

use crate::cell_data::GeneUmiHash;
use crate::feature_index::FeatureIndex;
use crate::{MatrixValueType, Scdata};

impl Scdata {
    /// Write a dense TSV matrix.
    ///
    /// Rows = cells  
    /// Columns = features
    pub fn write_dense<I: FeatureIndex>(
        &self,
        file_path: &PathBuf,
        index: &I,
    ) -> Result<(), String> {
        if file_path.exists() {
            let _ = fs::remove_file(file_path);
        }

        let file = File::create(file_path)
            .map_err(|err| format!("Error creating dense file {}: {err}", file_path.display()))?;

        let mut writer = BufWriter::new(file);

        // header
        let mut cols = Vec::new();
        cols.push("CellID".to_string());

        for fid in &self.feature_ids_with_data {
            cols.push(index.feature_name(*fid).to_string());
        }

        writeln!(writer, "{}", cols.join("\t"))
            .map_err(|e| format!("Dense header write failed: {e}"))?;

        // rows
        for cell in self.values() {
            if !cell.passing {
                continue;
            }

            let mut row = Vec::with_capacity(cols.len());

            row.push(IntToStr::u8_array_to_str(&cell.name.to_le_bytes()));

            for fid in &self.feature_ids_with_data {
                if let Some(value) = cell.total_reads.get( fid ){
                    row.push( value.to_string() );
                }
                
            }

            writeln!(writer, "{}", row.join("\t"))
                .map_err(|e| format!("Dense row write failed: {e}"))?;
        }

        Ok(())
    }

    /// Write a sparse MatrixMarket matrix (10x compatible layout).
    ///
    /// Produces:
    /// - matrix.mtx.gz
    /// - barcodes.tsv.gz
    /// - features.tsv.gz
    pub fn write_sparse<I: FeatureIndex>(
        &mut self,
        outdir: &PathBuf,
        index: &I,
    ) -> Result<String, String> {
        if !outdir.exists() {
            fs::create_dir_all(outdir)
                .map_err(|e| format!("Could not create output dir {}: {e}", outdir.display()))?;
        }

        let matrix_path = outdir.join("matrix.mtx.gz");
        let barcodes_path = outdir.join("barcodes.tsv.gz");
        let features_path = outdir.join("features.tsv.gz");

        let matrix_file = File::create(&matrix_path)
            .map_err(|e| format!("Could not create {}: {e}", matrix_path.display()))?;

        let barcodes_file = File::create(&barcodes_path)
            .map_err(|e| format!("Could not create {}: {e}", barcodes_path.display()))?;

        let features_file = File::create(&features_path)
            .map_err(|e| format!("Could not create {}: {e}", features_path.display()))?;

        let mut writer_m = BufWriter::new(GzEncoder::new(matrix_file, Compression::default()));
        let mut writer_b = BufWriter::new(GzEncoder::new(barcodes_file, Compression::default()));
        let mut writer_f = BufWriter::new(GzEncoder::new(features_file, Compression::default()));

        // write features.tsv
        for fid in &self.feature_ids_with_data {
            let line = index.to_10x_feature_line(*fid);
            writeln!(writer_f, "{line}")
                .map_err(|e| format!("Writing features.tsv.gz failed: {e}"))?;
        }

        self.keep_only_passing_cells();

        let mut entries = 0usize;
        let mut cell_id = 0usize;
        let value_type = self.value_type(None).clone();

        // header placeholder
        writeln!(
            writer_m,
            "%%MatrixMarket matrix coordinate {} general",
            value_type.as_matrix_market_str() 
            )
        .map_err(|e| format!("Matrix header write failed: {e}"))?;

        writeln!(
            writer_m,
            "{} {} {}",
            self.feature_ids_with_data.len(),
            self.passing_cells(),
            self.total_feature_data_entries
        )
        .map_err(|e| format!("Matrix dims write failed: {e}"))?;

        for cell in self.values() {

            cell_id += 1;

            let cell_name = IntToStr::u8_array_to_str(&cell.name.to_le_bytes());

            writeln!(writer_b, "{cell_name}")
                .map_err(|e| format!("Barcode write failed: {e}"))?;

			for (row_idx, feature_id) in self.feature_ids_with_data.iter().enumerate() {
			    if let Some(value) = cell.total_reads.get(feature_id) {
			        if *value > 0.0 {
			            writeln!(writer_m, "{} {} {}", row_idx + 1, cell_id, value_type.format_matrix_market_value(*value)? )
			                .map_err(|e| format!("Matrix entry write failed: {e}"))?;

			            entries += 1;
			        }
			    }
			}
        }

        assert_eq!(
		    entries,
		    self.total_feature_data_entries,
		    "Sparse export mismatch: wrote {} entries but expected {}",
		    entries,
		    self.total_feature_data_entries
		);

        Ok(format!(
            "sparse Matrix: {} cells, {} features, {} entries written to {}",
            cell_id,
            self.feature_ids_with_data.len(),
            entries,
            outdir.display()
        ))
    }

    /// Read a 10x-style MatrixMarket directory.
    ///
    /// Returns `(Scdata, FeatureIndex-like data)`
    /// for compatibility with previous pipelines.
    pub fn read_matrix_market<P: AsRef<Path>>(
        path: P,
    ) -> Result<Self, String> {
        let base = path.as_ref();

        let barcodes_path = base.join("barcodes.tsv.gz");
        let matrix_path = base.join("matrix.mtx.gz");

        let mut report = MappingInfo::new(None, 0.0, 0);

        let barcodes: Vec<u64> = {
            let f = File::open(&barcodes_path)
                .map_err(|e| format!("Cannot open barcodes: {e}"))?;

            let reader = BufReader::new(GzDecoder::new(BufReader::new(f)));

            reader
                .lines()
                .collect::<Result<Vec<_>, _>>()
                .map_err(|e| format!("Reading barcodes failed: {e}"))?
                .into_iter()
                .map(|barcode| IntToStr::new(&barcode).into_u64())
                .collect()
        };

        let file = File::open(&matrix_path)
            .map_err(|e| format!("Failed to open matrix file: {e}"))?;

        let reader = BufReader::new(GzDecoder::new(BufReader::new(file)));

        let mut lines = reader.lines();

        let mut scdata = Scdata::default();

        let mut header_seen = false;
        let mut dims_seen = false;

        for line in lines.by_ref() {
            let line = line.map_err(|e| format!("Matrix read error: {e}"))?;

            if line.trim().is_empty() {
                continue;
            }

            if !header_seen {
                header_seen = true;
                continue;
            }

            if line.starts_with('%') {
                continue;
            }

            if !dims_seen {
                dims_seen = true;
                continue;
            }

            let mut parts = line.split_whitespace();

            let row = parts
                .next()
                .ok_or("Missing row")?
                .parse::<usize>()
                .map_err(|e| e.to_string())?;

            let col = parts
                .next()
                .ok_or("Missing col")?
                .parse::<usize>()
                .map_err(|e| e.to_string())?;

            let val = parts
                .next()
                .ok_or("Missing value")?
                .parse::<f32>()
                .map_err(|e| e.to_string())?;

            for umi in 0..val as u64 {
                scdata.try_insert(
                    &barcodes[col - 1],
                    GeneUmiHash((row - 1) as u64 , umi),
                    0.0,
                    &mut report,
                );
            }
        }

        Ok(scdata)
    }
}