#[cfg(test)]
mod tests {
    use scdata::Scdata;
	use scdata::cell_data::GeneUmiHash;
    use scdata::IndexedGenes;
    //use rustody::singlecelldata::IndexedGenes;
    use mapping_info::MappingInfo;
    use std::collections::BTreeMap;
    use std::fs;
    use std::path::PathBuf;
    use flate2::bufread::GzDecoder;
    use std::io::BufReader;
    use std::io::BufRead;
    use scdata::scdata::MatrixValueType;

	static EMPTY_VEC: Vec<String> = Vec::new();


    #[test]
    fn singlecelldata_to_sparse() {
        let mut celldata = Scdata::new( 1, MatrixValueType::Integer ); // only one thread here

        let mut report = MappingInfo::new(None, 56.0, 10000 );

        // add one cell with two genes and each 20 umi counts

        for umi in 0..20 {
        	assert!( celldata.try_insert( &1_u64, GeneUmiHash(0, umi as u64), 0.0, &mut report), "I add Gene1 (0) umi {umi}" );
        }
        assert!( ! celldata.try_insert( &1_u64, GeneUmiHash(0, 0 as u64), 0.0, &mut report ),"I add Gene1 (0) umi 0 - again - and failed" );
        for umi in 20..30 {
            assert!( celldata.try_insert( &1_u64, GeneUmiHash(3, umi as u64), 0.0, &mut report ),"I add Gene4 (3) umi {umi} ");
        }

        // add one other cell with one gene and 20 umi counts

        for umi in 0..20 {
        	assert!( celldata.try_insert( &13452355_u64, GeneUmiHash(2, umi as u64), 0.0, &mut report), "I add Gene3 (2) umi {umi}" );
        	assert!( celldata.try_insert( &13452355_u64, GeneUmiHash(0, umi as u64), 0.0, &mut report), "I add Gene1 (0) umi {umi}" );
        }
        assert!( ! celldata.try_insert( &13452355_u64, GeneUmiHash(0, 0 as u64), 0.0, &mut report ),"I add Gene1 umi 0 - again - and failed" );

        let mut genes = BTreeMap::<String, usize>::new();
        genes.insert("Gene1".to_string(), 0 );
        genes.insert("Gene2".to_string(), 1 );
        genes.insert("Gene3".to_string(), 2 );
        genes.insert("Gene4".to_string(), 3 );

        let indexed_genes = IndexedGenes::new( &genes, 0);

        //to_str<'live>(&mut self, gene_info:&GeneIds, names: &Vec<String> ) 
        let  names= vec!("Gene1".to_string(), "Gene3".to_string(), "Gene4".to_string() );
        // this string counts: genes, cell, lines
        let  exp2:String = "3 2 4".to_string();
        celldata.update_genes_to_print( &indexed_genes , &names);
        let  val = celldata.mtx_counts( &indexed_genes, 1, 1 );

        assert_eq!( val,  exp2 );

        let out_dir = PathBuf::from("tests/out/integer");
        // Clean output dir before test
        if out_dir.exists() {
            fs::remove_dir_all(&out_dir).unwrap();
        }

        let result = celldata.write_sparse( &out_dir, &indexed_genes, 1 );
        
        // Paths to output files
        let matrix_path = out_dir.join("matrix.mtx.gz");
        let features_path = out_dir.join("features.tsv.gz");
        let barcodes_path = out_dir.join("barcodes.tsv.gz");

        // Validate files exist
        assert!(matrix_path.exists(), "matrix.mtx.gz not found");
        assert!(features_path.exists(), "features.tsv.gz not found");
        assert!(barcodes_path.exists(), "barcodes.tsv.gz not found");

        // Validate barcodes.tsv.gz
        let (scdata, mut indexed_genes) = Scdata::read_matrix_market( out_dir ).unwrap();

        eprintln!("the cell object: {}", scdata );
        eprintln!("the cell object cell ids: {:?}", scdata.keys() );

        assert_eq!( scdata.get( &13452355_u64 ).unwrap().n_umi_4_gene_id( &indexed_genes.get_gene_id("Gene3") ), 20);
        assert_eq!( scdata.get( &1_u64 ).unwrap().n_umi_4_gene_id( &indexed_genes.get_gene_id("Gene4") ), 10);

    }

    #[test]
    fn singlecelldata_to_sparse_real() {
        let mut celldata = Scdata::new( 1, MatrixValueType::Real ); // only one thread here

        let mut report = MappingInfo::new(None, 56.0, 10000 );

        // add one cell with two genes and each 20 umi counts

        for umi in 0..20 {
            assert!( celldata.try_insert( &1_u64, GeneUmiHash(0, umi as u64), 2.0+ umi as f32 / 20.0 , &mut report), "I add Gene1 (0) umi {umi}" );
        }
        assert!( ! celldata.try_insert( &1_u64, GeneUmiHash(0, 0 as u64), 2.0, &mut report ),"I add Gene1 (0) umi 0 - again - and failed" );
        for umi in 20..30 {
            assert!( celldata.try_insert( &1_u64, GeneUmiHash(3, umi as u64), 1.0+ umi as f32 / 20.0, &mut report ),"I add Gene4 (3) umi {umi} ");
        }

        // add one other cell with one gene and 20 umi counts

        for umi in 0..20 {
            assert!( celldata.try_insert( &13452355_u64, GeneUmiHash(2, umi as u64), 2.0+ umi as f32 / 20.0, &mut report), "I add Gene3 (2) umi {umi}" );
            assert!( celldata.try_insert( &13452355_u64, GeneUmiHash(0, umi as u64), 3.0+ umi as f32 / 20.0, &mut report), "I add Gene1 (0) umi {umi}" );
        }
        assert!( ! celldata.try_insert( &13452355_u64, GeneUmiHash(0, 0 as u64), 2.0 , &mut report ),"I add Gene1 umi 0 - again - and failed" );

        let mut genes = BTreeMap::<String, usize>::new();
        genes.insert("Gene1".to_string(), 0 );
        genes.insert("Gene2".to_string(), 1 );
        genes.insert("Gene3".to_string(), 2 );
        genes.insert("Gene4".to_string(), 3 );

        let indexed_genes = IndexedGenes::new( &genes, 0);

        //to_str<'live>(&mut self, gene_info:&GeneIds, names: &Vec<String> ) 
        let  names= vec!("Gene1".to_string(), "Gene3".to_string(), "Gene4".to_string() );
        // this string counts: genes, cell, lines
        let  exp2:String = "3 2 4".to_string();
        celldata.update_genes_to_print( &indexed_genes , &names);
        let  val = celldata.mtx_counts( &indexed_genes, 1, 1 );

        assert_eq!( val,  exp2 );

        let out_dir = PathBuf::from("tests/out/real");
        // Clean output dir before test
        if out_dir.exists() {
            fs::remove_dir_all(&out_dir).unwrap();
        }

        let result = celldata.write_sparse( &out_dir, &indexed_genes, 1 );
        
        // Paths to output files
        let matrix_path = out_dir.join("matrix.mtx.gz");
        let features_path = out_dir.join("features.tsv.gz");
        let barcodes_path = out_dir.join("barcodes.tsv.gz");

        // Validate files exist
        assert!(matrix_path.exists(), "matrix.mtx.gz not found");
        assert!(features_path.exists(), "features.tsv.gz not found");
        assert!(barcodes_path.exists(), "barcodes.tsv.gz not found");

        // Validate barcodes.tsv.gz
        let (scdata, mut indexed_genes) = Scdata::read_matrix_market( out_dir ).unwrap();

        eprintln!("the cell object: {}", scdata );
        eprintln!("the cell object cell ids: {:?}", scdata.keys() );

        assert_eq!( scdata.get( &13452355_u64 ).unwrap().n_umi_4_gene_id( &indexed_genes.get_gene_id("Gene3") ), 1);
        assert_eq!( scdata.get( &1_u64 ).unwrap().n_umi_4_gene_id( &indexed_genes.get_gene_id("Gene4") ), 1);

        assert_eq!( scdata.get( &13452355_u64 ).unwrap().get_mean_for_gene( &indexed_genes.get_gene_id("Gene3") ), Some(2.475));
        assert_eq!( scdata.get( &1_u64 ).unwrap().get_mean_for_gene( &indexed_genes.get_gene_id("Gene4") ), Some(2.225));

    }
        use std::collections::HashSet;

    #[test]
    fn passing_cell_set_by_umi_parallel_matches_expected() {
        let mut celldata = Scdata::new(1, MatrixValueType::Integer); // single thread here
        let mut report = MappingInfo::new(None, 56.0, 10000);

        // Cell 1: 30 UMIs total (20 for gene0 + 10 for gene3)
        for umi in 0..20 {
            assert!(celldata.try_insert(&1_u64, GeneUmiHash(0, umi as u64), 0.0, &mut report));
        }
        for umi in 20..30 {
            assert!(celldata.try_insert(&1_u64, GeneUmiHash(3, umi as u64), 0.0, &mut report));
        }

        // Cell 2: 40 UMIs total (20 for gene2 + 20 for gene0)
        for umi in 0..20 {
            assert!(celldata.try_insert(&13452355_u64, GeneUmiHash(2, umi as u64), 0.0, &mut report));
            assert!(celldata.try_insert(&13452355_u64, GeneUmiHash(0, umi as u64), 0.0, &mut report));
        }

        // Threshold 1: both pass
        let passing = celldata.passing_cell_set_by_umi(1);
        assert!(passing.contains(&1_u64));
        assert!(passing.contains(&13452355_u64));
        assert_eq!(passing.len(), 2);

        // Threshold 31: only cell2 passes (40 >= 31, 30 < 31)
        let passing = celldata.passing_cell_set_by_umi(31);
        assert!(!passing.contains(&1_u64));
        assert!(passing.contains(&13452355_u64));
        assert_eq!(passing.len(), 1);

        // Threshold 41: none pass
        let passing = celldata.passing_cell_set_by_umi(41);
        assert!(passing.is_empty());

        // IMPORTANT: this function must not mutate Scdata
        assert!(celldata.get(&1_u64).is_some());
        assert!(celldata.get(&13452355_u64).is_some());
    }

    #[test]
    fn restrict_to_cells_removes_others_and_invalidates_cached_state() {
        let mut celldata = Scdata::new(1, MatrixValueType::Integer);
        let mut report = MappingInfo::new(None, 56.0, 10000);

        // Populate the same two cells as in the other tests.
        for umi in 0..20 {
            assert!(celldata.try_insert(&1_u64, GeneUmiHash(0, umi as u64), 0.0, &mut report));
        }
        for umi in 20..30 {
            assert!(celldata.try_insert(&1_u64, GeneUmiHash(3, umi as u64), 0.0, &mut report));
        }

        for umi in 0..20 {
            assert!(celldata.try_insert(&13452355_u64, GeneUmiHash(2, umi as u64), 0.0, &mut report));
            assert!(celldata.try_insert(&13452355_u64, GeneUmiHash(0, umi as u64), 0.0, &mut report));
        }

        // Build whitelist with only cell 1
        let mut keep = HashSet::new();
        keep.insert(1_u64);

        celldata.restrict_to_cells(&keep);

        // Only whitelisted cell remains
        let keys = celldata.keys();
        assert_eq!(keys.len(), 1);
        assert_eq!(keys[0], 1_u64);

        assert!(celldata.get(&1_u64).is_some());
        assert!(celldata.get(&13452355_u64).is_none());

        // Data integrity: the remaining cell still has its original UMIs
        // gene0: 20, gene3: 10
        assert_eq!(celldata.get(&1_u64).unwrap().n_umi_4_gene_id(&0), 20);
        assert_eq!(celldata.get(&1_u64).unwrap().n_umi_4_gene_id(&3), 10);

        // And: after restriction, the passing set at different thresholds behaves as expected
        let passing_1 = celldata.passing_cell_set_by_umi(1);
        assert_eq!(passing_1.len(), 1);
        assert!(passing_1.contains(&1_u64));

        let passing_31 = celldata.passing_cell_set_by_umi(31);
        assert!(passing_31.is_empty()); // remaining cell has only 30 UMIs
    }
    
 }