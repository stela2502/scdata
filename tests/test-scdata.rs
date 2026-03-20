use scdata::{FeatureIndex, GeneUmiHash, MatrixValueType, Scdata};
use mapping_info::MappingInfo;

use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;

/// Minimal feature index for integration tests.
///
/// Feature IDs are just the index in `names` as `u64`.
struct TestFeatureIndex {
    names: Vec<String>,
    name_to_id: HashMap<String, u64>,
}

impl TestFeatureIndex {
    fn new(names: Vec<&str>) -> Self {
        let names: Vec<String> = names.into_iter().map(|s| s.to_string()).collect();
        let name_to_id = names
            .iter()
            .enumerate()
            .map(|(i, name)| (name.clone(), i as u64))
            .collect();

        Self { names, name_to_id }
    }

    fn feature_id(&self, name: &str) -> u64 {
        *self
            .name_to_id
            .get(name)
            .unwrap_or_else(|| panic!("unknown feature name: {name}"))
    }
}

impl FeatureIndex for TestFeatureIndex {
    fn feature_name(&self, feature_id: u64) -> &str {
        &self.names[feature_id as usize]
    }

    fn to_10x_feature_line(&self, feature_id: u64) -> String {
        let name = self.feature_name(feature_id);
        format!("{name}\t{name}\tGene Expression")
    }
}

fn test_out_dir(name: &str) -> PathBuf {
    PathBuf::from("tests").join("out").join(name)
}

#[test]
fn singlecelldata_to_sparse_integer_roundtrip() {
    let mut celldata = Scdata::new(1, MatrixValueType::Integer);
    let mut report = MappingInfo::new(None, 56.0, 10_000);

    let feature_index = TestFeatureIndex::new(vec!["Gene1", "Gene2", "Gene3", "Gene4"]);

    let gene1 = feature_index.feature_id("Gene1");
    let gene3 = feature_index.feature_id("Gene3");
    let gene4 = feature_index.feature_id("Gene4");

    // Cell 1: Gene1 x20, Gene4 x10
    for umi in 0..20 {
        assert!(
            celldata.try_insert(&1_u64, GeneUmiHash(gene1, umi as u64), 0.0, &mut report),
            "insert Gene1 umi {umi}"
        );
    }
    assert!(
        !celldata.try_insert(&1_u64, GeneUmiHash(gene1, 0), 0.0, &mut report),
        "duplicate Gene1 umi 0 should fail"
    );
    for umi in 20..30 {
        assert!(
            celldata.try_insert(&1_u64, GeneUmiHash(gene4, umi as u64), 0.0, &mut report),
            "insert Gene4 umi {umi}"
        );
    }

    // Cell 2: Gene3 x20, Gene1 x20
    for umi in 0..20 {
        assert!(
            celldata.try_insert(
                &13_452_355_u64,
                GeneUmiHash(gene3, umi as u64),
                0.0,
                &mut report
            ),
            "insert Gene3 umi {umi}"
        );
        assert!(
            celldata.try_insert(
                &13_452_355_u64,
                GeneUmiHash(gene1, umi as u64),
                0.0,
                &mut report
            ),
            "insert Gene1 umi {umi}"
        );
    }
    assert!(
        !celldata.try_insert(
            &13_452_355_u64,
            GeneUmiHash(gene1, 0),
            0.0,
            &mut report
        ),
        "duplicate Gene1 umi 0 should fail"
    );

    let out_dir = test_out_dir("integer");
    if out_dir.exists() {
        fs::remove_dir_all(&out_dir).unwrap();
    }

    let write_result = celldata.write_sparse(&out_dir, &feature_index);
    assert!(write_result.is_ok(), "write_sparse failed: {write_result:?}");

    let matrix_path = out_dir.join("matrix.mtx.gz");
    let features_path = out_dir.join("features.tsv.gz");
    let barcodes_path = out_dir.join("barcodes.tsv.gz");

    assert!(matrix_path.exists(), "matrix.mtx.gz not found");
    assert!(features_path.exists(), "features.tsv.gz not found");
    assert!(barcodes_path.exists(), "barcodes.tsv.gz not found");

    // Assumption:
    // read_matrix_market reconstructs feature IDs in the same order as written
    // from features.tsv.gz, so Gene1..Gene4 keep IDs 0..3.
    let scdata2 = Scdata::read_matrix_market(&out_dir).unwrap();

    assert_eq!(
        scdata2
            .get(&13_452_355_u64)
            .unwrap()
            .n_umi_4_gene_id(&gene3),
        20
    );
    assert_eq!(
        scdata2
            .get(&1_u64)
            .unwrap()
            .n_umi_4_gene_id(&gene4),
        10
    );
    assert_eq!(
        scdata2
            .get(&13_452_355_u64)
            .unwrap()
            .n_umi_4_gene_id(&gene1),
        20
    );
    assert_eq!(
        scdata2
            .get(&1_u64)
            .unwrap()
            .n_umi_4_gene_id(&gene1),
        20
    );
}

#[test]
fn singlecelldata_to_sparse_real_roundtrip() {
    let mut celldata = Scdata::new(1, MatrixValueType::Real);
    let mut report = MappingInfo::new(None, 56.0, 10_000);

    let feature_index = TestFeatureIndex::new(vec!["Gene1", "Gene2", "Gene3", "Gene4"]);

    let gene1 = feature_index.feature_id("Gene1");
    let gene3 = feature_index.feature_id("Gene3");
    let gene4 = feature_index.feature_id("Gene4");

    // Cell 1: Gene1 x20, Gene4 x10
    for umi in 0..20 {
        assert!(
            celldata.try_insert(
                &1_u64,
                GeneUmiHash(gene1, umi as u64),
                2.0 + umi as f32 / 20.0,
                &mut report
            ),
            "insert Gene1 umi {umi}"
        );
    }
    assert!(
        !celldata.try_insert(&1_u64, GeneUmiHash(gene1, 0), 2.0, &mut report),
        "duplicate Gene1 umi 0 should fail"
    );
    for umi in 20..30 {
        assert!(
            celldata.try_insert(
                &1_u64,
                GeneUmiHash(gene4, umi as u64),
                1.0 + umi as f32 / 20.0,
                &mut report
            ),
            "insert Gene4 umi {umi}"
        );
    }

    // Cell 2: Gene3 x20, Gene1 x20
    for umi in 0..20 {
        assert!(
            celldata.try_insert(
                &13_452_355_u64,
                GeneUmiHash(gene3, umi as u64),
                2.0 + umi as f32 / 20.0,
                &mut report
            ),
            "insert Gene3 umi {umi}"
        );
        assert!(
            celldata.try_insert(
                &13_452_355_u64,
                GeneUmiHash(gene1, umi as u64),
                3.0 + umi as f32 / 20.0,
                &mut report
            ),
            "insert Gene1 umi {umi}"
        );
    }
    assert!(
        !celldata.try_insert(
            &13_452_355_u64,
            GeneUmiHash(gene1, 0),
            2.0,
            &mut report
        ),
        "duplicate Gene1 umi 0 should fail"
    );

    let out_dir = test_out_dir("real");
    if out_dir.exists() {
        fs::remove_dir_all(&out_dir).unwrap();
    }

    let write_result = celldata.write_sparse(&out_dir, &feature_index);
    assert!(write_result.is_ok(), "write_sparse failed: {write_result:?}");

    let matrix_path = out_dir.join("matrix.mtx.gz");
    let features_path = out_dir.join("features.tsv.gz");
    let barcodes_path = out_dir.join("barcodes.tsv.gz");

    assert!(matrix_path.exists(), "matrix.mtx.gz not found");
    assert!(features_path.exists(), "features.tsv.gz not found");
    assert!(barcodes_path.exists(), "barcodes.tsv.gz not found");

    // Same assumption as above: feature order is preserved on read-back.
    let (scdata2, _read_features) = Scdata::read_matrix_market(&out_dir).unwrap();

    // In real mode, one entry per gene/cell should survive in sparse storage.
    assert_eq!(
        scdata2
            .get(&13_452_355_u64)
            .unwrap()
            .n_umi_4_gene_id(&gene3),
        1
    );
    assert_eq!(
        scdata2
            .get(&1_u64)
            .unwrap()
            .n_umi_4_gene_id(&gene4),
        1
    );

    assert_eq!(
        scdata2
            .get(&13_452_355_u64)
            .unwrap()
            .get_mean_for_gene(&gene3),
        Some(2.475)
    );
    assert_eq!(
        scdata2
            .get(&1_u64)
            .unwrap()
            .get_mean_for_gene(&gene4),
        Some(2.225)
    );
    assert_eq!(
        scdata2
            .get(&13_452_355_u64)
            .unwrap()
            .get_mean_for_gene(&gene1),
        Some(3.475)
    );
}