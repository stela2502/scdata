// src/lib.rs

pub mod sparse_matrix;
pub mod ambient_rna_detect;
pub mod cell_data;
pub mod feature_index;

pub use crate::sparse_matrix::{MatrixValueType, Scdata};
pub use crate::cell_data::GeneUmiHash;
pub use crate::feature_index::FeatureIndex;