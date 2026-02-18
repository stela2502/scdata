// /Scdata/mod.rs


pub mod indexed_genes;
pub mod scdata;
pub mod ambient_rna_detect;
pub mod cell_data;


pub use crate::scdata::Scdata as Scdata;
pub use crate::cell_data::CellData  as CellData;
pub use crate::ambient_rna_detect::AmbientRnaDetect  as AmbientRnaDetect;
pub use crate::indexed_genes::IndexedGenes as IndexedGenes;
pub use crate::scdata::MatrixValueType;
