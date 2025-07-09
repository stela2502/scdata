// /Scdata/mod.rs


pub mod indexed_genes;
pub mod Scdata;
pub mod ambient_rna_detect;
pub mod cell_data;


pub use Scdata::Scdata as Scdata;
pub use crate::scdata::cell_data::CellData  as CellData;
pub use crate::scdata::ambient_rna_detect::AmbientRnaDetect  as AmbientRnaDetect;
pub use crate::scdata::indexed_genes::IndexedGenes as IndexedGenes;