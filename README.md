# scdata

[![Rust](https://github.com/stela2502/scdata/actions/workflows/rust.yml/badge.svg?branch=main)](https://github.com/stela2502/scdata/actions/workflows/rust.yml)

A lightweight Rust library for handling **sparse single-cell expression
data**, with support for incremental insertion, merging, and Matrix
Market import/export.

This crate is designed for high-performance pipelines where you want: -
Fine-grained control over sparse data - Thread-aware data insertion -
Explicit gene indexing - Matrix Market compatibility - Clean Rust-native
APIs

------------------------------------------------------------------------

## Overview

`scdata` provides the following main components:

### `Scdata`

Core sparse container for single-cell count data.

### `IndexedGenes`

Gene name ↔ index mapping.

### `CellData`

Internal per-cell sparse representation.

### `AmbientRnaDetect`

Utilities for detecting ambient RNA contamination (experimental).

### `MatrixValueType`

Specifies the Matrix Market value type: - `Integer` - `Real` -
`Complex` - `Pattern` - `Unknown(String)`

------------------------------------------------------------------------

# Design Philosophy

- Sparse-first architecture
- Incremental insertion (`try_insert`)
- Explicit gene indexing
- Controlled merging
- Minimal dependencies
- Matrix Market interoperability

------------------------------------------------------------------------

# Installation

Add to your `Cargo.toml`:

```toml
scdata = { git = "https://github.com/stela2502/scdata" }
```

------------------------------------------------------------------------

# Basic Usage

## 1️⃣ Create Gene Index

```rust
use scdata::IndexedGenes;

let genes = IndexedGenes::from_names(vec![
    "GeneA",
    "GeneB",
    "GeneC"
]);
```

This creates a stable mapping between gene names and internal indices.

------------------------------------------------------------------------

## 2️⃣ Create a Sparse Dataset

```rust
use scdata::{Scdata, MatrixValueType};

let mut data = Scdata::new(
    4,                         // number of threads
    MatrixValueType::Integer   // matrix type
);
```

------------------------------------------------------------------------

## 3️⃣ Insert Counts

```rust
use scdata::Scdata;
use mapping_info::MappingInfo; // example external report struct

let mut report = MappingInfo::default();

let cell_id: u64 = 123456;
let gene_hash: GeneUmiHash = ...; // your gene/UMI hash structure

data.try_insert(
    &cell_id,
    gene_hash,
    1.0,
    &mut report
);
```

`try_insert`: - Inserts a value into the sparse structure - Returns
`true` if insertion succeeded - Updates external `MappingInfo`

------------------------------------------------------------------------

## 4️⃣ Merge Datasets

```rust
data.merge(&other_dataset);
```

This merges sparse cell structures efficiently.

------------------------------------------------------------------------

## 5️⃣ Export Data

### Write filtered matrix

```rust
use std::path::PathBuf;

let output = PathBuf::from("matrix.mtx");

data.write(
    &output,
    &genes,
    10    // min_count threshold
).unwrap();
```

### Write sparse format explicitly

```rust
data.write_sparse(
    &output,
    &genes,
    10
).unwrap();
```

------------------------------------------------------------------------

## 6️⃣ Read Matrix Market

```rust
let (data, genes) =
    Scdata::read_matrix_market("matrix.mtx").unwrap();
```

This loads both: - The sparse matrix - The gene index

------------------------------------------------------------------------

# Internal Structure

```text
Scdata
 ├── data: [BTreeMap<u64, CellData>; 255]
 ├── genes_with_data: HashSet<usize>
 ├── num_threads
 └── value_type
```

Key design features:

- Per-thread storage buckets
- BTreeMap for deterministic ordering
- Gene-level tracking via `genes_with_data`
- Lazy validation (`checked` flag)

------------------------------------------------------------------------

# Matrix Market Compatibility

Supports standard `.mtx` sparse format.

Typical 10x-style layout:

    matrix.mtx
    genes.tsv
    barcodes.tsv

You can integrate this library into pipelines that:

- Process CellRanger outputs
- Produce Scanpy-compatible matrices
- Export to downstream R/Python workflows

------------------------------------------------------------------------

# Threading Model

- `num_threads` controls internal parallelization
- Designed for deterministic merge behavior
- No implicit global locks

------------------------------------------------------------------------

# When To Use scdata

✔ Building custom RNA-seq pipelines in Rust\
✔ Integrating single-cell logic into Rust tools\
✔ Avoiding Python/R for production pipelines\
✔ High-performance sparse counting\
✔ Custom filtering logic before export

------------------------------------------------------------------------

# Roadmap Ideas

- Whitelist-based sample detection
- Improved ambient RNA modeling
- HDF5 export
- 10x-compatible directory export
- Compression support

------------------------------------------------------------------------

# License

MIT or Apache-2.0 (depending on repository settings)

------------------------------------------------------------------------

# Author

Stefan Lang\
Bioinformatics & Single-Cell Systems
