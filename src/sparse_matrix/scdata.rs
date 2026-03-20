use std::collections::{HashMap, HashSet};

use core::fmt;
use rayon::prelude::*;

use mapping_info::MappingInfo;

use crate::cell_data::GeneUmiHash;
use crate::cell_data::CellData;
use crate::MatrixValueType;

/// Sparse single-cell count store.
///
/// Cells are partitioned into 256 buckets by the top byte of the cell id.
/// Each bucket stores `CellData` objects keyed by the full cell id.
pub struct Scdata {
    /// 256 buckets, indexed by top 8 bits of the cell id.
    pub(crate) data: [HashMap<u64, CellData>; u8::MAX as usize + 1],

    /// Cached set of feature ids currently observed in retained cells.
    pub(crate) feature_ids_with_data: Vec<u64>,

    /// Cached number of stored feature entries across all cells
    pub(crate) total_feature_data_entries: usize,

    /// True once cell filtering has been evaluated.
    checked: bool,

    /// Number of currently passing cells.
    passing: usize,

    /// Preferred thread count for internal parallel helpers.
    pub num_threads: usize,

    /// Numeric matrix value type.
    pub(crate) value_type: MatrixValueType,
}

impl Default for Scdata {
    /// Create an integer-valued matrix with one default worker thread.
    fn default() -> Self {
        Self::new(1, MatrixValueType::Integer)
    }
}

impl fmt::Display for Scdata {
    /// Print a short structural summary of the current matrix state.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "Scdata Summary")?;
        writeln!(f, "===============")?;
        writeln!(f, "Checked: {}", self.checked)?;
        writeln!(f, "Passing Cells: {}", self.passing)?;
        writeln!(f, "Matrix Value Type: {:?}", self.value_type)?;
        writeln!(
            f,
            "Features with Data: [{} entries]",
            self.feature_ids_with_data.len()
        )?;

        let total_cells: usize = self.data.iter().map(|m| m.len()).sum();
        let non_empty_buckets = self.data.iter().filter(|m| !m.is_empty()).count();

        writeln!(f, "Total Cells: {}", total_cells)?;
        writeln!(f, "Non-empty Buckets: {}", non_empty_buckets)?;
        Ok(())
    }
}

impl Scdata {
    /// Create a new empty sparse matrix store.
    pub fn new(num_threads: usize, value_type: MatrixValueType) -> Self {
        let data = std::array::from_fn(|_| HashMap::<u64, CellData>::new());

        Self {
            data,
            feature_ids_with_data: Vec::new(),
            total_feature_data_entries: 0,
            checked: false,
            passing: 0,
            num_threads,
            value_type,
        }
    }

    /// Get the current matrix value type, or update it if a new one is provided.
    pub fn value_type(&mut self, value_type: Option<MatrixValueType>) -> &MatrixValueType {
        if let Some(val) = value_type {
            self.value_type = val;
        }
        &self.value_type
    }

    /// Map a cell id to one of the 256 storage buckets.
    #[inline]
    fn to_key(&self, name: &u64) -> usize {
        (*name >> 56) as usize
    }

    /// Invalidate cached filtering/export-derived state after matrix mutation.
    fn invalidate_export_cache(&mut self) {
        self.checked = false;
        self.passing = 0;
    }

    /// Iterate over all stored cells across all buckets.
    pub(crate) fn values(&self) -> impl Iterator<Item = &CellData> {
        self.data.iter().flat_map(|map| map.values())
    }

    /// Return all currently stored cell ids.
    pub fn keys(&self) -> Vec<u64> {
        self.data
            .iter()
            .flat_map(|map| map.keys())
            .copied()
            .collect()
    }

    /// Check whether no cells are currently stored.
    pub fn is_empty(&self) -> bool {
        self.data.iter().all(|m| m.is_empty())
    }

    /// Return the total number of stored cells.
    pub fn len(&self) -> usize {
        self.data.iter().map(|m| m.len()).sum()
    }

    /// Get mutable access to one cell by id.
    fn get_mut(&mut self, key: &u64) -> Option<&mut CellData> {
        let index = self.to_key(key);
        self.data[index].get_mut(key)
    }

    /// Get read-only access to one cell by id.
    pub fn get(&self, key: &u64) -> Option<&CellData> {
        let index = self.to_key(key);
        self.data[index].get(key)
    }

    /// Remove all cells currently marked as failing.
    pub fn keep_only_passing_cells(&mut self) {
        for map in &mut self.data {
            map.retain(|_, cell_data| cell_data.passing);
        }
        self.rebuild_feature_ids_with_data();
    }

    /// Merge another matrix into this one in parallel, bucket by bucket.
    pub fn merge(&mut self, other: &Scdata) {
        if other.is_empty() {
            return;
        }

        self.invalidate_export_cache();

        self.data
            .par_iter_mut()
            .enumerate()
            .for_each(|(index, self_bucket)| {
                if let Some(other_bucket) = other.data.get(index) {
                    for (cell_name, other_cell) in other_bucket {
                        match self_bucket.entry(*cell_name) {
                            std::collections::hash_map::Entry::Occupied(mut entry) => {
                                entry.get_mut().merge(other_cell);
                            }
                            std::collections::hash_map::Entry::Vacant(entry) => {
                                entry.insert(other_cell.clone());
                            }
                        }
                    }
                }
            });

        self.rebuild_feature_ids_with_data();
    }

    /// Merge another matrix into this one by consuming the source in a single thread.
    pub fn merge_single_thread(&mut self, mut other: Scdata) {
        if other.is_empty() {
            return;
        }

        self.invalidate_export_cache();

        for other_cell in other.data.iter_mut().flat_map(|map| map.values_mut()) {
            let index = self.to_key(&other_cell.name);
            match self.data[index].entry(other_cell.name) {
                std::collections::hash_map::Entry::Occupied(mut entry) => {
                    entry.get_mut().merge(other_cell);
                }
                std::collections::hash_map::Entry::Vacant(entry) => {
                    entry.insert(std::mem::take(other_cell));
                }
            }
        }

        self.rebuild_feature_ids_with_data();
    }

    /// Insert one unique feature/UMI observation for a cell.
    ///
    /// Returns `false` if the exact feature/UMI pair already existed.
    pub fn try_insert(
        &mut self,
        name: &u64,
        data: GeneUmiHash,
        _value: f32,
        report: &mut MappingInfo,
    ) -> bool {
        let index = self.to_key(name);
        self.invalidate_export_cache();

        let cell_info = self.data[index]
            .entry(*name)
            .or_insert_with(|| CellData::new(*name));

        report.ok_reads += 1;

        if !cell_info.add(data) {
            report.pcr_duplicates += 1;
            report.local_dup += 1;
            false
        } else{
        	true
        }
    }

    /// Insert or accumulate a numeric value for one feature/UMI observation.
    pub fn try_insert_value(
        &mut self,
        name: &u64,
        data: GeneUmiHash,
        value: f32,
        report: &mut MappingInfo,
    ) -> bool {
        let index = self.to_key(name);
        self.invalidate_export_cache();

        let cell_info = self.data[index]
            .entry(*name)
            .or_insert_with(|| CellData::new(*name));

        report.ok_reads += 1;

        cell_info.add_value(data, value)

    }

    /// Return the number of cells currently marked as passing.
    pub fn passing_cells(&self) -> usize {
        self.passing
    }

    /// Return all currently retained cell ids.
	pub fn cell_ids(&self) -> HashSet<u64> {
	    self.data
	        .iter()
	        .flat_map(|bucket| bucket.keys().copied())
	        .collect()
	}


    /// Return the set of cell ids with at least `min_count` UMIs.
    pub fn passing_cell_set_by_umi(&self, min_count: usize) -> HashSet<u64> {
        let keys = self.keys();
        if keys.is_empty() {
            return HashSet::new();
        }

        let n_threads = self.num_threads.max(1);
        let chunk_size = keys.len() / n_threads + 1;

        let passing: Vec<u64> = keys
            .par_chunks(chunk_size)
            .flat_map(|chunk| {
                let mut keep = Vec::<u64>::with_capacity(chunk.len());
                for key in chunk {
                    if let Some(cell) = self.get(key) {
                        if cell.n_umi() >= min_count {
                            keep.push(*key);
                        }
                    }
                }
                keep
            })
            .collect();

        passing.into_iter().collect()
    }

	/// Rebuild the cached set of feature ids observed in all cells
	/// and count the total number of stored sparse matrix entries.
	///
	/// The feature id cache is collected from the per-cell `total_reads` keys.
	/// The total number of exported sparse entries is the sum of
	/// `cell.total_reads.len()` across all cells.
	fn rebuild_feature_ids_with_data(&mut self) {
	    let (feature_ids, total_entries) = self
	        .data
	        .par_iter()
	        .map(|bucket| {
	            let mut local_ids = HashSet::with_capacity(64);
	            let mut local_entries = 0usize;

	            for cell in bucket.values() {
	                for feature_id in cell.total_reads.keys() {
	                    local_ids.insert(*feature_id);
	                }

	                local_entries += cell.total_reads.len();
	            }

	            (local_ids, local_entries)
	        })
	        .reduce(
	            || (HashSet::new(), 0usize),
	            |(mut ids_a, n_a), (ids_b, n_b)| {
	                ids_a.extend(ids_b);
	                (ids_a, n_a + n_b)
	            },
	        );

	    let mut feature_ids: Vec<u64> = feature_ids.into_iter().collect();

	    // Ensure deterministic export order.
	    feature_ids.sort_unstable();

	    self.feature_ids_with_data = feature_ids;
	    self.total_feature_data_entries = total_entries;
	}

    /// Restrict the matrix to cell witzh more than min_umis unique reads.
	pub fn retain_cells_with_min_umis(&mut self, min_umis: usize) {
	    for bucket in &mut self.data {
	        bucket.retain(|_, cell| cell.seen.len() >= min_umis);
	    }

	    // after filtering the export cache must be recomputed
	    self.rebuild_feature_ids_with_data();
	}

    /// Restrict the matrix to a predefined set of cell ids and mark them as passing.
    pub fn restrict_to_cells(&mut self, keep: &HashSet<u64>) {
        self.passing = 0;

        for bucket in &mut self.data {
            bucket.retain(|cell_id, _| keep.contains(cell_id));
        }

        for cell in self.data.iter_mut().flat_map(|m| m.values_mut()) {
            cell.passing = true;
            self.passing += 1;
        }

        self.rebuild_feature_ids_with_data();
        self.checked = true;
    }
}