use crate::gvcf_parser::GVcfRecordIterator;
use pyo3::prelude::*;
use pyo3::types::PyModule;
use std::collections::HashMap;

#[pyfunction]
pub fn collect_vars_positions(path: &str) -> PyResult<(Vec<u32>, Vec<u32>, Vec<u32>, Vec<String>)> {
    let mut chrom_map = HashMap::new();
    let mut id_to_chrom = Vec::new();
    let mut next_id = 0u32;

    let mut chrom_ids = Vec::new();
    let mut positions = Vec::new();
    let mut widths = Vec::new();

    let iter = GVcfRecordIterator::from_gzip_path(path)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;

    for r in iter {
        let rec = r.map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
        let (start, end) = rec
            .get_span()
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
        let chrom_id = match chrom_map.get(&rec.chrom) {
            Some(&id) => id,
            None => {
                chrom_map.insert(rec.chrom.clone(), next_id);
                id_to_chrom.push(rec.chrom.clone());
                next_id += 1;
                next_id - 1
            }
        };
        chrom_ids.push(chrom_id);
        positions.push(start);
        widths.push(end - start + 1);
    }

    Ok((chrom_ids, positions, widths, id_to_chrom))
}

#[pymodule]
fn gvcfparser(_py: Python<'_>, m: Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(collect_vars_positions, &m)?)?;
    Ok(())
}
