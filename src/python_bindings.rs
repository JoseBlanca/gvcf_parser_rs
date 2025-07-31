use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyBytes;

use arrow2::array::{Array, UInt32Array, Utf8Array};
use arrow2::chunk::Chunk;
use arrow2::datatypes::{DataType, Field, Schema};
use arrow2::io::ipc::write::FileWriter;
use arrow2::io::ipc::write::WriteOptions;
use std::io::Cursor;

use crate::errors::VcfParseError;
use crate::gvcf_parser::{GVcfRecord, GVcfRecordIterator, VcfResult};

pub fn collect_variant_coords_as_arrow<I>(iter: I) -> VcfResult<(Schema, Chunk<Box<dyn Array>>)>
where
    I: Iterator<Item = VcfResult<GVcfRecord>>,
{
    let mut chroms = Vec::new();
    let mut positions = Vec::new();
    let mut widths = Vec::new();

    for result in iter {
        match result {
            Ok(rec) => {
                let (start, end) = rec.get_span()?;
                chroms.push(rec.chrom);
                positions.push(start);
                widths.push(end - start + 1);
            }
            Err(VcfParseError::InvariantgVCFLine) => continue,
            Err(e) => return Err(e),
        }
    }

    let arrays: Vec<Box<dyn Array>> = vec![
        Box::new(Utf8Array::<i32>::from_slice(chroms)),
        Box::new(UInt32Array::from_slice(positions)),
        Box::new(UInt32Array::from_slice(widths)),
    ];

    let schema = Schema::from(vec![
        Field::new("chroms", DataType::Utf8, false),
        Field::new("positions", DataType::UInt32, false),
        Field::new("var_widths", DataType::UInt32, false),
    ]);

    Ok((schema, Chunk::new(arrays)))
}

#[pyfunction]
pub fn export_arrow_ipc(path: &str) -> PyResult<Py<PyBytes>> {
    let iter = GVcfRecordIterator::from_gzip_path(path)
        .map_err(|e| PyValueError::new_err(format!("Failed to open GVCF: {e}")))?;

    let (schema, chunk) = collect_variant_coords_as_arrow(iter)
        .map_err(|e| PyValueError::new_err(format!("Error collecting records: {e}")))?;

    let buffer = Cursor::new(Vec::new());
    let options = WriteOptions { compression: None };

    let mut writer = FileWriter::try_new(buffer, schema, None, options)
        .map_err(|e| PyValueError::new_err(format!("FileWriter init failed: {e}")))?;

    writer
        .write(&chunk, None)
        .map_err(|e| PyValueError::new_err(format!("IPC write failed: {e}")))?;

    writer
        .finish()
        .map_err(|e| PyValueError::new_err(format!("IPC finish failed: {e}")))?;

    let buffer = writer.into_inner();
    let bytes = buffer.into_inner();

    let pybytes = Python::with_gil(|py| PyBytes::new(py, &bytes).into());

    Ok(pybytes)
}

#[pymodule]
fn gvcfparser(_py: Python<'_>, m: Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(export_arrow_ipc, &m)?)?;
    Ok(())
}
