use crate::gvcf_parser::{GVcfRecord, GVcfRecordIterator};
use pyo3::prelude::*;
use pyo3::types::PyModule;
use pyo3::types::PyType;

#[pyclass]
pub struct PyGVcfRecord {
    inner: GVcfRecord,
}

#[pymethods]
impl PyGVcfRecord {
    #[getter]
    fn chrom(&self) -> &str {
        &self.inner.chrom
    }

    #[getter]
    fn pos(&self) -> u32 {
        self.inner.pos
    }

    #[getter]
    fn alleles(&self) -> Vec<String> {
        self.inner.alleles.clone()
    }

    pub fn get_span(&self) -> PyResult<(u32, u32)> {
        self.inner
            .get_span()
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(format!("get_span error: {}", e)))
    }

    pub fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}

#[pyclass]
pub struct GVcfGzipIterator {
    inner: GVcfRecordIterator<std::io::BufReader<flate2::read::MultiGzDecoder<std::fs::File>>>,
}

#[pymethods]
impl GVcfGzipIterator {
    #[classmethod]
    pub fn from_path(_cls: &Bound<'_, PyType>, path: &str) -> PyResult<Self> {
        GVcfRecordIterator::from_gzip_path(path)
            .map(|inner| Self { inner })
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("{}", e)))
    }

    fn __iter__(slf: PyRefMut<'_, Self>) -> PyRefMut<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<PyGVcfRecord> {
        while let Some(item) = slf.inner.next() {
            match item {
                Ok(record) => return Some(PyGVcfRecord { inner: record }),
                Err(_) => continue,
            }
        }
        None
    }
}

#[pymodule]
#[pyo3(name = "vcfparser")]
fn vcfparser(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<GVcfGzipIterator>()?;
    //m.add_class::<PyVcfRecordIterator>()?;
    //m.add_class::<PyVcfRecord>()?;
    Ok(())
}
