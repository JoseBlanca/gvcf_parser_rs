use crate::region_splitter::{GVcfRecord, GVcfRecordIterator};
use crate::vcf_iterator::{VcfRecord, VcfRecordIterator, VcfResult};
use pyo3::prelude::*;
use pyo3::types::PyModule;
use pyo3::types::PyType;
use rust_htslib::tpool::ThreadPool;

#[pyclass(name = "VcfRecord")]
#[derive(Debug, Clone)]
pub struct PyVcfRecord {
    #[pyo3(get)]
    pub chrom: String,
    #[pyo3(get)]
    pub pos: u32,
    #[pyo3(get)]
    pub alleles: Vec<String>,
    #[pyo3(get)]
    pub qual: f32,
    #[pyo3(get)]
    pub genotypes: Vec<i32>,
}

impl From<VcfRecord> for PyVcfRecord {
    fn from(rec: VcfRecord) -> Self {
        PyVcfRecord {
            chrom: rec.chrom,
            pos: rec.pos,
            alleles: rec.alleles,
            qual: rec.qual,
            genotypes: rec.genotypes,
        }
    }
}

#[pyclass(name = "VcfRecordIterator", unsendable)]
pub struct PyVcfRecordIterator {
    inner: Box<dyn Iterator<Item = VcfResult<VcfRecord>>>,
    _pool: Option<ThreadPool>,
}

#[pymethods]
impl PyVcfRecordIterator {
    #[new]
    fn new(path: String, n_threads: u32) -> PyResult<Self> {
        let (parser, pool) = VcfRecordIterator::from_gzipped_vcf_path(path, n_threads)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{:?}", e)))?;
        Ok(Self {
            inner: Box::new(parser),
            _pool: pool,
        })
    }

    fn __iter__(slf: PyRefMut<'_, Self>) -> Py<PyVcfRecordIterator> {
        slf.into()
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<PyResult<PyVcfRecord>> {
        slf.inner.next().map(|result| match result {
            Ok(rec) => Ok(PyVcfRecord::from(rec)),
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                "{:?}",
                e
            ))),
        })
    }
}

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
fn vcfparser(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyGVcfRecord>()?;
    m.add_class::<GVcfGzipIterator>()?;
    //m.add_class::<PyVcfRecordIterator>()?;
    //m.add_class::<PyVcfRecord>()?;
    Ok(())
}
