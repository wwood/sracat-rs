//! Safe-ish wrapper around the C shim (`shim.c`) over the ncbi-vdb cursor API.

use std::ffi::{c_char, c_int, CString};
use std::ptr;
use std::slice;

use anyhow::{bail, Result};

/// Opaque handle owned by the C shim.
#[repr(C)]
struct SracatRun {
    _private: [u8; 0],
}

extern "C" {
    fn sracat_open(
        path: *const c_char,
        with_quality: c_int,
        allow_aligned: c_int,
        out: *mut *mut SracatRun,
        errbuf: *mut c_char,
        errlen: usize,
    ) -> c_int;
    fn sracat_first_row(run: *const SracatRun) -> i64;
    fn sracat_row_count(run: *const SracatRun) -> u64;
    fn sracat_is_aligned(run: *const SracatRun) -> c_int;
    #[allow(clippy::too_many_arguments)]
    fn sracat_read_spot(
        run: *const SracatRun,
        row: i64,
        bases: *mut *const c_char,
        nbases: *mut u32,
        quals: *mut *const u8,
        read_len: *mut *const u32,
        read_type: *mut *const u8,
        nreads: *mut u32,
        errbuf: *mut c_char,
        errlen: usize,
    ) -> c_int;
    fn sracat_close(run: *mut SracatRun);
}

/// An opened SRA run. Closes the underlying cursor/table/manager on drop.
pub struct Run {
    ptr: *mut SracatRun,
}

/// Borrowed cell data for one spot. Valid only until the next `read_spot` call.
pub struct Spot<'a> {
    pub bases: &'a [u8],
    /// Per-base phred scores (same length as `bases`); `None` unless opened with
    /// quality.
    pub quals: Option<&'a [u8]>,
    pub read_len: &'a [u32],
    pub read_type: &'a [u8],
}

fn errstr(buf: &[u8]) -> String {
    let end = buf.iter().position(|&b| b == 0).unwrap_or(buf.len());
    String::from_utf8_lossy(&buf[..end]).into_owned()
}

impl Run {
    pub fn open(path: &str, with_quality: bool, allow_aligned: bool) -> Result<Run> {
        let c = CString::new(path)?;
        let mut out: *mut SracatRun = ptr::null_mut();
        let mut err = [0u8; 512];
        let rc = unsafe {
            sracat_open(
                c.as_ptr(),
                c_int::from(with_quality),
                c_int::from(allow_aligned),
                &mut out,
                err.as_mut_ptr() as *mut c_char,
                err.len(),
            )
        };
        if rc != 0 {
            bail!("{path}: {}", errstr(&err));
        }
        Ok(Run { ptr: out })
    }

    pub fn first_row(&self) -> i64 {
        unsafe { sracat_first_row(self.ptr) }
    }

    pub fn row_count(&self) -> u64 {
        unsafe { sracat_row_count(self.ptr) }
    }

    /// Whether the run is aligned (cSRA): `READ` is reconstructed from the
    /// alignment table via random access, which does not parallelise.
    pub fn is_aligned(&self) -> bool {
        unsafe { sracat_is_aligned(self.ptr) != 0 }
    }

    pub fn read_spot(&self, row: i64) -> Result<Spot<'_>> {
        let mut bases: *const c_char = ptr::null();
        let mut nbases: u32 = 0;
        let mut quals: *const u8 = ptr::null();
        let mut read_len: *const u32 = ptr::null();
        let mut read_type: *const u8 = ptr::null();
        let mut nreads: u32 = 0;
        let mut err = [0u8; 512];
        let rc = unsafe {
            sracat_read_spot(
                self.ptr,
                row,
                &mut bases,
                &mut nbases,
                &mut quals,
                &mut read_len,
                &mut read_type,
                &mut nreads,
                err.as_mut_ptr() as *mut c_char,
                err.len(),
            )
        };
        if rc != 0 {
            bail!("row {row}: {}", errstr(&err));
        }
        // SAFETY: on success the shim guarantees the pointers are non-null and
        // valid for the given lengths until the next cursor access.
        let spot = unsafe {
            Spot {
                bases: slice::from_raw_parts(bases as *const u8, nbases as usize),
                quals: if quals.is_null() {
                    None
                } else {
                    Some(slice::from_raw_parts(quals, nbases as usize))
                },
                read_len: slice::from_raw_parts(read_len, nreads as usize),
                read_type: slice::from_raw_parts(read_type, nreads as usize),
            }
        };
        Ok(spot)
    }
}

impl Drop for Run {
    fn drop(&mut self) {
        unsafe { sracat_close(self.ptr) }
    }
}
