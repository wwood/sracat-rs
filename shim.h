/* Minimal C shim over the ncbi-vdb VDB cursor API.
 *
 * Exposes just enough to iterate the SEQUENCE table of an SRA run in storage
 * (row) order and hand the raw READ / READ_LEN / READ_TYPE cell data back to
 * Rust. All rc_t / variadic handling lives here so the Rust side is plain FFI.
 */
#ifndef SRACAT_SHIM_H
#define SRACAT_SHIM_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SracatRun SracatRun;

/* Open an SRA run for reading.
 * If with_quality != 0, the QUALITY column is also opened so sracat_read_spot
 * returns per-base phred scores.
 * If allow_aligned == 0, aligned runs (a PRIMARY_ALIGNMENT table is present) are
 * refused, since READ is reconstructed from alignments rather than stored. If
 * allow_aligned != 0, such runs are opened anyway and the computed READ column
 * reconstructs each spot (correct and in spot order, but with random access into
 * the alignment table).
 * Returns 0 on success (and sets *out), nonzero on failure (message in errbuf). */
int sracat_open(const char *path, int with_quality, int allow_aligned,
                SracatRun **out, char *errbuf, size_t errlen);

/* First row id and number of rows (spots) in the SEQUENCE table. */
int64_t sracat_first_row(const SracatRun *run);
uint64_t sracat_row_count(const SracatRun *run);

/* Read one spot at absolute row id.
 * On success (returns 0):
 *   *bases     -> ASCII read bases for the whole spot, length *nbases
 *   *quals     -> per-base phred scores (length *nbases), or NULL if the run was
 *                 opened without quality
 *   *read_len  -> uint32 array of per-read lengths,  count *nreads
 *   *read_type -> uint8  array of per-read types,     count *nreads
 * The returned pointers are owned by the cursor and remain valid only until the
 * next sracat_read_spot / sracat_close call. */
int sracat_read_spot(const SracatRun *run, int64_t row,
                     const char **bases, uint32_t *nbases,
                     const uint8_t **quals,
                     const uint32_t **read_len,
                     const uint8_t **read_type, uint32_t *nreads,
                     char *errbuf, size_t errlen);

void sracat_close(SracatRun *run);

#ifdef __cplusplus
}
#endif

#endif /* SRACAT_SHIM_H */
