#include "shim.h"

#include <vdb/manager.h>
#include <vdb/database.h>
#include <vdb/table.h>
#include <vdb/cursor.h>
#include <klib/namelist.h>
#include <klib/rc.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct SracatRun {
    const VDBManager *mgr;
    const VDatabase *db; /* NULL for flat (legacy) table objects */
    const VTable *tbl;
    const VCursor *curs;
    uint32_t read_idx;
    uint32_t read_len_idx;
    uint32_t read_type_idx;
    uint32_t qual_idx; /* valid only if has_qual */
    int has_qual;
    int64_t first;
    uint64_t count;
};

static void seterr(char *errbuf, size_t errlen, const char *msg) {
    if (errbuf != NULL && errlen > 0) {
        snprintf(errbuf, errlen, "%s", msg);
    }
}

int sracat_open(const char *path, int with_quality, int allow_aligned,
                SracatRun **out, char *errbuf, size_t errlen) {
    *out = NULL;
    SracatRun *r = calloc(1, sizeof(*r));
    if (r == NULL) {
        seterr(errbuf, errlen, "out of memory");
        return 1;
    }
    r->has_qual = with_quality;

    if (VDBManagerMakeRead(&r->mgr, NULL) != 0) {
        seterr(errbuf, errlen, "VDBManagerMakeRead failed");
        goto fail;
    }

    /* Most SRA runs are databases; legacy runs are flat tables. */
    if (VDBManagerOpenDBRead(r->mgr, &r->db, NULL, "%s", path) == 0) {
        /* Refuse aligned runs unless explicitly allowed: a PRIMARY_ALIGNMENT
         * table means READ in the SEQUENCE table is reconstructed from
         * alignments, not stored. */
        KNamelist *tables = NULL;
        if (!allow_aligned && VDatabaseListTbl(r->db, &tables) == 0 && tables != NULL) {
            uint32_t n = 0;
            KNamelistCount(tables, &n);
            for (uint32_t i = 0; i < n; i++) {
                const char *nm = NULL;
                if (KNamelistGet(tables, i, &nm) == 0 && nm != NULL &&
                    strcmp(nm, "PRIMARY_ALIGNMENT") == 0) {
                    KNamelistRelease(tables);
                    seterr(errbuf, errlen,
                           "aligned run (PRIMARY_ALIGNMENT present): READ is "
                           "reconstructed from alignments, refusing");
                    goto fail;
                }
            }
            KNamelistRelease(tables);
        }
        if (VDatabaseOpenTableRead(r->db, &r->tbl, "SEQUENCE") != 0) {
            seterr(errbuf, errlen, "could not open SEQUENCE table");
            goto fail;
        }
    } else if (VDBManagerOpenTableRead(r->mgr, &r->tbl, NULL, "%s", path) != 0) {
        seterr(errbuf, errlen, "not a readable SRA database or table");
        goto fail;
    }

    if (VTableCreateCursorRead(r->tbl, &r->curs) != 0) {
        seterr(errbuf, errlen, "VTableCreateCursorRead failed");
        goto fail;
    }
    if (VCursorAddColumn(r->curs, &r->read_idx, "(INSDC:dna:text)READ") != 0) {
        seterr(errbuf, errlen, "no READ column (aligned or unsupported run)");
        goto fail;
    }
    if (VCursorAddColumn(r->curs, &r->read_len_idx, "(INSDC:coord:len)READ_LEN") != 0) {
        seterr(errbuf, errlen, "no READ_LEN column");
        goto fail;
    }
    if (VCursorAddColumn(r->curs, &r->read_type_idx, "(INSDC:SRA:read_type)READ_TYPE") != 0) {
        seterr(errbuf, errlen, "no READ_TYPE column");
        goto fail;
    }
    if (r->has_qual &&
        VCursorAddColumn(r->curs, &r->qual_idx, "(INSDC:quality:phred)QUALITY") != 0) {
        seterr(errbuf, errlen, "no QUALITY column");
        goto fail;
    }
    if (VCursorOpen(r->curs) != 0) {
        seterr(errbuf, errlen, "VCursorOpen failed");
        goto fail;
    }
    if (VCursorIdRange(r->curs, r->read_idx, &r->first, &r->count) != 0) {
        seterr(errbuf, errlen, "VCursorIdRange failed");
        goto fail;
    }

    *out = r;
    return 0;

fail:
    sracat_close(r);
    return 1;
}

int64_t sracat_first_row(const SracatRun *run) { return run->first; }
uint64_t sracat_row_count(const SracatRun *run) { return run->count; }

static int cell(const VCursor *curs, int64_t row, uint32_t idx,
                uint32_t want_bits, const void **base, uint32_t *len,
                char *errbuf, size_t errlen, const char *what) {
    uint32_t elem_bits = 0, boff = 0, row_len = 0;
    const void *b = NULL;
    if (VCursorCellDataDirect(curs, row, idx, &elem_bits, &b, &boff, &row_len) != 0) {
        seterr(errbuf, errlen, what);
        return 1;
    }
    if (elem_bits != want_bits) {
        seterr(errbuf, errlen, "unexpected column element size");
        return 1;
    }
    *base = b;
    *len = row_len;
    return 0;
}

int sracat_read_spot(const SracatRun *run, int64_t row,
                     const char **bases, uint32_t *nbases,
                     const uint8_t **quals,
                     const uint32_t **read_len,
                     const uint8_t **read_type, uint32_t *nreads,
                     char *errbuf, size_t errlen) {
    const void *b = NULL;
    uint32_t n = 0;

    if (cell(run->curs, row, run->read_idx, 8, &b, &n, errbuf, errlen,
             "reading READ failed"))
        return 1;
    *bases = (const char *)b;
    *nbases = n;

    *quals = NULL;
    if (run->has_qual) {
        const void *q = NULL;
        uint32_t nq = 0;
        if (cell(run->curs, row, run->qual_idx, 8, &q, &nq, errbuf, errlen,
                 "reading QUALITY failed"))
            return 1;
        if (nq != n) {
            seterr(errbuf, errlen, "READ / QUALITY length mismatch");
            return 1;
        }
        *quals = (const uint8_t *)q;
    }

    if (cell(run->curs, row, run->read_len_idx, 32, &b, &n, errbuf, errlen,
             "reading READ_LEN failed"))
        return 1;
    *read_len = (const uint32_t *)b;
    *nreads = n;

    uint32_t nt = 0;
    if (cell(run->curs, row, run->read_type_idx, 8, &b, &nt, errbuf, errlen,
             "reading READ_TYPE failed"))
        return 1;
    *read_type = (const uint8_t *)b;
    if (nt != n) {
        seterr(errbuf, errlen, "READ_LEN / READ_TYPE length mismatch");
        return 1;
    }
    return 0;
}

void sracat_close(SracatRun *run) {
    if (run == NULL)
        return;
    if (run->curs != NULL)
        VCursorRelease(run->curs);
    if (run->tbl != NULL)
        VTableRelease(run->tbl);
    if (run->db != NULL)
        VDatabaseRelease(run->db);
    if (run->mgr != NULL)
        VDBManagerRelease(run->mgr);
    free(run);
}
