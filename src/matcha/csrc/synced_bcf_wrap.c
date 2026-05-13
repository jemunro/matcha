/* Thin C wrappers for synced_bcf_reader macros (pure struct-field access,
   no external symbol references so no -lhts needed for this object). */
#include "htslib/synced_bcf_reader.h"

int        srs_has_line   (bcf_srs_t *sr, int i) { return bcf_sr_has_line(sr, i); }
bcf1_t    *srs_get_line   (bcf_srs_t *sr, int i) { return bcf_sr_get_line(sr, i); }
bcf_hdr_t *srs_get_header (bcf_srs_t *sr, int i) { return bcf_sr_get_header(sr, i); }
int        srs_nreaders   (bcf_srs_t *sr)         { return sr->nreaders; }
int        srs_errnum     (bcf_srs_t *sr)         { return (int)sr->errnum; }
