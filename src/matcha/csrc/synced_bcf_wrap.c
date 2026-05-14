/* Thin C wrappers for synced_bcf_reader macros (pure struct-field access,
   no external symbol references so no -lhts needed for this object). */
#include "htslib/synced_bcf_reader.h"
#include "htslib/hts.h"

int        srs_has_line   (bcf_srs_t *sr, int i) { return bcf_sr_has_line(sr, i); }
bcf1_t    *srs_get_line   (bcf_srs_t *sr, int i) { return bcf_sr_get_line(sr, i); }
bcf_hdr_t *srs_get_header (bcf_srs_t *sr, int i) { return bcf_sr_get_header(sr, i); }
int        srs_nreaders   (bcf_srs_t *sr)         { return sr->nreaders; }
int        srs_errnum     (bcf_srs_t *sr)         { return (int)sr->errnum; }

/* For htsThreadPool wiring: access each reader's underlying htsFile and surface
   the HTS_OPT_THREAD_POOL enum value (avoids hardcoding enum ordering). */
htsFile *srs_get_file              (bcf_srs_t *sr, int i) { return sr->readers[i].file; }
int      srs_hts_opt_thread_pool   (void)                 { return HTS_OPT_THREAD_POOL; }
