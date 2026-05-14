## FFI bindings for htslib synced_bcf_reader.
## synced_bcf_reader is not exposed by hts-nim; we bind it here via a thin
## C wrapper (synced_bcf_wrap.c) that converts macros to regular functions.

{.compile: "csrc/synced_bcf_wrap.c".}

import hts/private/hts_concat
from hts/private/hts_concat import libname

# Opaque handle — all interaction is via API functions / thin wrappers.
type bcf_srs_t* {.bycopy.} = object

# ---- Core API (real exported C functions) -----------------------------------

proc bcf_sr_init*(): ptr bcf_srs_t
    {.cdecl, importc: "bcf_sr_init", dynlib: libname.}
proc bcf_sr_destroy*(sr: ptr bcf_srs_t)
    {.cdecl, importc: "bcf_sr_destroy", dynlib: libname.}
proc bcf_sr_strerror*(errnum: cint): cstring
    {.cdecl, importc: "bcf_sr_strerror", dynlib: libname.}
proc bcf_sr_add_reader*(sr: ptr bcf_srs_t; fname: cstring): cint
    {.cdecl, importc: "bcf_sr_add_reader", dynlib: libname.}
proc bcf_sr_next_line*(sr: ptr bcf_srs_t): cint
    {.cdecl, importc: "bcf_sr_next_line", dynlib: libname.}
proc bcf_sr_seek*(sr: ptr bcf_srs_t; sq: cstring; pos: int64): cint
    {.cdecl, importc: "bcf_sr_seek", dynlib: libname.}
proc bcf_sr_set_threads*(sr: ptr bcf_srs_t; n: cint): cint
    {.cdecl, importc: "bcf_sr_set_threads", dynlib: libname.}

# ---- bcf_sr_set_opt with variadic support -----------------------------------

const
  BCF_SR_REQUIRE_IDX*     = 0.cint
  BCF_SR_PAIR_LOGIC*      = 1.cint
  BCF_SR_ALLOW_NO_IDX*    = 2.cint
  BCF_SR_REGIONS_OVERLAP* = 3.cint
  BCF_SR_TARGETS_OVERLAP* = 4.cint

proc bcf_sr_set_opt*(sr: ptr bcf_srs_t; opt: cint): cint
    {.varargs, cdecl, importc: "bcf_sr_set_opt", dynlib: libname.}

# ---- Thin wrappers from synced_bcf_wrap.c (pure struct-field access) --------

proc srs_has_line*(sr: ptr bcf_srs_t; i: cint): cint
    {.cdecl, importc: "srs_has_line".}
proc srs_get_line*(sr: ptr bcf_srs_t; i: cint): ptr bcf1_t
    {.cdecl, importc: "srs_get_line".}
proc srs_get_header*(sr: ptr bcf_srs_t; i: cint): ptr bcf_hdr_t
    {.cdecl, importc: "srs_get_header".}
proc srs_nreaders*(sr: ptr bcf_srs_t): cint
    {.cdecl, importc: "srs_nreaders".}
proc srs_errnum*(sr: ptr bcf_srs_t): cint
    {.cdecl, importc: "srs_errnum".}
proc srs_get_file*(sr: ptr bcf_srs_t; i: cint): ptr htsFile
    {.cdecl, importc: "srs_get_file".}
proc srs_hts_opt_thread_pool*(): cint
    {.cdecl, importc: "srs_hts_opt_thread_pool".}

# ---- bcf_translate (not in hts-nim) -----------------------------------------

proc bcf_translate*(dst_hdr, src_hdr: ptr bcf_hdr_t; line: ptr bcf1_t): cint
    {.cdecl, importc: "bcf_translate", dynlib: libname.}

# ---- htsThreadPool (shared BGZF decompression across readers/writers) -------

type htsThreadPool* {.bycopy.} = object
  pool*:  pointer
  qsize*: cint

proc hts_tpool_init*(n: cint): pointer
    {.cdecl, importc: "hts_tpool_init", dynlib: libname.}
proc hts_tpool_destroy*(p: pointer)
    {.cdecl, importc: "hts_tpool_destroy", dynlib: libname.}
proc hts_set_opt*(fp: ptr htsFile; opt: cint): cint
    {.varargs, cdecl, importc: "hts_set_opt", dynlib: libname.}

# ---- C stdlib free for bcf_get_info_values buffers --------------------------

proc c_free*(p: pointer) {.importc: "free", header: "<stdlib.h>".}

# ---- Non-owning hts-nim Variant view ---------------------------------------
#
# Construct a hts-nim Variant over a (hdr, rec) pair owned elsewhere (e.g. by
# bcf_sr_next_line). The constructor leaves Variant.p (private INFO scratch)
# at nil and Variant.own at false; no finalizer is registered, so the GC will
# not call bcf_destroy on the record or bcf_hdr_destroy on the header.
#
# Variant.p is allocated lazily by hts-nim's bcf_get_*_values; it is freed
# only by hts-nim's destroy_variant, which we cannot register from outside
# the module. Allocate ONE view per reader and re-point it at each record via
# setRecView — htslib reallocates p in place, so the buffer is reused. A
# single ~1 KB buffer leaks at process exit, which is acceptable.

import hts/vcf
export vcf.Variant

proc newVariantView*(): Variant =
  ## Allocate a reusable, non-owning Variant view. Target it at a record via
  ## setRecView. Caller is responsible for ensuring the header and record
  ## pointers outlive any reads from the view.
  Variant(vcf: VCF(header: Header()))

proc setRecView*(v: Variant; hdr: ptr bcf_hdr_t; rec: ptr bcf1_t) {.inline.} =
  ## Re-target an existing view at a new (header, record) pair.
  v.vcf.header.hdr = hdr
  v.c = rec
