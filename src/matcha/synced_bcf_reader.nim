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

# ---- bcf_translate (not in hts-nim) -----------------------------------------

proc bcf_translate*(dst_hdr, src_hdr: ptr bcf_hdr_t; line: ptr bcf1_t): cint
    {.cdecl, importc: "bcf_translate", dynlib: libname.}

# ---- C stdlib free for bcf_get_info_values buffers --------------------------

proc c_free*(p: pointer) {.importc: "free", header: "<stdlib.h>".}
