## mergecore.nim — shared clustering machinery for matcha collapse (and future merge).
##
## Exports:
##   Types: LinkageMethod, PriorityCriterion, CallerInput, FieldConflict,
##          FieldResolution, MergedHeader
##
##   resolveHeaders        — analyse N input VCF headers, produce MergedHeader
##   mergeSortSlimBcfs     — k-way merge + sort of per-caller slim BCFs
##   buildSimilarityMap    — MatchPair seq → (canonicalPair → similarity) table
##   buildComponents       — union-find over pairs → offset→componentId table
##   agglomerateComponent  — Lance-Williams agglomerative clustering
##   clusterAll            — full pipeline: components → agglomerate each
##   selectRepresentative  — priority-cascade representative selection

import std/[sequtils, sets, strutils, tables]
import hts
import hts/private/hts_concat
import preproc, matchcore, synced_bcf_reader

# ---------------------------------------------------------------------------
# Types
# ---------------------------------------------------------------------------

type
  LinkageMethod* = enum
    lmAverage  = "average"
    lmSingle   = "single"
    lmComplete = "complete"

  PriorityCriterion* = enum
    pcPass   = "PASS"
    pcQual   = "QUAL"
    pcCentre = "CENTRE"
    pcOrder  = "ORDER"

  CallerInput* = object
    name*: string   ## user-supplied prefix or basename from positional arg
    path*: string

  FieldConflict* = enum
    fcCompatible, fcIncompatibleInfo, fcIncompatibleFmt

  FieldResolution* = object
    ## How a field appearing in ≥2 callers was resolved.
    ## For fcCompatible: a single merged definition is used; no renames needed.
    ## For fcIncompatibleInfo: each caller that defines the field gets a renamed
    ##   version; renames[callerIdx] = renamed field name.
    ## For fcIncompatibleFmt: Number=. is emitted; no per-caller renaming.
    kind*:    FieldConflict
    renames*: Table[int, string]   ## callerIdx → new field name (Info conflicts only)

  MergedHeader* = object
    ## Result of resolving N input headers. Drives the merge pass and output.
    infoRes*: Table[string, FieldResolution]  ## orig INFO name → resolution
    fmtRes*:  Table[string, FieldResolution]  ## orig FORMAT name → resolution
    ## Synthesized header lines to add to the merged BCF header (in insertion order).
    ## Each entry is a complete ##INFO=<...> or ##FORMAT=<...> line.
    headerLines*: seq[string]
    warnings*: seq[string]

# Canonical pair key: always (min, max) so (a,b) and (b,a) collide.
proc pairKey*(a, b: int64): (int64, int64) {.inline.} =
  if a <= b: (a, b) else: (b, a)

# ---------------------------------------------------------------------------
# Header access helpers (low-level bcf_hrec_t traversal)
# ---------------------------------------------------------------------------

type HrecInfo = object
  id: string
  number: string
  typ: string
  desc: string
  line: string   # verbatim ##INFO=<...> or ##FORMAT=<...> line

proc parseHrec(hrec: ptr bcf_hrec_t): HrecInfo =
  let keys = cast[ptr UncheckedArray[cstring]](hrec.keys)
  let vals = cast[ptr UncheckedArray[cstring]](hrec.vals)
  result.line = if hrec.`type` == BCF_HEADER_TYPE.BCF_HL_INFO.cint:
                  "##INFO=<" else: "##FORMAT=<"
  for i in 0 ..< hrec.nkeys.int:
    let k = $keys[i]
    let v = $vals[i]
    if i > 0: result.line &= ","
    result.line &= k & "=" & v
    case k
    of "ID":     result.id   = v
    of "Number": result.number = v
    of "Type":   result.typ  = v
    of "Description":
      # htslib keeps surrounding quotes in vals for quoted fields; strip them
      # so makeInfoLine/makeFmtLine can re-add the quotes exactly once.
      result.desc = if v.len >= 2 and v[0] == '"' and v[^1] == '"':
                      v[1 ..< v.len - 1]
                    else: v
  result.line &= ">"

proc collectHrecs(h: ptr bcf_hdr_t;
                  hlType: cint): seq[HrecInfo] =
  let hrecs = cast[ptr UncheckedArray[ptr bcf_hrec_t]](h.hrec)
  for i in 0 ..< h.nhrec.int:
    let hr = hrecs[i]
    if hr.`type` != hlType: continue
    result.add(parseHrec(hr))

proc contigs(h: ptr bcf_hdr_t): seq[string] =
  var n: cint = 0
  let names = bcf_hdr_seqnames(h, n.addr)
  if names == nil: return
  for i in 0 ..< n.int:
    result.add($names[i])
  free(names)

# ---------------------------------------------------------------------------
# resolveHeaders
# ---------------------------------------------------------------------------

proc mergeDescriptions(descs: seq[string]): string =
  ## Concatenate unique descriptions with "; ".
  var seen: seq[string]
  for d in descs:
    if d notin seen: seen.add(d)
  seen.join("; ")

proc makeInfoLine(id, number, typ, desc: string): string =
  "##INFO=<ID=" & id & ",Number=" & number & ",Type=" & typ &
  ",Description=\"" & desc & "\">"

proc makeFmtLine(id, number, typ, desc: string): string =
  "##FORMAT=<ID=" & id & ",Number=" & number & ",Type=" & typ &
  ",Description=\"" & desc & "\">"

proc resolveHeaders*(callers: seq[CallerInput]): MergedHeader =
  ## Open each caller's header, collect INFO and FORMAT definitions, and
  ## produce a MergedHeader that drives the merge pass and output header.
  ##
  ## Compatible INFO/FORMAT (same Number+Type across all callers that define it):
  ##   → single merged def in output; no renames needed.
  ##
  ## Incompatible INFO (different Number or Type):
  ##   → rename every instance to FIELD_CALLERNAME (or FIELD_N if unnamed).
  ##     Warn once per conflict. The original field name is dropped.
  ##
  ## Incompatible FORMAT:
  ##   → emit Number=. with merged Type (best-effort) and warn.
  ##
  ## Fields that appear in only one caller are kept unchanged.

  # Collect INFO and FORMAT definitions per field, per caller.
  type FieldDef = object
    callerIdx: int
    number, typ, desc: string
    line: string

  var infoByField:  Table[string, seq[FieldDef]]
  var fmtByField:   Table[string, seq[FieldDef]]
  # Collect all contig lines (for the merged header; we use the first caller's order).
  var firstContigLines: seq[string]

  for ci, caller in callers.pairs:
    var vcf: VCF
    if not open(vcf, caller.path):
      raise newException(IOError, "cannot open: " & caller.path)
    let h = vcf.header.hdr

    # INFO
    for hi in collectHrecs(h, BCF_HEADER_TYPE.BCF_HL_INFO.cint):
      infoByField.mgetOrPut(hi.id, @[]).add(
        FieldDef(callerIdx: ci, number: hi.number, typ: hi.typ, desc: hi.desc, line: hi.line))

    # FORMAT
    for hi in collectHrecs(h, BCF_HEADER_TYPE.BCF_HL_FMT.cint):
      fmtByField.mgetOrPut(hi.id, @[]).add(
        FieldDef(callerIdx: ci, number: hi.number, typ: hi.typ, desc: hi.desc, line: hi.line))

    if ci == 0:
      for ctg in contigs(h):
        firstContigLines.add("##contig=<ID=" & ctg & ">")

    vcf.close()

  # Resolve INFO fields.
  for fld, defs in infoByField.pairs:
    let byCallers = defs  # one entry per caller that defines this field
    # Check compatibility: same Number+Type across all callers that define it.
    let refNum = byCallers[0].number
    let refTyp = byCallers[0].typ
    var compatible = true
    for d in byCallers:
      if d.number != refNum or d.typ != refTyp:
        compatible = false
        break
    if compatible or byCallers.len == 1:
      # Single definition (possibly from one caller only) or compatible → keep as-is.
      let desc = mergeDescriptions(byCallers.mapIt(it.desc))
      result.headerLines.add(makeInfoLine(fld, refNum, refTyp, desc))
      # No rename entry needed: resolution defaults to compatible.
    else:
      # Incompatible → rename each caller's instance.
      var res = FieldResolution(kind: fcIncompatibleInfo)
      var warnParts: seq[string]
      for d in byCallers:
        let suffix = callers[d.callerIdx].name
        let newName = fld & "_" & suffix
        res.renames[d.callerIdx] = newName
        result.headerLines.add(makeInfoLine(newName, d.number, d.typ, d.desc))
        warnParts.add(callers[d.callerIdx].name & ":" & d.number & "/" & d.typ)
      result.infoRes[fld] = res
      result.warnings.add("INFO conflict for " & fld & " (" & warnParts.join(", ") &
                          ") — renamed to " & res.renames.values.toSeq.join(", "))

  # Resolve FORMAT fields.
  for fld, defs in fmtByField.pairs:
    let refNum = defs[0].number
    let refTyp = defs[0].typ
    var compatible = true
    for d in defs:
      if d.number != refNum or d.typ != refTyp:
        compatible = false; break
    if compatible or defs.len == 1:
      let desc = mergeDescriptions(defs.mapIt(it.desc))
      result.headerLines.add(makeFmtLine(fld, refNum, refTyp, desc))
    else:
      # Incompatible FORMAT → Number=. with warning.
      var res = FieldResolution(kind: fcIncompatibleFmt)
      let desc = mergeDescriptions(defs.mapIt(it.desc))
      result.headerLines.add(makeFmtLine(fld, ".", refTyp, desc))
      result.fmtRes[fld] = res
      var warnParts: seq[string]
      for d in defs:
        warnParts.add(callers[d.callerIdx].name & ":" & d.number & "/" & d.typ)
      result.warnings.add("FORMAT conflict for " & fld & " (" & warnParts.join(", ") &
                          ") — emitting Number=.")

# ---------------------------------------------------------------------------
# Similarity map (from MatchPair seq)
# ---------------------------------------------------------------------------

proc buildSimilarityMap*(pairs: seq[MatchPair]): Table[(int64, int64), float64] =
  ## Build canonical-pair → similarity table from match pairs.
  ## Canonical key: (min(aOff,bOff), max(aOff,bOff)).
  for p in pairs:
    let key = pairKey(p.aOff, p.bOff)
    # Keep highest similarity if duplicate pairs arrive (shouldn't happen, but be safe).
    let cur = result.getOrDefault(key, 0.0)
    if p.sim > cur:
      result[key] = p.sim

# ---------------------------------------------------------------------------
# Union-find → connected components
# ---------------------------------------------------------------------------

type UnionFind = object
  parent: seq[int]
  rank:   seq[int]

proc initUnionFind(n: int): UnionFind =
  result.parent = newSeq[int](n)
  result.rank   = newSeq[int](n)
  for i in 0 ..< n:
    result.parent[i] = i

proc find(uf: var UnionFind; x: int): int =
  if uf.parent[x] != x:
    uf.parent[x] = find(uf, uf.parent[x])  # path compression
  uf.parent[x]

proc union(uf: var UnionFind; x, y: int) =
  let rx = find(uf, x)
  let ry = find(uf, y)
  if rx == ry: return
  if uf.rank[rx] < uf.rank[ry]:
    uf.parent[rx] = ry
  elif uf.rank[rx] > uf.rank[ry]:
    uf.parent[ry] = rx
  else:
    uf.parent[ry] = rx
    inc uf.rank[rx]

proc buildComponents*(simMap: Table[(int64, int64), float64];
                      allOffsets: seq[int64]): Table[int64, int] =
  ## Union-find over pairs → offset→componentId table.
  ## Singletons get their own component. ComponentIds are not contiguous.
  let n = allOffsets.len
  if n == 0: return
  var offIdx: Table[int64, int]
  for i, off in allOffsets.pairs:
    offIdx[off] = i
  var uf = initUnionFind(n)
  for (a, b) in simMap.keys:
    let ia = offIdx.getOrDefault(a, -1)
    let ib = offIdx.getOrDefault(b, -1)
    if ia >= 0 and ib >= 0:
      union(uf, ia, ib)
  for i, off in allOffsets.pairs:
    result[off] = find(uf, i)

# ---------------------------------------------------------------------------
# k-way merge-sort of per-caller slim BCFs
# ---------------------------------------------------------------------------

proc mergeSortSlimBcfs*(inputs: seq[string]; callerIdxs: seq[int];
                        outPath: string; chromOrder: seq[string]) =
  ## Merge N slim BCFs for the same (svtype, bin) into a coordinate-sorted
  ## output BCF with a CSI index using synced_bcf_reader for streaming.
  ## Rewrites MATCHA_BOFF as a composite value:
  ##   compositeOff = (callerIdx.int64 shl 48) or origOff
  ## Works for k=1 (degenerate single-input pass that still rewrites BOFF + indexes).
  ## Chromosome order follows the natural index order of the first reader.
  if inputs.len == 0: return

  let sr = bcf_sr_init()
  if sr == nil:
    raise newException(IOError, "bcf_sr_init failed")
  discard bcf_sr_set_opt(sr, BCF_SR_ALLOW_NO_IDX)

  for path in inputs:
    if bcf_sr_add_reader(sr, path.cstring) == 0:
      raise newException(IOError, "cannot open slim BCF for merge: " & path &
                         " (" & $bcf_sr_strerror(srs_errnum(sr)) & ")")

  let nreaders = srs_nreaders(sr).int

  # Build output header: dup first reader's slim header, add missing contigs.
  let outHdr = bcf_hdr_dup(srs_get_header(sr, 0))
  var knownCtgs: HashSet[string]
  block:
    var n: cint = 0
    let names = bcf_hdr_seqnames(outHdr, n.addr)
    if names != nil:
      for i in 0 ..< n.int: knownCtgs.incl($names[i])
      free(names)
  for ri in 1 ..< nreaders:
    var n: cint = 0
    let names = bcf_hdr_seqnames(srs_get_header(sr, ri.cint), n.addr)
    if names != nil:
      for i in 0 ..< n.int:
        let c = $names[i]
        if c notin knownCtgs:
          discard bcf_hdr_append(outHdr, ("##contig=<ID=" & c & ">").cstring)
          knownCtgs.incl(c)
      free(names)
  discard bcf_hdr_sync(outHdr)

  # Open writer: open the first input just to get a hts-nim header wrapper,
  # then replace the underlying hdr pointer with our merged outHdr.
  var tmpVcf: VCF
  if not open(tmpVcf, inputs[0]):
    raise newException(IOError, "cannot open for header init: " & inputs[0])
  var outVcf: VCF
  if not open(outVcf, outPath, mode = "wb"):
    raise newException(IOError, "cannot create merged slim BCF: " & outPath)
  outVcf.copy_header(tmpVcf.header)
  tmpVcf.close()
  bcf_hdr_destroy(outVcf.header.hdr)
  outVcf.header.hdr = outHdr
  discard outVcf.write_header()

  let headerOff = uint64(bgzf_tell(bgzfHandle(outVcf)))
  let idx = hts_idx_init(0, HTS_FMT_CSI.cint, headerOff, 14, 5)
  if idx == nil:
    raise newException(IOError, "cannot create CSI index: " & outPath)

  # Reusable buffers for bcf_get_info_values (reallocated by htslib as needed,
  # freed with c_free after the loop).
  var boffDst: pointer = nil
  var boffN:   cint    = 0
  var endDst:  pointer = nil
  var endN:    cint    = 0

  while bcf_sr_next_line(sr) > 0:
    for i in 0 ..< nreaders:
      if srs_has_line(sr, i.cint) == 0: continue
      let srcHdr  = srs_get_header(sr, i.cint)
      let rec     = bcf_dup(srs_get_line(sr, i.cint))

      # Read MATCHA_BOFF and END before modifying the record.
      let nboff = bcf_get_info_values(srcHdr, rec, "MATCHA_BOFF",
                                      boffDst.addr, boffN.addr, BCF_HT_INT.cint)
      let origOff =
        if nboff >= 2:
          let a = cast[ptr UncheckedArray[int32]](boffDst)
          decodeBoff(@[a[0], a[1]])
        else: 0'i64

      let nend = bcf_get_info_values(srcHdr, rec, "END",
                                     endDst.addr, endN.addr, BCF_HT_INT.cint)
      let endPos =
        if nend >= 1: int64(cast[ptr UncheckedArray[int32]](endDst)[0])
        else: int64(rec.pos) + 1

      # Rewrite MATCHA_BOFF as composite (callerIdx | origOff).
      var newBoff = encodeBoff((callerIdxs[i].int64 shl 48) or origOff)
      discard bcf_update_info(srcHdr, rec, "MATCHA_BOFF",
                              newBoff[0].addr, 2.cint, BCF_HT_INT.cint)

      # Translate record RIDs and tag IDs from srcHdr space to outHdr space.
      discard bcf_translate(outHdr, srcHdr, rec)

      discard bcf_write(vcfHtsFile(outVcf), outHdr, rec)
      if rec.rid >= 0:
        let woff = uint64(bgzf_tell(bgzfHandle(outVcf)))
        discard hts_idx_push(idx, rec.rid, int64(rec.pos), endPos, woff, 1)
      bcf_destroy(rec)

  if boffDst != nil: c_free(boffDst)
  if endDst  != nil: c_free(endDst)
  bcf_sr_destroy(sr)

  let finalOff = uint64(bgzf_tell(bgzfHandle(outVcf)))
  hts_idx_finish(idx, finalOff)
  outVcf.close()
  hts_idx_save(idx, outPath.cstring, HTS_FMT_CSI.cint)
  hts_idx_destroy(idx)

# ---------------------------------------------------------------------------
# Agglomerative clustering (Lance-Williams on similarity)
# ---------------------------------------------------------------------------

type
  ClusterState = object
    ## During agglomeration: each cluster is a set of member offsets.
    members: seq[int64]   ## offsets in this cluster
    size:    int

proc agglomerateComponent*(offsets: seq[int64];
                            simMap: Table[(int64, int64), float64];
                            linkage: LinkageMethod;
                            threshold: float64): seq[seq[int64]] =
  ## Agglomerative clustering of one connected component.
  ## Returns a seq of clusters (each cluster = seq[int64] of member offsets).
  ## Singletons that never exceed threshold are returned as single-element clusters.
  if offsets.len == 0: return
  if offsets.len == 1:
    return @[@[offsets[0]]]

  let n = offsets.len
  # idx→offset, offset→idx
  var offIdx: Table[int64, int]
  for i, off in offsets.pairs:
    offIdx[off] = i

  # clusterSim[i][j] = current similarity between active clusters i and j (i < j).
  # Initialise from simMap (missing entries → 0.0).
  var clusterSim = newSeq[seq[float64]](n)
  for i in 0 ..< n:
    clusterSim[i] = newSeq[float64](n)
  for i in 0 ..< n:
    for j in i + 1 ..< n:
      let key = pairKey(offsets[i], offsets[j])
      clusterSim[i][j] = simMap.getOrDefault(key, 0.0)

  # Track active cluster membership. Start: each record is its own cluster.
  var members = newSeq[seq[int64]](n)
  for i in 0 ..< n:
    members[i] = @[offsets[i]]
  var sizes = newSeq[int](n)
  for i in 0 ..< n: sizes[i] = 1
  var active = newSeq[bool](n)
  for i in 0 ..< n: active[i] = true

  while true:
    # Find highest-similarity pair among active clusters.
    var bestSim = 0.0
    var bestI = -1
    var bestJ = -1
    for i in 0 ..< n:
      if not active[i]: continue
      for j in i + 1 ..< n:
        if not active[j]: continue
        if clusterSim[i][j] > bestSim:
          bestSim  = clusterSim[i][j]
          bestI    = i
          bestJ    = j
    if bestSim < threshold or bestI < 0:
      break  # no more merges above threshold

    # Merge cluster bestJ into bestI (keep bestI, deactivate bestJ).
    let sA = float64(sizes[bestI])
    let sB = float64(sizes[bestJ])

    for x in 0 ..< n:
      if not active[x] or x == bestI or x == bestJ: continue
      let dAX = if bestI < x: clusterSim[bestI][x] else: clusterSim[x][bestI]
      let dBX = if bestJ < x: clusterSim[bestJ][x] else: clusterSim[x][bestJ]
      let merged =
        case linkage
        of lmSingle:   max(dAX, dBX)
        of lmComplete: min(dAX, dBX)
        of lmAverage:  (sA * dAX + sB * dBX) / (sA + sB)
      if bestI < x: clusterSim[bestI][x] = merged
      else:         clusterSim[x][bestI] = merged

    for off in members[bestJ]:
      members[bestI].add(off)
    sizes[bestI] += sizes[bestJ]
    active[bestJ] = false

  # Collect active clusters.
  for i in 0 ..< n:
    if active[i]:
      result.add(members[i])

# ---------------------------------------------------------------------------
# Representative selection
# ---------------------------------------------------------------------------

proc selectRepresentative*(cluster: seq[int64];
                            simMap: Table[(int64, int64), float64];
                            passQualMap: Table[int64, tuple[hasPASS: bool; qual: float32]];
                            priority: seq[PriorityCriterion]): int64 =
  ## Pick the representative record from a cluster using the priority cascade.
  ## callerIdx is decoded from the high 16 bits of each composite offset.
  ## Missing passQualMap entries default to (hasPASS=false, qual=0.0).
  if cluster.len == 1: return cluster[0]

  var candidates = cluster

  for crit in priority:
    if candidates.len == 1: break
    case crit
    of pcPass:
      let passOnes = candidates.filterIt(
        passQualMap.getOrDefault(it, (false, 0f32)).hasPASS)
      if passOnes.len > 0: candidates = passOnes
    of pcQual:
      var bestQ = passQualMap.getOrDefault(candidates[0], (false, 0f32)).qual
      for off in candidates:
        let q = passQualMap.getOrDefault(off, (false, 0f32)).qual
        if q > bestQ: bestQ = q
      candidates = candidates.filterIt(
        passQualMap.getOrDefault(it, (false, 0f32)).qual >= bestQ)
    of pcCentre:
      var bestMean = -1.0
      for off in candidates:
        var total = 0.0; var count = 0
        for other in cluster:
          if other == off: continue
          total += simMap.getOrDefault(pairKey(off, other), 0.0); inc count
        let mean = if count > 0: total / float64(count) else: 0.0
        if mean > bestMean: bestMean = mean
      candidates = candidates.filterIt(
        block:
          var t = 0.0; var c = 0
          for other in cluster:
            if other == it: continue
            t += simMap.getOrDefault(pairKey(it, other), 0.0); inc c
          let m = if c > 0: t / float64(c) else: 0.0
          m >= bestMean - 1e-12)
    of pcOrder:
      var bestIdx = high(int)
      for off in candidates:
        let ci = int(off shr 48)
        if ci < bestIdx: bestIdx = ci
      candidates = candidates.filterIt(int(it shr 48) == bestIdx)

  # Final tiebreak: lowest callerIdx, then lowest compositeOff for determinism.
  var best = candidates[0]
  for off in candidates:
    let mci = int(off shr 48); let bci = int(best shr 48)
    if mci < bci: best = off
    elif mci == bci and off < best: best = off
  best

# ---------------------------------------------------------------------------
# Cluster all records
# ---------------------------------------------------------------------------

proc clusterAll*(allOffsets: seq[int64];
                 simMap: Table[(int64, int64), float64];
                 linkage: LinkageMethod;
                 threshold: float64): seq[seq[int64]] =
  ## Build connected components, then agglomerate each independently.
  ## Missing pairs within a component are treated as similarity 0.
  let compId = buildComponents(simMap, allOffsets)
  var byComp: Table[int, seq[int64]]
  for off in allOffsets:
    byComp.mgetOrPut(compId.getOrDefault(off, -1), @[]).add(off)
  for offsets in byComp.values:
    for cl in agglomerateComponent(offsets, simMap, linkage, threshold):
      result.add(cl)
