## mergecore.nim — shared clustering machinery for matcha collapse (and future merge).
##
## Exports:
##   Types: LinkageMethod, PriorityCriterion, CallerInput, FieldConflict,
##          FieldResolution, MergedHeader
##
##   resolveHeaders        — analyse N input VCF headers, produce MergedHeader
##   buildSimilarityMap    — MatchPair seq → (canonicalPair → similarity) table
##   buildComponents       — union-find over pairs → offset→componentId table
##   agglomerateComponent  — Lance-Williams agglomerative clustering
##   clusterAll            — agglomerate every connected component
##   selectRepresentative  — priority-cascade representative selection

import std/[heapqueue, os, sequtils, sets, strutils, tables]
import hts
import hts/private/hts_concat
import utils, preproc, match, log, synced_bcf_reader

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
proc pairKey*(a, b: int32): (int32, int32) {.inline.} =
  if a <= b: (a, b) else: (b, a)

# ---------------------------------------------------------------------------
# BCF sentinel constants and low-level helpers (no nim binding for the C macros)
# ---------------------------------------------------------------------------

const
  bcfInt32Missing*      = low(int32)               ## INT32_MIN
  bcfInt32VectorEnd*    = low(int32) + 1'i32       ## INT32_MIN + 1
  bcfFloatMissingU32*   = 0x7F800001'u32           ## NaN-tagged missing
  bcfFloatVectorEndU32* = 0x7F800002'u32

proc bcfFloatMissing*(): float32 {.inline.} =
  cast[float32](bcfFloatMissingU32)

proc bcfFloatVectorEnd*(): float32 {.inline.} =
  cast[float32](bcfFloatVectorEndU32)

proc hdrInt2Id*(hdr: ptr bcf_hdr_t; typ: cint; id: cint): cstring {.inline.} =
  ## hts-nim keeps bcf_hdr_int2id private; replicate it here.
  let arr = cast[ptr UncheckedArray[bcf_idpair_t]](hdr.id[typ])
  arr[id].key

proc infoFieldIds*(h: ptr bcf_hdr_t): HashSet[string] =
  ## Set of all INFO IDs declared on `h`.
  let hrecs = cast[ptr UncheckedArray[ptr bcf_hrec_t]](h.hrec)
  for i in 0 ..< h.nhrec.int:
    let hr = hrecs[i]
    if hr.`type` != BCF_HEADER_TYPE.BCF_HL_INFO.cint: continue
    let keys = cast[ptr UncheckedArray[cstring]](hr.keys)
    let vals = cast[ptr UncheckedArray[cstring]](hr.vals)
    for j in 0 ..< hr.nkeys.int:
      if $keys[j] == "ID":
        result.incl($vals[j]); break

proc fmtFieldIds*(h: ptr bcf_hdr_t): HashSet[string] =
  ## Set of all FORMAT IDs declared on `h`.
  let hrecs = cast[ptr UncheckedArray[ptr bcf_hrec_t]](h.hrec)
  for i in 0 ..< h.nhrec.int:
    let hr = hrecs[i]
    if hr.`type` != BCF_HEADER_TYPE.BCF_HL_FMT.cint: continue
    let keys = cast[ptr UncheckedArray[cstring]](hr.keys)
    let vals = cast[ptr UncheckedArray[cstring]](hr.vals)
    for j in 0 ..< hr.nkeys.int:
      if $keys[j] == "ID":
        result.incl($vals[j]); break

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

proc resolveHeaders*(callers: seq[CallerInput];
                     infoFilter: seq[string] = @[];
                     fmtFilter:  seq[string] = @[]): MergedHeader =
  ## Open each caller's header, collect INFO and FORMAT definitions, and
  ## produce a MergedHeader that drives the merge pass and output header.
  ##
  ## Compatible INFO/FORMAT (same Number+Type across all callers that define it):
  ##   → single merged def in output; no renames needed.
  ##
  ## Incompatible INFO:
  ##   Number-only conflict (same Type): emit Number=. — silent, no rename.
  ##   Type conflict: rename every instance to FIELD_CALLERNAME; warn if fld
  ##     is in infoFilter (or infoFilter is empty, meaning no output filter).
  ##
  ## Incompatible FORMAT:
  ##   Number-only conflict: emit Number=. — warn if fld is in fmtFilter.
  ##   Type conflict: rename to FIELD_CALLERNAME — warn if fld is in fmtFilter.
  ##
  ## headerLines and infoRes/fmtRes are built unconditionally (slim BCFs may
  ## carry fields beyond the output filter); warnings are gated by the filters.
  ##
  ## Fields that appear in only one caller are kept unchanged.

  let infoFilterSet = toHashSet(infoFilter)
  let fmtFilterSet  = toHashSet(fmtFilter)

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
    let refNum = byCallers[0].number
    let refTyp = byCallers[0].typ
    var compatible  = true
    var typeConflict = false
    for d in byCallers:
      if d.number != refNum or d.typ != refTyp:
        compatible = false
      if d.typ != refTyp:
        typeConflict = true
    if compatible or byCallers.len == 1:
      let desc = mergeDescriptions(byCallers.mapIt(it.desc))
      result.headerLines.add(makeInfoLine(fld, refNum, refTyp, desc))
    elif not typeConflict:
      # Number-only conflict (same Type) → widen to Number=.; no rename, no warning.
      let desc = mergeDescriptions(byCallers.mapIt(it.desc))
      result.headerLines.add(makeInfoLine(fld, ".", refTyp, desc))
    else:
      # Type conflict → rename each caller's instance.
      var res = FieldResolution(kind: fcIncompatibleInfo)
      var warnParts: seq[string]
      for d in byCallers:
        let suffix = callers[d.callerIdx].name
        let newName = fld & "_" & suffix
        res.renames[d.callerIdx] = newName
        result.headerLines.add(makeInfoLine(newName, d.number, d.typ, d.desc))
        warnParts.add(callers[d.callerIdx].name & ":" & d.number & "/" & d.typ)
      result.infoRes[fld] = res
      if infoFilter.len == 0 or fld in infoFilterSet:
        result.warnings.add("INFO Type conflict for " & fld & " (" & warnParts.join(", ") &
                            ") — renamed to " & res.renames.values.toSeq.join(", "))

  # Resolve FORMAT fields.
  for fld, defs in fmtByField.pairs:
    let refNum = defs[0].number
    let refTyp = defs[0].typ
    var compatible = true
    var typeConflict = false
    for d in defs:
      if d.number != refNum or d.typ != refTyp:
        compatible = false
      if d.typ != refTyp:
        typeConflict = true
    if compatible or defs.len == 1:
      let desc = mergeDescriptions(defs.mapIt(it.desc))
      result.headerLines.add(makeFmtLine(fld, refNum, refTyp, desc))
    elif typeConflict:
      # Type conflict → rename each caller's instance (data semantics differ).
      var res = FieldResolution(kind: fcIncompatibleFmt)
      var warnParts: seq[string]
      for d in defs:
        let suffix = callers[d.callerIdx].name
        let newName = fld & "_" & suffix
        res.renames[d.callerIdx] = newName
        result.headerLines.add(makeFmtLine(newName, d.number, d.typ, d.desc))
        warnParts.add(callers[d.callerIdx].name & ":" & d.number & "/" & d.typ)
      result.fmtRes[fld] = res
      if fmtFilter.len == 0 or fld in fmtFilterSet:
        result.warnings.add("FORMAT Type conflict for " & fld & " (" & warnParts.join(", ") &
                            ") — renamed to " & res.renames.values.toSeq.join(", "))
    else:
      # Number-only conflict (compatible Type) → emit Number=. and warn.
      var res = FieldResolution(kind: fcIncompatibleFmt)
      let desc = mergeDescriptions(defs.mapIt(it.desc))
      result.headerLines.add(makeFmtLine(fld, ".", refTyp, desc))
      result.fmtRes[fld] = res
      var warnParts: seq[string]
      for d in defs:
        warnParts.add(callers[d.callerIdx].name & ":" & d.number & "/" & d.typ)
      if fmtFilter.len == 0 or fld in fmtFilterSet:
        result.warnings.add("FORMAT Number conflict for " & fld & " (" & warnParts.join(", ") &
                            ") — emitting Number=.")

# ---------------------------------------------------------------------------
# Similarity map (from MatchPair seq)
# ---------------------------------------------------------------------------

proc buildSimilarityMap*(pairs: seq[MatchPair]): Table[(int32, int32), float64] =
  ## Build canonical-pair → similarity table from match pairs.
  ## Canonical key: (min(srcIndexA, srcIndexB), max(...)).
  ## Singletons (srcIndexB == NO_MATCH) are skipped — they have no B partner.
  for p in pairs:
    if p.srcIndexB == NO_MATCH: continue
    let key = pairKey(p.srcIndexA, p.srcIndexB)
    # Keep highest similarity if duplicate pairs arrive (shouldn't happen, but be safe).
    let cur = result.getOrDefault(key, 0.0)
    if float64(p.sim) > cur:
      result[key] = float64(p.sim)

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

proc buildComponents*(simMap: Table[(int32, int32), float64];
                      allOffsets: seq[int32]): Table[int32, int] =
  ## Union-find over pairs → offset→componentId table.
  ## Singletons get their own component. ComponentIds are not contiguous.
  let n = allOffsets.len
  if n == 0: return
  var offIdx: Table[int32, int]
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
# Agglomerative clustering (Lance-Williams on similarity)
# ---------------------------------------------------------------------------

type
  ClusterState = object
    ## During agglomeration: each cluster is a set of member SRC_INDEX values.
    members: seq[int32]
    size:    int

const AggDenseThreshold* = 256
  ## Components at or below this size use the dense O(N³) reference impl.
  ## Above this, switch to the heap-based sparse impl which is O((N+E) log N).

proc agglomerateDense*(offsets: seq[int32];
                       simMap: Table[(int32, int32), float64];
                       linkage: LinkageMethod;
                       threshold: float64): seq[seq[int32]] =
  ## Reference O(N³) agglomerative clustering. Used for small components and
  ## as the trusted baseline against which agglomerateSparse is regression-tested.
  if offsets.len == 0: return
  if offsets.len == 1:
    return @[@[offsets[0]]]

  let n = offsets.len
  # idx→SRC_INDEX, SRC_INDEX→idx
  var offIdx: Table[int32, int]
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
  var members = newSeq[seq[int32]](n)
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

# Heap entry for sparse agglomeration. Tuple lex-order puts smallest negSim
# (= largest similarity) first, then smallest (i, j) — matching the dense
# scan's row-major tie-break. (vi, vj) ride along for lazy invalidation.
type SparseHeapEntry = tuple[negSim: float64; i, j, vi, vj: int32]

proc agglomerateSparse*(offsets: seq[int32];
                        simMap: Table[(int32, int32), float64];
                        linkage: LinkageMethod;
                        threshold: float64): seq[seq[int32]] =
  ## Heap-based sparse agglomeration with lazy invalidation. O((N + E) log N)
  ## for E = simMap edges within the component. Produces output identical to
  ## agglomerateDense for the same inputs (deterministic tie-break preserved).
  if offsets.len == 0: return
  if offsets.len == 1:
    return @[@[offsets[0]]]

  let n = offsets.len
  var offIdx: Table[int32, int32]
  for i, off in offsets.pairs:
    offIdx[off] = int32(i)

  # Sparse adjacency: neighbors[i][j] = current similarity between active
  # clusters i and j (stored symmetrically for cheap O(1) lookup during merges).
  var neighbors = newSeq[Table[int32, float64]](n)
  var heap: HeapQueue[SparseHeapEntry]
  for key, sim in simMap.pairs:
    let ia = offIdx.getOrDefault(key[0], -1'i32)
    let ib = offIdx.getOrDefault(key[1], -1'i32)
    if ia < 0 or ib < 0: continue  # endpoint not in this component
    neighbors[ia][ib] = sim
    neighbors[ib][ia] = sim
    if sim >= threshold:
      let lo = if ia < ib: ia else: ib
      let hi = if ia < ib: ib else: ia
      heap.push((-sim, lo, hi, 0'i32, 0'i32))

  var members = newSeq[seq[int32]](n)
  var sizes = newSeq[int](n)
  var active = newSeq[bool](n)
  var version = newSeq[int32](n)
  for i in 0 ..< n:
    members[i] = @[offsets[i]]
    sizes[i] = 1
    active[i] = true

  while heap.len > 0:
    let top = heap.pop()
    let i = top.i
    let j = top.j
    if not active[i] or not active[j]: continue
    if version[i] != top.vi or version[j] != top.vj: continue
    let sim = -top.negSim
    if sim < threshold: break  # no more merges above threshold

    # Merge cluster j into cluster i (keep i, deactivate j). Lance-Williams
    # update for every neighbor of either i or j; missing edges treated as 0.0
    # (matching the dense path's matrix-of-zeros init).
    let sA = float64(sizes[i])
    let sB = float64(sizes[j])
    var nbrs: HashSet[int32]
    for x in neighbors[i].keys: nbrs.incl(x)
    for x in neighbors[j].keys: nbrs.incl(x)
    nbrs.excl(i); nbrs.excl(j)

    # Invalidate any stale heap entries that referenced cluster i (its
    # similarity to every other cluster is about to change). version[j] does
    # not need bumping — j becomes inactive and is filtered on pop.
    inc version[i]

    for x in nbrs:
      let dAX = neighbors[i].getOrDefault(x, 0.0)
      let dBX = neighbors[j].getOrDefault(x, 0.0)
      let merged =
        case linkage
        of lmSingle:   max(dAX, dBX)
        of lmComplete: min(dAX, dBX)
        of lmAverage:  (sA * dAX + sB * dBX) / (sA + sB)
      if merged > 0:
        neighbors[i][x] = merged
        neighbors[x][i] = merged
      else:
        neighbors[i].del(x)
        neighbors[x].del(i)
      neighbors[x].del(j)
      if merged >= threshold:
        let lo = if i < x: i else: x
        let hi = if i < x: x else: i
        heap.push((-merged, lo, hi, version[lo], version[hi]))

    for off in members[j]:
      members[i].add(off)
    sizes[i] += sizes[j]
    active[j] = false
    neighbors[j].clear()  # free per-cluster adjacency once inactive

  for i in 0 ..< n:
    if active[i]:
      result.add(members[i])

proc agglomerateComponent*(offsets: seq[int32];
                            simMap: Table[(int32, int32), float64];
                            linkage: LinkageMethod;
                            threshold: float64): seq[seq[int32]] =
  ## Agglomerative clustering of one connected component.
  ## Dispatches to a dense O(N³) reference impl for small components and a
  ## heap-based sparse impl for large ones. Output is identical either way.
  if offsets.len <= AggDenseThreshold:
    agglomerateDense(offsets, simMap, linkage, threshold)
  else:
    agglomerateSparse(offsets, simMap, linkage, threshold)

# ---------------------------------------------------------------------------
# Representative selection
# ---------------------------------------------------------------------------

proc selectRepresentative*(cluster: seq[int32];
                            simMap: Table[(int32, int32), float64];
                            passQualMap: Table[int32, tuple[hasPASS: bool; qual: uint16; callerIdx: int32]];
                            priority: seq[PriorityCriterion]): int32 =
  ## Pick the representative record from a cluster using the priority cascade.
  ## callerIdx comes from passQualMap (populated from MatchPair.callerIdxA).
  ## Missing passQualMap entries default to (hasPASS=false, qual=0, callerIdx=high(int32)).
  let missing = (hasPASS: false, qual: 0'u16, callerIdx: high(int32))
  if cluster.len == 1: return cluster[0]

  var candidates = cluster

  for crit in priority:
    if candidates.len == 1: break
    case crit
    of pcPass:
      let passOnes = candidates.filterIt(
        passQualMap.getOrDefault(it, missing).hasPASS)
      if passOnes.len > 0: candidates = passOnes
    of pcQual:
      var bestQ = passQualMap.getOrDefault(candidates[0], missing).qual
      for idx in candidates:
        let q = passQualMap.getOrDefault(idx, missing).qual
        if q > bestQ: bestQ = q
      candidates = candidates.filterIt(
        passQualMap.getOrDefault(it, missing).qual >= bestQ)
    of pcCentre:
      var bestMean = -1.0
      for idx in candidates:
        var total = 0.0; var count = 0
        for other in cluster:
          if other == idx: continue
          total += simMap.getOrDefault(pairKey(idx, other), 0.0); inc count
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
      var bestCi = high(int32)
      for idx in candidates:
        let ci = passQualMap.getOrDefault(idx, missing).callerIdx
        if ci < bestCi: bestCi = ci
      let bestCiFinal = bestCi
      candidates = candidates.filterIt(
        passQualMap.getOrDefault(it, missing).callerIdx == bestCiFinal)

  # Final tiebreak: lowest callerIdx, then lowest SRC_INDEX for determinism.
  var best = candidates[0]
  for idx in candidates:
    let mci = passQualMap.getOrDefault(idx, missing).callerIdx
    let bci = passQualMap.getOrDefault(best, missing).callerIdx
    if mci < bci: best = idx
    elif mci == bci and idx < best: best = idx
  best

# ---------------------------------------------------------------------------
# Cluster all records
# ---------------------------------------------------------------------------

const LargeComponentWarn* = 500
  ## Components at or above this size trigger a logWarn flagging the chrom,
  ## svtype, and dominant caller. Surface for caller pathology like delly's
  ## known chr2 LowQual SV flood on some samples.

proc clusterAll*(byComp: Table[int, seq[int32]];
                 simMap: Table[(int32, int32), float64];
                 linkage: LinkageMethod;
                 threshold: float64): seq[seq[int32]] =
  ## Agglomerate each connected component independently.
  ## Missing pairs within a component are treated as similarity 0.
  for offsets in byComp.values:
    for cl in agglomerateComponent(offsets, simMap, linkage, threshold):
      result.add(cl)

# ---------------------------------------------------------------------------
# selfMatchAndCluster — full post-integratedMerge pipeline.
#
# Both `matcha collapse` and `matcha merge` share this body: self-match the
# merged slim BCFs with singleton emission, build the similarity map and
# location lookup, cluster, fetch PASS/QUAL/CALLER_IDX per cluster member
# via grouped CSI queries, and select representatives.
# ---------------------------------------------------------------------------

type
  ClusterPipelineResult* = object
    finalClusters*: seq[seq[int32]]  ## each cluster: representative first
    locByIdx*:      Table[int32, tuple[chromIdx: int16; pos: int32; fileIdx: int16]]
    passQualMap*:   Table[int32, tuple[hasPASS: bool; qual: uint16; callerIdx: int32]]
    fileList*:      seq[string]      ## file index → slim BCF path

proc selfMatchAndCluster*(mergedPreproc: PreprocOutput;
                          matchCfg: MatchConfig;
                          linkage: LinkageMethod;
                          threshold: float64;
                          priority: seq[PriorityCriterion];
                          modeTag = "self-match"): ClusterPipelineResult =
  let (jobs, fileList) = buildWorkQueue(mergedPreproc, mergedPreproc, matchCfg)
  logVerbose(modeTag & ": " & $jobs.len & " job(s)")
  var allPairs: seq[MatchPair]
  for jrs in runMatchPairJobsWithPool(jobs, matchCfg):
    for mp in jrs: allPairs.add(mp)
  logVerbose(modeTag & ": " & $allPairs.len & " pair(s)")

  let simMap = buildSimilarityMap(allPairs)
  var seenOffsets: HashSet[int32]
  var allOffsets: seq[int32]
  result.locByIdx = initTable[int32,
                              tuple[chromIdx: int16; pos: int32; fileIdx: int16]]()
  for p in allPairs:
    if p.srcIndexA notin seenOffsets:
      seenOffsets.incl(p.srcIndexA)
      allOffsets.add(p.srcIndexA)
      result.locByIdx[p.srcIndexA] = (p.chromIdx, p.posA, p.fileIdxA)
    if p.srcIndexB != NO_MATCH and p.srcIndexB notin seenOffsets:
      seenOffsets.incl(p.srcIndexB)
      allOffsets.add(p.srcIndexB)
      result.locByIdx[p.srcIndexB] = (p.chromIdx, p.posB, p.fileIdxB)
  logVerbose(modeTag & ": " & $allOffsets.len & " unique record(s)")

  # Build connected components and agglomerate each. For components at or
  # above LargeComponentWarn, emit a logWarn (lazily building the per-srcIndex
  # caller metadata the first time a large component is encountered).
  let compId = buildComponents(simMap, allOffsets)
  var byComp: Table[int, seq[int32]]
  for off in allOffsets:
    byComp.mgetOrPut(compId.getOrDefault(off, -1), @[]).add(off)

  type Meta = tuple[chromIdx: int16; svtype: int8; callerIdx: int16]
  var meta: Table[int32, Meta]
  var metaBuilt = false
  var clusters: seq[seq[int32]]
  for offsets in byComp.values:
    if offsets.len >= LargeComponentWarn:
      if not metaBuilt:
        for p in allPairs:
          if p.srcIndexA notin meta:
            meta[p.srcIndexA] = (p.chromIdx, p.svtype, p.callerIdxA)
        metaBuilt = true
      let repMeta = meta.getOrDefault(offsets[0])
      var callerCounts: Table[int16, int]
      for srcIdx in offsets:
        let ci = meta.getOrDefault(srcIdx).callerIdx
        callerCounts[ci] = callerCounts.getOrDefault(ci, 0) + 1
      var callerK: int16 = 0
      var callerBest = -1
      for k, c in callerCounts.pairs:
        if c > callerBest: callerBest = c; callerK = k
      let chromName = if repMeta.chromIdx.int < mergedPreproc.chromOrder.len:
                        mergedPreproc.chromOrder[repMeta.chromIdx.int]
                      else: "?"
      logWarn(modeTag & ": large cluster component: " & chromName & "/" &
              $SvType(repMeta.svtype) & " N=" & $offsets.len &
              " dominant=caller" & $callerK & ":" &
              formatFloat(100.0 * float(callerBest) / float(offsets.len),
                          ffDecimal, 2) & "% — possible caller artifact")
    for cl in agglomerateComponent(offsets, simMap, linkage, threshold):
      clusters.add(cl)

  # Build passQualMap from MatchPair metadata stamped by matchcore.
  # Every offset appears as srcIndexA in at least one pair (emitSingletons=true
  # guarantees a singleton entry for unmatched records), so A-side fields suffice.
  for p in allPairs:
    if p.srcIndexA notin result.passQualMap:
      result.passQualMap[p.srcIndexA] = (hasPASS: p.passA,
                                          qual: p.qualQ,
                                          callerIdx: int32(p.callerIdxA))

  for cl in clusters:
    let rep = selectRepresentative(cl, simMap, result.passQualMap, priority)
    var ordered = @[rep]
    for idx in cl:
      if idx != rep: ordered.add(idx)
    result.finalClusters.add(ordered)
  logInfo(modeTag & ": " & $allOffsets.len & " record(s) -> " & $result.finalClusters.len & " cluster(s)")
  result.fileList = fileList

# ---------------------------------------------------------------------------
# integratedMerge — fused preproc + merge in one synced_bcf_reader pass
# ---------------------------------------------------------------------------

type
  MergeStreamConfig* = object
    ## Configuration for integratedMerge; decoupled from CollapseConfig so the
    ## same streaming kernel can be reused by matcha merge.
    formatFields*:      seq[string]   ## FORMAT fields to carry; empty = no FORMAT
    nThreads*:          int
    tmpDir*:            string
    stampSID*:          bool          ## merge mode: write FORMAT/SID per record
    sampleIdByCaller*:  seq[string]   ## merge mode: SID value for caller ci
    preserveBndAlt*:    bool          ## merge mode: keep source BND ALT verbatim
    preserveInsAlt*:    bool          ## keep source INS REF/ALT verbatim (sequence-resolved)
    keptChrs*:          HashSet[string] ## --chrs filter; empty = no filter

# Header utility: collect ##FILTER=<...> lines from a header.
proc collectFilterLines*(h: ptr bcf_hdr_t): seq[string] =
  let hrecs = cast[ptr UncheckedArray[ptr bcf_hrec_t]](h.hrec)
  for i in 0 ..< h.nhrec.int:
    let hr = hrecs[i]
    if hr.`type` != BCF_HEADER_TYPE.BCF_HL_FLT.cint: continue
    let keys = cast[ptr UncheckedArray[cstring]](hr.keys)
    let vals = cast[ptr UncheckedArray[cstring]](hr.vals)
    var line = "##FILTER=<"
    for j in 0 ..< hr.nkeys.int:
      if j > 0: line &= ","
      line &= $keys[j] & "=" & $vals[j]
    line &= ">"
    result.add(line)

const AlwaysKeepInMerged* = ["SVTYPE", "SVLEN", "END", "CHR2", "POS2",
                              "CALLERS", "N_CALLERS", "N_MERGED",
                              "SRC_INDEX", "CALLER_IDX"]

# ---------------------------------------------------------------------------
# Header skeleton builders — shared by collapse.buildFinalHdr,
# merge.buildSlimHdr, and merge.buildOutputHdr. Each callsite handles its
# mode-specific additions (sample columns, cohort/provenance INFO defs,
# `##source` / `##matcha_cmdline` lines) on top.
# ---------------------------------------------------------------------------

proc addContigsUnion*(hdr: ptr bcf_hdr_t; callers: seq[CallerInput];
                     keptChrs: HashSet[string] = initHashSet[string]()) =
  ## Append `##contig=<ID=...>` for every contig seen across callers; first
  ## caller wins for order. When `keptChrs` is non-empty, contigs not in the
  ## set are skipped.
  var seen: HashSet[string]
  for caller in callers:
    var vcf: VCF
    if not open(vcf, caller.path): continue
    var n: cint = 0
    let names = bcf_hdr_seqnames(vcf.header.hdr, n.addr)
    if names != nil:
      for i in 0 ..< n.int:
        let c = $names[i]
        if c notin seen:
          seen.incl(c)
          if keptChrs.len > 0 and c notin keptChrs: continue
          discard bcf_hdr_append(hdr, ("##contig=<ID=" & c & ">").cstring)
      free(names)
    vcf.close()

proc addFiltersUnion*(hdr: ptr bcf_hdr_t; callers: seq[CallerInput]) =
  ## Append the union of `##FILTER=<...>` lines across callers (first occurrence
  ## wins; duplicates from later callers are dropped).
  var seen: HashSet[string]
  for caller in callers:
    var vcf: VCF
    if not open(vcf, caller.path): continue
    for line in collectFilterLines(vcf.header.hdr):
      if line notin seen:
        seen.incl(line)
        discard bcf_hdr_append(hdr, line.cstring)
    vcf.close()

proc addStandardSvInfoDefs*(hdr: ptr bcf_hdr_t) =
  ## Append SVTYPE/SVLEN/END/CHR2/POS2 INFO defs that aren't already declared.
  let have = infoFieldIds(hdr)
  if "SVTYPE" notin have:
    discard bcf_hdr_append(hdr,
      "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">".cstring)
  if "SVLEN" notin have:
    discard bcf_hdr_append(hdr,
      "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the SV (absolute value)\">".cstring)
  if "END" notin have:
    discard bcf_hdr_append(hdr,
      "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the SV (1-based, inclusive)\">".cstring)
  if "CHR2" notin have:
    discard bcf_hdr_append(hdr,
      "##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome of mate breakend\">".cstring)
  if "POS2" notin have:
    discard bcf_hdr_append(hdr,
      "##INFO=<ID=POS2,Number=1,Type=Integer,Description=\"Position of mate breakend\">".cstring)

proc addHeaderLinesFiltered*(hdr: ptr bcf_hdr_t; mh: MergedHeader;
                              keepInfo: proc (id: string): bool {.closure.};
                              keepFmt:  proc (id: string): bool {.closure.}) =
  ## Walk `mh.headerLines` (synthesized ##INFO=<...> / ##FORMAT=<...>) and
  ## append those whose ID satisfies the relevant predicate. Non-ID lines
  ## are appended unconditionally.
  for line in mh.headerLines:
    let idStart = line.find("ID=")
    if idStart < 0:
      discard bcf_hdr_append(hdr, line.cstring)
      continue
    let idEnd = line.find(',', idStart + 3)
    let fieldId =
      if idEnd > 0: line[idStart + 3 ..< idEnd]
      else:         line[idStart + 3 ..< line.len - 1]
    if line.startsWith("##FORMAT"):
      if keepFmt(fieldId):
        discard bcf_hdr_append(hdr, line.cstring)
    else:
      if keepInfo(fieldId):
        discard bcf_hdr_append(hdr, line.cstring)

proc keepInfoForMerged*(name: string; infoFilter: seq[string]): bool =
  ## Keep set for INFO fields in the merged slim BCFs.
  ## - Always keeps matcha-internal + matchcore-required fields.
  ## - With infoFilter empty: keep all non-internal fields.
  ## - With infoFilter set: keep listed fields (or their *_<caller> renames).
  for n in AlwaysKeepInMerged:
    if name == n: return true
  if infoFilter.len == 0: return true
  for tok in infoFilter:
    if name == tok or name.startsWith(tok & "_"): return true
  false

proc btToHt(bt: cint): cint =
  if bt == BCF_BT_FLOAT: BCF_HT_REAL.cint
  elif bt == BCF_BT_CHAR: BCF_HT_STR.cint
  else: BCF_HT_INT.cint

proc applyInfoRename*(srcHdr: ptr bcf_hdr_t; rec: ptr bcf1_t;
                      origName, newName: string;
                      buf: var pointer; bufN: var cint) =
  let info = bcf_get_info(srcHdr, rec, origName.cstring)
  if info == nil: return
  let htType = btToHt(info.`type`)
  let n = bcf_get_info_values(srcHdr, rec, origName.cstring,
                              buf.addr, bufN.addr, htType)
  if n > 0:
    discard bcf_update_info(srcHdr, rec, newName.cstring, buf, n.cint, htType)
  discard bcf_update_info(srcHdr, rec, origName.cstring, nil, 0.cint, htType)

proc applyFmtRename*(srcHdr: ptr bcf_hdr_t; rec: ptr bcf1_t;
                     origName, newName: string;
                     buf: var pointer; bufN: var cint) =
  let fmt = bcf_get_fmt(srcHdr, rec, origName.cstring)
  if fmt == nil: return
  let htType = btToHt(fmt.`type`)
  let n = bcf_get_format_values(srcHdr, rec, origName.cstring,
                                buf.addr, bufN.addr, htType)
  if n > 0:
    discard bcf_update_format(srcHdr, rec, newName.cstring, buf, n.cint, htType)
  discard bcf_update_format(srcHdr, rec, origName.cstring, nil, 0.cint, htType)

proc augmentSrcHdrForRenames*(srcHdr: ptr bcf_hdr_t; ci: int;
                               mh: MergedHeader) =
  ## For each renamed INFO/FORMAT field this caller writes under, append the
  ## renamed def to srcHdr so subsequent bcf_update_info/_format calls have
  ## the field's Number/Type available. Must be done before any rename write.
  for origName, res in mh.infoRes.pairs:
    if res.kind != fcIncompatibleInfo: continue
    if ci notin res.renames: continue
    let newName = res.renames[ci]
    for hi in collectHrecs(srcHdr, BCF_HEADER_TYPE.BCF_HL_INFO.cint):
      if hi.id == origName:
        discard bcf_hdr_append(srcHdr,
          ("##INFO=<ID=" & newName & ",Number=" & hi.number & ",Type=" & hi.typ &
           ",Description=\"" & hi.desc & "\">").cstring)
        break
  for origName, res in mh.fmtRes.pairs:
    if res.kind != fcIncompatibleFmt: continue
    if ci notin res.renames: continue
    let newName = res.renames[ci]
    for hi in collectHrecs(srcHdr, BCF_HEADER_TYPE.BCF_HL_FMT.cint):
      if hi.id == origName:
        discard bcf_hdr_append(srcHdr,
          ("##FORMAT=<ID=" & newName & ",Number=" & hi.number & ",Type=" & hi.typ &
           ",Description=\"" & hi.desc & "\">").cstring)
        break
  # Also ensure SRC_INDEX / CALLER_IDX / SV defs exist for the
  # authoritative writes issued against srcHdr in the per-record loop.
  proc ensureInfo(h: ptr bcf_hdr_t; name, num, typ, desc: string) =
    for hi in collectHrecs(h, BCF_HEADER_TYPE.BCF_HL_INFO.cint):
      if hi.id == name: return
    discard bcf_hdr_append(h,
      ("##INFO=<ID=" & name & ",Number=" & num & ",Type=" & typ &
       ",Description=\"" & desc & "\">").cstring)
  ensureInfo(srcHdr, "SRC_INDEX",  "1", "Integer", "matcha-internal: sequential record index")
  ensureInfo(srcHdr, "CALLER_IDX", "1", "Integer", "matcha-internal: caller index (0-based)")
  ensureInfo(srcHdr, "SVTYPE",      "1", "String",  "Type of structural variant")
  ensureInfo(srcHdr, "SVLEN",       "1", "Integer", "Length of the SV")
  ensureInfo(srcHdr, "END",         "1", "Integer", "End position of the SV")
  ensureInfo(srcHdr, "CHR2",        "1", "String",  "Chromosome of mate breakend")
  ensureInfo(srcHdr, "POS2",        "1", "Integer", "Position of mate breakend")
  discard bcf_hdr_sync(srcHdr)

type IntegratedMergeResult* = object
  paths*:     Table[SvtypeBin, string]
  populated*: Table[SvtypeBin, HashSet[string]]

proc integratedMerge*(callers: seq[CallerInput]; mh: MergedHeader;
                      finalHdr: ptr bcf_hdr_t;
                      cfg: MergeStreamConfig;
                      chromOrder: seq[string]): IntegratedMergeResult =
  ## Stream all N caller VCFs via one synced_bcf_reader, normalize each
  ## record, filter INFO/FORMAT to the user-selected fields, write per-
  ## (svtype, bin) merged slim BCFs. finalHdr is pre-built by the caller
  ## (collapse builds it with buildFinalHdr; merge will supply its own).

  # 1. Init synced_bcf_reader + thread pool.
  let sr = bcf_sr_init()
  if sr == nil:
    raise newException(IOError, "bcf_sr_init failed")
  discard bcf_sr_set_opt(sr, BCF_SR_REQUIRE_IDX)

  let nThr = max(1, cfg.nThreads)
  var tpool = htsThreadPool(pool: nil)
  if nThr >= 2:
    tpool.pool = hts_tpool_init(nThr.cint)

  for caller in callers:
    if bcf_sr_add_reader(sr, caller.path.cstring) == 0:
      let msg = $bcf_sr_strerror(srs_errnum(sr))
      bcf_sr_destroy(sr)
      if tpool.pool != nil: hts_tpool_destroy(tpool.pool)
      raise newException(IOError, "cannot open: " & caller.path & " (" & msg & ")")

  let nreaders = srs_nreaders(sr).int

  # Attach thread pool to each reader's underlying file.
  if tpool.pool != nil:
    for i in 0 ..< nreaders:
      discard hts_set_opt(srs_get_file(sr, i.cint),
                          srs_hts_opt_thread_pool(), tpool.addr)

  # Subset samples on each reader's header (one sample per caller for
  # collapse, or none).
  for ci in 0 ..< nreaders:
    let srcHdr = srs_get_header(sr, ci.cint)
    if cfg.formatFields.len == 0:
      discard bcf_hdr_set_samples(srcHdr, cstring(nil), 0.cint)
    else:
      let nsamp = bcf_hdr_nsamples(srcHdr).int
      if nsamp > 0:
        let firstSample = $cast[cstringArray](srcHdr.samples)[0]
        discard bcf_hdr_set_samples(srcHdr, firstSample.cstring, 0.cint)
      else:
        discard bcf_hdr_set_samples(srcHdr, cstring(nil), 0.cint)
    augmentSrcHdrForRenames(srcHdr, ci, mh)
    if cfg.stampSID:
      # SID FORMAT def must exist on srcHdr so bcf_update_format can resolve
      # the key; also on finalHdr (which is the slimHdr in merge mode).
      var sidPresent = false
      for hi in collectHrecs(srcHdr, BCF_HEADER_TYPE.BCF_HL_FMT.cint):
        if hi.id == "SID": sidPresent = true; break
      if not sidPresent:
        discard bcf_hdr_append(srcHdr,
          ("##FORMAT=<ID=SID,Number=1,Type=String," &
           "Description=\"matcha-internal: source sample ID\">").cstring)
        discard bcf_hdr_sync(srcHdr)

  # Per-caller WarnState + lineno counter (for synthetic ID generation).
  # globalSrcIndex is a single counter across all callers for SRC_INDEX writes.
  var wsList: seq[WarnState]
  var perCallerLineno = newSeq[int](nreaders)
  var globalSrcIndex: int32 = 0
  for i in 0 ..< nreaders:
    wsList.add(initWarnState(callers[i].name))

  # Keep sets derived from finalHdr: post-filter records contain exactly the
  # fields finalHdr defines → bcf_translate is clean.
  var infoKeepSet, fmtKeepSet: HashSet[string]
  for hi in collectHrecs(finalHdr, BCF_HEADER_TYPE.BCF_HL_INFO.cint):
    infoKeepSet.incl(hi.id)
  for hi in collectHrecs(finalHdr, BCF_HEADER_TYPE.BCF_HL_FMT.cint):
    fmtKeepSet.incl(hi.id)

  # Reusable Variant view over the synced reader's records.
  let view = newVariantView()
  var svtypeBuf: string
  var endBuf, svlenBuf, inslenBuf: seq[int32]
  var leftSeqBuf, rightSeqBuf: string
  var renameBuf: pointer = nil
  var renameN:   cint    = 0

  # Lazy-opened writers + CSI indexes per (svtype, bin).
  var writers:    Table[SvtypeBin, VCF]
  var indexes:    Table[SvtypeBin, ptr hts_idx_t]
  var writerHdrs: Table[SvtypeBin, ptr bcf_hdr_t]

  try:
    # 2. Stream records.
    while bcf_sr_next_line(sr) > 0:
      for ci in 0 ..< nreaders:
        if srs_has_line(sr, ci.cint) == 0: continue
        let srcHdr = srs_get_header(sr, ci.cint)
        let rawRec = srs_get_line(sr, ci.cint)

        if cfg.keptChrs.len > 0 and
           getChromName(srcHdr, rawRec.rid) notin cfg.keptChrs:
          inc wsList[ci].nRead
          inc wsList[ci].skipped[skChromFiltered]
          continue

        let rec = bcf_dup(rawRec)
        inc perCallerLineno[ci]

        let nr = normalizeRecord(srcHdr, rec, perCallerLineno[ci],
                                 wsList[ci], view, svtypeBuf,
                                 endBuf, svlenBuf, inslenBuf,
                                 leftSeqBuf, rightSeqBuf,
                                 callers[ci].path)
        if not nr.ok:
          bcf_destroy(rec)
          continue
        if cfg.keptChrs.len > 0 and nr.svt == svBND and
           nr.bndChr2 notin cfg.keptChrs:
          inc wsList[ci].skipped[skChromFiltered]
          # normalizeRecord already incremented nKept; back it out.
          dec wsList[ci].nKept
          bcf_destroy(rec)
          continue

        # 2a. Apply INFO renames.
        for origName, res in mh.infoRes.pairs:
          if res.kind != fcIncompatibleInfo: continue
          if ci notin res.renames: continue
          applyInfoRename(srcHdr, rec, origName, res.renames[ci],
                          renameBuf, renameN)

        # 2b. Apply FORMAT renames.
        for origName, res in mh.fmtRes.pairs:
          if res.kind != fcIncompatibleFmt: continue
          if ci notin res.renames: continue
          applyFmtRename(srcHdr, rec, origName, res.renames[ci],
                         renameBuf, renameN)

        # 2c+2d. Filter INFO and FORMAT fields.
        discard bcf_unpack(rec, BCF_UN_ALL.cint)
        let idPairs = cast[ptr UncheckedArray[bcf_idpair_t]](srcHdr.id[BCF_DT_ID])
        var infoToDel: seq[string]
        let nInfo = rec.n_info.int
        if nInfo > 0:
          let infoArr = cast[ptr UncheckedArray[bcf_info_t]](rec.d.info)
          for i in 0 ..< nInfo:
            let key = infoArr[i].key
            if key < 0: continue
            let nm = $idPairs[key].key
            if nm notin infoKeepSet: infoToDel.add(nm)
        for nm in infoToDel:
          discard bcf_update_info(srcHdr, rec, nm.cstring,
                                  nil, 0.cint, BCF_HT_INT.cint)
        var fmtToDel: seq[string]
        let nFmt = rec.n_fmt.int
        if nFmt > 0:
          let fmtArr = cast[ptr UncheckedArray[bcf_fmt_t]](rec.d.fmt)
          for i in 0 ..< nFmt:
            let id = fmtArr[i].id
            if id < 0: continue
            let nm = $idPairs[id].key
            if nm notin fmtKeepSet: fmtToDel.add(nm)
        for nm in fmtToDel:
          discard bcf_update_format(srcHdr, rec, nm.cstring,
                                    nil, 0.cint, BCF_HT_INT.cint)

        # 2e. Authoritative writes for END / CHR2 / POS2 / SVTYPE / SVLEN.
        var svtStr = $nr.svt
        discard bcf_update_info(srcHdr, rec, "SVTYPE".cstring,
                                svtStr[0].addr, svtStr.len.cint, BCF_HT_STR.cint)
        var svlenVal = nr.svlen.int32
        discard bcf_update_info(srcHdr, rec, "SVLEN".cstring,
                                svlenVal.addr, 1.cint, BCF_HT_INT.cint)
        if nr.svt == svBND:
          var chr2Str = nr.bndChr2
          if chr2Str.len > 0:
            discard bcf_update_info(srcHdr, rec, "CHR2".cstring,
                                    chr2Str[0].addr, chr2Str.len.cint, BCF_HT_STR.cint)
          var pos2Val = nr.bndPos2.int32
          discard bcf_update_info(srcHdr, rec, "POS2".cstring,
                                  pos2Val.addr, 1.cint, BCF_HT_INT.cint)
        else:
          var endVal = nr.endPos.int32
          discard bcf_update_info(srcHdr, rec, "END".cstring,
                                  endVal.addr, 1.cint, BCF_HT_INT.cint)

        # 2f. Write SRC_INDEX (global sequential) and CALLER_IDX (which caller).
        var srcIdxVal = globalSrcIndex
        inc globalSrcIndex
        var callerIdxVal = ci.int32
        discard bcf_update_info(srcHdr, rec, "SRC_INDEX".cstring,
                                srcIdxVal.addr, 1.cint, BCF_HT_INT.cint)
        discard bcf_update_info(srcHdr, rec, "CALLER_IDX".cstring,
                                callerIdxVal.addr, 1.cint, BCF_HT_INT.cint)

        # 2g. Write FORMAT/SID (merge mode only) — source sample identity for
        # the dummy single-sample slim record.
        if cfg.stampSID and ci < cfg.sampleIdByCaller.len:
          var sidArr: array[1, cstring]
          sidArr[0] = cfg.sampleIdByCaller[ci].cstring
          discard bcf_update_format_string(srcHdr, rec, "SID".cstring,
                                           cast[cstringArray](sidArr[0].addr),
                                           1.cint)

        # 2h. REF/ALT trim to keep records small. For BND in merge mode we
        # preserve the source ALT verbatim (strand orientation lives in the
        # bracket form and is not derivable from CHR2/POS2 at output time).
        if (cfg.preserveBndAlt and nr.svt == svBND) or
           (cfg.preserveInsAlt and nr.svt == svINS):
          discard  # leave REF/ALT as the source had them
        else:
          discard bcf_update_alleles_str(srcHdr, rec, "N,.".cstring)

        # 2i. Track metadata + lazy-open writer.
        let mKey: SvtypeBin = (nr.svt, nr.binIdx)
        result.populated.mgetOrPut(mKey, initHashSet[string]()).incl(
          getChromName(srcHdr, rec.rid))

        let key: SvtypeBin = (nr.svt, nr.binIdx)
        if key notin writers:
          let outPath = cfg.tmpDir / "matcha_" & $getCurrentProcessId() &
                        "_M_" & $nr.svt & "_b" & $nr.binIdx & ".bcf"
          var wtr: VCF
          if not open(wtr, outPath, mode = "wb"):
            raise newException(IOError, "cannot create merged BCF: " & outPath)
          # Init wtr.header via a dummy source, then replace with dup of finalHdr
          # (each writer owns its dup; finalHdr remains valid for writeOutput).
          var dummy: VCF
          discard open(dummy, callers[0].path)
          wtr.copy_header(dummy.header)
          dummy.close()
          bcf_hdr_destroy(wtr.header.hdr)
          let wHdr = bcf_hdr_dup(finalHdr)
          wtr.header.hdr = wHdr
          writerHdrs[key] = wHdr
          discard wtr.write_header()
          writers[key] = wtr
          result.paths[key] = outPath
          let headerOff = uint64(bgzf_tell(bgzfHandle(wtr)))
          let idx = hts_idx_init(0.cint, HTS_FMT_CSI.cint, headerOff,
                                 14.cint, 5.cint)
          if idx == nil:
            raise newException(IOError, "cannot create CSI index: " & outPath)
          indexes[key] = idx

        # 2j. Translate field IDs to finalHdr space, write record, push to index.
        discard bcf_translate(writerHdrs[key], srcHdr, rec)
        discard bcf_write(vcfHtsFile(writers[key]), writerHdrs[key], rec)
        let woff = uint64(bgzf_tell(bgzfHandle(writers[key])))
        discard hts_idx_push(indexes[key], rec.rid, int64(rec.pos),
                             int64(rec.pos) + int64(rec.rlen), woff, 1.cint)

        bcf_destroy(rec)

  finally:
    # Teardown order: destroy synced reader first (it owns reader file handles),
    # then close writers (each holds its own dup of finalHdr), then destroy
    # the shared thread pool last (after every htsFile* that referenced it).
    bcf_sr_destroy(sr)
    for key, wtr in writers.mpairs:
      let finalOff = uint64(bgzf_tell(bgzfHandle(wtr)))
      hts_idx_finish(indexes[key], finalOff)
      wtr.close()  # bcf_hdr_destroy on wtr.header.hdr (the dup) — finalHdr is safe
      let path = result.paths[key]
      hts_idx_save(indexes[key], path.cstring, HTS_FMT_CSI.cint)
      hts_idx_destroy(indexes[key])
    if tpool.pool != nil: hts_tpool_destroy(tpool.pool)

    if renameBuf != nil: c_free(renameBuf)

  # Emit per-caller summaries.
  for ws in wsList: emitSummary(ws)
