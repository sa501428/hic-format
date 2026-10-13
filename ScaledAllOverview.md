# Stored All overview conventions for V10

These are application conventions using existing V10 chromosome records,
attributes, resolutions, matrices, blocks and directories. They introduce no
binary fields, flags, enum values or file-version changes. The normative wire
format remains [HiCFormatV10.md](HiCFormatV10.md).

Native producers now write convention **version 2**, using the V9-style genome
span in kilobases and unpadded chromosome boundaries. Readers also support the
earlier compact convention **version 1**, and legacy All without these
attributes. Unknown convention versions retain ordinary core V10 decoding;
applications must not assume either convention's display semantics.

## Version 2: V9-style span, finest-resolution storage

Let `L_i` be each real chromosome length in header order, skipping All, and
`T = sum(L_i)`. The coordinate scale `C` is normally **1000**, exactly as in V9.
For unusually fragmented or enormous genomes it is increased to
`max(1000, ceil(T / max(L_i)), ceil(T / INT32_MAX))` to preserve expected-vector
lengths and desktop boundary-array compatibility. The stored All chromosome
length is `K = max(1, floor(T / C))`. For the hg38 chromosome table used in the
sample, `K = 3,088,286`, matching V9. Matching requires the same real chromosome
names, order and lengths in both files.

The useful All matrix remains stored at the **smallest BP bin size already
advertised**, `F`. The ordinary global resolution list is unchanged. These
finest-resolution entries are sparse overview measurements, not fine-resolution
whole-genome measurements.

For a requested axis budget `B` (default 500), use
`H = F * max(1, ceil(max(1, floor(K / B)) / F))`. `H` is the display-bin width in
pseudo-coordinate units, and `D = H / F` is the storage stride. If needed, raise
`H` to a multiple of `F` so `N = ceil(K / H) <= B + 1`. This follows V9's
`floor(All length / 500)` binning on ordinary large genomes, with rounding when
that interval is not representable at the finest advertised resolution.

For example, hg38 with `F=1` gives `H=6176` and 501 display bins at the default
budget. The same genome with `F=100` gives `H=6200` and 499 display bins. Both
have **All length 3,088,286**. Their effective physical bin widths are 6,176,000
and 6,200,000 bp respectively. `--all-bins 1250` requests a grid closer to 2.5 Mb
while retaining the same V9-style chromosome length.

Concatenate actual chromosome lengths with **unpadded** physical offsets
`O_i = sum(L_j, j < i)`. A raw cell at source resolution `S` maps as follows:

```
p_i = min(K - 1, floor((O_i + x*S) / C))
p_j = min(K - 1, floor((O_j + y*S) / C))
a = floor(p_i / H)
b = floor(p_j / H)
AllFinest[min(a,b)*D, max(a,b)*D] += value
```

Raw queries use ordinary All pseudo coordinates and bin size `F`. Returned bin
indices are sparse multiples of `D`; BP coordinates are those indices times
`F`. Viewers query the finest stored level, divide indices by `D`, and display
an `N`-bin overview at virtual bin size `H`. This virtual size is not inserted
into the global resolution list. Merely changing the chromosome length without
remapping records and viewer geometry would produce a mostly blank map.

Every other BP descriptor follows the existing global storage/derivation rules.
Materialized levels contain valid aggregates of these sparse finest coordinates;
derived levels are exact SUM aggregations of their materialized source. Coarser
raw levels may remain useful, but the viewer selects the finest stored level.

The default preferred source is **500,000 bp**, giving smaller positional errors
than 2.5 Mb near overview boundaries. Select the materialized anchor nearest to
that preference by multiplicative distance, preferring finer on ties; a coarser
anchor is eligible only if it does not exceed `H*C`. This avoids accidentally
scanning a 1 bp matrix when 500 kb is absent but 2.5 Mb is available. Preflight
reports the actual source. If no coarse anchor is eligible, use the finest
anchor. Direct and normal large writers use prepared coarse cells; root-only
large assembly uses its existing root-cell pass.

Repair selects from anchors populated by real chromosomes, so a resolution
retained only for a legacy V9 All does not accidentally produce an empty
overview. If a real pair lacks data at the selected anchor but has data at
another materialized level, it uses a populated fallback for that pair. The
source-resolution attribute records the default anchor; Slurm preflight's
`source_resolution_pairs` reports actual per-pair source choices.

Only canonical raw cis and trans cells are added, without mirroring or
normalization. Checked uint64 accumulation preserves the total of every supplied
count cell; score sources use double accumulation and float32 storage. The
original V9 All may use float precision. Repair cannot recover positions within
coarse cells: contacts near an overview-bin boundary can land in a neighboring
bin, although the total raw count remains exact. Chromosome-end padding no
longer shifts every subsequent chromosome.

### Version-two attributes

Integer attributes are positive decimal strings. `chromosomeBins` retains
source-bin counts for validation but does **not** define padded boundaries.

| Attribute | Meaning |
| --- | --- |
| `hicAll.version` | `2` |
| `hicAll.finestResolution` | `F`, first advertised BP bin size |
| `hicAll.sourceResolution` | `S`, actual real-chromosome source |
| `hicAll.groupFactor` | `1`, retained for API compatibility; not a physical scale |
| `hicAll.binCount` | `N`, compact display-axis bins |
| `hicAll.chromosomeBins` | Comma-separated `ceil(L_i/S)` in real chromosome order |
| `hicAll.coordinateScale` | `C`, physical bp per All pseudo-coordinate unit |
| `hicAll.storageStride` | `D`, finest stored bins per compact display bin |

The metadata APIs expose `version`, `displayResolution=H`, `storageStride=D`,
`coordinateScale=C`, `syntheticLength=K` where applicable, and
`effectiveResolution=H*C`. JS `layout` contains unpadded physical offsets and
lengths; `displayLength` includes the partial terminal display bin so the whole
compact grid fits in the viewport. Java/C++/Rust `sourceOffsets` contains physical bp in version 2;
`chromosomeLengths` contains actual physical lengths. In version 1,
`sourceOffsets` retains its earlier source-bin units.

## Version 1: compact synthetic span (backward compatibility)

Earlier files concatenate `P_i = ceil(L_i/S)` source bins, with source-bin offsets
`O_i = sum(P_j, j < i)`. A grouping factor `G` gives `N = ceil(sum(P_i)/G)`.
Their stored All length is `N*F`, so a 1 bp synthetic level with 1249 overview
bins has All length **1249**, irrespective of the physical genome span.
Cells map to `floor((O_i+x)/G)` and `floor((O_j+y)/G)`, with raw finest indices
already equal to compact display indices. Physical display bins have width
`S*G`, and chromosome boundaries use source-bin padding.

Its attributes are `hicAll.version=1`, `finestResolution=F`,
`sourceResolution=S`, `groupFactor=G`, `binCount=N` and `chromosomeBins=P_i`.
The new readers retain this interpretation; existing files are not silently
reinterpreted. `--rebuild-all` migrates them from their real-chromosome data.

## Shared constraints and repair

All is identified by its complete case-insensitive name, not index zero.
Fresh builds append it, preserving real chromosome IDs. All-to-real matrix
pairs are omitted and therefore zero. Real interchromosomal contacts are pixels
inside All/All. All has no normalization vectors or expected-value scale factors
and is excluded from normalization. FRAG files give it an empty site list and
zero FRAG matrices. Existing EVI/NEVI maximum lengths remain valid because the
new pseudo-length cannot exceed the longest real chromosome.

`hic_v10 add-all INPUT OUTPUT` creates All if absent and preserves usable
existing All. `--rebuild-all` explicitly reconstructs/replaces it from the real
matrices and removes active old All/real pairs. It preserves chromosome IDs,
real counts, resolution policies, normalization and expected data. Rebuilding
an All with normalization vectors is rejected rather than silently retaining
vectors of the wrong length. `--in-place` requests atomic replacement.

The repair copies the existing file body in parallel, shifts absolute pointers
by the header-growth offset, patches indexes in bounded batches and appends the
new All and active footer. Obsolete data remain inert in the copied body. A
padding attribute is permitted when rebuilding would otherwise shrink the
header. The output is synced and validated before publication. It requires
space and I/O for a complete replacement and is not resumable. See
[Slurm repair](../hictools-c/v10/SlurmAddAll.md).

Updated metadata APIs are JS `allOverview`, Java `V10Header.allOverview`,
C++ `File::allOverview(AllOverview&)`, and Rust `HicFile::all_overview()`.
Ordinary V10 readers still decode the sparse matrices; convention-aware viewers
provide the compact overview, correct boundaries and chromosome navigation.
