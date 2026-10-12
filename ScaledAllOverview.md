# Scaled All overview convention for V10

This is an optional application convention using the existing V10 chromosome
records, attributes, resolutions, matrix descriptors, blocks, and directories.
It introduces no binary fields, enumeration values, reserved flags, or version
change. The normative wire format remains [HiCFormatV10.md](HiCFormatV10.md).

## Stored chromosome and resolutions

The native V10 producers include an `All` chromosome by default. Its useful
matrix is stored at the **smallest BP bin size already advertised** by the file.
That bin size describes synthetic coordinates, not actual fine-resolution
whole-genome measurements. A file advertising 100 bp uses a synthetic 100 bp
All level; it does not acquire a new global 1 bp resolution.

Let `F` be that smallest bin size, `S` the coarse real-chromosome source bin size,
`L_i` each real chromosome length, and `P_i = ceil(L_i / S)`. Concatenate the
chromosomes in header order, skipping All, with source-bin offsets
`O_i = sum(P_j, j < i)`. Choose a positive integer grouping factor `G` and let
`N = ceil(sum(P_i) / G)`. The All chromosome length is `N * F`. A source cell
`(x, y, value)` in pair `(i, j)` maps to:

```
a = floor((O_i + x) / G)
b = floor((O_j + y) / G)
All[min(a,b), max(a,b)] += value
```

Only canonical raw source cells are added; neither cis mirroring nor
normalization is applied. Raw count accumulation uses checked uint64 addition.
Score sources are accumulated in double precision and stored as float32. The
source-bin grid pads each chromosome to a multiple of `S`; the last All bin can
also contain padding. Boundary labels use this layout instead of unpadded
chromosome-length proportions. An All pixel may span a chromosome boundary.
Rebinning coarse cells preserves count totals, but cannot recover positions
within those cells.

Every All matrix descriptor follows the global resolution mode/source rules.
Materialized BP levels contain aggregates of the base overview, and derived
levels are exact SUM aggregations of their declared materialized source. Most
coarser levels collapse to a few bins or one bin. Readers still support querying
these levels; viewers normally select only the useful finest All level.
The All/All matrix uses the ordinary rotated-cis grid. All-to-real pairs are
omitted: as specified in the existing format, absent pairs are zero matrices at
every resolution. Interchromosomal contacts between real chromosomes are pixels
inside All/All and do not require All/chrN matrices.

New All is appended to the chromosome table so existing chromosome indices stay
unchanged. Consumers identify it by the complete case-insensitive name `all`,
not index zero or a substring match. Existing legacy All data without the
attributes below retains its legacy coordinate interpretation.

## Attributes for convention version 1

All values are ordinary NUL-terminated attribute strings. Integer values are
positive decimal integers. Unknown convention versions do not change core V10
decoding; applications must not assume version-one display semantics for them.

| Attribute | Meaning |
| --- | --- |
| `hicAll.version` | `1` |
| `hicAll.finestResolution` | `F`, equal to the first advertised BP bin size |
| `hicAll.sourceResolution` | `S`, the real-chromosome source bin size |
| `hicAll.groupFactor` | `G`, consecutive source bins per base overview bin |
| `hicAll.binCount` | `N`, base overview axis width |
| `hicAll.chromosomeBins` | Comma-separated `P_i` values for real chromosomes in header order |

The physical scale per overview bin is `S * G`, while raw All queries use
synthetic coordinates with bin size `F`. Applications validate the metadata
against the chromosome lengths and header resolutions. The JS metadata API
exposes `allOverview`; Java exposes `V10Header.allOverview`; native C++ exposes
`File::allOverview(AllOverview&)`; Rust exposes `HicFile::all_overview()`.

## Writer defaults and vectors

The default bin budget is 500; the grouping factor is chosen by integer ceiling
so the result may be smaller than 500. The preferred source is 2,500,000 bp.
Repair and normal large assembly select the coarsest stored materialized BP
anchor no coarser than that preference, falling back to the finest anchor.
Direct preprocessing can use its already prepared coarse matrices; root-only
large assembly accumulates during its existing root-cell pass.

The bin budget also limits the synthetic span to signed 32-bit coordinates
when the finest bin permits it, for desktop boundary-array compatibility.
It is capped so the synthetic All does not increase any
maximum chromosome bin count. Consequently, existing EVI/NEVI vector lengths
remain valid during repair, including for small genomes. All has no NVI vectors
or expected-value scale factors; normalization computations exclude it. For
files advertising FRAG, All has an empty cut-site list, one fragment, and zero
FRAG matrices with otherwise valid descriptors. Its overview is BP-only.

## Repair and compatibility

`hic_v10 add-all INPUT OUTPUT` reconstructs the overview from the existing map;
it does not require HBS or original read pairs. `--in-place` explicitly requests
atomic replacement. Existing usable All is preserved, and repeated in-place
calls leave it unchanged. An existing All header without an All/All matrix is
reported as an error rather than mistaken for a successful repair.

The repacker copies original compressed matrix blocks and vector chunks,
relocates their absolute pointers, appends the new matrix, and rebuilds the
header and footer. It validates the new file before publishing it. Repair needs
space and I/O for a replacement file; aggregation time and full-file copy time
are separate costs. Native fresh builds insert All before final publication and
avoid that extra full-file copy. Conversion of a V9 file lacking All currently
uses the repacker after conversion; existing V9 All is preserved.

Ordinary V10 parsers can decode these files without a format update. Viewers
need convention support to select the useful level and interpret genome-wide
boundaries. Real-chromosome resolutions, counts, chromosome IDs, normalizations,
and expected values remain unchanged by add-All repair.
