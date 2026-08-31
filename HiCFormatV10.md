# `.hic` V10: Lossless, Query-Optimized Format

## Status

**Final consolidated V10 specification**

V10 substantially reduces `.hic` file size for deeply sequenced, high-resolution Hi-C maps while preserving:

- exact raw contact counts and scores;
- interactive local and remote querying;
- all advertised normalization vectors;
- all advertised expected-value vectors;
- exact existing resolution semantics;
- random access to chromosome pairs and genomic windows.

V10 is **not** a cold-storage format and does not introduce lossy approximations.

The key words **MUST**, **MUST NOT**, **SHOULD**, **SHOULD NOT**, **MAY**, and **RECOMMENDED** in this document describe conformance requirements and implementation guidance.

The central design principle is:

> **Do not store information that can be reconstructed exactly from another stored representation, but do precompute and store information that cannot be reconstructed cheaply or exactly at query time.**

The major V10 changes are:

1. Store only selected matrix-resolution anchors; the 20, 50, 200, 500, and 2,000 bp matrices are always derived from raw anchor counts.
2. Store normalization and expected-value vectors for every advertised logical resolution, including derived resolutions.
3. Replace ZLib with Zstandard.
4. Compress every non-empty logical block independently for true block-level random access.
5. Replace the sparse list-of-rows representation with flat delta-coded cell positions.
6. Separate position and value streams.
7. Exploit the count≈1 regime with implicit/default values and exception streams.
8. Add adaptive sparse, bitmap, and dense block representations.
9. Introduce an explicit integer count type rather than falling back from `short` to floating point.
10. Compress normalization and expected-value vectors losslessly in independently accessible chunks.
11. Use an exact binary block index with one locator per independently compressed logical block.

---

## Part I — Normative Binary Specification

This part defines the V10 byte stream. A reader or writer MUST follow this part. Later sections define required semantics, the standard resolution hierarchy, implementation guidance, and rationale. If an example or rationale conflicts with this part, this part controls.

### A. Primitive Encodings and Conventions

#### A.1 Byte order and scalar types

All multibyte scalars are stored in **little-endian** byte order. Records are packed without implicit alignment or padding. Reserved fields are part of the byte stream and MUST be written as zero. Readers MUST ignore a reserved field after validating any requirement explicitly stated for it.

| Type | Bytes | Encoding |
|---|---:|---|
| `u8` | 1 | Unsigned integer, 0 through 255 |
| `u16` | 2 | Unsigned little-endian integer |
| `u32` | 4 | Unsigned little-endian integer |
| `u64` | 8 | Unsigned little-endian integer |
| `f32` | 4 | IEEE 754 binary32; raw bits stored little-endian |
| `f64` | 8 | IEEE 754 binary64; raw bits stored little-endian |
| `cstr` | variable | UTF-8 bytes followed by one zero byte |
| `uleb128` | 1–10 | Canonical unsigned LEB128 representing a `u64` |

The legacy names `byte`, `short`, `int`, `long`, `float`, and `double` correspond to `u8`, 16-bit integer, 32-bit integer, 64-bit integer, `f32`, and `f64` respectively. V10 uses explicit signedness in its normative tables.

A `cstr` does not include a byte length, MUST be valid UTF-8, and MUST NOT contain an embedded zero byte. For example, `HIC` is the four bytes `48 49 43 00` in hexadecimal. Readers SHOULD impose a configurable maximum string length and MUST fail cleanly if the terminator is not present within the containing record.

#### A.2 Canonical unsigned LEB128

For a value `v`, a writer emits seven payload bits per byte, least-significant group first. Bit 7 is one when another byte follows and zero on the final byte. Zero is encoded as the single byte `00`.

A conforming encoding:

- uses the minimum possible number of bytes;
- uses no more than 10 bytes;
- does not set payload bits above bit 63;
- does not overflow `u64` during decoding.

Readers MUST reject non-canonical, unterminated, or overflowing encodings. Signed LEB128 and zig-zag encoding are not used in V10.

#### A.3 Offsets, lengths, and intervals

All file positions are absolute `u64` byte offsets from the first byte of the file. All byte lengths are unsigned and describe exactly the stored bytes of the referenced record. The interval described by `(position, length)` is half-open:

```text
[position, position + length)
```

The addition MUST be checked for overflow and the resulting interval MUST lie within the file. A zero position and zero length mean that an optional structure is absent. A locator pair with only one zero value is invalid. No valid present structure has a zero length.

Genomic coordinates and bin intervals are also zero-based and half-open. For base-pair resolution `R`:

```text
binIndex = floor(position / R)
binStart = binIndex * R
binEnd   = min((binIndex + 1) * R, chromosomeLength)
nBins    = ceil(chromosomeLength / R)
```

#### A.4 Floating-point requirements

Contact scores use `f32`. Normalization values, expected values, chromosome scale factors, standard deviations, and percentiles also use `f32`. A score-resolution sum uses `f64`, as defined in the matrix descriptor.

Writers MUST preserve the exact `f32` bit patterns of score, normalization, and expected-value data. All IEEE 754 bit patterns, including signed zero, infinities, and NaNs, are representable. A reader MUST NOT treat a NaN score as an absent cell; presence is encoded separately. Vector consumers MAY treat NaN normalization values as unusable bins, but the storage layer MUST reproduce their bits exactly.

#### A.5 Enumerations

Unknown enumeration values are errors unless a later revision explicitly defines an extension mechanism for the containing record.

| Name | Value | Meaning |
|---|---:|---|
| `Unit.BP` | 0 | Base-pair bins |
| `Unit.FRAG` | 1 | Restriction-fragment bins |
| `StorageMode.MATERIALIZED` | 0 | Independently compressed logical blocks are stored |
| `StorageMode.DERIVED` | 1 | Raw matrix is derived from a finer materialized resolution |
| `Aggregation.NONE` | 0 | Resolution cannot be derived |
| `Aggregation.SUM` | 1 | Target cells are sums of source cells |
| `ValueType.COUNT_UINT` | 0 | Non-negative integer contacts, up to `u64` |
| `ValueType.SCORE_FLOAT32` | 1 | Exact IEEE 754 binary32 scores |
| `GridType.RECTANGULAR` | 0 | Ordinary row-major block grid |
| `GridType.ROTATED_CIS` | 1 | V9-compatible rotated cis grid |
| `BlockCodec.ZSTD` | 1 | Zstandard frame without a preset dictionary |
| `BlockRepresentation.SPARSE_DELTA` | 0 | Sorted occupied positions encoded as deltas |
| `BlockRepresentation.BITMAP` | 1 | One presence bit per cell |
| `BlockRepresentation.DENSE` | 2 | Every cell has a value slot |
| `ValueMode.ALL_DEFAULT` | 0 | One value applies to every encoded value slot |
| `ValueMode.DEFAULT_EXCEPTIONS` | 1 | One default plus ordinal exceptions |
| `ValueMode.DIRECT` | 2 | One value per encoded value slot |
| `VectorTransform.RAW` | 0 | Unmodified little-endian `f32` words |
| `VectorTransform.BYTE_SHUFFLE` | 1 | Four byte lanes, defined below |
| `VectorTransform.XOR32` | 2 | Adjacent `f32` words XORed bitwise |

#### A.6 Global invariants

A conforming file MUST satisfy all of the following:

- the first four bytes are `HIC\0` and the version is 10;
- every referenced interval is in bounds;
- chromosome indices and resolution indices refer to header entries;
- chromosome-pair keys are canonical with `chr1Index <= chr2Index`;
- matrix and vector directory keys are unique;
- every logical matrix cell is encoded at most once;
- logical block numbers are strictly increasing within each resolution block index;
- every indexed matrix block has one independently compressed stored-block record;
- sparse positions and exception ordinals are strictly increasing;
- counts never overflow `u64` while decoding or deriving;
- the exact full-matrix count sum fits `u64`;
- every chromosome's bin count at every advertised resolution fits `u32`;
- every computed logical block number fits `u32`;
- all advertised derived resolutions point directly to a finer materialized resolution in the same unit;
- advertised BP resolutions 20, 50, 200, 500, and 2,000 use the mandatory sources defined in Section C.3;
- an advertised 500,000 bp resolution is materialized;
- cis matrix descriptors use `ROTATED_CIS` and trans descriptors use `RECTANGULAR`;
- all unused flag bits and reserved bytes are zero.

### B. Top-Level File Organization

The physical order produced by the reference writer is:

```text
Header
Matrix metadata records
Independently compressed matrix blocks
Materialized-resolution block indexes
Normalization-vector chunks
Expected-value chunks
Normalized-expected-value chunks
Normalization-vector index
Expected-value index
Normalized-expected-value index
Footer / matrix master index
```

The offsets are authoritative, so a writer MAY choose another order. Writers SHOULD store the blocks of one materialized matrix resolution contiguously in increasing block-number order so readers can coalesce adjacent local or HTTP range reads, but readers MUST use the exact index locators and MUST NOT depend on physical contiguity.

Seekable writers normally reserve and write the header, emit the remaining structures, then backpatch header offsets and lengths. A streaming producer MUST stage enough metadata to emit a correct final seekable file; V10 does not define an unpatched streaming variant.

Fixed-size record summary:

| Record | Bytes |
|---|---:|
| Fixed header prefix | 88 |
| Footer header | 24 |
| Footer matrix entry | 24 |
| Matrix record header | 24 |
| Matrix resolution descriptor | 76 |
| Block-index header | 24 |
| Block-index entry | 16 |
| Stored block-record header | 16 |
| Logical block header | 40 |
| Vector-index header | 16 |
| Vector chunk descriptor | 32 |
| Stored vector chunk header | 16 |

Variable records have explicit lengths or formulas in their defining sections. This summary is informative; the field tables are authoritative.

### C. Header

#### C.1 Fixed header prefix

The file starts with this fixed 88-byte prefix:

| Offset | Field | Type | Required value or meaning |
|---:|---|---|---|
| 0 | `magic` | 4 bytes | `48 49 43 00` (`HIC\0`) |
| 4 | `version` | `u32` | 10 |
| 8 | `headerByteLength` | `u64` | Entire header, fixed prefix through final variable field |
| 16 | `footerPosition` | `u64` | Footer/master-index position |
| 24 | `footerByteLength` | `u64` | Footer/master-index stored length |
| 32 | `normVectorIndexPosition` | `u64` | Normalization-vector index, or zero if absent |
| 40 | `normVectorIndexLength` | `u64` | Normalization-vector index length, or zero |
| 48 | `expectedValueIndexPosition` | `u64` | Raw expected-value index, or zero if absent |
| 56 | `expectedValueIndexLength` | `u64` | Raw expected-value index length, or zero |
| 64 | `normExpectedValueIndexPosition` | `u64` | Normalized expected-value index, or zero if absent |
| 72 | `normExpectedValueIndexLength` | `u64` | Normalized expected-value index length, or zero |
| 80 | `fileFlags` | `u32` | Zero in this specification |
| 84 | `reserved` | `u32` | Zero |

The three direct vector-index locators preserve inexpensive remote access without first downloading the footer. `headerByteLength` MUST be at least 88 and MUST equal the reader's position immediately after parsing the variable header.

The footer locator is mandatory and both of its fields MUST be nonzero. Each optional vector-index locator MUST be either a valid nonzero pair or the absent pair `(0, 0)`.

#### C.2 Variable header

The following fields immediately follow the fixed prefix, with no alignment padding:

| Field | Type | Cardinality and meaning |
|---|---|---|
| `genomeId` | `cstr` | Genome or assembly identifier; may be empty |
| `nAttributes` | `u32` | Number of attribute pairs |
| `attributeKey`, `attributeValue` | `cstr`, `cstr` | Repeated `nAttributes` times in stored order |
| `nChromosomes` | `u32` | Number of chromosome records |
| `chromosomeName`, `chromosomeLength` | `cstr`, `u64` | Repeated `nChromosomes` times |
| `nBpResolutions` | `u32` | Number of BP resolution records |
| `bpResolution` | 12-byte record | Repeated `nBpResolutions` times |
| `nFragResolutions` | `u32` | Number of FRAG resolution records |
| `fragResolution` | 12-byte record | Repeated `nFragResolutions` times |
| fragment-site lists | variable | Present only when `nFragResolutions > 0` |
| `nNormalizationTypes` | `u32` | Number of normalization type strings |
| `normalizationTypeName` | `cstr` | Repeated `nNormalizationTypes` times |

`nChromosomes` MUST be positive. Chromosome index is the zero-based array index. Names MUST be non-empty and unique, and chromosome lengths MUST be positive. `chromosomeLength` is stored as `u64`; however, the bin count at every advertised BP or FRAG resolution MUST fit `u32`, because block coordinates are `u32`. Implementations interoperating with signed-language coordinate APIs SHOULD also reject chromosome lengths above `2^63 - 1`.

The normalization type ID is the zero-based index into `normalizationTypeName`. Normalization type names MUST be non-empty and unique. Names such as `VC`, `VC_SQRT`, `KR`, `INTER_KR`, `INTER_VC`, `GW_KR`, and `GW_VC` have their established meanings, but the dictionary is extensible. `NONE` denotes the absence of normalization and MUST NOT be entered as a normalization type.

Attributes are an ordered application-metadata list and do not alter core decoding unless this specification assigns a meaning to a key. Duplicate keys are permitted and MUST NOT be collapsed by a repacker. Common keys include `statistics`, `graphs`, and `software`. A reader MUST preserve unknown attributes, order, and duplicates when repacking a file.

#### C.3 Logical-resolution record

Each BP and FRAG resolution list contains fixed 12-byte records:

| Field | Type | Meaning |
|---|---|---|
| `binSize` | `u32` | Positive bin size in the list's unit |
| `storageMode` | `u8` | `MATERIALIZED` or `DERIVED` |
| `aggregation` | `u8` | `NONE` or `SUM` |
| `reserved` | `u16` | Zero |
| `sourceResolutionIndex` | `u32` | Index in this same unit list, or `0xFFFFFFFF` |

Resolution lists MUST be strictly increasing by `binSize`. A materialized record MUST use `sourceResolutionIndex = 0xFFFFFFFF`. A derived record MUST use `aggregation = SUM`, and its source MUST be materialized, finer, in the same list, and divide the target bin size exactly. Chained derivation is forbidden.

The following BP storage policy is mandatory whenever the target resolution is advertised:

| Target BP resolution | Required storage | Required source |
|---:|---|---:|
| 20 | `DERIVED` | 10 |
| 50 | `DERIVED` | 10 |
| 200 | `DERIVED` | 100 |
| 500 | `DERIVED` | 100 |
| 2,000 | `DERIVED` | 1,000 |
| 500,000 | `MATERIALIZED` | — |

Advertising one of the five mandatory derived targets therefore also requires advertising its stated materialized source. A writer MUST NOT materialize those targets or select another source. A writer MUST NOT derive 500,000 bp. These rules apply to ordinary production and V9 conversion; a conversion that cannot reproduce a mandatory target exactly from its required source MUST fail rather than emit a nonconforming V10 file.

The `(unit, resolutionIndex)` pair is the canonical resolution identifier. `binSize` is duplicated in later records for validation; all copies MUST match the header.

#### C.4 Fragment-site lists and FRAG coordinates

When at least one FRAG resolution is advertised, one site list follows for every chromosome in chromosome order:

| Field | Type | Meaning |
|---|---|---|
| `nSites` | `u32` | Number of restriction cut sites |
| `sitePosition` | `u64` | Repeated `nSites` times |

Site positions MUST be strictly increasing, greater than zero, and less than the chromosome length. They are zero-based cut coordinates. Chromosome coordinate `p` belongs to fragment:

```text
fragmentIndex = number of sitePosition values <= p
```

Thus a chromosome with `nSites` has `nSites + 1` fragments. At FRAG bin size `R`, `binIndex = floor(fragmentIndex / R)` and `nBins = ceil((nSites + 1) / R)`.

### D. Footer and Matrix Master Index

`footerPosition` addresses the footer, which has this exact layout:

| Field | Type | Meaning |
|---|---|---|
| `footerMagic` | 4 bytes | ASCII `H10F` |
| `footerVersion` | `u32` | 1 |
| `footerByteLength` | `u64` | Must match the header |
| `matrixEntryCount` | `u32` | Number of matrix-directory entries |
| `reserved` | `u32` | Zero |
| matrix entries | 24 bytes each | Repeated `matrixEntryCount` times |

Each matrix entry is:

| Field | Type | Meaning |
|---|---|---|
| `chr1Index` | `u32` | First chromosome index |
| `chr2Index` | `u32` | Second chromosome index; `chr1Index <= chr2Index` |
| `matrixPosition` | `u64` | Position of the matrix metadata record |
| `matrixByteLength` | `u64` | Exact matrix metadata record length |

Entries MUST be sorted by `(chr1Index, chr2Index)` and unique. Unlike V9, V10 does not store textual keys such as `0_1`; the numeric pair is the key. A chromosome pair with no directory entry is an all-zero raw matrix at every logical resolution. The footer contains metadata locators only; stored block bytes are not included in `matrixByteLength`.

The footer length is exactly `24 + 24 * matrixEntryCount` bytes.

### E. Matrix Metadata Record

#### E.1 Matrix record header

The record selected by a footer entry begins:

| Field | Type | Meaning |
|---|---|---|
| `matrixMagic` | 4 bytes | ASCII `H10M` |
| `matrixRecordVersion` | `u32` | 1 |
| `chr1Index` | `u32` | Must match the directory key |
| `chr2Index` | `u32` | Must match the directory key |
| `resolutionCount` | `u32` | Must equal `nBpResolutions + nFragResolutions` |
| `reserved` | `u32` | Zero |
| resolution descriptors | 76 bytes each | BP descriptors followed by FRAG descriptors |

The matrix record length is exactly `24 + 76 * resolutionCount` bytes and MUST equal the footer entry's `matrixByteLength`.

#### E.2 Resolution descriptor

Every logical resolution, including a derived resolution, has one fixed 76-byte descriptor:

| Field | Type | Meaning |
|---|---|---|
| `unit` | `u8` | `BP` or `FRAG` |
| `storageMode` | `u8` | Must match the header resolution |
| `aggregation` | `u8` | Must match the header resolution |
| `valueType` | `u8` | `COUNT_UINT` or `SCORE_FLOAT32` |
| `resolutionIndex` | `u32` | Index in the corresponding header list |
| `binSize` | `u32` | Must match the header record |
| `sourceResolutionIndex` | `u32` | Must match the header record |
| `gridType` | `u8` | `RECTANGULAR` or `ROTATED_CIS` |
| `reserved0` | 3 bytes | Zero |
| `sumCountsOrScores` | 8 bytes | `u64` for counts; `f64` for scores |
| `occupiedCellCount` | `u64` | Number of non-absent cells in the full logical matrix |
| `stdDev` | `f32` | Population standard deviation among occupied values; NaN if not computed |
| `percent95` | `f32` | Estimated 95th percentile among occupied values; NaN if not computed |
| `blockBinCount` | `u32` | Positive logical block scale in bins |
| `blockColumnCount` | `u32` | Positive number of columns along the grid's primary axis |
| `blockIndexPosition` | `u64` | Exact block-index position, or zero when absent |
| `blockIndexLength` | `u64` | Exact block-index length, or zero when absent |
| `logicalBlockCount` | `u32` | Total physically stored non-empty blocks; zero for a derived or empty resolution |
| `reserved1` | `u32` | Zero |

For `COUNT_UINT`, `sumCountsOrScores` is the exact unsigned sum over the stored canonical matrix cells. For `SCORE_FLOAT32`, scores are visited by increasing global bin row and then bin column, converted exactly to binary64, and added using round-to-nearest, ties-to-even; the resulting binary64 bits are stored, including an IEEE result such as NaN or infinity. `occupiedCellCount` is an integer, not a float. Zero-valued count cells are absent; a present score may have any `f32` bit pattern, including `+0`, `-0`, or NaN.

`stdDev` and `percent95` are descriptive metadata and do not affect decoding. Writers that do not compute either field MUST store the canonical quiet-NaN bits `0x7fc00000`, not zero.

A non-empty materialized descriptor MUST address a block index and have `logicalBlockCount > 0`. An all-empty materialized resolution and every derived resolution MUST use `(blockIndexPosition, blockIndexLength, logicalBlockCount, reserved1) = (0, 0, 0, 0)`. The block-index entry count MUST equal `logicalBlockCount`. A derived descriptor's `valueType`, grid metadata, statistics, and auxiliary vectors describe the target resolution, not the source.

Every descriptor in a cis matrix record, including a derived or empty descriptor, MUST use `ROTATED_CIS`. Every descriptor in a trans matrix record MUST use `RECTANGULAR`. This invariant lets readers reject geometry ambiguity before fetching blocks.

Within one matrix record, a derived descriptor and its declared source descriptor MUST have the same `valueType`. A missing chromosome-pair footer entry denotes an all-zero `COUNT_UINT` matrix. Therefore a file using `SCORE_FLOAT32` for any omitted chromosome pair MUST instead emit an explicit matrix record for that pair, even when it has no occupied cells.

### F. Logical Block Geometry

#### F.1 Canonical chromosome-pair orientation

The matrix key always uses `chr1Index <= chr2Index`. In canonical orientation, `binColumn` belongs to chromosome 1 and `binRow` belongs to chromosome 2. A query with chromosomes in the opposite order is transposed into canonical order before block lookup and transposed back in the result.

For cis matrices, writers store only cells satisfying `binRow >= binColumn`. The upper triangle is reconstructed by transposition. A diagonal cell is emitted once. Counts and statistics refer to this canonical stored triangle, not to a duplicated symmetric square.

#### F.2 Rectangular trans grid

`GridType.RECTANGULAR` is required when `chr1Index != chr2Index` and forbidden for cis matrices. With `B = blockBinCount`:

```text
blockColumn = floor(binColumn / B)
blockRow    = floor(binRow / B)
blockNumber = blockRow * blockColumnCount + blockColumn
```

`blockColumnCount` MUST equal `ceil(nBins(chr1) / B)`. The number of possible block rows is `ceil(nBins(chr2) / B)`.

#### F.3 Rotated cis grid

`GridType.ROTATED_CIS` is required when `chr1Index == chr2Index` and forbidden for trans matrices. It preserves the V9 diagonal/anti-diagonal organization. First canonicalize `binRow >= binColumn`. Let `B = blockBinCount` and `d = binRow - binColumn`. Then:

```text
alongDiagonal = floor((binColumn + binRow) / (2 * B))

alongAntiDiagonal = largest integer a >= 0 such that
    d^2 >= 2 * B^2 * (2^a - 1)^2

blockNumber       = alongAntiDiagonal * blockColumnCount + alongDiagonal
```

The integer inequality is the normative definition and is exactly equivalent to `floor(log2(1 + d / (sqrt(2) * B)))`. It avoids platform-dependent behavior near a floating-point boundary. Implementations MUST use checked sufficiently wide integer arithmetic (or an algebraically equivalent exact comparison) and MUST verify that the final block number fits `u32`.

For a rotated grid, `blockColumnCount` MUST equal `ceil(nBins(chr1) / B)` just as it does along the primary axis of a rectangular grid. Multiplication by this value separates distance bands without collisions.

All cells assigned the same number form one logical block, even though the region is rotated in genomic coordinates. A stored block therefore carries explicit genomic bin offsets, width, and height; readers MUST use those fields when reconstructing cells.

For a cis query rectangle, a reader MAY reproduce the V9 permissive block-candidate calculation, but it MUST filter decoded cells against the requested bin intervals. An implementation may instead enumerate the query's canonical bins or maintain an equivalent inverse index. Query optimization MUST NOT change cell membership.

#### F.4 Standard adaptive block scale

V10 restores the established V9 adaptive block-sizing policy so ultra-fine matrices do not fragment into hundreds of millions of nearly empty logical blocks. For an ordinary BP matrix, let `N = max(nBins(chr1), nBins(chr2))`, `R = binSize`, `K = 1000`, and let `C = 500` for cis or `C = 5000` for trans. Compute:

```text
if R < C:
    targetColumns = floor(N * R / (K * C)) + 1
else:
    targetColumns = floor(N / K) + 1

adaptiveB = floor(N / targetColumns) + 1
blockBinCount = max(requestedMinimumBlockBinCount, adaptiveB)
```

All products MUST use checked sufficiently wide arithmetic. A writer MUST further increase `blockBinCount` if needed to keep every possible block number in `u32`. The standard requested minimum is 256 bins. A writer option may increase that minimum, but MUST NOT reduce the adaptive minimum. This policy deliberately makes blocks much larger as BP resolution becomes finer while retaining bounded logical grids and V9-compatible cis distance bands.

### G. Resolution Block Index

Every non-empty materialized resolution has one uncompressed exact block index. Its length is `blockIndexLength`:

| Field | Type | Meaning |
|---|---|---|
| `indexMagic` | 4 bytes | ASCII `H10I` |
| `indexVersion` | `u32` | 2 |
| `indexByteLength` | `u64` | Exact length of this index record |
| `blockCount` | `u32` | Must equal the matrix descriptor's `logicalBlockCount` |
| `reserved` | `u32` | Zero |
| block entries | 16 bytes each | Repeated `blockCount` times |

The block-index length is exactly `24 + 16 * blockCount` bytes and MUST equal both `indexByteLength` and `blockIndexLength`.

Each block entry is:

| Field | Type | Meaning |
|---|---|---|
| `blockNumber` | `u32` | Exact logical block number |
| `storedByteLength` | `u32` | Stored-block header plus one Zstandard frame; greater than 16 |
| `blockPosition` | `u64` | Absolute file position of the stored-block record |

Entries MUST be sorted by strictly increasing `blockNumber`. Every indexed interval MUST be in bounds, MUST contain exactly one stored-block record defined in Section H, and MUST NOT overlap any other indexed block interval. An absent block number denotes an empty logical block.

The exact locators permit direct binary search and at most one direct range read per requested block. Readers MAY sort and coalesce contiguous requested intervals into a single local or HTTP range read, but every block remains independently decompressible and no query may require bytes from an unrequested block.

`indexVersion = 2` intentionally distinguishes this exact block index from the grouped-block `H10I` version 1 used by pre-final V10 drafts. A conforming reader MUST reject version 1 rather than interpreting it as this layout.

### H. Stored Matrix Block

Each non-empty logical block is independently decompressible and starts with this 16-byte uncompressed stored-record header:

| Field | Type | Meaning |
|---|---|---|
| `recordMagic` | 4 bytes | ASCII `H10B` |
| `codec` | `u8` | `ZSTD` (1) |
| `recordVersion` | `u8` | 1 |
| `recordFlags` | `u16` | Zero |
| `uncompressedBytes` | `u32` | Exact decompressed logical-block length |
| `blockNumber` | `u32` | Must match the selected block-index entry |
| `compressedPayload` | bytes | One complete Zstandard frame containing one logical block |

The Zstandard frame consumes exactly `storedByteLength - 16` bytes and MUST decompress to exactly `uncompressedBytes`. The decompressed bytes are exactly one logical block encoded by Section I; `uncompressedBytes` MUST equal `40 + positionStreamBytes + valueStreamBytes`. No block directory or additional block payload may appear in the frame.

The payload MUST be one ordinary Zstandard data frame as defined by [RFC 8878](https://www.rfc-editor.org/rfc/rfc8878.html), not a skippable frame or concatenation of frames. Preset dictionaries are forbidden. The frame checksum flag MAY be used and MUST be verified by readers when present. Compression level is a writer choice and is not stored because it does not affect decoding.

### I. Logical Block Encoding

#### I.1 Fixed block header

Every logical block begins with this fixed 40-byte header:

| Field | Type | Meaning |
|---|---|---|
| `blockVersion` | `u8` | 1 |
| `representation` | `u8` | `SPARSE_DELTA`, `BITMAP`, or `DENSE` |
| `valueMode` | `u8` | `ALL_DEFAULT`, `DEFAULT_EXCEPTIONS`, or `DIRECT` |
| `valueType` | `u8` | Must match the matrix descriptor |
| `blockFlags` | `u8` | Bit 0 is `EXPLICIT_PRESENCE`; all other bits zero |
| `reserved` | 3 bytes | Zero |
| `binColumnOffset` | `u32` | Global bin column represented by local column zero |
| `binRowOffset` | `u32` | Global bin row represented by local row zero |
| `blockWidth` | `u32` | Positive number of local columns, `W` |
| `blockHeight` | `u32` | Positive number of local rows, `H` |
| `occupiedCellCount` | `u64` | Number of present cells, `N` |
| `positionStreamBytes` | `u32` | Exact position-stream length |
| `valueStreamBytes` | `u32` | Exact value-stream length |

The two streams immediately follow the header, position stream first. Their lengths plus 40 MUST equal the enclosing stored record's `uncompressedBytes`. Local position `p` maps to:

```text
localRow    = floor(p / W)
localColumn = p % W
binRow      = binRowOffset + localRow
binColumn   = binColumnOffset + localColumn
```

`W * H` MUST be computed with checked `u64` arithmetic. Every reconstructed bin must lie within the chromosomes' bin counts, belong to the block number under the declared grid, and obey the canonical cis triangle.

#### I.2 Sparse-delta position stream

For `SPARSE_DELTA`, the position stream contains exactly `N` `uleb128` values. The first is absolute `p0`; every later value is `p[i] - p[i-1]`. Positions MUST be strictly increasing and less than `W * H`. `EXPLICIT_PRESENCE` MUST be zero because the position list itself defines presence.

#### I.3 Bitmap position stream

For `BITMAP`, the position stream is exactly `ceil(W * H / 8)` bytes. Bit `p` is bit `(p % 8)`—least significant bit first—of byte `floor(p / 8)`. One means present. Padding bits above `W * H` MUST be zero. The bitmap population count MUST equal `N`. `EXPLICIT_PRESENCE` MUST be one.

#### I.4 Dense position stream

For a dense `COUNT_UINT` block, the position stream is empty, `EXPLICIT_PRESENCE` is zero, and every one of the `W * H` cells has a value slot. A zero count means absent; `N` MUST equal the number of nonzero decoded values.

For a dense `SCORE_FLOAT32` block, `EXPLICIT_PRESENCE` MUST be one and the position stream is the same bitmap defined above. Every cell still has a value slot. A writer MUST store positive-zero bits (`0x00000000`) in every absent slot; a reader MUST validate this but return only slots whose presence bit is one. Present values may use any `f32` bit pattern, including zero or NaN, and therefore remain distinguishable from absence. `N` is the bitmap population count.

Dense blocks MUST use `ValueMode.DIRECT`.

#### I.5 Value scalar encoding

A count scalar is canonical `uleb128`. Sparse and bitmap count values MUST be positive; dense count values MAY be zero. A score scalar is exactly four bytes containing its original little-endian `f32` bits.

For sparse and bitmap blocks, the number of value slots `S` is `N`. For dense blocks, `S = W * H`.

#### I.6 `ALL_DEFAULT`

The value stream contains exactly one scalar, `defaultValue`, and no other bytes. The decoded value of every one of the `S` slots is `defaultValue`. This mode is forbidden for `DENSE` and requires `S > 0`.

#### I.7 `DEFAULT_EXCEPTIONS`

The value stream is:

```text
defaultValue                 scalar
exceptionCount               uleb128
exceptionOrdinalDeltas       exceptionCount × uleb128
exceptionValues              exceptionCount × scalar
```

`exceptionCount` MUST be greater than zero and less than `S`. The first exception ordinal is absolute. Each later ordinal is the previous ordinal plus a strictly positive delta. Every ordinal is less than `S`. An exception value MUST differ from the default (`u64` comparison for counts, bitwise comparison for `f32`). This mode is forbidden for `DENSE`.

#### I.8 `DIRECT`

The value stream contains exactly `S` scalars in position or dense row-major order and no other bytes.

#### I.9 Representation selection

The writer MAY choose any valid representation and value mode per block. Choice affects size and speed only, never logical values. Empty blocks (`N = 0`) MUST NOT be stored. A count block in which all occupied values equal one is represented as `ALL_DEFAULT` with default count 1; no legacy `allCountsOne` flag exists because this mode carries the same semantics explicitly.

V9's `short` count limit, float fallback for large integer counts, `-32768` dense sentinel, list-of-rows hierarchy, and per-row position-width flags are not used. Their information is represented without loss by `COUNT_UINT`, explicit presence, flat positions, and stream lengths.

### J. Vector Indexes and Chunks

#### J.1 Shared chunk descriptor

Normalization, raw expected, and normalized expected arrays are `f32` arrays divided into independently addressable chunks. Every index entry contains fixed 32-byte chunk descriptors:

| Field | Type | Meaning |
|---|---|---|
| `firstValueIndex` | `u64` | First array index in the chunk |
| `valueCount` | `u32` | Positive number of `f32` values |
| `transform` | `u8` | `RAW`, `BYTE_SHUFFLE`, or `XOR32` |
| `codec` | `u8` | `ZSTD` (1) |
| `flags` | `u16` | Zero |
| `filePosition` | `u64` | Position of the stored vector chunk |
| `storedByteLength` | `u32` | 16-byte header plus Zstandard frame |
| `uncompressedByteLength` | `u32` | Must equal `4 * valueCount` |

Descriptors MUST cover `[0, valueCountOfVector)` contiguously with no overlaps or gaps. Each descriptor's `valueCount` MUST be no greater than `floor(0xFFFFFFFF / 4)`, and `storedByteLength` MUST be greater than 16. For every vector entry, a positive vector length requires positive `nominalChunkValueCount` and `chunkCount`; a zero vector length requires `chunkCount = 0`. A nominal chunk size of 65,536 values is recommended; the last chunk may be shorter.

#### J.2 Stored vector chunk

At `filePosition` is:

| Field | Type | Meaning |
|---|---|---|
| `chunkMagic` | 4 bytes | ASCII `H10V` |
| `codec` | `u8` | Must match the descriptor |
| `transform` | `u8` | Must match the descriptor |
| `flags` | `u16` | Zero |
| `uncompressedByteLength` | `u32` | Must match the descriptor |
| `valueCount` | `u32` | Must match the descriptor |
| `compressedPayload` | bytes | One complete Zstandard frame |

The payload MUST be one ordinary Zstandard data frame, not a skippable frame or a concatenation of frames. It MUST decompress to exactly `4 * valueCount` transformed bytes. Preset dictionaries are forbidden, and a frame checksum MUST be verified when present.

Transform definitions operate on the original `u32` words `b[i]`, where each word is the exact bit pattern of an `f32` value:

- `RAW`: emit each `b[i]` as four little-endian bytes.
- `BYTE_SHUFFLE`: emit byte 0 of every word, then byte 1 of every word, then byte 2, then byte 3. Within a lane, values remain in array order.
- `XOR32`: emit `b[0]`, then emit `b[i] XOR b[i-1]` for each `i > 0`, each as a little-endian `u32`.

Inverse transformation MUST reproduce every original `f32` bit exactly.

#### J.3 Index container header

All three indexes begin with:

| Field | Type | Meaning |
|---|---|---|
| `indexMagic` | 4 bytes | `NVI0`, `EVI0`, or `NEVI` |
| `indexVersion` | `u32` | 1 |
| `entryCount` | `u32` | Number of entries |
| `reserved` | `u32` | Zero |
| entries | variable | Exactly `entryCount` length-delimited entries |

Every entry starts with `entryByteLength u32`, which includes that length field and permits a reader to skip an entry after validating it remains within the index interval. No padding occurs between entries. Keys MUST be unique and entries MUST be lexicographically sorted by the keys specified below.

For each index, `16 + sum(entryByteLength)` MUST equal the index length recorded in the header. Trailing bytes are invalid.

Index presence is the capability advertisement: an `NVI0` entry advertises that normalization for one chromosome and resolution; an expected capability is advertised only by the matching `EVI0` or `NEVI` entry. A reader MUST NOT infer missing capabilities from another resolution or normalization type.

#### J.4 Normalization-vector index (`NVI0`)

The key is `(normalizationTypeId, chrIndex, unit, resolutionIndex)`. Each entry is:

| Field | Type | Meaning |
|---|---|---|
| `entryByteLength` | `u32` | Entire entry |
| `normalizationTypeId` | `u32` | Header dictionary index |
| `chrIndex` | `u32` | Chromosome index |
| `unit` | `u8` | `BP` or `FRAG` |
| `reserved0` | 3 bytes | Zero |
| `resolutionIndex` | `u32` | Header resolution index |
| `binSize` | `u32` | Validation copy |
| `vectorValueCount` | `u64` | Number of `f32` values |
| `nominalChunkValueCount` | `u32` | Writer target; informational |
| `chunkCount` | `u32` | Number of descriptors |
| chunks | 32 bytes each | Shared chunk descriptors |

`vectorValueCount` MUST equal the chromosome's `nBins` at that resolution. If it is positive, `nominalChunkValueCount` and `chunkCount` MUST be positive; if it is zero, `chunkCount` MUST be zero. Each advertised normalization capability, including one at a derived resolution, requires its own target-resolution entry.

The entry length is exactly `40 + 32 * chunkCount` bytes.

#### J.5 Raw expected-value index (`EVI0`)

The key is `(unit, resolutionIndex)`. Expected array index `d` is genomic bin distance `abs(binRow - binColumn)`. Each entry is:

| Field | Type | Meaning |
|---|---|---|
| `entryByteLength` | `u32` | Entire entry |
| `unit` | `u8` | `BP` or `FRAG` |
| `reserved0` | 3 bytes | Zero |
| `resolutionIndex` | `u32` | Header resolution index |
| `binSize` | `u32` | Validation copy |
| `vectorValueCount` | `u64` | Number of expected `f32` values |
| `nominalChunkValueCount` | `u32` | Writer target; informational |
| `chunkCount` | `u32` | Number of chunk descriptors |
| `scaleFactorCount` | `u32` | Number of chromosome scale factors |
| `reserved1` | `u32` | Zero |
| scale factors | 8 bytes each | `(chrIndex u32, scaleFactor f32)` |
| chunks | 32 bytes each | Shared chunk descriptors |

`vectorValueCount` MUST equal the maximum chromosome `nBins` for this unit and resolution, covering every possible cis distance from zero through `nBins - 1`. Scale factors MUST be sorted by chromosome index and unique. A missing factor means 1.0. Raw expected vectors apply only to cis matrices.

The entry length is exactly `40 + 8 * scaleFactorCount + 32 * chunkCount` bytes.

#### J.6 Normalized expected-value index (`NEVI`)

The key is `(normalizationTypeId, unit, resolutionIndex)`. Its entry is identical to an `EVI0` entry except that `normalizationTypeId u32` appears immediately after `entryByteLength`. The expected array and scale factors are those of the named target-resolution normalization, and `vectorValueCount` has the same full-distance requirement as `EVI0`.

The entry length is exactly `44 + 8 * scaleFactorCount + 32 * chunkCount` bytes.

An advertised O/E capability at a derived resolution requires an expected entry at that target resolution; readers MUST NOT infer it from the source resolution.

#### J.7 Normalization and observed/expected arithmetic

For a present raw contact `C(i,j)` and normalization-vector values `N1[i]` and `N2[j]`:

```text
normalized(i,j) = C(i,j) / (N1[i] * N2[j])
```

If either normalization value is zero, NaN, infinite, absent, or outside its vector, the normalized result is unavailable. Storage is `f32`; applications MAY perform arithmetic in `f64`.

For cis bin distance `d` and chromosome scale factor `s` (default 1.0):

```text
effectiveExpected(d) = expectedVector[d] / s
observedExpected(i,j) = observedOrNormalized(i,j) / effectiveExpected(d)
```

The raw expected index is used with `NONE`; the normalized expected index matching the requested normalization type is used otherwise. An out-of-range, zero, NaN, or infinite expected value makes O/E unavailable for that cell.

### K. Derived Resolution Semantics

A derived raw matrix is computed only from its declared materialized source. Let target bin size be `Rt`, source bin size be `Rs`, and `factor = Rt / Rs`.

For each source contact `(x, y, value)`:

```text
targetX = floor(x / factor)
targetY = floor(y / factor)
target[targetX, targetY] += value
```

Equivalently, target bin `t` covers source-bin interval `[t * factor, min((t + 1) * factor, sourceNBins))`. Readers MUST use this clipped half-open interval at chromosome ends.

Count accumulation uses checked `u64` arithmetic and overflow is a format or query error. For `SCORE_FLOAT32`, derivation is permitted only when `aggregation = SUM` and every contributing source score is finite. Values are visited in increasing source row and then source column, accumulated in IEEE 754 binary64 using round-to-nearest, ties-to-even, and rounded once with the same mode to binary32 for the target cell. A writer MUST materialize a score target if a source value is non-finite, the binary64 accumulation overflows, or this rule does not reproduce a target matrix being preserved.

Aggregation occurs before normalization or expected-value division. Target-resolution vectors are then applied as specified above. Source vectors MUST NOT be aggregated or reused.

### L. Required Reader Procedure

A minimal random-access reader performs these steps:

1. Read the 88-byte header prefix; validate magic, version, lengths, flags, and in-bounds top-level intervals.
2. Read through `headerByteLength`; build chromosome, resolution, fragment-site, attribute, and normalization-type tables.
3. Read the footer and binary-search its numeric matrix entries for the canonical chromosome pair.
4. Read the selected matrix metadata record and select the descriptor by `(unit, resolutionIndex)`.
5. For a materialized resolution, read its exact block index, map candidate logical block numbers to entries, fetch and independently decompress only the indexed blocks that exist, and decode their block streams.
6. For a derived resolution, query the declared materialized source over the expanded source-bin rectangle and aggregate raw values using Section K.
7. Filter decoded or aggregated cells to the exact requested half-open interval and transpose if the original chromosome order was noncanonical.
8. If requested, range-read the target-resolution normalization chunks and apply normalization.
9. If requested, range-read the matching target-resolution expected chunks and apply O/E.

A remote reader does not need to download structures unrelated to the requested chromosome pair, resolution, logical blocks, vector ranges, or expected-distance ranges.

### M. Required Writer Procedure

A conforming writer:

1. canonicalizes chromosome pairs and cis triangles;
2. constructs the logical resolution lists and validates every derived source;
3. stores exact raw counts as `COUNT_UINT` and genuine scores as `SCORE_FLOAT32`;
4. assigns every materialized contact to exactly one logical block;
5. sorts cells within each block, chooses a valid adaptive representation and value mode, and emits canonical varints;
6. compresses every non-empty logical block as its own independent Zstandard stored-block record;
7. builds an exact block index containing one locator for every stored logical block;
8. computes or preserves every advertised target-resolution normalization and expected vector, chunks it, applies the smallest supported exact transform, and writes independent Zstandard chunks;
9. emits matrix and vector indexes, then the footer;
10. backpatches every header locator and verifies every referenced interval and duplicated field.

Writers SHOULD target 65,536 `f32` values per vector chunk. Matrix block boundaries are determined exclusively by the logical block geometry; writers MUST NOT combine logical blocks merely to improve compression ratio.

### N. Validation and Defensive Reading

Before allocation, a reader MUST validate counts and products against the containing record length and its configured resource limits. In particular it MUST check:

- `headerByteLength`, footer length, index lengths, and every `entryByteLength`;
- `W * H`, bitmap byte count, value-slot count, and scalar decode count;
- block and vector decompressed sizes before invoking Zstandard;
- block-index length, strict block-number order, stored lengths, and non-overlapping indexed intervals;
- ULEB128 overflow and maximum length;
- derived-resolution factors and `u64` sum overflow;
- vector chunk coverage and `4 * valueCount` overflow;
- chromosome, resolution, normalization-type, and block-coordinate bounds.

A malformed optional matrix or vector is not silently treated as absent. Readers MUST report corruption rather than returning a plausible partial result.

### O. V9 Field Migration and Preserved Information

V10 is a new wire format selected by `version = 10`; a V9 reader cannot parse a V10 byte stream. The information carried by the established V9 fields is retained or replaced as follows. Conversion tools MUST use this mapping rather than dropping fields.

| V9 information | V10 representation |
|---|---|
| `magic`, `version` | Fixed header prefix |
| `footerPosition` | Fixed header `footerPosition` and `footerByteLength` |
| `genomeId` | Variable header `genomeId` `cstr` |
| normalization-index position and length | Fixed header `normVectorIndexPosition` and `normVectorIndexLength` |
| attributes | Counted key/value `cstr` pairs; unknown pairs preserved |
| chromosome name and length | `cstr` name plus `u64` length |
| BP and FRAG bin-size lists | Typed logical-resolution records, including storage/source semantics |
| per-chromosome fragment sites | Counted, strictly ordered `u64` site positions |
| `chr1Idx`, `chr2Idx`, `nResolutions` | Matrix header and footer numeric key |
| `unit`, `resIdx`, `binSize` | Matrix resolution descriptor and validation copy |
| `sumCounts` | Exact `u64` for counts or summary `f64` for scores |
| `occupiedCellCount` | Exact `u64` |
| `stdDev`, `percent95` | Explicit `f32`; unavailable values are NaN |
| `blockSize`, `blockColumnCount` | `blockBinCount`, `blockColumnCount`, and `gridType` |
| `blockCount` and per-block position/size index | Resolution `logicalBlockCount` and exact binary block index with `u64` positions |
| `nRecords` | Block `occupiedCellCount` |
| `binColumnOffset`, `binRowOffset` | Fixed block header offsets |
| short/int position flags | Canonical `uleb128` position streams |
| short/float contact flag | Explicit `COUNT_UINT` or `SCORE_FLOAT32` |
| list-of-rows/dense selector | `SPARSE_DELTA`, `BITMAP`, or `DENSE` |
| `rowCount`, `rowNumber`, `recordCount`, `binColumn` | Flat sorted local positions plus `occupiedCellCount` |
| dense width `w` | Explicit `blockWidth` and `blockHeight` |
| `allCountsOne` | `ALL_DEFAULT` with count value 1 |
| delta-column flag | Flat row-major position deltas |
| `-32768` or NaN absence sentinel | Zero count semantics or an explicit presence bitmap |
| master-index text key and matrix locator | Numeric chromosome-pair footer entry with `u64` locator and length |
| `nBytesV5` | Exact `footerByteLength`; the version-specific legacy name is retired |
| raw and normalized expected arrays | `EVI0` and `NEVI` entries plus exact chunked `f32` storage |
| chromosome scale factors | Sorted `(chrIndex u32, scaleFactor f32)` pairs |
| normalization-vector index and arrays | `NVI0` entries plus exact chunked `f32` storage |

V9 files already store chromosome lengths and many offsets as 64-bit values and V9 expected/normalization values as `f32`; a repacker MUST preserve those values exactly. Older source formats that contain `f64` vectors require an explicit, separately documented conversion policy and cannot claim bitwise lossless conversion to V10 `f32`.

Because some V9 writers emitted placeholder or low-precision matrix statistics, a repacker MUST compute the exact V10 `sumCountsOrScores` and `occupiedCellCount` from decoded logical cells. It MUST recompute `stdDev` and `percent95` under the V10 definitions or write the required canonical NaN; it MUST NOT reinterpret a legacy placeholder zero as a measured statistic.

---

## Part II — Resolution Policy, Semantics, and Rationale

### 1. Resolution Pyramid

#### 1.1 Logical resolutions versus materialized resolutions

V10 separates the concept of a **queryable resolution** from a **physically materialized matrix resolution**.

Every advertised resolution is one of:

- `MATERIALIZED`: contact matrix blocks are physically stored.
- `DERIVED`: raw matrix data is reconstructed by exact summation from a finer materialized resolution.

Derived resolutions remain first-class resolutions from the reader's perspective. For every normalization or expected-value capability the file advertises at that resolution, they retain their own target-resolution:

- bin size;
- expected-value vectors;
- normalized expected-value vectors;
- normalization vectors;
- chromosome scale factors;
- resolution metadata.

Only the redundant contact matrix blocks are omitted.

#### 1.2 Required high-resolution pyramid and standard resolution set

For a file extending from 10 bp through 2.5 Mb, the standard layout is:

| Resolution | Matrix storage | Source |
|---:|---|---:|
| 10 bp | **Materialized** | — |
| 20 bp | Derived | 10 bp |
| 50 bp | Derived | 10 bp |
| 100 bp | **Materialized** | — |
| 200 bp | Derived | 100 bp |
| 500 bp | Derived | 100 bp |
| 1 kb | **Materialized** | — |
| 2 kb | Derived | 1 kb |
| 5 kb | **Materialized** | — |
| 10 kb | **Materialized** | — |
| 25 kb | **Materialized** | — |
| 50 kb | **Materialized** | — |
| 100 kb | **Materialized** | — |
| 250 kb | **Materialized** | — |
| 500 kb | **Materialized** | — |
| 1 Mb | **Materialized** | — |
| 2.5 Mb | **Materialized** | — |

The five high-resolution virtual levels below are mandatory whenever advertised because duplicated matrix storage is most expensive there:

```text
10 bp
 ├── 20 bp   DERIVED
 └── 50 bp   DERIVED

100 bp
 ├── 200 bp  DERIVED
 └── 500 bp  DERIVED

1 kb
 └── 2 kb    DERIVED
```

Two especially important working resolutions remain directly stored:

- **5 kb**
- **50 kb**

These are kept because they are common interactive and analysis resolutions.

At resolutions finer than 1 kb, the 5× resolutions remain derived:

- **50 bp is derived from 10 bp**
- **500 bp is derived from 100 bp**

They are not materialized simply because they are 5× levels.

The established coarse-resolution hierarchy remains deliberately conservative and fully materialized:

```text
5 kb       MATERIALIZED
10 kb      MATERIALIZED
25 kb      MATERIALIZED
50 kb      MATERIALIZED
100 kb     MATERIALIZED
250 kb     MATERIALIZED
500 kb     MATERIALIZED
1 Mb       MATERIALIZED
2.5 Mb     MATERIALIZED
```

The 25 kb, 250 kb, 500 kb, and 2.5 Mb matrices remain materialized for compatibility with the established `.hic` hierarchy. They are comparatively inexpensive because the data are already highly aggregated, so derivation would add query overhead for little storage benefit. In particular, the small savings from deriving 500 kb do not justify its query-time fan-out.

The standard pyramid intentionally omits 20 kb and 200 kb. Although both can be produced exactly as 2× aggregations, their practical value does not justify the additional resolution metadata, normalization and expected-value vectors, reader complexity, and conformance testing. A specialized file MAY advertise them or other additional logical resolutions.

The format permits additional materialized or derived resolutions when a specialized file requires them, but no extension may override the five mandatory derivations or derive 500 kb. The table above defines the standard V10 hierarchy.

#### 1.3 Derived resolution semantics

Derivation always operates on the **raw unnormalized matrix**.

For target resolution `R_t` derived from source resolution `R_s`:

```text
factor = R_t / R_s
```

where the factor is normally `2` or `5`.

Each target contact is the exact sum of all source contacts whose bins fall within that target pixel.

For example, a 500 bp pixel derived from 100 bp data is:

```text
C500(i,j) =
    sum of C100(x,y)
    for the corresponding 5 × 5 region
```

This summation MUST occur before normalization.

##### Correct query order

A normalized 500 bp query therefore performs:

```text
100 bp raw contacts
        ↓
exact 5×5 summation
        ↓
500 bp raw contacts
        ↓
500 bp normalization vector
        ↓
normalized 500 bp contacts
```

For observed/expected:

```text
100 bp raw contacts
        ↓
exact aggregation to 500 bp
        ↓
500 bp normalization, if requested
        ↓
500 bp expected-value vector
        ↓
500 bp O/E result
```

Readers MUST NOT:

- sum normalized fine-resolution pixels;
- derive normalization vectors from finer resolutions;
- derive expected-value vectors from finer resolutions.

Those operations are not equivalent.

#### 1.4 Required auxiliary data for derived resolutions

If a V10 file advertises a normalization or expected-value capability at a derived resolution, the corresponding target-resolution data MUST already be present in the file.

For example, a file supporting KR-normalized 500 bp queries MUST contain the actual 500 bp KR normalization vectors.

Likewise, 500 bp O/E requires the actual 500 bp expected-value information.

Matrix storage may therefore be virtual while normalization and expected-value storage remains materialized.

#### 1.5 Aggregation semantics

Matrix derivation is valid only when the matrix has additive semantics.

V10 adds an explicit matrix aggregation field:

```text
aggregation = SUM | NONE
```

`SUM` permits derived resolutions.

`NONE` requires the resolution to be materialized.

Ordinary raw contact-count matrices use `SUM`.

Score matrices MUST NOT be implicitly aggregated unless the writer explicitly declares their aggregation semantics to be `SUM` and satisfies the finite-value and deterministic-arithmetic requirements in Section K.

#### 1.6 Source-resolution rule

A derived resolution MUST point directly to a materialized source resolution.

For the five fixed high-resolution targets, the direct source is not a writer choice: `20 -> 10`, `50 -> 10`, `200 -> 100`, `500 -> 100`, and `2000 -> 1000` BP. The 500,000 bp level is never derived.

For example:

```text
50 bp  -> 10 bp
500 bp -> 100 bp
2 kb   -> 1 kb
```

Readers should not need to execute chains such as:

```text
10 bp -> 20 bp -> 100 bp -> 500 bp
```

This keeps query cost predictable and prevents cascading decode operations.

---

### 2. Independently Compressed Logical Blocks

#### 2.1 The logical block is the random-access unit

V10 retains the concept of a **logical block** for genomic addressing and makes that same block the fundamental compression and random-access unit:

```text
Matrix
  └── Resolution
       ├── Independently compressed logical block
       ├── Independently compressed logical block
       └── ...
```

Every non-empty logical block has:

- one exact block-index entry;
- one stored-block header;
- one independent Zstandard frame;
- one logical-block payload.

A reader requesting one block never needs to read, decompress, or decode another block. This property controls read amplification for small local and remote queries.

#### 2.2 Physical ordering and read coalescing

Writers SHOULD store blocks from one matrix resolution contiguously in increasing block-number order. This layout lets readers combine adjacent indexed intervals into one disk or HTTP range read for broad queries without changing the independent compression boundary.

Physical ordering is an optimization only. Exact `u64` locators in the block index remain authoritative, so readers do not reconstruct positions from cumulative sizes and do not depend on block adjacency.

#### 2.3 Cis ordering

V10 retains the useful V9 principle of organizing cis blocks according to diagonal/anti-diagonal locality. The required rotated cis geometry and adaptive block scale in Section F.4 preserve distance locality and prevent ultra-fine cis maps from degenerating into enormous collections of nearly empty rectangular blocks.

Readers MUST compute both the nearest and farthest candidate distance bands for a cis query. They MUST NOT search every band from zero through the farthest band unless the requested rectangle actually crosses the diagonal.

---

### 3. Zstandard Compression

All stored matrix blocks and compressed vector chunks use **Zstandard** rather than ZLib. Section H fixes the interoperable matrix-block frame requirements and uses the framing defined by [RFC 8878](https://www.rfc-editor.org/rfc/rfc8878.html).

Each logical matrix block is independently decompressible.

The compression level is a writer policy, not part of the logical data model.

A reasonable default implementation should begin around:

```text
zstd level 6
```

and allow applications to select faster or denser writer profiles.

For example:

```text
--compression fast
--compression default
--compression compact
```

All profiles remain fully interoperable because decompression does not depend on the writer's selected compression level.

The reader MUST never require decompression of another matrix block to decode the requested block.

---

### 4. V10 Block Representations

Each logical block independently selects one of three representations:

```text
SPARSE_DELTA
BITMAP
DENSE
```

The writer chooses the smallest appropriate representation for that block.

This allows the format to adapt naturally across:

- extremely sparse fine-resolution blocks;
- intermediate-density blocks;
- dense near-diagonal blocks.

---

### 5. Sparse Delta Representation

#### 5.1 Remove list-of-rows overhead

V10 removes the V9-style sparse row hierarchy as the primary sparse representation.

Records are instead sorted in row-major order.

For a block of width `W`:

```text
localIndex = localRow * W + localColumn
```

The occupied positions become a monotonically increasing integer sequence:

```text
p0, p1, p2, p3, ...
```

V10 stores:

```text
p0
p1 - p0
p2 - p1
p3 - p2
...
```

using an unsigned variable-length integer representation.

This removes repeated:

- row numbers;
- row record counts;
- fixed-width column positions.

The sparse block becomes conceptually:

```text
nRecords
positionDeltaStream
valueStream
```

rather than:

```text
row
  record count
  column, value
  column, value

row
  record count
  column, value
...
```

#### 5.2 Position and value separation

Position data and contact values MUST be separate streams inside the block representation.

This allows each stream to expose its own statistical regularity before the logical block is compressed.

The compressor therefore sees runs of:

- small coordinate deltas;
- repeated or small count values;

rather than an interleaved coordinate/value structure.

---

### 6. Bitmap Representation

Intermediate-density blocks may use:

```text
BITMAP
```

The position stream is one bit per possible block cell:

```text
0 = absent
1 = occupied
```

Values are stored only for occupied cells, in bitmap traversal order.

The same value codecs used for sparse blocks are reused here.

A bitmap becomes preferable when the sparse position deltas cost more than approximately one bit per possible matrix position after considering block density.

The exact selection threshold is a writer implementation decision.

---

### 7. Dense Representation

Dense representation remains available.

V10 does **not** assume that fine-resolution data can never become dense.

For sufficiently dense blocks, storing no position information at all can remain optimal.

For count matrices:

```text
0 = no contact
positive integer = observed count
```

For scored matrices, the dense encoding always uses the explicit presence bitmap defined in Section I.4, so absence remains distinguishable from a numerical zero or NaN value.

The writer compares applicable representations and chooses the smallest.

---

### 8. Explicit Count Data Type

V10 introduces an explicit non-negative integer contact-count type.

Conceptually:

```text
COUNT_UINT
SCORE_FLOAT32
```

Raw contact counts MUST NOT be converted to floating-point merely because they exceed the range of a signed 16-bit integer.

`COUNT_UINT` is encoded as a canonical unsigned variable-length integer over the full `u64` range.

This gives V10:

- exact integer semantics;
- compact encoding of small counts;
- no `short` overflow limitation;
- no unnecessary float representation of integer counts.

`SCORE_FLOAT32` preserves exact 32-bit floating-point values for genuinely floating-point matrices.

---

### 9. Value-Stream Encoding

Fine-resolution contact maps frequently contain blocks in which nearly every occupied cell has the same count, usually `1`.

V10 therefore treats values independently from positions.

Each sparse or bitmap block selects a value mode. Dense blocks use `DIRECT`, as required by Section I.4.

#### 9.1 `ALL_DEFAULT`

For a sparse or bitmap block, if every occupied cell has the same value:

```text
defaultValue
```

is stored once.

No per-record values are written.

The common fine-resolution case becomes:

```text
defaultValue = 1
```

with zero bytes of value payload per occupied contact.

#### 9.2 `DEFAULT_EXCEPTIONS`

For a sparse or bitmap block, if most occupied cells have the default value:

```text
defaultValue
```

is stored once.

Only exceptional records are stored.

Exception locations refer to the **ordinal in the position stream**, not the genomic coordinate.

Example:

```text
record 0     value 1
record 1     value 1
...
record 137   value 2
...
record 921   value 3
```

may be represented approximately as:

```text
default = 1

exception ordinal deltas:
137, 784

exception values:
2, 3
```

The coordinate of the exceptional record is not duplicated.

#### 9.3 `DIRECT`

When the distribution is not dominated by one value, one value is stored per value slot. A sparse or bitmap value slot is an occupied record; a dense value slot is every cell in its `W × H` rectangle.

For count matrices these values are unsigned variable-length integers.

For score matrices they are exact `float32` values.

The writer chooses among these modes per block.

---

### 10. Stored Block Records

Each non-empty logical block is wrapped in one stored record:

```text
H10B stored-record header
one Zstandard frame
    └── one encoded logical block
```

The stored header supplies the codec, decompressed length, and logical block number. The independently compressed payload contains the Section I block header followed by its position and value streams. No stored record contains more than one logical block.

This one-to-one relationship makes read amplification explicit and bounded: fetching a requested block never entails decoding unrelated blocks. Readers may still coalesce adjacent stored-record byte ranges when a query needs several blocks.

---

### 11. Exact Block Index

V10 stores one fixed-width locator for every non-empty logical block:

```text
blockNumber
storedByteLength
absoluteFilePosition
```

The entries are strictly ordered by block number, so a reader can use binary search without decoding preceding entries or reconstructing cumulative positions. Empty logical blocks have no entry.

The modest index-space cost is deliberate. Exact `u64` positions avoid the read amplification and lookup ambiguity caused by grouping blocks behind range-only descriptors. Implementations may memory-map or cache an index and may merge adjacent block intervals into larger reads.

#### 11.1 Numeric matrix keys

Chromosome-pair matrix keys MUST use chromosome indices directly rather than constructing textual keys such as chromosome-index strings.

Units MUST use the `Unit` enum and normalization types MUST use identifiers into the header dictionary.

---

### 12. Lossless Compression of Normalization Vectors

Normalization vectors remain fully materialized for every logical resolution at which the corresponding normalization capability is advertised.

They are no longer stored as one large raw float array.

Each vector is divided into independently retrievable chunks.

A recommended starting chunk size is approximately:

```text
64K float values
```

with larger chunks permitted when appropriate.

Each chunk is:

1. transformed losslessly;
2. independently Zstandard-compressed;
3. independently indexed.

A viewport query therefore does not need to fetch an entire chromosome-scale normalization vector merely to normalize a small genomic interval.

#### 12.1 Exact float transforms

A vector chunk may choose between exact transforms such as:

```text
RAW
BYTE_SHUFFLE
XOR32
```

##### `BYTE_SHUFFLE`

The four bytes of each `float32` are transposed into four homogeneous byte streams before compression.

This groups similar exponent and high-mantissa bytes together.

##### `XOR32`

Each `float32` is treated as its exact 32-bit representation:

```text
encoded[i] = bits[i] XOR bits[i - 1]
```

No numerical conversion is performed.

Smooth or repeated vector values therefore produce many zero bits while preserving the original bit pattern exactly.

The writer may try the supported transformations and retain the smallest encoded chunk.

There is:

- no float16 conversion;
- no rounding;
- no fixed-point approximation.

---

### 13. Lossless Compression of Expected-Value Vectors

Expected-value vectors and normalized expected-value vectors use the same chunked compression framework.

They remain explicitly stored at every logical resolution at which they are advertised.

An expected-value index SHOULD allow readers to directly locate:

```text
unit
binSize
normalization type, if applicable
distance chunk
```

without reading unrelated expected vectors.

Expected vectors are often smooth and should therefore benefit strongly from exact XOR and byte-shuffle transforms before Zstandard compression.

Chromosome scale factors remain stored exactly.

---

### 14. Derived-Resolution Query Path

A query to a derived resolution proceeds as follows.

Example:

```text
500 bp KR query
```

##### Step 1 — identify source

```text
500 bp -> source 100 bp
```

##### Step 2 — identify source blocks

The requested genomic rectangle is converted to the corresponding 100 bp source region.

The exact block index identifies the stored source blocks overlapping that region.

##### Step 3 — decode blocks

Only those independently compressed blocks are fetched and decompressed.

##### Step 4 — reconstruct source contacts

Sparse, bitmap, or dense blocks are decoded into their raw contacts or scores with the declared `valueType`.

##### Step 5 — aggregate

100 bp contacts are summed into the corresponding 500 bp pixels.

The result is the exact raw 500 bp matrix for the requested region.

##### Step 6 — normalize

The stored **500 bp KR normalization vector** is applied.

##### Step 7 — expected values, if requested

The stored **500 bp expected or KR-expected vector** is applied.

No normalization or expected-value quantity is inferred from the 100 bp level.

---

### 15. Resolution-Family Alignment

Derived resolution families SHOULD be designed so that 2× and 5× aggregation maps cleanly onto source storage.

For example:

```text
10 bp   -> 20 bp, 50 bp
100 bp  -> 200 bp, 500 bp
1 kb    -> 2 kb
```

Source block dimensions SHOULD, where practical, be divisible by both `2` and `5`.

The rotated or rectangular block geometry SHOULD preserve genomic locality. This minimizes source-block fan-out when serving derived queries and improves reuse between zoom levels.

The optimization does not alter the mathematical bin boundaries: all target bins remain aligned to chromosome coordinate zero and are exact integer aggregations of the source bins.

---

### 16. Reader Caching

Caching is not part of the persistent binary semantics, but V10 readers SHOULD maintain two independent caches.

#### 16.1 Decompressed block cache

Keyed approximately by:

```text
chromosome pair
materialized resolution
logical block number
```

This avoids repeated network reads and Zstandard decompression while panning.

#### 16.2 Derived tile cache

Keyed approximately by:

```text
chromosome pair
target resolution
genomic tile
```

This stores already aggregated 20 bp, 50 bp, 200 bp, 500 bp, or 2 kb results.

This is particularly useful when:

- zooming in and out;
- repeatedly drawing the same region;
- applying different normalization modes to the same raw derived tile.

The raw derived tile can be cached before normalization and reused for KR, VC, NONE, and O/E requests.

---

### 17. Resolution Metadata

Each logical resolution carries the following core metadata. Section E.2 defines the complete 76-byte wire record and is authoritative.

```text
unit
storageMode
aggregation
valueType
resolutionIndex
binSize
sourceResolutionIndex
gridType
sumCountsOrScores
occupiedCellCount
stdDev
percent95
blockBinCount
blockColumnCount
blockIndexPosition
blockIndexLength
logicalBlockCount
reserved1
```

For example:

```text
binSize       = 500
storageMode   = DERIVED
sourceResolutionIndex = index of 100 bp
aggregation   = SUM
```

A materialized resolution has:

```text
storageMode   = MATERIALIZED
sourceResolutionIndex = 0xFFFFFFFF
```

Derived resolutions do not have their own matrix block index because no matrix blocks exist for that resolution.

They still participate normally in normalization and expected-value indexes.

---

### 18. V10 File Organization

Conceptually:

```text
HEADER
    magic
    version = 10
    genome
    attributes
    chromosomes
    logical resolutions
    fragment sites
    normalization type dictionary

MATRIX METADATA RECORDS
    chromosome pair
        resolution descriptors

MATERIALIZED MATRIX DATA
    exact block index
    independently compressed matrix blocks
    exact block index
    independently compressed matrix blocks
    ...

NORMALIZATION-VECTOR CHUNKS
NORMALIZATION-VECTOR INDEX

EXPECTED-VALUE CHUNKS
EXPECTED-VALUE INDEX

NORMALIZED-EXPECTED CHUNKS
NORMALIZED-EXPECTED INDEX

FOOTER / MATRIX MASTER INDEX
```

Part I defines the exact byte ordering. The conceptual lookup path is:

```text
logical resolution
        ↓
materialized or derived
        ↓
logical block
        ↓
position stream + value stream
        ↓
Zstandard
```

---

### 19. V9 → V10 Lossless Repacking

A `hic repack` command should be a first-class part of V10 deployment.

Its purpose is to convert existing large V9 files without rerunning alignment or contact generation.

#### 19.1 Matrix conversion

For each matrix and materialized resolution:

1. decompress V9 blocks;
2. reconstruct raw records;
3. sort records into V10 logical-block order;
4. select sparse/bitmap/dense representation;
5. select the best value mode;
6. wrap and Zstandard-compress every non-empty logical block independently;
7. record each block's exact position and stored length;
8. construct the exact binary block index.

#### 19.2 Virtual-resolution verification

Before deleting an existing V9 resolution that is intended to become derived, the repacker MUST verify it.

For example, before discarding the V9 500 bp matrix:

```text
aggregate existing 100 bp raw matrix to 500 bp
```

and compare the result exactly against the existing V9 500 bp raw matrix.

Only if they are identical may the 500 bp matrix blocks be removed.

This can be exposed as:

```text
hic repack --verify-derived-resolutions
```

If verification fails:

- conversion of a mandatory target (20, 50, 200, 500, or 2,000 bp) fails;
- a nonstandard optional derived target may remain materialized, unless strict mode requests failure.

There is no silent approximation.

#### 19.3 Preserve vectors

Existing normalization and expected-value vectors are copied exactly into the V10 vector encoding.

The repacker changes their physical compression representation, not their numerical values.

This makes V9 → V10 repacking a lossless transformation.

---

### 20. Correctness Requirements

A conforming V10 implementation MUST satisfy the following.

##### Raw matrix invariance

For every materialized or derived logical resolution:

```text
V10 raw query == corresponding exact V9/raw matrix
```

where the source data represents an additive contact-count matrix.

##### Derived matrix invariance

For a derived resolution:

```text
derived target matrix ==
exact integer aggregation of source matrix
```

##### Vector invariance

After decompression:

```text
normalization float32 bits ==
original normalization float32 bits
```

and:

```text
expected float32 bits ==
original expected float32 bits
```

##### No lossy transforms

V10 core matrix and vector storage performs no:

- count thresholding;
- distance truncation;
- count quantization;
- float quantization;
- precision reduction;
- stochastic sampling.

##### Bounded random access

No ordinary regional query should require decompressing an entire chromosome matrix, entire resolution, or entire file.

Every matrix block and vector chunk is independently addressable.

---

### 21. Benchmarking and Validation Requirements

V10 implementations SHOULD be benchmarked using representative deeply sequenced `.hic` files, including the largest available sub-kilobase datasets.

At minimum, measure:

##### File-size components

```text
stored matrix blocks
matrix indexes
normalization vectors
expected vectors
other metadata
total file
```

##### Query latency

For both local disk and HTTP range access:

```text
p50 latency
p95 latency
p99 latency
bytes fetched per query
blocks fetched per query
blocks decompressed per query
unrequested blocks decompressed per query
```

Test separately:

- materialized resolution query;
- 2× derived resolution query;
- 5× derived resolution query;
- near-diagonal cis;
- far-cis;
- trans;
- cold cache;
- warm cache.

Particular attention should be paid to:

```text
50 bp  from 10 bp
500 bp from 100 bp
2 kb   from 1 kb
```

because these exercise the virtual-resolution design.

##### Decode throughput

Measure:

```text
SPARSE_DELTA
BITMAP
DENSE
ALL_DEFAULT
DEFAULT_EXCEPTIONS
DIRECT
```

independently.

The format should optimize total interactive query cost, not compression ratio in isolation.

---

### 22. Required V10 Baseline

The required V10 baseline is therefore:

##### Resolution storage

**Materialized by default:**

```text
10 bp
100 bp
1 kb
5 kb
10 kb
25 kb
50 kb
100 kb
250 kb
500 kb
1 Mb
2.5 Mb
```

**Derived by exact raw summation:**

```text
20 bp
50 bp
200 bp
500 bp
2 kb
```

Every advertised normalization or expected-value capability uses its own target-resolution data, including at derived matrix resolutions.

##### Matrix compression

```text
logical blocks
    ↓
SPARSE_DELTA / BITMAP / DENSE
    ↓
separate positions and values
    ↓
ALL_DEFAULT / DEFAULT_EXCEPTIONS / DIRECT
    ↓
one independent Zstandard frame per logical block
```

##### Count representation

```text
non-negative integer counts
variable-length encoded
64-bit-capable
```

rather than `short` with float fallback.

##### Vector compression

```text
normalization / expected vectors
    ↓
independent chunks
    ↓
RAW / BYTE_SHUFFLE / XOR32
    ↓
Zstandard
```

with exact reconstruction of every original float32 bit.

##### Indexing

```text
exact block-level indexes
strictly ordered block numbers
absolute `u64` positions
explicit stored byte lengths
```

rather than a range-only index over groups of compressed logical blocks.

---

### 23. Summary

V10 should achieve its compression gains primarily by eliminating **structural redundancy**, not by discarding biological information.

The largest change is the resolution pyramid:

> Store the raw matrix at selected anchor resolutions, while treating intermediate 2× and 5× matrices as exact virtual views over those anchors.

The important working resolutions **5 kb and 50 kb** remain materialized. The established coarse levels **25 kb, 250 kb, 500 kb, 1 Mb, and 2.5 Mb** are also retained physically for compatibility and inexpensive direct access. The five mandatory virtual levels are **20 bp, 50 bp, 200 bp, 500 bp, and 2 kb**. The standard hierarchy does not advertise 20 kb or 200 kb.

The second major change is to retain V9-compatible rotated cis distance bands and adaptive, increasingly large blocks at fine resolutions, while redesigning the contents of those blocks around what the data actually look like:

> a sorted sparse set of occupied positions whose values are overwhelmingly small integers and frequently equal to one.

Flat position deltas, implicit values, bitmap/dense alternatives, independent block compression, and Zstandard exploit those properties without changing a single contact count.

Finally, auxiliary arrays and indexes become compressed and range-addressable rather than remaining disproportionately expensive metadata at ultra-fine resolutions.

The resulting V10 format remains an **interactive `.hic` file**, not an archive:

- exact;
- randomly queryable;
- normalization-aware;
- O/E-aware;
- remote-range friendly;
- substantially less redundant than V9.
