# `.hic` V10: Lossless, Query-Optimized Format Proposal

## Status

**Proposed V10 strategy**

V10 is intended to substantially reduce `.hic` file size for deeply sequenced, high-resolution Hi-C maps while preserving:

- exact raw contact counts and scores;
- interactive local and remote querying;
- all advertised normalization vectors;
- all advertised expected-value vectors;
- exact existing resolution semantics;
- random access to chromosome pairs and genomic windows.

V10 is **not** a cold-storage format and does not introduce lossy approximations.

The central design principle is:

> **Do not store information that can be reconstructed exactly from another stored representation, but do precompute and store information that cannot be reconstructed cheaply or exactly at query time.**

The major V10 changes are:

1. Store only selected matrix-resolution anchors and derive redundant 2×/5× resolutions from raw counts.
2. Store normalization and expected-value vectors for every advertised logical resolution, including derived resolutions.
3. Replace ZLib with Zstandard.
4. Compress neighboring logical blocks together into bounded random-access pages.
5. Replace the sparse list-of-rows representation with flat delta-coded cell positions.
6. Separate position and value streams.
7. Exploit the count≈1 regime with implicit/default values and exception streams.
8. Add adaptive sparse, bitmap, and dense block representations.
9. Introduce an explicit integer count type rather than falling back from `short` to floating point.
10. Compress normalization and expected-value vectors losslessly in independently accessible chunks.
11. Replace the per-block absolute-position index with a compact page-oriented index.

---

# 1. Resolution Pyramid

## 1.1 Logical resolutions versus materialized resolutions

V10 separates the concept of a **queryable resolution** from a **physically materialized matrix resolution**.

Every advertised resolution is one of:

- `MATERIALIZED`: contact matrix blocks are physically stored.
- `DERIVED`: raw matrix data is reconstructed by exact summation from a finer materialized resolution.

Derived resolutions remain first-class resolutions from the reader's perspective. They retain their own:

- bin size;
- expected-value vectors;
- normalized expected-value vectors;
- normalization vectors;
- chromosome scale factors;
- resolution metadata.

Only the redundant contact matrix blocks are omitted.

## 1.2 Recommended V10 resolution set

For a file extending from 10 bp through 1 Mb, the recommended layout is:

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
| 20 kb | Derived | 10 kb |
| 50 kb | **Materialized** | — |
| 100 kb | **Materialized** | — |
| 200 kb | Derived | 100 kb |
| 500 kb | Derived | 100 kb |
| 1 Mb | **Materialized** | — |

Thus the normal anchors are powers of ten, with two intentionally materialized 5× exceptions:

- **5 kb**
- **50 kb**

These are kept because they are common interactive and analysis resolutions.

At resolutions finer than 1 kb, the 5× resolutions remain derived:

- **50 bp is derived from 10 bp**
- **500 bp is derived from 100 bp**

They are not materialized simply because they are 5× levels.

The format itself MUST permit additional materialized resolutions if a specialized file requires them, but the table above is the recommended default.

## 1.3 Derived resolution semantics

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

### Correct query order

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

## 1.4 Required auxiliary data for derived resolutions

If a V10 file advertises a normalization or expected-value capability at a derived resolution, the corresponding target-resolution data MUST already be present in the file.

For example, a file supporting KR-normalized 500 bp queries MUST contain the actual 500 bp KR normalization vectors.

Likewise, 500 bp O/E requires the actual 500 bp expected-value information.

Matrix storage may therefore be virtual while normalization and expected-value storage remains materialized.

## 1.5 Aggregation semantics

Matrix derivation is valid only when the matrix has additive semantics.

V10 adds an explicit matrix aggregation field:

```text
aggregation = SUM | NONE
```

`SUM` permits derived resolutions.

`NONE` requires the resolution to be materialized.

Ordinary raw contact-count matrices use `SUM`.

Arbitrary score matrices MUST NOT be implicitly aggregated unless the writer explicitly declares their aggregation semantics to be `SUM`.

## 1.6 Source-resolution rule

A derived resolution SHOULD point directly to a materialized source resolution.

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

# 2. Logical Blocks and Compression Pages

## 2.1 Separate spatial blocks from compression units

V10 retains the concept of a **logical block** for genomic addressing.

However, logical blocks are no longer individually compressed.

Instead, several spatially adjacent logical blocks are grouped into a **compression page**.

```text
Matrix
  └── Resolution
       └── Compression Page
            ├── Logical Block
            ├── Logical Block
            ├── Logical Block
            └── ...
```

The page is the fundamental:

- compressed unit;
- disk-read unit;
- HTTP range-read unit;
- decompression unit;
- cache unit.

Logical blocks remain the spatial indexing unit inside the page.

## 2.2 Page sizing

The recommended target is approximately:

```text
128 KiB compressed per page
```

with a practical target range of approximately:

```text
64–256 KiB
```

Writers MAY adjust page size according to local density, but SHOULD impose a bounded maximum page size.

Pages MUST NOT span unrelated chromosome-pair matrices.

Pages SHOULD contain blocks that are adjacent in block-number/genomic order.

This provides enough data for the compressor to exploit statistical redundancy without turning a small viewport query into a multi-megabyte read.

## 2.3 Cis ordering

V10 retains the useful V9 principle of organizing cis blocks according to diagonal/anti-diagonal locality.

Blocks grouped into the same page SHOULD come from the same or neighboring distance bands.

This improves both:

- compression similarity;
- probability that blocks are requested together.

V10 does not require replacing the existing rotated cis geometry solely for compression purposes.

---

# 3. Zstandard Compression

All matrix pages and compressed vector chunks use **Zstandard** rather than ZLib.

Each page is independently decompressible.

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

The reader MUST never require decompression of data preceding the requested page.

---

# 4. V10 Block Representations

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

# 5. Sparse Delta Representation

## 5.1 Remove list-of-rows overhead

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

## 5.2 Position and value separation

Position data and contact values MUST be separate streams inside the block representation.

This allows each stream to expose its own statistical regularity before page compression.

The compressor therefore sees runs of:

- small coordinate deltas;
- repeated or small count values;

rather than an interleaved coordinate/value structure.

---

# 6. Bitmap Representation

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

# 7. Dense Representation

Dense representation remains available.

V10 does **not** assume that fine-resolution data can never become dense.

For sufficiently dense blocks, storing no position information at all can remain optimal.

For count matrices:

```text
0 = no contact
positive integer = observed count
```

For generic scored matrices where absence must be distinguishable from a numerical zero or NaN value, the dense encoding MUST retain an explicit presence representation when required.

The writer compares applicable representations and chooses the smallest.

---

# 8. Explicit Count Data Type

V10 introduces an explicit non-negative integer contact-count type.

Conceptually:

```text
COUNT_UINT
SCORE_FLOAT32
```

Raw contact counts MUST NOT be converted to floating-point merely because they exceed the range of a signed 16-bit integer.

`COUNT_UINT` is encoded as an unsigned variable-length integer with support for at least 64-bit count values.

This gives V10:

- exact integer semantics;
- compact encoding of small counts;
- no `short` overflow limitation;
- no unnecessary float representation of integer counts.

`SCORE_FLOAT32` preserves exact 32-bit floating-point values for genuinely floating-point matrices.

---

# 9. Value-Stream Encoding

Fine-resolution contact maps frequently contain blocks in which nearly every occupied cell has the same count, usually `1`.

V10 therefore treats values independently from positions.

Each block selects a value mode.

## 9.1 `ALL_DEFAULT`

If every occupied cell has the same value:

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

## 9.2 `DEFAULT_EXCEPTIONS`

If most cells have the default value:

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

## 9.3 `DIRECT`

When the distribution is not dominated by one value, one value is stored per occupied record.

For count matrices these values are unsigned variable-length integers.

For score matrices they are exact `float32` values.

The writer chooses among these modes per block.

---

# 10. Page Contents

A matrix page contains:

```text
Page Header
Block Directory
Encoded Logical Blocks
```

The page header includes at minimum:

```text
codec
uncompressedSize
blockCount
```

The internal block directory identifies the logical block numbers contained in the page and their locations within the decompressed page.

Once a page has been decompressed, multiple nearby block queries therefore require no additional I/O or decompression.

---

# 11. Compact Page-Oriented Index

V10 replaces the large per-block:

```text
blockNumber
absoluteFilePosition
compressedSize
```

index with a page-oriented index.

## 11.1 Resolution-level page directory

For each materialized matrix resolution, the metadata records pages in physical file order.

Page descriptors contain enough information to identify:

- the range/set of block numbers in the page;
- compressed page size;
- page location.

Because page data is written sequentially, file positions SHOULD be represented primarily by:

```text
baseOffset + cumulative compressed sizes
```

rather than storing a complete 64-bit absolute file offset for every page.

## 11.2 Delta coding

Monotonically increasing quantities SHOULD use unsigned deltas:

```text
page/block identifiers
file positions
page lengths where beneficial
```

These deltas are variable-length encoded.

## 11.3 Checkpoints

To avoid scanning the entire size list to locate an arbitrary page, the index includes periodic absolute checkpoints.

A recommended initial interval is:

```text
one checkpoint every 64 pages
```

A lookup therefore performs:

1. locate the nearest checkpoint;
2. decode at most a bounded number of descriptors;
3. issue the corresponding range read.

This keeps the index both compact and random-access friendly.

## 11.4 Numeric matrix keys

Chromosome-pair matrix keys SHOULD use chromosome indices directly rather than constructing textual keys such as chromosome-index strings.

Repeated strings such as units and normalization types SHOULD similarly use compact enum/dictionary identifiers in V10 metadata.

---

# 12. Lossless Compression of Normalization Vectors

Normalization vectors remain fully materialized for every supported logical resolution.

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

## 12.1 Exact float transforms

A vector chunk may choose between exact transforms such as:

```text
RAW
BYTE_SHUFFLE
XOR32
```

### `BYTE_SHUFFLE`

The four bytes of each `float32` are transposed into four homogeneous byte streams before compression.

This groups similar exponent and high-mantissa bytes together.

### `XOR32`

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

# 13. Lossless Compression of Expected-Value Vectors

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

# 14. Derived-Resolution Query Path

A query to a derived resolution proceeds as follows.

Example:

```text
500 bp KR query
```

### Step 1 — identify source

```text
500 bp -> source 100 bp
```

### Step 2 — identify source pages

The requested genomic rectangle is converted to the corresponding 100 bp source region.

The page index identifies the source pages overlapping that region.

### Step 3 — decode pages

Only those pages are fetched and decompressed.

### Step 4 — reconstruct source contacts

Sparse, bitmap, or dense blocks are decoded into their raw integer contacts.

### Step 5 — aggregate

100 bp contacts are summed into the corresponding 500 bp pixels.

The result is the exact raw 500 bp matrix for the requested region.

### Step 6 — normalize

The stored **500 bp KR normalization vector** is applied.

### Step 7 — expected values, if requested

The stored **500 bp expected or KR-expected vector** is applied.

No normalization or expected-value quantity is inferred from the 100 bp level.

---

# 15. Resolution-Family Alignment

Derived resolution families SHOULD be designed so that 2× and 5× aggregation maps cleanly onto source storage.

For example:

```text
10 bp   -> 20 bp, 50 bp
100 bp  -> 200 bp, 500 bp
1 kb    -> 2 kb
10 kb   -> 20 kb
100 kb  -> 200 kb, 500 kb
```

Source block dimensions SHOULD, where practical, be divisible by both `2` and `5`.

Pages SHOULD preserve genomic locality.

This minimizes source-page fan-out when serving derived queries and improves reuse between zoom levels.

The optimization does not alter the mathematical bin boundaries: all target bins remain aligned to chromosome coordinate zero and are exact integer aggregations of the source bins.

---

# 16. Reader Caching

Caching is not part of the persistent binary semantics, but V10 readers SHOULD maintain two independent caches.

## 16.1 Decompressed page cache

Keyed approximately by:

```text
chromosome pair
materialized resolution
page ID
```

This avoids repeated network reads and Zstandard decompression while panning.

## 16.2 Derived tile cache

Keyed approximately by:

```text
chromosome pair
target resolution
genomic tile
```

This stores already aggregated 20/50/200/500/etc. results.

This is particularly useful when:

- zooming in and out;
- repeatedly drawing the same region;
- applying different normalization modes to the same raw derived tile.

The raw derived tile can be cached before normalization and reused for KR, VC, NONE, and O/E requests.

---

# 17. Proposed Resolution Metadata

Each logical resolution should include metadata equivalent to:

```text
unit
binSize
storageMode
sourceBinSize
aggregation
sumCounts
occupiedCellCount
```

For example:

```text
binSize       = 500
storageMode   = DERIVED
sourceBinSize = 100
aggregation   = SUM
```

A materialized resolution has:

```text
storageMode   = MATERIALIZED
sourceBinSize = 0
```

Derived resolutions do not have their own matrix page index because no matrix pages exist for that resolution.

They still participate normally in normalization and expected-value indexes.

---

# 18. Proposed V10 File Organization

Conceptually:

```text
HEADER
    magic
    version = 10
    genome
    chromosomes
    logical resolutions
    attributes

MATRIX DIRECTORY
    chromosome pair
        resolution descriptors

MATERIALIZED MATRIX DATA
    page index
    compressed matrix pages
    page index
    compressed matrix pages
    ...

EXPECTED-VALUE INDEX
EXPECTED-VALUE CHUNKS

NORMALIZED-EXPECTED INDEX
NORMALIZED-EXPECTED CHUNKS

NORMALIZATION-VECTOR INDEX
NORMALIZATION-VECTOR CHUNKS

MASTER INDEX
```

Exact byte ordering may remain close to V9 where useful, but V10 readers should reason in terms of:

```text
logical resolution
        ↓
materialized or derived
        ↓
logical block
        ↓
compression page
        ↓
position stream + value stream
        ↓
Zstandard
```

---

# 19. V9 → V10 Lossless Repacking

A `hic repack` command should be a first-class part of V10 deployment.

Its purpose is to convert existing large V9 files without rerunning alignment or contact generation.

## 19.1 Matrix conversion

For each matrix and materialized resolution:

1. decompress V9 blocks;
2. reconstruct raw records;
3. sort records into V10 logical-block order;
4. select sparse/bitmap/dense representation;
5. select the best value mode;
6. group blocks into pages;
7. Zstandard-compress pages;
8. construct the compact page index.

## 19.2 Virtual-resolution verification

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

- the resolution remains materialized; or
- the conversion fails in strict mode.

There is no silent approximation.

## 19.3 Preserve vectors

Existing normalization and expected-value vectors are copied exactly into the V10 vector encoding.

The repacker changes their physical compression representation, not their numerical values.

This makes V9 → V10 repacking a lossless transformation.

---

# 20. Correctness Requirements

A conforming V10 implementation MUST satisfy the following.

### Raw matrix invariance

For every materialized or derived logical resolution:

```text
V10 raw query == corresponding exact V9/raw matrix
```

where the source data represents an additive contact-count matrix.

### Derived matrix invariance

For a derived resolution:

```text
derived target matrix ==
exact integer aggregation of source matrix
```

### Vector invariance

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

### No lossy transforms

V10 core matrix and vector storage performs no:

- count thresholding;
- distance truncation;
- count quantization;
- float quantization;
- precision reduction;
- stochastic sampling.

### Bounded random access

No ordinary regional query should require decompressing an entire chromosome matrix, entire resolution, or entire file.

Every matrix page and vector chunk is independently addressable.

---

# 21. Benchmarking Requirements Before Format Freeze

V10 should be benchmarked using representative deeply sequenced `.hic` files, including the largest available sub-kilobase datasets.

At minimum, measure:

### File-size components

```text
matrix pages
matrix indexes
normalization vectors
expected vectors
other metadata
total file
```

### Query latency

For both local disk and HTTP range access:

```text
p50 latency
p95 latency
p99 latency
bytes fetched per query
pages decompressed per query
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
20 kb  from 10 kb
500 kb from 100 kb
```

because these exercise the virtual-resolution design.

### Decode throughput

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

# 22. Recommended V10 Baseline

The recommended V10 implementation is therefore:

### Resolution storage

**Materialized by default:**

```text
10 bp
100 bp
1 kb
5 kb
10 kb
50 kb
100 kb
1 Mb
```

**Derived by exact raw summation:**

```text
20 bp
50 bp
200 bp
500 bp
2 kb
20 kb
200 kb
500 kb
```

All advertised resolutions retain their own normalization and expected-value data.

### Matrix compression

```text
logical blocks
    ↓
SPARSE_DELTA / BITMAP / DENSE
    ↓
separate positions and values
    ↓
ALL_DEFAULT / DEFAULT_EXCEPTIONS / DIRECT
    ↓
bounded multi-block pages
    ↓
Zstandard
```

### Count representation

```text
non-negative integer counts
variable-length encoded
64-bit-capable
```

rather than `short` with float fallback.

### Vector compression

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

### Indexing

```text
page-level indexes
delta-coded identifiers
delta/cumulative positions
periodic absolute checkpoints
```

rather than an absolute position and size for every compressed logical block.

---

# 23. Summary

V10 should achieve its compression gains primarily by eliminating **structural redundancy**, not by discarding biological information.

The largest change is the resolution pyramid:

> Store the raw matrix at selected anchor resolutions, while treating intermediate 2× and 5× matrices as exact virtual views over those anchors.

The important exceptions are **5 kb and 50 kb**, which remain materialized because they are common working resolutions. The sub-kilobase 5× resolutions, **50 bp and 500 bp**, remain derived.

The second major change is to redesign fine-resolution blocks around what the data actually look like:

> a sorted sparse set of occupied positions whose values are overwhelmingly small integers and frequently equal to one.

Flat position deltas, implicit values, bitmap/dense alternatives, page compression, and Zstandard exploit those properties without changing a single contact count.

Finally, auxiliary arrays and indexes become compressed and range-addressable rather than remaining disproportionately expensive metadata at ultra-fine resolutions.

The resulting V10 format remains an **interactive `.hic` file**, not an archive:

- exact;
- randomly queryable;
- normalization-aware;
- O/E-aware;
- remote-range friendly;
- substantially less redundant than V9.