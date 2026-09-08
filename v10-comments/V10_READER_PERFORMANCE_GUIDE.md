# HIC V10 Reader Performance Guide

## Purpose

HIC V10 can provide substantially smaller files without making materialized-resolution queries
slower than V9. Achieving that result depends on the reader preserving the format's random-access
properties all the way through block selection, I/O, decompression, and contact decoding.

This document describes the implementation requirements and optimizations needed to obtain that
performance. It does not change the V10 wire format. The normative format remains defined in
`HiCFormatV10.md`.

The most important rule is:

> A regional query must identify and decode only the logical blocks that can intersect the requested
> region.

Violating that rule can make a valid V10 file appear 10–20 times slower than V9, especially for
fine-resolution far-cis queries. The compact block codec can also be made unnecessarily expensive
if a reader allocates intermediate arrays or performs general-purpose validation work for every
decoded byte.

## Performance model

The cost of a materialized-resolution query should be approximately:

```text
index lookup
+ bytes read for intersecting non-empty blocks
+ Zstandard decompression of those blocks
+ one streaming pass over their encoded contacts
+ output callback cost
```

It should not include:

- scanning every cis distance band between the diagonal and the query;
- reading or decompressing non-intersecting blocks;
- linearly testing every occupied block in the resolution when the exact index is available;
- reconstructing full position and value arrays before emitting contacts;
- creating temporary diagnostic strings for successful bounds checks;
- constructing a new Zstandard decompression context for every block;
- repeated parsing of matrix descriptors or block indexes for consecutive queries;
- per-contact type-erased callback dispatch inside multiple internal reader layers.

Derived resolutions have an additional aggregation cost and must be benchmarked separately from
materialized resolutions. A derived query being slower does not justify a slowdown at a physically
stored resolution.

## 1. Cache the random-access metadata

Opening a V10 file should parse and validate the fixed header and footer once. A long-lived reader
should then cache, on first use:

- the numeric chromosome-pair matrix locator;
- the matrix's resolution descriptors;
- each materialized resolution's exact block index;
- normalization and expected-vector index metadata where applicable.

The block index is ordered by `blockNumber`. Keep it in an ordered contiguous representation such
as a vector. It is small relative to the matrix payload and is designed for binary lookup. Do not
read and revalidate it for every query window.

An implementation may additionally build a hash map from block number to entry, but this is not
required. Binary search in the sorted vector avoids hash-table overhead and supports efficient
iteration over contiguous block-number ranges.

## 2. Compute exact candidate block-number ranges

### 2.1 Do not filter the complete block index

The reader should transform the query rectangle into one or more inclusive block-number ranges.
For each range, use `lower_bound`, an equivalent binary search, or a keyed lookup to find the first
stored block. Iterate only until the end of that range.

The following pattern is intentionally avoided:

```cpp
for (const BlockEntry& entry : completeBlockIndex) {
    if (mightIntersect(entry.blockNumber, query)) {
        readAndDecode(entry);
    }
}
```

Even if the predicate is correct, this pays a linear index scan for every query. If the predicate
is overly conservative, it also causes expensive read and decompression amplification.

The desired pattern is:

```cpp
for (BlockRange range : candidateRanges(query)) {
    auto it = lower_bound(index.begin(), index.end(), range.first);
    while (it != index.end() && it->blockNumber <= range.last) {
        readAndDecode(*it++);
    }
}
```

Empty logical blocks have no index entries and are skipped naturally.

### 2.2 Rectangular trans grid

For a half-open bin rectangle `[x0, x1) × [y0, y1)` and block dimension `B`:

```text
firstColumn = floor(x0 / B)
lastColumn  = floor((x1 - 1) / B)
firstRow    = floor(y0 / B)
lastRow     = floor((y1 - 1) / B)
```

For each row in the inclusive row range, request:

```text
[row * blockColumnCount + firstColumn,
 row * blockColumnCount + lastColumn]
```

These ranges are already ordered and non-overlapping.

### 2.3 Rotated cis grid

V10 retains V9's rotated cis geometry. One axis is the position along the diagonal and the other is
the logarithmic distance band away from the diagonal.

For a half-open query rectangle, first compute inclusive endpoints:

```text
lastX = x1 - 1
lastY = y1 - 1
```

The enclosing diagonal-position range is:

```text
firstPad = floor((x0 + y0) / (2 * B))
lastPad  = floor((lastX + lastY) / (2 * B))
```

The nearest absolute distance represented by the rectangle is:

```text
if x1 <= y0:
    nearest = y0 - lastX
else if y1 <= x0:
    nearest = x0 - lastY
else:
    nearest = 0
```

The farthest distance is the maximum absolute difference at the two opposite corners:

```text
farthest = max(abs(lastY - x0), abs(lastX - y0))
```

Map both distances through the same distance-to-depth function used by block-number assignment.
For every depth from `firstDepth` through `lastDepth`, request the inclusive range:

```text
[depth * blockColumnCount + firstPad,
 depth * blockColumnCount + lastPad]
```

The reader must not request depths zero through `lastDepth` unless the query rectangle actually
reaches or crosses the diagonal. Doing so is the primary cause of the previously observed 10–20×
far-cis slowdown: dense near-diagonal blocks are read and decoded for a sparse far-cis query and
their contacts are discarded only afterward.

The distance-to-depth calculation must be identical to the block geometry calculation. An exact
integer implementation avoids floating-point disagreement at band boundaries. One equivalent
form is:

```text
lhs   = distance²
scale = 2 * B²
depth = 0

while depth < 32:
    threshold = 2^(depth + 1) - 1
    if threshold² > floor(lhs / scale):
        break
    depth += 1
```

Use sufficiently wide intermediate arithmetic for the squared products.

## 3. Fetch exactly one independent frame per selected block

Each exact block-index entry supplies:

- `blockNumber`;
- `storedByteLength`;
- absolute `blockPosition`.

The reader should use those fields directly. It should not reconstruct positions from cumulative
lengths or read adjacent blocks merely to locate the requested block.

For local files, read the indexed byte interval. For remote files, issue an exact HTTP range request.
When a broad query needs multiple physically adjacent entries, a reader may coalesce adjacent byte
ranges to reduce system-call or HTTP overhead, but each logical block remains an independent
decompression frame. Range coalescing must not cause unrelated blocks to be decompressed.

Useful per-query diagnostics include:

```text
candidate block ranges
stored blocks found
compressed bytes fetched
blocks decompressed
contacts decoded
contacts emitted
```

For a correct materialized query, `unrequested blocks decompressed` should be zero.

## 4. Reuse the Zstandard decompression context

Create one `ZSTD_DCtx` per reader instance, thread, or equivalent synchronization domain and reuse
it for matrix blocks and vector chunks. Use `ZSTD_decompressDCtx` or the corresponding streaming API.

Do not use an API pattern that allocates and initializes a new decompression context for every
small block. Context setup becomes visible for queries that touch many small independent frames.

The reader must still validate:

- the Zstandard frame magic;
- that dictionaries are not required when forbidden by the format;
- that the indexed record contains exactly one frame;
- the declared decompressed length;
- decompression errors and checksum failures.

Context reuse changes only the implementation cost, not validation behavior.

## 5. Keep successful binary reads allocation-free

The V10 decoder performs many bounds checks and ULEB128 byte reads. These checks must remain, but
their successful path should consist only of comparisons and pointer/index advancement.

A subtle but severe anti-pattern is a validation helper that accepts only `const std::string&`:

```cpp
void require(bool condition, const std::string& message);
```

Calling that function with a string literal may construct a temporary `std::string` before every
successful check. In a ULEB128 loop this can happen once or several times per contact value and
position. This was a major source of V10 decoder overhead.

Provide a literal overload whose successful path constructs nothing:

```cpp
inline void require(bool condition, const char* message) {
    if (!condition)
        throw FormatError(message);
}
```

A separate `std::string` overload can handle dynamically constructed messages. Construct the final
exception text only on failure.

Similarly, specialize the common fixed-width operations:

- `byte()` performs one bounds check and one indexed load;
- `u32()` performs one bounds check and four little-endian loads;
- `u64()` performs one bounds check and eight little-endian loads;
- `uleb128()` advances directly over bytes without routing each byte through a generic N-byte loop.

The implementation must preserve overflow, truncation, canonical-ULEB128, and trailing-byte checks.

## 6. Decode block streams without full intermediate arrays

The V10 block payload deliberately separates positions and values, but the reader does not need to
materialize both streams into full arrays before producing contacts.

Avoid this multi-pass design:

```text
decode every occupied position into vector A
decode every value slot into vector B
iterate A and B to create records
```

It adds allocations, memory writes, cache pressure, and another complete pass through each block.

### 6.1 `SPARSE_DELTA`

Decode the next position delta, calculate the absolute local position, obtain the value for the same
ordinal, validate it, and emit immediately. Retain only the previous position.

### 6.2 `BITMAP`

Walk the bitmap in local-position order. For each set bit, obtain the value for the next occupied
ordinal and emit immediately. Count set bits while walking and validate that the population equals
`occupiedCellCount`.

Efficient implementations may process a byte or machine word at a time and use bit-scan/popcount
operations, but a single allocation-free pass is already preferable to building an occupied-position
array.

### 6.3 `DENSE`

Decode each value slot in local-position order. For count blocks, zero means absent. For score blocks,
use the presence bitmap and require the canonical positive-zero representation for absent values.
Emit only present cells.

### 6.4 Value modes

- `ALL_DEFAULT`: decode the default once and reuse it for every occupied ordinal.
- `DIRECT`: decode one value when its position is visited.
- `DEFAULT_EXCEPTIONS`: decode and retain only the ordered exception ordinals, then consume exception
  values as those ordinals are reached. Do not allocate a value array sized to every block slot.

The exception-ordinal list is the only normally necessary temporary vector in the streaming decoder.
Its allocation must remain bounded by the containing record and the reader's resource limit.

## 7. Avoid repeated type-erased callbacks inside the hot loop

The public API may expose `std::function` or another type-erased callback, but internal layers should
use templates, function objects, or another inlineable mechanism where practical.

A contact should not pass through separate type-erased callbacks for:

```text
block decoder
→ materialized-resolution reader
→ coordinate filter
→ raw API
→ public API
```

Template the internal block, materialized, and raw emitters so the compiler can inline filters and
statistics accounting. Leave only the unavoidable public callback boundary type-erased.

This is especially useful for coarse, dense queries that emit hundreds of thousands of contacts.

## 8. Preserve validation without putting unnecessary work on every byte

Performance is not a reason to weaken malformed-file handling. Readers must retain:

- containing-record bounds checks;
- allocation limits before allocation;
- checked integer addition and multiplication;
- canonical ULEB128 validation;
- bitmap length, population, and padding validation;
- ordered and unique sparse positions;
- ordered and unique exception ordinals;
- value-slot and occupied-cell counts;
- chromosome coordinate bounds;
- cis upper-triangle and block-geometry validation;
- stored-record/index block-number agreement;
- exact decompressed-length and trailing-byte checks;
- full-matrix sum and occupied-count validation when the complete matrix is read.

Optimize the successful path rather than removing the checks. Examples include literal error-message
overloads, one bounds check per fixed-width scalar, streaming validation, and validating invariant
block/header properties once rather than reconstructing diagnostic objects per contact.

## 9. Derived-resolution query path

Derived resolutions aggregate exact source contacts and therefore have unavoidable work not present
in a materialized query. They should still use the same optimized source-block selection and streaming
decoder.

Recommended practices include:

- point directly to the declared finer materialized source;
- multiply the target rectangle by the exact aggregation factor with checked arithmetic;
- fetch only source blocks intersecting that expanded rectangle;
- filter source contacts before inserting them into an accumulator;
- use packed integer coordinate keys or an efficient pair hash;
- reserve hash-table capacity when a safe estimate is available;
- accumulate integer counts directly with checked `u64` addition;
- preserve deterministic source order for floating-point score aggregation when exact behavior
  requires it;
- cache decompressed source blocks or raw derived tiles in interactive applications with repeated
  overlapping queries.

Materialized and derived timings must be reported separately. Readers may optionally cache derived
tiles, but the on-disk source data remains authoritative.

## 10. Vector access

Normalization, expected, and normalized-expected vectors are independently chunked and indexed.
To preserve their random-access advantage:

- locate only chunks intersecting the requested bin or distance range;
- cache the parsed vector index;
- reuse the Zstandard context;
- reverse `RAW`, `BYTE_SHUFFLE`, or `XOR32` directly into the destination range;
- avoid decoding a chromosome-scale vector for a small viewport;
- preserve every original `float32` bit pattern before conversion to the public numeric type.

## 11. Reader lifetime and concurrency

For repeated queries, use a persistent file/reader object. It should retain header, footer, matrix,
block-index, and vector-index caches.

A reader containing a file stream and reusable `ZSTD_DCtx` should either:

- be confined to one thread;
- protect mutable state with appropriate synchronization; or
- provide one lightweight query context per thread while sharing immutable parsed indexes.

Do not add a global lock around all queries merely to reuse one decoder. Per-thread decoder contexts
are normally inexpensive relative to serialized query execution.

## 12. Benchmarking and acceptance criteria

Benchmark the reader with identical logical queries against a V9 file and its V10 conversion. At
minimum, stratify by:

- materialized versus derived resolution;
- fine, medium, and coarse bin sizes;
- near-diagonal cis;
- middle-distance cis;
- far-cis;
- inter-chromosomal;
- small and large query windows;
- cold and warm index/block caches;
- local files and HTTP range access.

Record p50, p95, p99, bytes fetched, selected blocks, decompressed blocks, decoded contacts, and
emitted contacts. Alternate file order during a comparison to reduce page-cache bias.

For materialized resolutions, the target is:

```text
p50 and p95 no more than 5–10% slower than V9,
with no systematic resolution or distance stratum exceeding that range.
```

Faster V10 results are expected when smaller frames reduce I/O and the hot decoder is implemented as
described above.

During the Straw C++ optimization that motivated this guide, the same 25-per-stratum benchmark over
5,100 timed reads per file changed from broad V10 slowdowns to:

| Metric | V9 | Optimized V10 | Result |
|---|---:|---:|---:|
| Total materialized and derived region-read time | 19.651 s | 11.347 s | V10 1.73× faster |
| p50 latency | 0.829 ms | 0.358 ms | V10 2.31× faster |
| p95 latency | 18.345 ms | 10.607 ms | V10 1.73× faster |

The benchmark reported identical query results. Representative far-cis medians changed to:

| Resolution | V9 | Optimized V10 |
|---:|---:|---:|
| 10 bp | 0.354 ms | 0.159 ms |
| 100 bp | 0.432 ms | 0.213 ms |
| 1 kb | 0.626 ms | 0.259 ms |
| 5 kb | 1.679 ms | 0.932 ms |

These figures are implementation evidence, not normative performance guarantees. Hardware, cache
state, file layout, storage, and query mix affect absolute results.

## 13. Common regression checklist

When a V10 reader is unexpectedly slower, check these in order:

1. Does a far-cis query begin at its nearest distance band, or incorrectly at depth zero?
2. Are exact block numbers located by binary/keyed lookup, or is the complete index scanned?
3. How many blocks and compressed bytes are fetched versus the mathematically required set?
4. Are matrix descriptors and block indexes cached across windows?
5. Is the Zstandard decompression context reused?
6. Does each successful validation check allocate or construct a string?
7. Does ULEB128 decoding call a generic scalar reader for every byte?
8. Are full position and decoded-value vectors allocated for every block?
9. How many type-erased callback boundaries are crossed per emitted contact?
10. Are derived and materialized timings being combined in a way that hides their different costs?
11. Are local and remote tests measuring the same cache state and query semantics?
12. Do correctness and corruption tests still pass after optimization?

## 14. Required verification after optimization

Every optimized reader should pass:

- independent byte-level fixtures for every block representation and value mode;
- malformed/truncated record tests;
- ULEB128 overflow and non-canonical encoding tests;
- bitmap population and padding tests;
- block geometry and chromosome-bound tests;
- exact `u64` count tests, including values above `2^53`;
- derived-resolution aggregation and overflow tests;
- vector transform reconstruction tests;
- local and HTTP exact-range tests;
- V9-to-V10 region comparisons across all query strata;
- memory and undefined-behavior sanitizer runs where supported.

Performance optimizations are complete only when they preserve identical results and the full defensive
reading behavior required by the V10 specification.
