# V10 Independent-Block Migration Plan

## Purpose and scope

The V10 specification no longer groups matrix blocks into compression pages. Every non-empty logical matrix block is now an independently compressed, independently indexed random-access record.

This document describes the implementation work required outside `hic-format`. It is a handoff plan only: the corresponding changes to `hictools-c`, `straw`, build scripts, converters, and downstream readers have not yet been made.

Vector chunks are unchanged. The removal applies only to matrix pages.

## Wire-format changes

### Matrix resolution descriptor

The descriptor remains 76 bytes. Bytes 52 through 75 change as follows:

| Offset | Previous draft field | Revised field |
|---:|---|---|
| 52 | `pageIndexPosition u64` | `blockIndexPosition u64` |
| 60 | `pageIndexLength u64` | `blockIndexLength u64` |
| 68 | `pageCount u32` | `logicalBlockCount u32` |
| 72 | `logicalBlockCount u32` | `reserved1 u32`, required to be zero |

A non-empty materialized resolution has a block index and a positive logical-block count. A derived or empty resolution stores zero in all four fields.

### Exact block index

The new index retains magic `H10I` but uses `indexVersion = 2`. Its size is exactly:

```text
24 + 16 * blockCount
```

Header:

```text
indexMagic       4 bytes  "H10I"
indexVersion     u32      2
indexByteLength  u64
blockCount       u32
reserved         u32      0
```

Each sorted entry is:

```text
blockNumber       u32
storedByteLength  u32
blockPosition     u64
```

The absolute locator is intentional. The index cost is small compared with matrix data and restores direct binary-search lookup with no cumulative-position reconstruction.

### Independently stored block

Every indexed block points to one `H10B` record:

```text
recordMagic        4 bytes  "H10B"
codec              u8       1 (Zstandard)
recordVersion      u8       1
recordFlags        u16      0
uncompressedBytes  u32
blockNumber        u32
compressedPayload  one complete Zstandard frame
```

The frame decompresses to exactly one existing 40-byte logical-block header plus its position and value streams. There is no page header, internal block directory, or concatenation of logical blocks.

### Compatibility rule

Page-based draft V10 files use `H10I` version 1. Revised readers must reject index version 1. They must not guess which layout is present.

All existing draft V10 files, including benchmark fixtures, must be rebuilt or explicitly converted. V9 files are unaffected.

## `hictools-c` changes

Repository: `aidenlab/hictools-c`, V10 implementation under `v10/`.

### Shared model and command-line options

In `v10/format.h`:

- Remove `Options::pageBytes`.
- Remove page-specific structs, constants, and comments.
- Add a block-index entry model containing block number, absolute position, and stored length.
- Keep the existing logical block geometry and block payload encoding unchanged.
- Keep the compression-level option because each block still uses Zstandard.

In `v10/main.cpp` and `v10/README.md`:

- Remove `--page-bytes` from usage, parsing, validation, and documentation.
- Document that every logical block is independently compressed.
- Document that broad readers may coalesce adjacent block byte ranges without changing compression boundaries.

### V10 writer and direct preprocessor

In `v10/writer.cpp`:

1. Keep `block(...)` as the encoder for one logical block payload.
2. Delete `EncodedPage`, `encode_page(...)`, pending-page assembly, page flush thresholds, and page compression jobs.
3. For each non-empty logical block:
   - produce its encoded logical-block payload;
   - Zstandard-compress that payload independently;
   - prepend the 16-byte `H10B` header;
   - write the record;
   - retain `{blockNumber, storedByteLength, blockPosition}`.
4. Emit one `H10I` version-2 block index per non-empty materialized matrix resolution.
5. Patch the resolution descriptor with block-index position, block-index length, logical-block count, and zero `reserved1`.
6. Preserve increasing block-number physical order when practical so large reads can be coalesced.
7. Continue parallel compression, but schedule independent blocks rather than assembled pages. Preserve output order when consuming worker results.

The direct V10 preprocessing path in `v10/pre.cpp` should require no change to contact aggregation or logical-block assignment. Its output path must use the revised writer.

### V9-to-V10 conversion and repacking

In `v10/convert.cpp` and `v10/repack.cpp`:

- Preserve existing V9 block decoding, raw-cell reconstruction, derived-resolution verification, and vector preservation.
- Feed each encoded V10 logical block directly to independent block compression.
- Build an exact version-2 block index instead of grouping encoded blocks.
- Reject page-based draft V10 input unless an explicit one-time draft migration command is implemented.
- If a draft migration command is implemented, it must decompress every old page, recover every logical block, and rewrite each block as its own `H10B` record. It must not copy old compressed page bytes.

### Addnorm and relocation

`v10/addnorm.cpp` and the relocation helpers in `v10/reader.cpp` rewrite file sections and therefore need special care:

- Rename descriptor relocation fields from page-index to block-index terminology.
- When matrix data moves, update every `blockPosition` stored inside every exact block index.
- Continue updating the descriptor's `blockIndexPosition` when the index itself moves.
- Validate every rewritten block interval and require `reserved1 == 0`.

### Internal V10 reader

In `v10/reader.cpp` and `v10/reader.h`:

- Parse the revised descriptor fields.
- Parse only `H10I` version 2 using the fixed 24-byte header and 16-byte entries.
- Validate strict block-number order, exact index length, in-bounds locators, non-overlapping block intervals, `H10B` header fields, matching block number, exact Zstandard frame length, and exact decompressed length.
- For full-matrix operations, visit index entries in order and independently decompress each block.
- Remove page directory parsing and page-count validation.

### Tests

Update at least:

- `v10/tests/inspect_v10.py`
- `v10/tests/test_writer.py`
- `v10/tests/test_geometry.cpp`
- conversion and addnorm integration tests

Tests should assert:

- `H10I` version 2;
- exact `24 + 16 * blockCount` index length;
- one `H10B` record per index entry;
- one logical block per Zstandard frame;
- no overlapping indexed intervals;
- strict block-number order;
- empty and derived descriptor storage fields are all zero;
- page-based version-1 indexes are rejected;
- direct preprocessing and V9 conversion remain bitwise or numerically lossless under the existing conformance rules.

## `straw` changes

Repository: `aidenlab/straw`, primary implementation under `C++/`.

### V10 binary reader

In `C++/straw_v10.cpp` and `C++/straw_v10.h`:

- Replace `Page` with an exact block-index entry type.
- Replace `pages(const Zoom&)` with block-index loading and binary-search lookup.
- Parse the revised resolution descriptor fields and require `reserved1 == 0`.
- Accept only `H10I` version 2.
- Remove page candidate testing, page directory decoding, and the loop that decodes every block in a selected page.
- Compute the candidate logical block-number set first.
- For every candidate number, binary-search the block index and skip it if absent.
- Range-read only the corresponding `H10B` record, validate its block number, independently decompress it, and decode the single logical block.
- For a derived resolution, expand the source rectangle, enumerate source block numbers, and independently fetch only existing source blocks before aggregation.

The cis candidate calculation must compute both nearest and farthest rotated distance bands. It must not enumerate every band from zero to the farthest band unless the query crosses the diagonal.

### Local and HTTP read optimization

Independent compression does not require one operating-system or HTTP request per block:

- Sort selected entries by `blockPosition`.
- Merge contiguous requested stored-record intervals into bounded range reads.
- Slice the combined response into the indexed `H10B` records.
- Decompress each requested block independently.

This optimization must never fetch a nonadjacent block merely because it has a nearby logical block number, and it must never require decompression of an unrequested block.

### Caching and API lifetime

- Cache parsed block indexes per `(chromosome pair, unit, resolution)` in a persistent `File`.
- Cache decompressed blocks by file identity, chromosome pair, resolution, and block number where appropriate.
- Continue allowing stateless calls, but avoid a separate format-probe open when the caller already holds a V10 `File`.

### Tests and fixtures

Update:

- `C++/tests/test_v10.py`, whose fixtures currently construct `H10P` pages and a version-1 page index;
- `C++/tests/v10_api_probe.cpp` as needed;
- local and loopback-HTTP range tests;
- compare/timing tests.

Add assertions that a one-block query:

- reads only the exact indexed block record;
- decompresses exactly one logical block;
- performs no unrelated block decoding;
- returns the same records locally and over HTTP;
- rejects old grouped-block index version 1.

Any language binding that contains an independent V10 decoder must make the same wire-format changes. Bindings that call the C++ core inherit the revised behavior after rebuilding.

## Build and benchmark scripts

Scripts that invoke `hic_v10 pre` or `hic_v10 convert` generally retain the same positional arguments and resolution options. Required updates are:

- remove any `--page-bytes` argument;
- rebuild the `hic_v10` executable before regenerating files;
- delete or archive old draft V10 outputs so they are not mistaken for revised V10 files;
- regenerate V10 files from the original HBS/pairs input or from V9 input;
- record the writer commit and block-index version in benchmark metadata.

For `/Users/muhammad/Desktop/hic/subsample_9_10/build_v9_v10.sh`, the current command does not pass `--page-bytes`, so its command line should remain valid after rebuilding `hic_v10`. Its existing `hct_1B.v10.hic` must be regenerated because it contains the old grouped-block layout.

## Recommended implementation order

1. Land this specification change.
2. Update the `hictools-c` shared format model, writer, and exact block index.
3. Update `hictools-c` inspection tests and generate a small revised fixture.
4. Update the straw V10 reader against that fixture, including local and HTTP tests.
5. Update `hictools-c` conversion, addnorm, repacking, and relocation paths.
6. Regenerate the large benchmark V10 file.
7. Run exhaustive V9-versus-V10 content comparison.
8. Run the stratified latency benchmark and file-size comparison.
9. Update remaining bindings and public documentation.

## Acceptance criteria

Correctness:

- Revised V10 and V9 raw matrices compare exactly at every advertised resolution under the existing conversion rules.
- Normalization and expected-value vectors preserve their existing exactness guarantees.
- Every non-empty materialized logical block has exactly one index entry and one `H10B` record.
- No matrix Zstandard frame contains more than one logical block.
- Old index version 1 fails with an explicit unsupported-draft error.

Performance:

- A materialized-resolution query decompresses no block outside its candidate block-number set.
- The benchmark reports requested blocks, existing blocks fetched, blocks decompressed, compressed bytes fetched, and records decoded.
- Warm materialized-resolution p50 and p95 latency are measured against V9 at every stored resolution.
- Near-, mid-, far-cis, and trans queries show no distance-dependent over-selection of cis bands.

Space:

- Report matrix block bytes, exact block-index bytes, vectors, other metadata, and total file size separately.
- Attribute any size increase specifically to independent block framing, loss of cross-block compression, and the exact locator index.
- Do not reintroduce grouped matrix compression solely to recover size without a new latency and read-amplification review.
