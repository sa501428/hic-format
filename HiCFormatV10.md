# hic file format

## Structure

* Header
* Body
    * Matrix
    * Block
* Footer
    * Master index
    * Expected value vectors


## V10 changes from V9

V10 introduces three changes to the block encoding, all confined to the **Blocks** section of the Body.
No other section of the format is changed.

1. **Block compression algorithm** — blocks may now be compressed with [Zstandard (zstd)](https://github.com/facebook/zstd)
   in addition to zlib.  A new `compressionType` byte in the block header selects the algorithm.  The default value
   `0` means zlib, preserving backward compatibility with v9 readers that encounter a v10 file written with zlib.
   Value `1` means zstd.

2. **All-counts-one flag** — a new `allCountsOne` byte in the block header signals that every contact value in the
   block is exactly `1`.  When set, the `value` field is omitted from every contact record, halving per-record
   storage for sparse, far-from-diagonal blocks.  This flag is only meaningful for raw observed (non-normalized)
   data; writers MUST leave it `0` for normalized data.

3. **Delta-encoded column positions** — a new `useDeltaColumn` byte in the block header signals that `binColumn`
   values in the list-of-rows format are stored as successive deltas rather than absolute offsets.  Within each
   row the first `binColumn` value is absolute (relative to `binColumnOffset` as usual); each subsequent value is
   the difference from the preceding one.  Because columns within a row are sorted in ascending order, all deltas
   are non-negative.  This flag has no effect when `matrixRepresentation == 2` (dense).



## Header

|Field | Description |	Type | Value | V10 change |
|------|------------|------|-------|------|
|Magic|HiC magic string|String|HIC||
|Version|Version number|int|10||
|footerPosition|File position of footer|long|||
|genomeId| Genome identifier (e.g. hg19, mm9, etc)|String|||
|normVectorIndexPosition|  File position for normalization vector index|long|||
|normVectorIndexLength|  Length to read for normalization vector index|long|||

#### Attribute list
*List of key-value pair attributes.  See notes on common attributes below.*

|Field | Description |	Type | Value | V10 change |
|------|------------|------|-------|------|
|nAttributes	|Number of key-value pair attributes|	int|||
||
|*Repeat for each attribute (n = nAttributes)*|||
|key	|Attribute key|	String	|||
|value|Attribute value|		String|||

#### Chromosome list

*List of chromosome name and lengths*

|Field | Description |	Type | Value | V10 change |
|------|------------|------|-------|------|
|nChrs|	Number of chromosomes|int|||
||
|*Repeat for each chromosome (n = nChrs)*|
|chrName	|Chromosome name	|String|||
|chrLength|	Chromosome length |	long	|||

#### Base-pair resolution list

*List of base-pair resolutions*

|Field | Description |	Type | Value | V10 change |
|------|------------|------|-------|------|
|nBpResolutions	|Number of base pair resolutions|	int|||
||
|*Repeat for each resolution (n = nBpResolutions)*|||
|resBP	|Bin size in base pairs	|int|||

#### Fragment resolution list

*List of bin sizes for fragment resolution levels*

|Field | Description |	Type | Value | V10 change |
|------|------------|------|-------|------|
|nFragResolutions	|Number of fragment resolutions	|int|||
||
|*Repeat for each resolution (n = nFragResolutions)*|
|resFrag	|Bin size in fragment units (1, 2, 5, etc)|	int|||

#### Fragment site positions list

*List of fragment site positions per chromosome, in same order as chromosome list above (n = nChrs).  This section absent if nFragResolutions = 0.*

|Field | Description |	Type | Value | V10 change |
|------|------------|------|-------|------|
|nSites|	Number of sites for this chromosome|	int|||
||
|*Repeat for each site (n = nSites)*|
|sitePosition|	Site position in base pairs|	int|||



## Body

The **Header** section is followed immediately by the **Body**, which contains the contact map data for each
chromosome-chromosome pairing and each resolution.



### Matrix metadata

This section contains metadata for the contact matrices.  It is repeated for each chromosome-chromosome pair.
The master index contains an entry for each combination and is used to randomly access a specific
matrix as needed.  The metadata in this section includes an index for data blocks which contain the actual
contact data.


|Field	|Description|	Type|	Value| V10 Change |
|------|------------|------|-------|--------|
|chr1Idx| Index for chromosome 1.  This is the index into the array of chromosomes defined in the header above.  The first chromosome has index **0**.|	int|||
|chr2Idx| Index for chromosome 2. |	int	|||
|nResolutions	|Total number of resolutions for this chromosome-chromosome pair, including base pair and fragment resolutions.	|int|||

#### Resolution (zoom level) metadata

*The section below is repeated for each resolution (n = nResolutions)*

|Field	|Description|	Type|	Value| V10 Change |
|------|------------|------|-------|--------|
|unit|	Distance unit, base-pairs or fragments	|String	|BP or FRAG||
|resIdx	|Index number for this resolution level, an Array index into the bin size list of the header, first element is **0**. |	int|||
|sumCounts|	Sum of all counts (or scores) across all bins at current resolution.|	float|||
|occupiedCellCount|	Total count of cells that are occupied.  **Not currently used**|float|0||
|stdDev|	Standard deviation of counts among occupied bins. **Not currently used**|float|0||
|percent95|	Estimate of 95th percentile of counts among occupied bins. **Not currently used**|float|0||
|binSize|	The bin size in base-pairs or fragments	|int|||
|blockSize			|Dimension of each block in bins.  In v9 and v10 interchromosomal blocks are square, so the total number of bins is ```blockSize^2```. But intrachromosomal blocks are rotated and not necessarily square. In this case, blockSize specifies the dimension of the block along the diagonal axis.  See description of grid structure below|int|||
|blockColumnCount|The number of columns in the grid of blocks. For intrachromosomal block structure, this specifies the number of columns in the grid of blocks along the diagonal. |int|||
|blockCount|The number of blocks.  Note empty blocks are not stored.|int|||
||
|*repeat for each block (n = blockCount)*|
|blockNumber	|Numeric id for block.  This is the linear position of the block in the grid when counted in row-major order.   ```blockNumber = row * blockColumnCount + column``` where first row and column are **0**. **IMPORTANT: block index entries must be ordered by blockNumber**	|int||
|blockPosition|	File position of the start of the block |	long||
|blockSizeBytes	|Size of block in bytes| int||

***End of Matrix metadata section***



### Blocks

A block represents a square sub-matrix of a contact map.

***Note: Blocks are individually compressed. The compression algorithm is indicated by the `compressionType` field
in the block header below. Supported algorithms are zlib (type 0) and zstd (type 1).***

|Field	|Description|	Type|	Value| V10 Change|
|------|------------|------|-------|---------|
|nRecords	|Number of contact records in this block|	int	||
|binColumnOffset | Column offset for the contact records in this block.  The `binColumn` value below is relative to this offset.| int ||
|binRowOffset | Row offset for the contact records in this block.  The `rowNumber` value below is relative to this offset.| int ||
|useFloatContact | Flag indicating the `value` field in contact records for this block are recorded with data type `float`.  If == 1 a `float` is used, otherwise type is `short`| byte ||
|useIntXPos | Flag indicating the `recordCount` and `binColumn` fields in contact records for this block are recorded with data type `int`. If == 1 an `int` is used, otherwise type is `short` | byte ||
|useIntYPos | Flag indicating the `rowCount` and `rowNumber` fields in contact records for this block are recorded with data type `int`. If == 1 an `int` is used, otherwise type is `short` | byte ||
|matrixRepresentation | Representation of matrix used for the contact records.  If == 1 the representation is a `list of rows`, if == 2 `dense`. | byte ||
|compressionType | Compression algorithm used for this block. `0` = zlib, `1` = zstd. | byte || (ADDED FROM v9)|
|allCountsOne | If == 1, all contact values in this block are exactly `1` and the `value` field is omitted from every contact record.  The implicit value is `1` (as float: `1.0f`).  Writers MUST set this to `0` for normalized data. | byte || (ADDED FROM v9)|
|useDeltaColumn | If == 1, `binColumn` values in list-of-rows format are delta-encoded within each row (first value is absolute, subsequent values are differences from the preceding value).  Has no effect when `matrixRepresentation == 2`. | byte || (ADDED FROM v9)|
|blockData| The block matrix data.  See descriptions below.| | |

##### Block data - list of rows

|Field	|Description|	Type|	Value| V10 Change|
|------|------------|------|-------|--------|
|rowCount | Number of rows. The data type is determined by the `useIntYPos` flag above. | int : short ||
||
|*repeat for each row (n = rowCount)*|
|rowNumber | Matrix row number, relative to `binRowOffset`. First row is `0`. The data type is determined by the `useIntYPos` flag above. | int : short ||
|recordCount | Number of records for this row. Row is sparse, zeroes are not recorded. The data type is determined by the `useIntXPos` flag above. | int : short ||
||
|*repeat for each contact record (n = recordCount)*|
|binColumn	|Column index relative to `binColumnOffset`. If `useDeltaColumn == 1`, this is a delta from the previous `binColumn` in this row (first record is absolute). The data type is determined by the `useIntXPos` flag above. |	int : short||
|value	|Value (counts or score). The data type is determined by the `useFloatContact` flag above.  **Omitted for every record in this block if `allCountsOne == 1`.**|	float : short|| (CHANGED FROM v9)|
||
|*End of loop through contact records (n = recordCount)*|
||
|*End of loop through rows (n = rowCount)*|


##### Block data - dense

|Field	|Description|	Type|	Value|
|------|------------|------|-------|
|nRecords | Number of contact records in this block.  | int ||
|w | Width of the dense block.  This can be < the blockSize if the edge columns on either side are zeroes.  See discussion on block representation below | short ||
||
|*repeat for each contact record (n = nRecords)*||
|value	|Value (counts or score). The data type is determined by the `useFloatContact` flag above.  ***Note:  no value is flagged by the value -32768 if data type is short, NaN if data type is float***|	float : short||

### Footer

| Field |	Description|	Type |	Value | V10 change |
|------|------------|------|-------|-------|
|nBytesV5|	Number of bytes for the "version 5" footer, that is everything up to the normalized expected vectors	|long||

#### Master index

| Field |	Description|	Type |	Value |
|------|------------|------|-------|
|nEntries|	Number of index entries|	int||
||
||*List of index entries (n = nEntries)*||
|key|	A key constructed from the indices of the two chromosomes for this matrix.  The indices are defined by the list of chromosomes in the header section with the first chromosome occupying index **0**|String||
|position	|Position of the start of the chromosome-chromosome matrix record in bytes	|long||
|size	|Size of the chromosome-chromosome matrix record in bytes.  This does not include the **Block** data.| int||



#### Expected value vectors

| Field |	Description|	Type |	Value | V10 Change|
|------|------------|------|-------|--------|
|nExpectedValueVectors|	Number of expected value vectors to follow.  These are expected values from the non-normalized observed matrix.| int||
||
|*List of expected value vectors (n = nExpectedValueVectors)*||
|unit|	Bin units either FRAG or BP.	|String	|FRAG : BP||
|binSize	|Bin (grid) size for this calculation	|int|||
|nValues	|Size of the vector|	long||
||
||*List of expected values (n = nValues)*|
|value	|Expected value|	float||
|nChrScaleFactors| Number of chromosome normalization factors| int|||
||
||*List of normalization factors (n = nChrScaleFactors)*|||
|chrIndex|	Chromosome index|	int|||
|chrScaleFactor|	Chromosome scale factor	|float||



#### Normalized expected value vectors

| Field |	Description|	Type |	Value | V10 Change|
|------|------------|------|-------|---------|
|nNormExpectedValueVectors|	Number of normalized expected value vectors to follow	|int|||
||
|*List of normalized vectors (n = nNormExpectedValueVectors)*||
|type|	Indicates type of normalization	|String|	VC:KR:INTER_KR:INTER_VC:GW_KR:GW_VC||
|unit	|Bin units either FRAG or BP.	|String|	FRAG : BP||
|binSize|	Bin (grid) size for this calculation	|int|||
|nValues|	Size of the vector	|long	||
||
|*List of expected values (n = nValues)*|
|value	|Expected value	|float||
|nChrScaleFactors|Number of normalization factors for this vector| int|||
||
|*List of normalization factors (n = nChrScaleFactors)*|
|chrIndex|	Chromosome index	|int	|||
|chrScaleFactor|	Chromosome scale factor	|float||



#### Normalization vector index

| Field |	Description|	Type |	Value | V10 change |
|------|------------|------|-------|---------|
|nNormVectors|	Number of normalization vectors |	int|||
||
|*Repeat for each norm vector (n = nNormalizationVectors)*|
|type	|Indicates type of normalization	|String|	VC:KR:INTER_KR:INTER_VC:GW_KR:GW_VC||
|chrIdx|	Chromosome index	|int|	||
|unit|	Bin units either FRAG or BP.|	String|	FRAG : BP||
|binSize	|Resolution 	|int|||
|position|	File position of value array, described below|	long	|||
|nBytes|	Size in bytes of value array	| long	||

#### Normalization vector arrays, 1 per normalization vector.

| Field |	Description|	Type |	Value | V10 change |
|------|------------|------|-------|---------|
|nValues|	Number of values in array|	long||
||
|*Normalization vector values (n = nValues)*|
| value | Norm value | float ||



#### Notes

##### Data types

* Strings are null (0) terminated.  So for example the string "HIC" is represented by 4 bytes [48 49 43 0]
* Other data types are Java
    * short - 16 bit integer
    * int - 32 bit integer
    * long  -  64 bit integer
    * float - 32 bit floating point
    * double - 64 bit floating point

##### Attributes

The attributes table in the header can contain an arbitrary number of key-value string pairs.  The **Juicer** tool
inserts one or more of the following attributes.
* "statistics":
* "graphs":
* "software":
* "nviIndex":  reserved for future use
* "nviLength":  reserved for future use

##### Writer guidance for ultra-high-resolution data

At resolutions of 500 BP and below, the following practices are strongly recommended to minimise file size.

* Use `compressionType = 1` (zstd) for all blocks.  Zstd consistently achieves 30–50% better compression ratios
  than zlib on sparse integer block data, with faster decompression.

* Set `allCountsOne = 1` for any block in which every contact value is exactly `1`.  At ultra-high resolution,
  blocks at depth ≥ 2 (far from the diagonal) are typically composed almost entirely of singleton contacts.
  Setting this flag for those blocks eliminates the value field entirely, reducing per-record storage by 2–4 bytes
  before compression.  Near-diagonal blocks (depth 0–1) accumulate repeated contacts and will generally not
  qualify; writers simply inspect the block before writing and set the flag accordingly.

* Set `useDeltaColumn = 1` when contacts within a row are spread across a large column range.  Delta encoding
  replaces absolute column offsets with small positive differences, which reduces the magnitude of values fed to
  the compressor and improves compression ratios further.

* Use a larger `blockBinCount` (e.g. 2000–5000) at resolutions ≤ 500 BP.  This reduces the total number of blocks
  and amortizes per-block header and index overhead across more contact records.

#### Grid structure

Each chr-chr matrix at a given resolution is subdivided into a grid structure of square **blocks**.
Each block consists of NxN bins, where N is referred to as **blockSize**.  In older versions of the spec,
and in code, this parameter is referred to as **blockBinCount**.

For intra chromosome matrices (chr1 == chr2) only the lower diagonal is stored (row >= column).  The upper diagonal
can be inferred upon reading by transposition.


#### Intrachromosomal Block matrix representation

For intrachromosomal matrices, blocks are stored in a rotated manner, with the axes defined along the diagonal and perpendicular to the diagonal. A visual example of this is included at [https://bcm.box.com/v/hic-file-version-9](https://bcm.box.com/v/hic-file-version-9)

Furthermore, the block size increases by a factor of 2 along the anti-diagonal axis, as the number of contacts also decrease further from the diagonal. This allows for a natural and dynamic block size to decrease overall file size.

The spatial unit for a block is still a ```bin```, which can be computed from a genomic position with the formula

```bin = floor(position / binSize)```.

The origin of a block is then

```binX = floor(x / binSize), binY = floor(y / binSize)```

where x and y are genomic positions in either base pairs or fragment number.

To identify the block number data is stored in, we calculate

```
position_along_diagonal = (binX + binY) / 2 / blockBinCount;
position_along_anti_diagonal = log2(1 + Math.abs(binX - binY) / Math.sqrt(2) / blockBinCount);
block_number = position_along_anti_diagonal * blockColumnCount + positionAlongDiagonal
```

Because the 2D heatmap viewers are often at a 45 degree rotation from the representation of the block, it is necessary to identify all the blocks that overlap this region. For a rectangular region spanning binX1 to binX2 and binY1 to binY2, the rotation along the diagonal and antidiagonal correspond to:

```
// pad = position along diagonal
padMin = (binX1 + binY1) / 2 / blockBinCount;
padMax = (binX2 + binY2) / 2 / blockBinCount + 1;

// anti = position along anti diagonal
// UR = upper right corner, LL = lower left corner
antiUR = log2(1 + Math.abs(binX1 - binY2) / Math.sqrt(2) / blockBinCount);
antiLL = log2(1 + Math.abs(binX2 - binY1) / Math.sqrt(2) / blockBinCount);
```
We determine the appropriate boundaries for the anti-diagonal axis.
```
antiMin = Math.min(antiLL, antiUR);
antiMax = Math.max(antiLL, antiUR) + 1;
```
If the diagonal is contained in the viewer, we explicitly set the lower bound along the anti-diagonal axis.
```
if ((binX1 > binY2 && binX2 < binY1) || (binX2 > binY1 && binX1 < binY2)) {
   antiMin = 0;
}
```
This calculates a permissive region for the viewer to ensure all data is captured for the region,
resulting in block numbers defined by the (inclusive) boundaries:
```
for each p in [padMin, padMax]
   for each a in [antiMin, antiMax]
      block_number = a * blockColumnCount + p
```

#### Block matrix representation

The spatial unit for a block is a ```bin```, which can be computed from a genomic position with the formula

```bin = floor(position / binSize)```.

The origin of a block is then

```floor(x / binSize), floor(y / binSize)```

where x and y are genomic positions in either base pairs or fragment number, depending on the unit (BP or FRAG).

* List of rows

The list of rows is a sparse matrix format.  Each row is represented as follows

```rowNumber rowSize [binX1 value1, binX2 value2, ...]```

The first row in the matrix has ```rowNumber = 0```.  The highest row number possible is ```blockSize - 1```.

When `useDeltaColumn == 1`, the binX values above are replaced by deltas:

```rowNumber rowSize [binX1 (delta2) (delta3) ...]```

where `delta_i = binX_i - binX_{i-1}` for i > 1, and all deltas are non-negative since columns are sorted.

* Dense

In dense matrix format all values including zero are output in row major order.  Allowance is made however for the
possibility that only a sub-matrix of the block is populated, specifically that leading or trailing columns of
the block might have no contacts (value = 0).   To account for this possibility the maximum column number within the block
which has at least 1 non-zero value is determined, which we will call ```binXMax```.   The width of the block can
then be determined and used to obtain the x and y coordinates in bin units for each value as follows.

         w = (binXMax - binXOffset + 1);
         row = floor(i / w);
         col = i - row * w;
         binX = binXOffset + col;
         binY = binYOffset + row;
