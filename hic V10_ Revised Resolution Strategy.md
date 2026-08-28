# `.hic` V10: Revised Resolution Strategy

The V10 proposal is unchanged except for the resolution pyramid policy below.

## Resolution Pyramid

V10 distinguishes between:

- **Materialized resolutions** — raw matrix blocks are physically stored.
- **Derived resolutions** — raw matrices are reconstructed exactly by summing a finer materialized resolution.

Normalization vectors and expected-value vectors remain materialized for **every advertised resolution**, including derived resolutions.

### Recommended V10 resolution set

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
| 500 kb | Derived | 100 kb |
| 1 Mb | **Materialized** | — |
| 2.5 Mb | **Materialized** | — |

### Resolutions intentionally omitted

V10 does **not** introduce 20 kb or 200 kb as standard derived resolutions.

Although they could be produced exactly as 2× aggregations of 10 kb and 100 kb respectively, they provide relatively little practical value compared with the additional resolution metadata, normalization vectors, expected vectors, reader complexity, and testing surface they would introduce.

The standard pyramid therefore favors resolutions that are either:

1. established `.hic` ecosystem resolutions;
2. frequently used analysis/display resolutions; or
3. especially valuable at sub-kilobase scales, where avoiding duplicated matrices provides substantial storage savings.

### Why 25 kb, 250 kb, and 2.5 Mb remain materialized

These resolutions are retained as physical matrices primarily for backward compatibility with the established `.hic` resolution hierarchy.

They are also relatively inexpensive to store because the matrices are already highly aggregated.

There is consequently little reason to introduce derivation overhead at these coarse resolutions merely to save a comparatively small amount of space.

The coarse-resolution part of the hierarchy therefore remains deliberately conservative:

```text
5 kb       MATERIALIZED
10 kb      MATERIALIZED
25 kb      MATERIALIZED
50 kb      MATERIALIZED
100 kb     MATERIALIZED
250 kb     MATERIALIZED
500 kb     DERIVED from 100 kb
1 Mb       MATERIALIZED
2.5 Mb     MATERIALIZED
```

The main use of virtual resolutions remains concentrated where duplicate matrix storage is expensive:

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

while particularly important working resolutions remain directly stored:

```text
5 kb
50 kb
```

### Final recommended materialized set

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
1 Mb
2.5 Mb
```

### Final recommended derived set

```text
20 bp    <- 10 bp
50 bp    <- 10 bp

200 bp   <- 100 bp
500 bp   <- 100 bp

2 kb     <- 1 kb

500 kb   <- 100 kb
```

The format MAY support other resolutions, either materialized or derived, but these constitute the recommended standard V10 hierarchy.

All other aspects of the V10 strategy remain unchanged:

- exact aggregation of **raw** contacts before normalization;
- target-resolution normalization and expected vectors stored explicitly;
- Zstandard compression;
- bounded multi-block compression pages;
- sparse-delta / bitmap / dense adaptive block encoding;
- separate position and value streams;
- implicit/default count values with exception streams;
- explicit unsigned integer count storage;
- losslessly compressed normalization and expected-vector chunks;
- compact page-oriented indexing;
- exact V9 → V10 repacking and verification.