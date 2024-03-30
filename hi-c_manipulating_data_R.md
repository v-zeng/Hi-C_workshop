Manipulating Hi-C data in R
================
Vinson Zeng
March 29, 2024

The objective of this notebook is to develop working knowledge
of:<br> 1. Modifying information associated with an existing
HiCExperiment object<br> 2. Subsetting a HiCExperiment object<br> 3.
Coercing a HiCExperiment object in a base data structure

#### Load packages and objects

``` r
setwd("~/repos/Hi-C_workshop")
.libPaths("~/rlibrary")

library(ggplot2)
library(GenomicRanges)
```

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
    ##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    ##     Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
    ##     table, tapply, union, unique, unsplit, which.max, which.min

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:utils':
    ## 
    ##     findMatches

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## Loading required package: GenomeInfoDb

``` r
library(InteractionSet)
```

    ## Loading required package: SummarizedExperiment

    ## Loading required package: MatrixGenerics

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'MatrixGenerics'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
    ##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
    ##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
    ##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
    ##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
    ##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
    ##     colWeightedMeans, colWeightedMedians, colWeightedSds,
    ##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
    ##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
    ##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
    ##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
    ##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
    ##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
    ##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
    ##     rowWeightedSds, rowWeightedVars

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## 
    ## Attaching package: 'Biobase'

    ## The following object is masked from 'package:MatrixGenerics':
    ## 
    ##     rowMedians

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     anyMissing, rowMedians

``` r
library(HiCExperiment)
```

    ## Consider using the `HiContacts` package to perform advanced genomic operations 
    ## on `HiCExperiment` objects.
    ## 
    ## Read "Orchestrating Hi-C analysis with Bioconductor" online book to learn more:
    ## https://js2264.github.io/OHCA/

    ## 
    ## Attaching package: 'HiCExperiment'

    ## The following object is masked from 'package:SummarizedExperiment':
    ## 
    ##     metadata<-

    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     metadata<-

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     resolution

``` r
library(HiContactsData)
```

    ## Loading required package: ExperimentHub

    ## Loading required package: AnnotationHub

    ## Loading required package: BiocFileCache

    ## Loading required package: dbplyr

    ## 
    ## Attaching package: 'AnnotationHub'

    ## The following object is masked from 'package:Biobase':
    ## 
    ##     cache

#### Generate an example `hic` object

The HiCExperiment object will be created from an example `.cool` file in
the `HiContactsData` package.

``` r
library(HiCExperiment)
library(HiContactsData)

#=== download example `.mcool` file and cache it locally
coolf <- HiContactsData('yeast_wt', 'mcool')
```

    ## see ?HiContactsData and browseVignettes('HiContactsData') for documentation

    ## loading from cache

``` r
#=== create connection to disk-stored `.mcool` file
cf <- CoolFile(coolf)
cf
```

    ## CoolFile object
    ## .mcool file: /home/scatsac/.cache/R/ExperimentHub/283eb5e62932f_7752 
    ## resolution: 1000 
    ## pairs file: 
    ## metadata(0):

``` r
#=== import contacts from the long arm of chromosome `II`, at resolution `2000`
hic <- import(cf, focus = 'II:300001-813184', resolution = 2000)
hic
```

    ## `HiCExperiment` object with 306,212 contacts over 257 regions 
    ## -------
    ## fileName: "/home/scatsac/.cache/R/ExperimentHub/283eb5e62932f_7752" 
    ## focus: "II:300,001-813,184" 
    ## resolutions(5): 1000 2000 4000 8000 16000
    ## active resolution: 2000 
    ## interactions: 18513 
    ## scores(2): count balanced 
    ## topologicalFeatures: compartments(0) borders(0) loops(0) viewpoints(0) 
    ## pairsFile: N/A 
    ## metadata(0):

#### Subsetting a contact matrix

The two approaches to subset a Hi-C contact matrix include: <br> 1.
Subsetting before importing - leverages random access disk-stored
contact matrix to only import interactions overlapping with a genomic
locus of interest. 2. Subsetting after importing - parse entire contact
matrix in memory, subset interactions overlapping with a genomic locus
of interest.

#### Subset before import with `focus`:

``` r
import(cf, focus = 'II:300001-800000', resolution = 2000)
```

    ## `HiCExperiment` object with 301,018 contacts over 250 regions 
    ## -------
    ## fileName: "/home/scatsac/.cache/R/ExperimentHub/283eb5e62932f_7752" 
    ## focus: "II:300,001-800,000" 
    ## resolutions(5): 1000 2000 4000 8000 16000
    ## active resolution: 2000 
    ## interactions: 17974 
    ## scores(2): count balanced 
    ## topologicalFeatures: compartments(0) borders(0) loops(0) viewpoints(0) 
    ## pairsFile: N/A 
    ## metadata(0):

Subset to off-diagonal genomic location using pairs of coordinates
query:

``` r
import(cf, focus = 'II:300001-400000|II:600001-700000', resolution = 2000)
```

    ## `HiCExperiment` object with 402 contacts over 100 regions 
    ## -------
    ## fileName: "/home/scatsac/.cache/R/ExperimentHub/283eb5e62932f_7752" 
    ## focus: "II:300001-400000|II:600001-700000" 
    ## resolutions(5): 1000 2000 4000 8000 16000
    ## active resolution: 2000 
    ## interactions: 357 
    ## scores(2): count balanced 
    ## topologicalFeatures: compartments(0) borders(0) loops(0) viewpoints(0) 
    ## pairsFile: N/A 
    ## metadata(0):

Subset interactions constrained within single chromosome:

``` r
import(cf, focus = 'II', resolution = 2000)
```

    ## `HiCExperiment` object with 471,364 contacts over 407 regions 
    ## -------
    ## fileName: "/home/scatsac/.cache/R/ExperimentHub/283eb5e62932f_7752" 
    ## focus: "II" 
    ## resolutions(5): 1000 2000 4000 8000 16000
    ## active resolution: 2000 
    ## interactions: 34063 
    ## scores(2): count balanced 
    ## topologicalFeatures: compartments(0) borders(0) loops(0) viewpoints(0) 
    ## pairsFile: N/A 
    ## metadata(0):

Subset interactions to retain those between two chromosomes:

``` r
import(cf, focus = 'II|III', resolution = 2000)
```

    ## `HiCExperiment` object with 9,092 contacts over 566 regions 
    ## -------
    ## fileName: "/home/scatsac/.cache/R/ExperimentHub/283eb5e62932f_7752" 
    ## focus: "II|III" 
    ## resolutions(5): 1000 2000 4000 8000 16000
    ## active resolution: 2000 
    ## interactions: 7438 
    ## scores(2): count balanced 
    ## topologicalFeatures: compartments(0) borders(0) loops(0) viewpoints(0) 
    ## pairsFile: N/A 
    ## metadata(0):

Subset interactions to retain those between parts of two chromosomes:

``` r
import(cf, focus = 'II:300001-800000|V:1-500000', resolution = 2000)
```

    ## `HiCExperiment` object with 7,147 contacts over 500 regions 
    ## -------
    ## fileName: "/home/scatsac/.cache/R/ExperimentHub/283eb5e62932f_7752" 
    ## focus: "II:300001-800000|V:1-500000" 
    ## resolutions(5): 1000 2000 4000 8000 16000
    ## active resolution: 2000 
    ## interactions: 6523 
    ## scores(2): count balanced 
    ## topologicalFeatures: compartments(0) borders(0) loops(0) viewpoints(0) 
    ## pairsFile: N/A 
    ## metadata(0):

#### Subsetting after import with `subsetByOverlaps` or `[`

``` r
telomere <- GRanges("II:700001-813184")
subsetByOverlaps(hic, telomere) |> interactions()
```

    ## GInteractions object with 1540 interactions and 4 metadata columns:
    ##          seqnames1       ranges1     seqnames2       ranges2 |   bin_id1
    ##              <Rle>     <IRanges>         <Rle>     <IRanges> | <numeric>
    ##      [1]        II 700001-702000 ---        II 700001-702000 |       466
    ##      [2]        II 700001-702000 ---        II 702001-704000 |       466
    ##      [3]        II 700001-702000 ---        II 704001-706000 |       466
    ##      [4]        II 700001-702000 ---        II 706001-708000 |       466
    ##      [5]        II 700001-702000 ---        II 708001-710000 |       466
    ##      ...       ...           ... ...       ...           ... .       ...
    ##   [1536]        II 804001-806000 ---        II 810001-812000 |       518
    ##   [1537]        II 806001-808000 ---        II 806001-808000 |       519
    ##   [1538]        II 806001-808000 ---        II 808001-810000 |       519
    ##   [1539]        II 806001-808000 ---        II 810001-812000 |       519
    ##   [1540]        II 808001-810000 ---        II 808001-810000 |       520
    ##            bin_id2     count  balanced
    ##          <numeric> <numeric> <numeric>
    ##      [1]       466        30 0.0283618
    ##      [2]       467       145 0.0709380
    ##      [3]       468       124 0.0704979
    ##      [4]       469        59 0.0510221
    ##      [5]       470        59 0.0384004
    ##      ...       ...       ...       ...
    ##   [1536]       521         1       NaN
    ##   [1537]       519        15 0.0560633
    ##   [1538]       520        25       NaN
    ##   [1539]       521         1       NaN
    ##   [1540]       520        10       NaN
    ##   -------
    ##   regions: 57 ranges and 4 metadata columns
    ##   seqinfo: 16 sequences from an unspecified genome

Subsetting with `[`:

``` r
hic["II:800001-813184"]
```

    ## `HiCExperiment` object with 1,040 contacts over 6 regions 
    ## -------
    ## fileName: "/home/scatsac/.cache/R/ExperimentHub/283eb5e62932f_7752" 
    ## focus: "II:800,001-813,184" 
    ## resolutions(5): 1000 2000 4000 8000 16000
    ## active resolution: 2000 
    ## interactions: 19 
    ## scores(2): count balanced 
    ## topologicalFeatures: compartments(0) borders(0) loops(0) viewpoints(0) 
    ## pairsFile: N/A 
    ## metadata(0):

Subset off-diagonal genomic location using pairs of coordinates query:

``` r
hic["II:300001-320000|II:800001-813184"]
```

    ## `HiCExperiment` object with 3 contacts over 6 regions 
    ## -------
    ## fileName: "/home/scatsac/.cache/R/ExperimentHub/283eb5e62932f_7752" 
    ## focus: "II:300001-320000|II:800001-813184" 
    ## resolutions(5): 1000 2000 4000 8000 16000
    ## active resolution: 2000 
    ## interactions: 3 
    ## scores(2): count balanced 
    ## topologicalFeatures: compartments(0) borders(0) loops(0) viewpoints(0) 
    ## pairsFile: N/A 
    ## metadata(0):

Subsetting for vectors of several chromosomes cannot be performed with
`focus` from disk-stored data. This scenario requires `[`-based
in-memory subsetting of pre-imported data.

``` r
hic[c('II', 'III', 'IV')]
```

    ## `HiCExperiment` object with 306,212 contacts over 257 regions 
    ## -------
    ## fileName: "/home/scatsac/.cache/R/ExperimentHub/283eb5e62932f_7752" 
    ## focus: "II, III, IV" 
    ## resolutions(5): 1000 2000 4000 8000 16000
    ## active resolution: 2000 
    ## interactions: 18513 
    ## scores(2): count balanced 
    ## topologicalFeatures: compartments(0) borders(0) loops(0) viewpoints(0) 
    ## pairsFile: N/A 
    ## metadata(0):

#### Zooming in a `HiCExperiment`

The term “zooming” refers to dynamically changing the resolution of a
HiCExperiment. The `zoom` function only works for multi-resolution
contact matrices (e.g. `.mcool` or `.hic`) and it does not change the
`focus`. It only affects the `resolution` and `interactions`.

``` r
hic
```

    ## `HiCExperiment` object with 306,212 contacts over 257 regions 
    ## -------
    ## fileName: "/home/scatsac/.cache/R/ExperimentHub/283eb5e62932f_7752" 
    ## focus: "II:300,001-813,184" 
    ## resolutions(5): 1000 2000 4000 8000 16000
    ## active resolution: 2000 
    ## interactions: 18513 
    ## scores(2): count balanced 
    ## topologicalFeatures: compartments(0) borders(0) loops(0) viewpoints(0) 
    ## pairsFile: N/A 
    ## metadata(0):

``` r
# `HiCExperiment` object with 306,212 contacts over 257 regions 
# -------
# fileName: "/home/scatsac/.cache/R/ExperimentHub/283eb5e62932f_7752" 
# focus: "II:300,001-813,184" 
# resolutions(5): 1000 2000 4000 8000 16000
# active resolution: 2000 
# interactions: 18513 
# scores(2): count balanced 
# topologicalFeatures: compartments(0) borders(0) loops(0) viewpoints(0) 
# pairsFile: N/A 
# metadata(0):

zoom(hic, 4000)
```

    ## `HiCExperiment` object with 306,212 contacts over 129 regions 
    ## -------
    ## fileName: "/home/scatsac/.cache/R/ExperimentHub/283eb5e62932f_7752" 
    ## focus: "II:300,001-813,184" 
    ## resolutions(5): 1000 2000 4000 8000 16000
    ## active resolution: 4000 
    ## interactions: 6800 
    ## scores(2): count balanced 
    ## topologicalFeatures: compartments(0) borders(0) loops(0) viewpoints(0) 
    ## pairsFile: N/A 
    ## metadata(0):

``` r
# `HiCExperiment` object with 306,212 contacts over 129 regions 
# -------
# fileName: "/home/scatsac/.cache/R/ExperimentHub/283eb5e62932f_7752" 
# focus: "II:300,001-813,184" 
# resolutions(5): 1000 2000 4000 8000 16000
# active resolution: 4000 
# interactions: 6800 
# scores(2): count balanced 
# topologicalFeatures: compartments(0) borders(0) loops(0) viewpoints(0) 
# pairsFile: N/A 
# metadata(0):
```

#### Updating a `HiCExperiment` object

Immutable `HiCExperiment` slots:<br>
<ul>
<li>
`fileName(hic)`
</li>
<li>
`focus(hic)`
</li>
<li>
`resolutions(hic)`
</li>
<li>
`interactions(hic)`
</li>
</ul>
Mutable `HiCExperiment` slots:<br>
<ul>
<li>
`scores(hic)`
</li>
<li>
`topologicalFeatures(hic)`
</li>
<li>
`pairsFile(hic)`
</li>
<li>
`metadata(hic)`
</li>
</ul>

``` r
# example of mutable slots for a HiCExperiment object

#=== create additional topologicalFeatures or modify existing ones with topologicalFeatures()<-
topologicalFeatures(hic, 'CTCF') <- GRanges(c(
    "II:340-352", 
    "II:3520-3532", 
    "II:7980-7992", 
    "II:9240-9252" 
))
topologicalFeatures(hic, 'CTCF')
```

    ## GRanges object with 4 ranges and 0 metadata columns:
    ##       seqnames    ranges strand
    ##          <Rle> <IRanges>  <Rle>
    ##   [1]       II   340-352      *
    ##   [2]       II 3520-3532      *
    ##   [3]       II 7980-7992      *
    ##   [4]       II 9240-9252      *
    ##   -------
    ##   seqinfo: 1 sequence from an unspecified genome; no seqlengths

``` r
# GRanges object with 4 ranges and 0 metadata columns:
#       seqnames    ranges strand
#          <Rle> <IRanges>  <Rle>
#   [1]       II   340-352      *
#   [2]       II 3520-3532      *
#   [3]       II 7980-7992      *
#   [4]       II 9240-9252      *
#   -------
#   seqinfo: 1 sequence from an unspecified genome; no seqlengths
# 

topologicalFeatures(hic, 'loops') <- GInteractions(
    topologicalFeatures(hic, 'CTCF')[rep(1:3, each = 3)],
    topologicalFeatures(hic, 'CTCF')[rep(1:3, 3)]
)
topologicalFeatures(hic, 'loops')
```

    ## GInteractions object with 9 interactions and 0 metadata columns:
    ##       seqnames1   ranges1     seqnames2   ranges2
    ##           <Rle> <IRanges>         <Rle> <IRanges>
    ##   [1]        II   340-352 ---        II   340-352
    ##   [2]        II   340-352 ---        II 3520-3532
    ##   [3]        II   340-352 ---        II 7980-7992
    ##   [4]        II 3520-3532 ---        II   340-352
    ##   [5]        II 3520-3532 ---        II 3520-3532
    ##   [6]        II 3520-3532 ---        II 7980-7992
    ##   [7]        II 7980-7992 ---        II   340-352
    ##   [8]        II 7980-7992 ---        II 3520-3532
    ##   [9]        II 7980-7992 ---        II 7980-7992
    ##   -------
    ##   regions: 3 ranges and 0 metadata columns
    ##   seqinfo: 1 sequence from an unspecified genome; no seqlengths

``` r
# GInteractions object with 9 interactions and 0 metadata columns:
#       seqnames1   ranges1     seqnames2   ranges2
#           <Rle> <IRanges>         <Rle> <IRanges>
#   [1]        II   340-352 ---        II   340-352
#   [2]        II   340-352 ---        II 3520-3532
#   [3]        II   340-352 ---        II 7980-7992
#   [4]        II 3520-3532 ---        II   340-352
#   [5]        II 3520-3532 ---        II 3520-3532
#   [6]        II 3520-3532 ---        II 7980-7992
#   [7]        II 7980-7992 ---        II   340-352
#   [8]        II 7980-7992 ---        II 3520-3532
#   [9]        II 7980-7992 ---        II 7980-7992
#   -------
#   regions: 3 ranges and 0 metadata columns
#   seqinfo: 1 sequence from an unspecified genome; no seqlengths

#=== adding pairsFile after importing ContactFile into a HiCExperiment object
pairsf <- HiContactsData('yeast_wt', 'pairs.gz')
```

    ## see ?HiContactsData and browseVignettes('HiContactsData') for documentation

    ## loading from cache

``` r
pairsFile(hic) <- pairsf
# hic
# `HiCExperiment` object with 306,212 contacts over 257 regions 
# -------
# fileName: "/home/scatsac/.cache/R/ExperimentHub/283eb5e62932f_7752" 
# focus: "II:300,001-813,184" 
# resolutions(5): 1000 2000 4000 8000 16000
# active resolution: 2000 
# interactions: 18513 
# scores(2): count balanced 
# topologicalFeatures: compartments(0) borders(0) loops(9) viewpoints(0) CTCF(4) 
# pairsFile: /home/scatsac/.cache/R/ExperimentHub/283eb379aa82f_7753 
# metadata(0):

#=== update metadata for a HiCExperiment
metadata(hic) <- list(
    info = "HiCExperiment created from an example .mcool file from `HiContactsData`", 
    date = date()
)
metadata(hic)
```

    ## $info
    ## [1] "HiCExperiment created from an example .mcool file from `HiContactsData`"
    ## 
    ## $date
    ## [1] "Fri Mar 29 21:10:45 2024"

``` r
# $info
# [1] "HiCExperiment created from an example .mcool file from `HiContactsData`"
# 
# $date
# [1] "Fri Mar 29 20:58:31 2024"
```

#### Coercing `HiCExperiment` objects

Coercing functions allow for transformation of data stored as a
`HiCExperiment` into another class. These include:.<br>
<ul>
<li>
`as.matrix()`
</li>
<li>
`as.data.frame()`
</li>
</ul>

``` r
#=== `as.matrix` coerces into dense matrix by default
as.matrix(hic) |> class()
```

    ## [1] "dgTMatrix"
    ## attr(,"package")
    ## [1] "Matrix"

``` r
# [1] "dgTMatrix"
# attr(,"package")
# [1] "Matrix"

as.matrix(hic) |> dim()
```

    ## [1] 257 257

``` r
# [1] 257 257

#=== specify scores to be used when coercing into a matrix
as.matrix(hic, use.scores = "balanced")[1:5, 1:5]
```

    ## 5 x 5 sparse Matrix of class "dgTMatrix"
    ##                                                             
    ## [1,] 0.009657438 0.07662234 0.05410199 0.04294051 0.04090521
    ## [2,] 0.076622340 0.05128277 0.09841564 0.06926737 0.05263611
    ## [3,] 0.054101992 0.09841564 0.05657589 0.08723160 0.07316890
    ## [4,] 0.042940512 0.06926737 0.08723160 0.03699543 0.08403496
    ## [5,] 0.040905212 0.05263611 0.07316890 0.08403496 0.04787415

``` r
# 5 x 5 sparse Matrix of class "dgTMatrix"
#                                                             
# [1,] 0.009657438 0.07662234 0.05410199 0.04294051 0.04090521
# [2,] 0.076622340 0.05128277 0.09841564 0.06926737 0.05263611
# [3,] 0.054101992 0.09841564 0.05657589 0.08723160 0.07316890
# [4,] 0.042940512 0.06926737 0.08723160 0.03699543 0.08403496
# [5,] 0.040905212 0.05263611 0.07316890 0.08403496 0.04787415

as.matrix(hic, use.scores = "count")[1:5, 1:5]
```

    ## 5 x 5 sparse Matrix of class "dgTMatrix"
    ##                        
    ## [1,]  7  92  75  61  38
    ## [2,] 92 102 226 163  81
    ## [3,] 75 226 150 237 130
    ## [4,] 61 163 237 103 153
    ## [5,] 38  81 130 153  57

``` r
# 5 x 5 sparse Matrix of class "dgTMatrix"
#                        
# [1,]  7  92  75  61  38
# [2,] 92 102 226 163  81
# [3,] 75 226 150 237 130
# [4,] 61 163 237 103 153
# [5,] 38  81 130 153  57

#=== for a spare matrix, include sparse = TRUE
as.matrix(hic, use.scores = "count", sparse = TRUE)[1:5, 1:5]
```

    ## 5 x 5 sparse Matrix of class "dgTMatrix"
    ##                        
    ## [1,]  7  92  75  61  38
    ## [2,] 92 102 226 163  81
    ## [3,] 75 226 150 237 130
    ## [4,] 61 163 237 103 153
    ## [5,] 38  81 130 153  57

``` r
# 5 x 5 sparse Matrix of class "dgTMatrix"
#                        
# [1,]  7  92  75  61  38
# [2,] 92 102 226 163  81
# [3,] 75 226 150 237 130
# [4,] 61 163 237 103 153
# [5,] 38  81 130 153  57
```

Coercing `interactions` into a data frame (rectangular).

``` r
as.data.frame(hic) |> head()
```

    ##   seqnames1 start1   end1 width1 strand1 bin_id1    weight1 center1 seqnames2
    ## 1        II 300001 302000   2000       *     266 0.03714342  301000        II
    ## 2        II 300001 302000   2000       *     266 0.03714342  301000        II
    ## 3        II 300001 302000   2000       *     266 0.03714342  301000        II
    ## 4        II 300001 302000   2000       *     266 0.03714342  301000        II
    ## 5        II 300001 302000   2000       *     266 0.03714342  301000        II
    ## 6        II 300001 302000   2000       *     266 0.03714342  301000        II
    ##   start2   end2 width2 strand2 bin_id2    weight2 center2 count    balanced
    ## 1 300001 302000   2000       *     266 0.03714342  301000     7 0.009657438
    ## 2 302001 304000   2000       *     267 0.02242258  303000    92 0.076622340
    ## 3 304001 306000   2000       *     268 0.01942093  305000    75 0.054101992
    ## 4 306001 308000   2000       *     269 0.01895202  307000    61 0.042940512
    ## 5 308001 310000   2000       *     270 0.02898098  309000    38 0.040905212
    ## 6 310001 312000   2000       *     271 0.01834118  311000    43 0.029293930

``` r
#    seqnames1 start1   end1 width1 strand1 bin_id1    weight1 center1
#  1        II 300001 302000   2000       *     266 0.03714342  301000
#  2        II 300001 302000   2000       *     266 0.03714342  301000
#  3        II 300001 302000   2000       *     266 0.03714342  301000
#  4        II 300001 302000   2000       *     266 0.03714342  301000
#  5        II 300001 302000   2000       *     266 0.03714342  301000
#  6        II 300001 302000   2000       *     266 0.03714342  301000
#    seqnames2 start2   end2 width2 strand2 bin_id2    weight2 center2 count
#  1        II 300001 302000   2000       *     266 0.03714342  301000     7
#  2        II 302001 304000   2000       *     267 0.02242258  303000    92
#  3        II 304001 306000   2000       *     268 0.01942093  305000    75
#  4        II 306001 308000   2000       *     269 0.01895202  307000    61
#  5        II 308001 310000   2000       *     270 0.02898098  309000    38
#  6        II 310001 312000   2000       *     271 0.01834118  311000    43
#       balanced    random
#  1 0.009657438 0.7227776
#  2 0.076622340 0.4908796
#  3 0.054101992 0.8696518
#  4 0.042940512 0.1918134
#  5 0.040905212 0.5860512
#  6 0.029293930 0.5564170
```
