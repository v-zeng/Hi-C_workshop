Hi-C Data Visualization
================
Vinson Zeng
March 29, 2024

The objective of this notebook is to review the visualization tools in R
offered by:<br> 1) `HiContacts`<br> 2) `HiCExperiment`

#### Load packages and objects

``` r
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

``` r
library(HiContacts)
```

    ## Registered S3 methods overwritten by 'readr':
    ##   method                    from 
    ##   as.data.frame.spec_tbl_df vroom
    ##   as_tibble.spec_tbl_df     vroom
    ##   format.col_spec           vroom
    ##   print.col_spec            vroom
    ##   print.collector           vroom
    ##   print.date_names          vroom
    ##   print.locale              vroom
    ##   str.col_spec              vroom

``` r
library(rtracklayer)
```

    ## 
    ## Attaching package: 'rtracklayer'

    ## The following object is masked from 'package:AnnotationHub':
    ## 
    ##     hubUrl

### Visualizing Hi-C contact maps

#### Single map

Disk-stored Hi-C contact matrices can be visualized by importing the
interactions of interest over a genomic location into a `HiCExperiment`
object, followed by using `plotMatrix` from `HiContacts`.

``` r
#=== generate example hic object
# create HiCExperiment object from example .cool file in HiContactsData package

# download `.mcool` file and cache locally
coolf <- HiContactsData('yeast_wt', 'mcool')
```

    ## see ?HiContactsData and browseVignettes('HiContactsData') for documentation

    ## loading from cache

``` r
# create connection to disk-stored `.mcool` file
cf <- CoolFile(coolf)
cf
```

    ## CoolFile object
    ## .mcool file: /home/scatsac/.cache/R/ExperimentHub/283eb5e62932f_7752 
    ## resolution: 1000 
    ## pairs file: 
    ## metadata(0):

``` r
# import contacts from chromosome V at a resolution of 2000
hic <- import(cf, focus = 'V', resolution = 2000)

# plot the matrix
plotMatrix(hic)
```

![](hi-c_data_visualization_R_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

#### Horizontal map

A horizontal style Hi-C map can be generated when providing a
`maxDistance` argument to `plotMatrix`.

``` r
plotMatrix(hic, maxDistance = 200000)
```

![](hi-c_data_visualization_R_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

#### Side-by-side maps

Two Hi-C samples can be plotted side by side to compare interaction
landscopes over the same genomic locus.

``` r
# import a second `.mcool` file from Hi-C experiment performed in a eco1 yeast mutant
hic2 <- import(
    CoolFile(HiContactsData('yeast_eco1', 'mcool')), 
    focus = 'V', 
    resolution = 2000
)
```

    ## see ?HiContactsData and browseVignettes('HiContactsData') for documentation

    ## loading from cache

``` r
#=== plot 2 matrices side by side - top right will be first matrix, bottom left will be second matrix
plotMatrix(hic, compare.to = hic2)
```

![](hi-c_data_visualization_R_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

#### Plot multiple chromosomes

A Hi-C heatmap can be used to visualize interactions from multiple
chromosomes through the following steps:<br> 1. Parse entire contact
matrix in `R` <br> 2. Subset interactions over chromosomes of interest
with `[` <br> 3. Use `plotMatrix` to generate the multi-chromosome plot

Step 1.

``` r
# parse entire contact matrix
full_hic <- import(cf, resolution = 4000)
plotMatrix(full_hic)
```

![](hi-c_data_visualization_R_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

Steps 2 and 3.

``` r
# subset interactions over chromosomes of interest, plotMatrix
hic_subset <- full_hic[c("II", "III", "IV")]
plotMatrix(hic_subset)
```

![](hi-c_data_visualization_R_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

### Hi-C maps customization options

The `plotMatrix` function has several customization options:
<ul>
<li>
pick scores of interest to represent in a Hi-C heatmap
</li>
<li>
change numeric scale and boundaries
</li>
<li>
change the color map
</li>
<li>
extra customization options
</li>
</ul>

``` r
# plot count scores (un-normalized raw contact counts obtained from binning a .pairs file)
plotMatrix(hic, use.scores = 'count')
```

![](hi-c_data_visualization_R_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
# choose a scale by providing the `limits` argument
plotMatrix(hic, limits = c(-3.5, -1))
```

![](hi-c_data_visualization_R_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
# choose a colour map using `cmap`, ?HiContacts::palettes will show list of available colour maps
afmhotrColors() # this colour map is provided in the `HiContacts` package
```

    ## [1] "#ffffff" "#f8f5c3" "#f4ee8d" "#f6be35" "#ee7d32" "#c44228" "#821d19"
    ## [8] "#381211" "#050606"

``` r
plotMatrix(
    hic, 
    use.scores = 'balanced',
    limits = c(-4, -1),
    cmap = afmhotrColors()
)
```

![](hi-c_data_visualization_R_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->
\### Advanced visualization

#### Overlaying topological features

Topological features such as chromatin loops, domain borders, and A/B
compartments can be displayed over a Hi-C heatmap.

``` r
# import pre-computed chromatin loops identified using`chromosight` on the contact matrix previously used to import interactions from

library(rtracklayer)
library(InteractionSet)
loops <- system.file('extdata', 'S288C-loops.bedpe', package = 'HiCExperiment') |> 
    import() |> 
    makeGInteractionsFromGRangesPairs()
loops
```

    ## GInteractions object with 162 interactions and 0 metadata columns:
    ##         seqnames1       ranges1     seqnames2       ranges2
    ##             <Rle>     <IRanges>         <Rle>     <IRanges>
    ##     [1]         I     3001-4000 ---         I   29001-30000
    ##     [2]         I   29001-30000 ---         I   50001-51000
    ##     [3]         I   95001-96000 ---         I 128001-129000
    ##     [4]         I 133001-134000 ---         I 157001-158000
    ##     [5]        II     8001-9000 ---        II   46001-47000
    ##     ...       ...           ... ...       ...           ...
    ##   [158]       XVI 773001-774000 ---       XVI 803001-804000
    ##   [159]       XVI 834001-835000 ---       XVI 859001-860000
    ##   [160]       XVI 860001-861000 ---       XVI 884001-885000
    ##   [161]       XVI 901001-902000 ---       XVI 940001-941000
    ##   [162]       XVI 917001-918000 ---       XVI 939001-940000
    ##   -------
    ##   regions: 316 ranges and 0 metadata columns
    ##   seqinfo: 16 sequences from an unspecified genome; no seqlengths

Borders have also been mapped with `chromosight` and can be imported in
`R`.

``` r
borders <- system.file('extdata', 'S288C-borders.bed', package = 'HiCExperiment') |> 
    import()
borders
```

    ## GRanges object with 814 ranges and 0 metadata columns:
    ##         seqnames        ranges strand
    ##            <Rle>     <IRanges>  <Rle>
    ##     [1]        I   73001-74000      *
    ##     [2]        I 108001-109000      *
    ##     [3]        I 181001-182000      *
    ##     [4]       II   90001-91000      *
    ##     [5]       II 119001-120000      *
    ##     ...      ...           ...    ...
    ##   [810]      XVI 777001-778000      *
    ##   [811]      XVI 796001-797000      *
    ##   [812]      XVI 811001-812000      *
    ##   [813]      XVI 890001-891000      *
    ##   [814]      XVI 933001-934000      *
    ##   -------
    ##   seqinfo: 16 sequences from an unspecified genome; no seqlengths

`GInteractions` stores chromatin loops, displayed as off-diagonal
circles. Borders are in `GRanges` and are displayed as on-diagonal
diamonds in the Hi-C heatmap.

``` r
plotMatrix(hic, loops = loops, borders = borders)
```

![](hi-c_data_visualization_R_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

#### Aggregated Hi-C maps

“Snippets” (i.e. extracts) are typically aggregated together to display
an average signal. This is sometimes referred to as “Aggregated Plot
Analysis” (APA).<br> The `aggregate` function can compute aggregated
Hi-C maps over a collection of `targets`. Targets can be `GRanges` for
on-diagonal snippets or `GInteractions` for off-diagonal snippets.
Additionally, `flankingBins` specifies the number of matrix bins
extracted on each side of the `targets` of interest.

``` r
# compute aggregated Hi-C snippets of +/- 15kb around each chromatin loop listed in `loops`
hic <- zoom(hic, 1000)
aggr_loops <- aggregate(hic, targets = loops, flankingBins = 15)
```

    ## Going through preflight checklist...

    ## Parsing the entire contact matrice as a sparse matrix...

    ## Modeling distance decay...

    ## Filtering for contacts within provided targets...

``` r
aggr_loops
```

    ## `AggrHiCExperiment` object over 148 targets 
    ## -------
    ## fileName: "/home/scatsac/.cache/R/ExperimentHub/283eb5e62932f_7752" 
    ## focus: 148 targets 
    ## resolutions(5): 1000 2000 4000 8000 16000
    ## active resolution: 1000 
    ## interactions: 961 
    ## scores(4): count balanced expected detrended 
    ## slices(4): count balanced expected detrended 
    ## topologicalFeatures: targets(148) compartments(0) borders(0) loops(0) viewpoints(0) 
    ## pairsFile: N/A 
    ## metadata(0):

The `aggregate` function generates a `AggrHiCExperiment` object, which
has an extra `slices` slot for arrays.

``` r
slices(aggr_loops)
```

    ## List of length 4
    ## names(4): count balanced expected detrended

``` r
dim(slices(aggr_loops, 'count'))
```

    ## [1]  31  31 148

``` r
topologicalFeatures(aggr_loops, 'targets')
```

    ## Pairs object with 148 pairs and 0 metadata columns:
    ##                     first            second
    ##                 <GRanges>         <GRanges>
    ##     [1]     I:14501-44500     I:35501-65500
    ##     [2]    I:80501-110500   I:113501-143500
    ##     [3]   I:118501-148500   I:142501-172500
    ##     [4]    II:33501-63500    II:63501-93500
    ##     [5]  II:134501-164500  II:159501-189500
    ##     ...               ...               ...
    ##   [144] XVI:586501-616500 XVI:606501-636500
    ##   [145] XVI:733501-763500 XVI:754501-784500
    ##   [146] XVI:758501-788500 XVI:788501-818500
    ##   [147] XVI:819501-849500 XVI:844501-874500
    ##   [148] XVI:845501-875500 XVI:869501-899500

The `AggrHiCExperiment` can be plotted with `plotMatrix`.

``` r
plotMatrix(
    aggr_loops, 
    use.scores = 'detrended', 
    scale = 'linear', 
    limits = c(-1, 1), 
    cmap = bgrColors()
)
```

![](hi-c_data_visualization_R_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->
