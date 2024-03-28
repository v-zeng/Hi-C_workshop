Data Gateways: Accessing Public Hi-C Data Portals
================
Vinson Zeng
March 28, 2024

The objective of this notebook is to develop working knowledge of
accessing Hi-C datasets from:<br> 1) 4DN Consortium<br> 2) DNA Zoo
Project

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
library(fourDNData)
library(DNAZooData)
library(rtracklayer)
```

## 4DN Data Portal: 4D Nucleome Data Coordination and Integration Center (DCIC)

Use fourDNData() function to access gateway to 4DN-hosted Hi-C files.

``` r
library(fourDNData)
head(fourDNData())
```

    ##   experimentSetAccession     fileType     size organism experimentType details
    ## 1           4DNES18BMU79        pairs 10151.53    mouse   in situ Hi-C   DpnII
    ## 3           4DNES18BMU79          hic  5285.82    mouse   in situ Hi-C   DpnII
    ## 4           4DNES18BMU79        mcool  6110.75    mouse   in situ Hi-C   DpnII
    ## 5           4DNES18BMU79   boundaries     0.12    mouse   in situ Hi-C   DpnII
    ## 6           4DNES18BMU79   insulation     7.18    mouse   in situ Hi-C   DpnII
    ## 7           4DNES18BMU79 compartments     0.18    mouse   in situ Hi-C   DpnII
    ##                                dataset
    ## 1 Hi-C on Mouse Olfactory System cells
    ## 3 Hi-C on Mouse Olfactory System cells
    ## 4 Hi-C on Mouse Olfactory System cells
    ## 5 Hi-C on Mouse Olfactory System cells
    ## 6 Hi-C on Mouse Olfactory System cells
    ## 7 Hi-C on Mouse Olfactory System cells
    ##                                                         condition
    ## 1 Mature olfactory sensory neurons with conditional Ldb1 knockout
    ## 3 Mature olfactory sensory neurons with conditional Ldb1 knockout
    ## 4 Mature olfactory sensory neurons with conditional Ldb1 knockout
    ## 5 Mature olfactory sensory neurons with conditional Ldb1 knockout
    ## 6 Mature olfactory sensory neurons with conditional Ldb1 knockout
    ## 7 Mature olfactory sensory neurons with conditional Ldb1 knockout
    ##                 biosource biosourceType             publication
    ## 1 olfactory receptor cell  primary cell Monahan K et al. (2019)
    ## 3 olfactory receptor cell  primary cell Monahan K et al. (2019)
    ## 4 olfactory receptor cell  primary cell Monahan K et al. (2019)
    ## 5 olfactory receptor cell  primary cell Monahan K et al. (2019)
    ## 6 olfactory receptor cell  primary cell Monahan K et al. (2019)
    ## 7 olfactory receptor cell  primary cell Monahan K et al. (2019)
    ##                                                                                                                                   URL
    ## 1 https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/49504f97-904e-48c1-8c20-1033680b66da/4DNFIC5AHBPV.pairs.gz
    ## 3      https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/6cd4378a-8f51-4e65-99eb-15f5c80abf8d/4DNFIT4I5C6Z.hic
    ## 4    https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/01fb704f-2fd7-48c6-91af-c5f4584529ed/4DNFIVPAXJO8.mcool
    ## 5   https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/5c07cdee-53e2-43e0-8853-cfe5f057b3f1/4DNFIR3XCIMA.bed.gz
    ## 6       https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/d1f4beb9-701f-4188-abe2-6271fe658770/4DNFIXKKNMS7.bw
    ## 7       https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/3d429647-51c8-4e3a-a18b-eec0b1480905/4DNFIN13N8C1.bw

#### Query individual files with specific file parameters.

``` r
cf <- fourDNData(experimentSetAccession = '4DNESJNPEKZD', type = 'mcool')
```

This will download and cache any queried files locally.

``` r
cf
```

    ## [1] "/home/scatsac/.cache/R/fourDNData/283eb612fa19e_4DNFIZL8OZE1.mcool"

``` r
availableChromosomes(cf)
```

    ## Seqinfo object with 24 sequences from an unspecified genome:
    ##   seqnames seqlengths isCircular genome
    ##   chr1      248956422       <NA>   <NA>
    ##   chr2      242193529       <NA>   <NA>
    ##   chr3      198295559       <NA>   <NA>
    ##   chr4      190214555       <NA>   <NA>
    ##   chr5      181538259       <NA>   <NA>
    ##   ...             ...        ...    ...
    ##   chr20      64444167       <NA>   <NA>
    ##   chr21      46709983       <NA>   <NA>
    ##   chr22      50818468       <NA>   <NA>
    ##   chrX      156040895       <NA>   <NA>
    ##   chrY       57227415       <NA>   <NA>

``` r
availableResolutions(cf)
```

    ## resolutions(13): 1000 2000 ... 5e+06 1e+07

    ## 

``` r
import(cf, focus = "chr4:10000001-20000000", resolution = 5000)
```

    ## `HiCExperiment` object with 656 contacts over 2,000 regions 
    ## -------
    ## fileName: "/home/scatsac/.cache/R/fourDNData/283eb612fa19e_4DNFIZL8OZE1.mcool" 
    ## focus: "chr4:10,000,001-20,000,000" 
    ## resolutions(13): 1000 2000 ... 5000000 10000000
    ## active resolution: 5000 
    ## interactions: 614 
    ## scores(2): count balanced 
    ## topologicalFeatures: compartments(0) borders(0) loops(0) viewpoints(0) 
    ## pairsFile: N/A 
    ## metadata(0):

To fetch a specific Hi-C file type, such as ‘pairs’, use the ‘type’
argument.

``` r
pairs_f <- fourDNData(experimentSetAccession = '4DNESJNPEKZD', type = 'pairs') 
print(pairs_f)
```

    ## [1] "/home/scatsac/.cache/R/fourDNData/353733cf5c326_4DNFIOZ7D1OQ.pairs.gz"

``` r
import(pairs_f)
```

    ## GInteractions object with 453301 interactions and 3 metadata columns:
    ##            seqnames1   ranges1     seqnames2   ranges2 |       frag1     frag2
    ##                <Rle> <IRanges>         <Rle> <IRanges> | <character> <numeric>
    ##        [1]      chr1    618317 ---      chr1  10432996 |          UU        21
    ##        [2]      chr1    634055 ---      chr1 109762561 |          UU        10
    ##        [3]      chr1    676764 ---      chr1 207054459 |          UU         3
    ##        [4]      chr1    783553 ---      chr1 211368282 |          UU        29
    ##        [5]      chr1    791439 ---      chr1    791877 |          UU        31
    ##        ...       ...       ... ...       ...       ... .         ...       ...
    ##   [453297]      chrY  56836692 ---      chrY  56869006 |          UU         3
    ##   [453298]      chrY  56842643 ---      chrY  56843054 |          UU        60
    ##   [453299]      chrY  56850437 ---      chrY  56850710 |          UU        40
    ##   [453300]      chrY  56858287 ---      chrY  56858888 |          UU        60
    ##   [453301]      chrY  56885246 ---      chrY  56885709 |          UU        15
    ##             distance
    ##            <integer>
    ##        [1]   9814679
    ##        [2] 109128506
    ##        [3] 206377695
    ##        [4] 210584729
    ##        [5]       438
    ##        ...       ...
    ##   [453297]     32314
    ##   [453298]       411
    ##   [453299]       273
    ##   [453300]       601
    ##   [453301]       463
    ##   -------
    ##   regions: 905085 ranges and 0 metadata columns
    ##   seqinfo: 24 sequences from an unspecified genome; no seqlengths

For type = ‘insulation’, this fetches a .bigwig track file. This track
type corresponds to a genome-wide insulation score. Once fetched from
the 4DN portal, local files can be imported in R.

``` r
library(rtracklayer)
fourDNData(experimentSetAccession = '4DNES25ABNZ1', type = 'insulation') |> 
    import(as = 'Rle')
```

    ## RleList of length 21
    ## $chr1
    ## numeric-Rle of length 195471971 with 38145 runs
    ##   Lengths:      3065000         5000         5000 ...         5000       171971
    ##   Values :  0.00000e+00  1.01441e-01  7.41053e-02 ...     0.807009     0.000000
    ## 
    ## $chr10
    ## numeric-Rle of length 130694993 with 25100 runs
    ##   Lengths:     3175000        5000        5000 ...        5000      169993
    ##   Values :  0.00000000  0.37584546  0.33597839 ...    0.628601    0.000000
    ## 
    ## $chr11
    ## numeric-Rle of length 122082543 with 23536 runs
    ##   Lengths:    3165000       5000       5000 ...       5000     162543
    ##   Values :  0.0000000 -0.7906257 -0.7930040 ...   0.515919   0.000000
    ## 
    ## $chr12
    ## numeric-Rle of length 120129022 with 22578 runs
    ##   Lengths:   3075000      5000      5000 ...      5000      5000    164022
    ##   Values :  0.000000  0.411216  0.400357 ... 0.1650951 0.2175749 0.0000000
    ## 
    ## $chr13
    ## numeric-Rle of length 120421639 with 22807 runs
    ##   Lengths:     3080000        5000        5000 ...        5000      171639
    ##   Values :  0.00000000  0.17005745  0.10652249 ...  1.14856148  0.00000000
    ## 
    ## ...
    ## <16 more elements>

For type = ‘boundaries’, it fetches a .bed file which lists annotated
borders between topological domains. This file type can also be imported
in R and will generate a GRanges object.

``` r
fourDNData(experimentSetAccession = '4DNES25ABNZ1', type = 'boundaries') |> 
    import()
```

    ## GRanges object with 6103 ranges and 2 metadata columns:
    ##          seqnames            ranges strand |        name     score
    ##             <Rle>         <IRanges>  <Rle> | <character> <numeric>
    ##      [1]     chr1   4380001-4385000      * |      Strong  0.695274
    ##      [2]     chr1   4760001-4765000      * |        Weak  0.444476
    ##      [3]     chr1   4910001-4915000      * |        Weak  0.353184
    ##      [4]     chr1   5180001-5185000      * |      Strong  0.565763
    ##      [5]     chr1   6170001-6175000      * |      Strong  1.644911
    ##      ...      ...               ...    ... .         ...       ...
    ##   [6099]     chrY 89725001-89730000      * |        Weak  0.258094
    ##   [6100]     chrY 89790001-89795000      * |        Weak  0.442186
    ##   [6101]     chrY 89895001-89900000      * |        Weak  0.279879
    ##   [6102]     chrY 90025001-90030000      * |      Strong  0.660699
    ##   [6103]     chrY 90410001-90415000      * |      Strong  1.160018
    ##   -------
    ##   seqinfo: 21 sequences from an unspecified genome; no seqlengths

For type = ‘compartments’, this will fetch a .bigwig track file. This
track corresponds to a selected genome-wide eigenvector computed by
cooltools, which represents A/B compartments.

``` r
fourDNData(experimentSetAccession = '4DNES25ABNZ1', type = 'compartments') |> 
    import()
```

    ## GRanges object with 10911 ranges and 1 metadata column:
    ##           seqnames            ranges strand |     score
    ##              <Rle>         <IRanges>  <Rle> | <numeric>
    ##       [1]     chr1          1-250000      * |       NaN
    ##       [2]     chr1     250001-500000      * |       NaN
    ##       [3]     chr1     500001-750000      * |       NaN
    ##       [4]     chr1    750001-1000000      * |       NaN
    ##       [5]     chr1   1000001-1250000      * |       NaN
    ##       ...      ...               ...    ... .       ...
    ##   [10907]     chrY 90500001-90750000      * | 0.0237907
    ##   [10908]     chrY 90750001-91000000      * |       NaN
    ##   [10909]     chrY 91000001-91250000      * |       NaN
    ##   [10910]     chrY 91250001-91500000      * |       NaN
    ##   [10911]     chrY 91500001-91744698      * |       NaN
    ##   -------
    ##   seqinfo: 21 sequences from an unspecified genome

#### Query complete experiment datasets

Instead of importing multiple files individually for each single
experimentSet accession ID, we can import all available files
simultaneously using fourDNHiCEperiment().

``` r
hic <- fourDNHiCExperiment('4DNESJNPEKZD')
```

    ## Fetching local Hi-C contact map from Bioc cache

    ## Fetching local compartments bigwig file from Bioc cache

    ## Insulation not found for the provided experimentSet accession.

    ## Borders not found for the provided experimentSet accession.

    ## Importing contacts in memory

This aggregates all the components into a single HiCExperiment object
with topologicalFeatures and metadata slots.

``` r
metadata(hic)
```

    ## $`4DN_info`
    ##      experimentSetAccession     fileType   size organism experimentType details
    ## 1376           4DNESJNPEKZD        pairs   6.67    human   in situ Hi-C    MboI
    ## 1378           4DNESJNPEKZD          hic 179.51    human   in situ Hi-C    MboI
    ## 1379           4DNESJNPEKZD        mcool  30.17    human   in situ Hi-C    MboI
    ## 1380           4DNESJNPEKZD compartments   0.21    human   in situ Hi-C    MboI
    ##                                          dataset
    ## 1376 Hi-C on GM12878 cells - protocol variations
    ## 1378 Hi-C on GM12878 cells - protocol variations
    ## 1379 Hi-C on GM12878 cells - protocol variations
    ## 1380 Hi-C on GM12878 cells - protocol variations
    ##                                                      condition biosource
    ## 1376 in situ Hi-C on GM12878 with MboI and bio-dUTP (Tri-Link)   GM12878
    ## 1378 in situ Hi-C on GM12878 with MboI and bio-dUTP (Tri-Link)   GM12878
    ## 1379 in situ Hi-C on GM12878 with MboI and bio-dUTP (Tri-Link)   GM12878
    ## 1380 in situ Hi-C on GM12878 with MboI and bio-dUTP (Tri-Link)   GM12878
    ##               biosourceType          publication
    ## 1376 immortalized cell line Rao SS et al. (2014)
    ## 1378 immortalized cell line Rao SS et al. (2014)
    ## 1379 immortalized cell line Rao SS et al. (2014)
    ## 1380 immortalized cell line Rao SS et al. (2014)
    ##                                                                                                                                      URL
    ## 1376 https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/0bdd4745-7203-49d0-adf6-291cef1a96b7/4DNFIOZ7D1OQ.pairs.gz
    ## 1378      https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/1201682a-a223-482d-913d-3c3972b8eb65/4DNFIIRIHBR2.hic
    ## 1379    https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/356fab42-5562-4cfd-a3f8-592aa060b992/4DNFIZL8OZE1.mcool
    ## 1380       https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/333aabfd-b747-447c-b93a-8138f9488fad/4DNFIO9V5G93.bw
    ## 
    ## $eigens
    ## GRanges object with 11280 ranges and 2 metadata columns:
    ##           seqnames              ranges strand |       score       eigen
    ##              <Rle>           <IRanges>  <Rle> |   <numeric>   <numeric>
    ##       [1]     chr1      750001-1000000      * |   1.6911879   1.6911879
    ##       [2]     chr1     1000001-1250000      * |   0.0809129   0.0809129
    ##       [3]     chr1     1250001-1500000      * |   0.0690173   0.0690173
    ##       [4]     chr1     1500001-1750000      * |  -0.1903324  -0.1903324
    ##       [5]     chr1     1750001-2000000      * |   0.3283633   0.3283633
    ##       ...      ...                 ...    ... .         ...         ...
    ##   [11276]     chrX 154750001-155000000      * | -0.10909061 -0.10909061
    ##   [11277]     chrX 155000001-155250000      * | -1.39655280 -1.39655280
    ##   [11278]     chrX 155250001-155500000      * |  0.00264734  0.00264734
    ##   [11279]     chrX 155500001-155750000      * | -0.15279847 -0.15279847
    ##   [11280]     chrX 155750001-156000000      * | -1.41699576 -1.41699576
    ##   -------
    ##   seqinfo: 24 sequences from an unspecified genome

## DNA Zoo

The DNA Zoo Consortium is a data portal maintained by a group who aims
to correct and refine genome assemblies using Hi-C approaches. The
DNAZooData() function from the similarly titled package provides a
gateway to their Hi-C files.

``` r
library(DNAZooData)
head(DNAZooData())
```

    ##                   species                              readme
    ## 1        Acinonyx_jubatus        Acinonyx_jubatus/README.json
    ## 2      Acropora_millepora      Acropora_millepora/README.json
    ## 3     Addax_nasomaculatus     Addax_nasomaculatus/README.json
    ## 4           Aedes_aegypti           Aedes_aegypti/README.json
    ## 5   Aedes_aegypti__AaegL4   Aedes_aegypti__AaegL4/README.json
    ## 6 Aedes_aegypti__AaegL5.0 Aedes_aegypti__AaegL5.0/README.json
    ##                                                           readme_link
    ## 1        https://dnazoo.s3.wasabisys.com/Acinonyx_jubatus/README.json
    ## 2      https://dnazoo.s3.wasabisys.com/Acropora_millepora/README.json
    ## 3     https://dnazoo.s3.wasabisys.com/Addax_nasomaculatus/README.json
    ## 4           https://dnazoo.s3.wasabisys.com/Aedes_aegypti/README.json
    ## 5   https://dnazoo.s3.wasabisys.com/Aedes_aegypti__AaegL4/README.json
    ## 6 https://dnazoo.s3.wasabisys.com/Aedes_aegypti__AaegL5.0/README.json
    ##   original_assembly     new_assembly
    ## 1           aciJub1      aciJub1_HiC
    ## 2       amil_sf_1.1  amil_sf_1.1_HiC
    ## 3      ASM1959352v1 ASM1959352v1_HiC
    ## 4        AGWG.draft         AaegL5.0
    ## 5            AaegL3           AaegL4
    ## 6        AGWG.draft         AaegL5.0
    ##                                                               new_assembly_link
    ## 1         https://dnazoo.s3.wasabisys.com/Acinonyx_jubatus/aciJub1_HiC.fasta.gz
    ## 2   https://dnazoo.s3.wasabisys.com/Acropora_millepora/amil_sf_1.1_HiC.fasta.gz
    ## 3 https://dnazoo.s3.wasabisys.com/Addax_nasomaculatus/ASM1959352v1_HiC.fasta.gz
    ## 4               https://dnazoo.s3.wasabisys.com/Aedes_aegypti/AaegL5.0.fasta.gz
    ## 5         https://dnazoo.s3.wasabisys.com/Aedes_aegypti__AaegL4/AaegL4.fasta.gz
    ## 6     https://dnazoo.s3.wasabisys.com/Aedes_aegypti__AaegL5.0/AaegL5.0.fasta.gz
    ##   new_assembly_link_status
    ## 1                      200
    ## 2                      200
    ## 3                      200
    ## 4                      404
    ## 5                      200
    ## 6                      200
    ##                                                                   hic_link
    ## 1    https://dnazoo.s3.wasabisys.com/Acinonyx_jubatus/aciJub1.rawchrom.hic
    ## 2   https://dnazoo.s3.wasabisys.com/Acropora_millepora/amil_sf_1.1_HiC.hic
    ## 3 https://dnazoo.s3.wasabisys.com/Addax_nasomaculatus/ASM1959352v1_HiC.hic
    ## 4                                                                     <NA>
    ## 5         https://dnazoo.s3.wasabisys.com/Aedes_aegypti__AaegL4/AaegL4.hic
    ## 6     https://dnazoo.s3.wasabisys.com/Aedes_aegypti__AaegL5.0/AaegL5.0.hic

Fetching a Hi-C dataset generated from a tardigrade sample.

``` r
hicfile <- DNAZooData(species = 'Hypsibius_dujardini')
hicfile
```

    ## HicFile object
    ## .hic file: /home/scatsac/.cache/R/DNAZooData/353737c133efc_nHd_3.1_HiC.hic 
    ## resolution: 5000 
    ## pairs file: 
    ## metadata(6): organism draftAssembly ... credits assemblyURL

We can check the metadata in the HicFile parsed from the DNA Zoo portal.

``` r
metadata(hicfile)$organism
```

    ## $vernacular
    ## [1] "Tardigrade"
    ## 
    ## $binomial
    ## [1] "Hypsibius dujardini"
    ## 
    ## $funFact
    ## [1] "<i>Hypsibius dujardini</i> is a species of tardigrade, a tiny microscopic organism. They are also commonly called water bears. This species is found world-wide!"
    ## 
    ## $extraInfo
    ## [1] "on BioKIDS website"
    ## 
    ## $extraInfoLink
    ## [1] "http://www.biokids.umich.edu/critters/Hypsibius_dujardini/"
    ## 
    ## $image
    ## [1] "https://static.wixstatic.com/media/2b9330_82db39c219f24b20a75cb38943aad1fb~mv2.jpg"
    ## 
    ## $imageCredit
    ## [1] "By Willow Gabriel, Goldstein Lab - https://www.flickr.com/photos/waterbears/1614095719/ Template:Uploader Transferred from en.wikipedia to Commons., CC BY-SA 2.5, https://commons.wikimedia.org/w/index.php?curid=2261992"
    ## 
    ## $isChromognomes
    ## [1] "FALSE"
    ## 
    ## $taxonomy
    ## [1] "Species:202423-914154-914155-914158-155166-155362-710171-710179-710192-155390-155420"
