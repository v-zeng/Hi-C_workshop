### Workshop: Orchestrating Hi-C analysis with Bioconductor
## Part 1: Preprocessing Hi-C Data
## Part 2: Investigating Hi-C data structures in R
#======================================================
# The purpose of this notebook is to investigate the four of the main classes
# leveraged by Bioconductor in Hi-C analyses, their structures and how to
# interact with them.
# These classes include:
#   1. GRanges
#   2. GInteractions
#   3. ContactFile
#   4. HiCExperiment
#=======================================================


# set libPaths
.libPaths("/home/scatsac/rlibrary")
setwd("/home/scatsac/rlibrary")
### Part 1. Hi-C Pre-Processing=================================================

# install BiocManager and packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")

BiocManager::install("HiCExperiment", ask = FALSE)
BiocManager::install("HiCool", ask = FALSE)
BiocManager::install("HiContacts", ask = FALSE)
BiocManager::install("HiContactsData", ask = FALSE)
BiocManager::install("fourDNData", ask = FALSE)
BiocManager::install("DNAZooData", ask = FALSE)
BiocManager::install("ggbio", ask = FALSE)

# processing example fastq files with HiCool
library(HiContactsData)
r1 <- HiContactsData(sample = 'yeast_wt', format = 'fastq_R1')
r2 <- HiContactsData(sample = 'yeast_wt', format = 'fastq_R2')

# load HiCool
library(HiCool)

# check documentation for functions
?HiCool

# process paired-end HiC sequencing files with HiCool
HiCool(
  r1, # path to fastq file (R1 read)
  r2, # path to fastq file (R2 read)
  restriction = 'DpnII,HinfI', # REase
  resolutions = c(4000, 8000, 16000), # resolution to bin mcool file
  genome = 'R64-1-1', # genome used to map reads on
  output = './HiCool/' # output folder
)

# check generated output in HiCool directory in tree structure format
fs::dir_tree('HiCool/')

# Below the *.pairs and *.mcool files are pairs and contact matrix files, resp.
# *.html file is a report summarizing pairs numbers, filtering, etc.
# *.log contains all output, error messages, and commands used for preprocessing
# *.pdf graphic files are visual representation of distribution of pairs

# HiCool/
#   ├── 4f6ff9f4b51_7833^mapped-R64-1-1^BYD5IS.html
# ├── logs
# │   └── 4f6ff9f4b51_7833^mapped-R64-1-1^BYD5IS.log
# ├── matrices
# │   └── 4f6ff9f4b51_7833^mapped-R64-1-1^BYD5IS.mcool
# ├── pairs
# │   └── 4f6ff9f4b51_7833^mapped-R64-1-1^BYD5IS.pairs
# └── plots
# ├── 4f6ff9f4b51_7833^mapped-R64-1-1^BYD5IS_event_distance.pdf
# └── 4f6ff9f4b51_7833^mapped-R64-1-1^BYD5IS_event_distribution.pdf

### Part 2. Hi-C Data Structures in R ==========================================

## GRanges: GenomicRanges class-------------------------------------
# e.g., promoters, SNPs, chromatin loop anchors, etc.
library(GenomicRanges)
?GenomicRanges
?GenomicRanges::`intra-range-methods`

# can generate GRanges object with UCSC format from genomic coordinates vector
gr <- GRanges(c(
  "chr2:2004-7853:+", 
  "chr4:4482-9873:-", 
  "chr5:1943-4203:+", 
  "chr5:4103-5004:+"  
))

gr
# GRanges object with 4 ranges and 0 metadata columns:
#   seqnames    ranges strand
# <Rle> <IRanges>  <Rle>
#   [1]     chr2 2004-7853      +
#   [2]     chr4 4482-9873      -
#   [3]     chr5 1943-4203      +
#   [4]     chr5 4103-5004      +
#   -------
#   seqinfo: 3 sequences from an unspecified genome; no seqlengths

# navigate GRanges object ranges using a subset
gr[1]
# GRanges object with 1 range and 0 metadata columns:
#   seqnames    ranges strand
# <Rle> <IRanges>  <Rle>
#   [1]     chr2 2004-7853      +
#   -------
#   seqinfo: 3 sequences from an unspecified genome; no seqlengths

gr[2:3]
# GRanges object with 2 ranges and 0 metadata columns:
#   seqnames    ranges strand
# <Rle> <IRanges>  <Rle>
#   [1]     chr4 4482-9873      -
#   [2]     chr5 1943-4203      +
#   -------
#   seqinfo: 3 sequences from an unspecified genome; no seqlengths


# GRanges objects provide the following:
# seqnames (i.e., chromosome names), start(), end(), strand()

seqnames(gr)
# factor-Rle of length 4 with 3 runs
# Lengths:    1    1    2
# Values : chr2 chr4 chr5
# Levels(3): chr2 chr4 chr5

start(gr)
# [1] 2004 4482 1943 4103

end(gr)
# [1] 7853 9873 4203 5004

strand(gr)
# factor-Rle of length 4 with 3 runs
# Lengths: 1 1 2
# Values : + - +
# Levels(3): + - *

## GRanges can have optional metadata stored in a DataFrame
mcols(gr)
# DataFrame with 4 rows and 0 columns
mcols(gr)$GC <- c(0.45, 0.43, 0.44, 0.42)
mcols(gr)$annotation <- factor(c(NA, 'promoter', 'enhancer', 'centromere'))
mcols(gr)$extended.info <- c(
  list(c(NA)), 
  list(c(date = 2023, source = 'manual')), 
  list(c(date = 2021, source = 'manual')), 
  list(c(date = 2019, source = 'homology'))
)

mcols(gr)
# DataFrame with 4 rows and 3 columns
#       GC annotation extended.info
#     <numeric>   <factor>        <list>
# 1      0.45 NA                    NA
# 2      0.43 promoter     2023,manual
# 3      0.44 enhancer     2021,manual
# 4      0.42 centromere 2019,homology

## Genomic arithmetic on GRanges objects

gr
# GRanges object with 4 ranges and 3 metadata columns:
#   seqnames    ranges strand |        GC annotation extended.info
#         <Rle> <IRanges>  <Rle> | <numeric>   <factor>        <list>
# [1]     chr2 2004-7853      + |      0.45 NA                  <NA>
# [2]     chr4 4482-9873      - |      0.43 promoter     2023,manual
# [3]     chr5 1943-4203      + |      0.44 enhancer     2021,manual
# [4]     chr5 4103-5004      + |      0.42 centromere 2019,homology
# -------
#   seqinfo: 3 sequences from an unspecified genome; no seqlengths

#=== shift all genomic ranges towards the "right" (downstream in '+' strand) by 1000bp
shift(gr, 1000)
# GRanges object with 4 ranges and 3 metadata columns:
#   seqnames     ranges strand |        GC annotation extended.info
#         <Rle>  <IRanges>  <Rle> | <numeric>   <factor>        <list>
# [1]     chr2  3004-8853      + |      0.45 NA                  <NA>
# [2]     chr4 5482-10873      - |      0.43 promoter     2023,manual
# [3]     chr5  2943-5203      + |      0.44 enhancer     2021,manual
# [4]     chr5  5103-6004      + |      0.42 centromere 2019,homology
# -------
#   seqinfo: 3 sequences from an unspecified genome; no seqlengths

#=== now shift "left" (upstream in '+' strand) by 1000bp
shift(gr, -1000)
# GRanges object with 4 ranges and 3 metadata columns:
#   seqnames    ranges strand |        GC annotation extended.info
#         <Rle> <IRanges>  <Rle> | <numeric>   <factor>        <list>
# [1]     chr2 1004-6853      + |      0.45 NA                  <NA>
# [2]     chr4 3482-8873      - |      0.43 promoter     2023,manual
# [3]     chr5  943-3203      + |      0.44 enhancer     2021,manual
# [4]     chr5 3103-4004      + |      0.42 centromere 2019,homology
# -------
#   seqinfo: 3 sequences from an unspecified genome; no seqlengths

# narrow each genomic range in GRanges object by x number of bases
narrow(gr, start = 21, end = 40)
# GRanges object with 4 ranges and 3 metadata columns:
#   seqnames    ranges strand |        GC annotation extended.info
# <Rle> <IRanges>  <Rle> | <numeric>   <factor>        <list>
# [1]     chr2 2024-2043      + |      0.45 NA                  <NA>
# [2]     chr4 4502-4521      - |      0.43 promoter     2023,manual
# [3]     chr5 1963-1982      + |      0.44 enhancer     2021,manual
# [4]     chr5 4123-4142      + |      0.42 centromere 2019,homology
# -------
#   seqinfo: 3 sequences from an unspecified genome; no seqlengths

width(narrow(gr, start = 21, end = 40))
# [1] 20 20 20 20

# resize each genomic range to certain number of bases
#=== resize 'gr' entires to 100, fixed at start of each range
resize(gr, 100, fix = "start")
# GRanges object with 4 ranges and 3 metadata columns:
#   seqnames    ranges strand |        GC annotation extended.info
#         <Rle> <IRanges>  <Rle> | <numeric>   <factor>        <list>
# [1]     chr2 2004-2103      + |      0.45 NA                  <NA>
# [2]     chr4 9774-9873      - |      0.43 promoter     2023,manual
# [3]     chr5 1943-2042      + |      0.44 enhancer     2021,manual
# [4]     chr5 4103-4202      + |      0.42 centromere 2019,homology
# -------
#   seqinfo: 3 sequences from an unspecified genome; no seqlengths

#=== resize 'gr' to 100, fixed at start and ignoring strand information
resize(gr, 100, fix = "start", ignore.strand = TRUE)
# GRanges object with 4 ranges and 3 metadata columns:
#   seqnames    ranges strand |        GC annotation extended.info
#         <Rle> <IRanges>  <Rle> | <numeric>   <factor>        <list>
# [1]     chr2 2004-2103      + |      0.45 NA                  <NA>
# [2]     chr4 4482-4581*      - |      0.43 promoter     2023,manual
# [3]     chr5 1943-2042      + |      0.44 enhancer     2021,manual
# [4]     chr5 4103-4202      + |      0.42 centromere 2019,homology
# -------
#   seqinfo: 3 sequences from an unspecified genome; no seqlengths

#=== Resize `gr` entries to 1 bp, fixed at the center of each range
resize(gr, 1, fix = "center")
# GRanges object with 4 ranges and 3 metadata columns:
#   seqnames    ranges strand |        GC annotation extended.info
#        <Rle> <IRanges>  <Rle> | <numeric>   <factor>        <list>
# [1]     chr2      4928      + |      0.45 NA                  <NA>
# [2]     chr4      7177      - |      0.43 promoter     2023,manual
# [3]     chr5      3073      + |      0.44 enhancer     2021,manual
# [4]     chr5      4553      + |      0.42 centromere 2019,homology
# -------
#   seqinfo: 3 sequences from an unspecified genome; no seqlengths

# extracting flanking coordinates for each entry in 'gr'

# === Extract 100bp UPSTREAM of each genomic range, according to range strandness
flank(gr, 100, start = TRUE)
# GRanges object with 4 ranges and 3 metadata columns:
#   seqnames    ranges strand |        GC annotation extended.info
#         <Rle> <IRanges>  <Rle> | <numeric>   <factor>        <list>
# [1]     chr2 1904-2003      + |      0.45 NA                  <NA>
# [2]     chr4 9874-9973      - |      0.43 promoter     2023,manual
# [3]     chr5 1843-1942      + |      0.44 enhancer     2021,manual
# [4]     chr5 4003-4102      + |      0.42 centromere 2019,homology
# -------
#   seqinfo: 3 sequences from an unspecified genome; no seqlengths

#=== Extract 1bp DOWNSTREAM of each genomic range, according to range strandness
flank(gr, 1, start = FALSE)
# GRanges object with 4 ranges and 3 metadata columns:
#   seqnames    ranges strand |        GC annotation extended.info
#         <Rle> <IRanges>  <Rle> | <numeric>   <factor>        <list>
# [1]     chr2      7854      + |      0.45 NA                  <NA>
# [2]     chr4      4481      - |      0.43 promoter     2023,manual
# [3]     chr5      4204      + |      0.44 enhancer     2021,manual
# [4]     chr5      5005      + |      0.42 centromere 2019,homology
# -------
#   seqinfo: 3 sequences from an unspecified genome; no seqlengths

# arithmetic operators on GRanges
gr + 500 #=== Extend each side of the `GRanges` by a given number of bases
# GRanges object with 4 ranges and 3 metadata columns:
#   seqnames     ranges strand |        GC annotation extended.info
# <Rle>  <IRanges>  <Rle> | <numeric>   <factor>        <list>
# [1]     chr2  1504-8353      + |      0.45 NA                  <NA>
# [2]     chr4 3982-10373      - |      0.43 promoter     2023,manual
# [3]     chr5  1443-4703      + |      0.44 enhancer     2021,manual
# [4]     chr5  3603-5504      + |      0.42 centromere 2019,homology
# -------
#   seqinfo: 3 sequences from an unspecified genome; no seqlengths

gr - 200 #=== Shrink each side of the `GRanges` by a given number of bases 
# GRanges object with 4 ranges and 3 metadata columns:
#   seqnames    ranges strand |        GC annotation extended.info
#         <Rle> <IRanges>  <Rle> | <numeric>   <factor>        <list>
# [1]     chr2 2204-7653      + |      0.45 NA                  <NA>
# [2]     chr4 4682-9673      - |      0.43 promoter     2023,manual
# [3]     chr5 2143-4003      + |      0.44 enhancer     2021,manual
# [4]     chr5 4303-4804      + |      0.42 centromere 2019,homology
# -------
#   seqinfo: 3 sequences from an unspecified genome; no seqlengths

gr * 1000 #=== Zoom in by a given factor (decreases `GRanges` width by the same factor)
# GRanges object with 4 ranges and 3 metadata columns:
#   seqnames    ranges strand |        GC annotation extended.info
#         <Rle> <IRanges>  <Rle> | <numeric>   <factor>        <list>
# [1]     chr2 4926-4930      + |      0.45 NA                  <NA>
# [2]     chr4 7175-7179      - |      0.43 promoter     2023,manual
# [3]     chr5 3072-3073      + |      0.44 enhancer     2021,manual
# [4]     chr5 4554-4553      + |      0.42 centromere 2019,homology
# -------
#   seqinfo: 3 sequences from an unspecified genome; no seqlengths

## Inter-range methods
# result of each function in this section depends on entire set of ranges in object

# compute "inverse" genomic ranges (i.e., ranges in-bewteen input ranges)
gr
# GRanges object with 4 ranges and 3 metadata columns:
#   seqnames    ranges strand |        GC annotation extended.info
#        <Rle> <IRanges>  <Rle> | <numeric>   <factor>        <list>
# [1]     chr2 2004-7853      + |      0.45 NA                  <NA>
# [2]     chr4 4482-9873      - |      0.43 promoter     2023,manual
# [3]     chr5 1943-4203      + |      0.44 enhancer     2021,manual
# [4]     chr5 4103-5004      + |      0.42 centromere 2019,homology
# -------
#   seqinfo: 3 sequences from an unspecified genome; no seqlengths

gaps(gr)
# GRanges object with 3 ranges and 0 metadata columns:
#   seqnames    ranges strand
# <Rle> <IRanges>  <Rle>
#   [1]     chr2    1-2003      +
#   [2]     chr4    1-4481      -
#   [3]     chr5    1-1942      +
#   -------
#   seqinfo: 3 sequences from an unspecified genome; no seqlengths

# find index of preceding/following/nearest genomic range for each entry
precede(gr)
# [1] NA NA NA NA
follow(gr)
# [1] NA NA NA NA
nearest(gr)
# [1] NA NA  4  3 #=== for [3] nearest range is [4] and vice versa

# computing a coverage over a genome
coverage(gr, weight = 'GC')
# RleList of length 3
# $chr2
# numeric-Rle of length 7853 with 2 runs
# Lengths: 2003 5850
# Values : 0.00 0.45
# 
# $chr4
# numeric-Rle of length 9873 with 2 runs
# Lengths: 4481 5392
# Values : 0.00 0.43
# 
# $chr5
# numeric-Rle of length 5004 with 4 runs
# Lengths: 1942 2160  101  801
# Values : 0.00 0.44 0.86 0.42

## Comparing multiple GRanges objects
# peaks represents dummy 8 Chromatin Immunoprecipitation Sequencing (ChIP-seq) peaks
peaks <- GRanges(c(
  'chr1:320-418',
  'chr1:512-567',
  'chr1:843-892',
  'chr1:1221-1317', 
  'chr1:1329-1372', 
  'chr1:1852-1909', 
  'chr1:2489-2532', 
  'chr1:2746-2790'
))
peaks
# GRanges object with 8 ranges and 0 metadata columns:
#   seqnames    ranges strand
#           <Rle> <IRanges>  <Rle>
#   [1]     chr1   320-418      *
#   [2]     chr1   512-567      *
#   [3]     chr1   843-892      *
#   [4]     chr1 1221-1317      *
#   [5]     chr1 1329-1372      *
#   [6]     chr1 1852-1909      *
#   [7]     chr1 2489-2532      *
#   [8]     chr1 2746-2790      *
#   -------
#   seqinfo: 1 sequence from an unspecified genome; no seqlengths

# TSSs represents dummy 3 gene promoters (+/- 10bp around the TSS)
genes <- GRanges(c(
  'chr1:358-1292:+',
  'chr1:1324-2343:+', 
  'chr1:2732-2751:+'
))
TSSs <- resize(genes, width = 1, fix = 'start') + 10
TSSs
# GRanges object with 3 ranges and 0 metadata columns:
#   seqnames    ranges strand
#          <Rle> <IRanges>  <Rle>
#   [1]     chr1   348-368      +
#   [2]     chr1 1314-1334      +
#   [3]     chr1 2722-2742      +
#   -------
#   seqinfo: 1 sequence from an unspecified genome; no seqlengths

# visually check which ChIP-seq peaks overlap with a TSS
library(ggplot2)
peaks$type <- 'peaks'
TSSs$type <- 'TSSs'
ggplot() + 
  ggbio::geom_rect(c(peaks, TSSs), aes(fill = type), facets = type~.) + 
  ggbio::theme_alignment() + 
  coord_fixed(ratio = 300)

## find overlap between two GRanges sets; between query and subject
?IRanges::`findOverlaps-methods`

ov <- findOverlaps(query = peaks, subject = TSSs)
ov
# Hits object with 3 hits and 0 metadata columns:
#    queryHits subjectHits
#     <integer>   <integer>
# [1]         1           1
# [2]         4           2
# [3]         5           2
# -------
#   queryLength: 8 / subjectLength: 3
#=== Hits output above indicates query (peak) overlaps with subject (TSS)
#=== query (peak) #1 overlaps with subject (TSS) #1
#=== query (peak) #4 overlaps with subject (TSS) #2
#=== query (peak) #5 overlaps with subject (TSS) #2

# subset overlaps between query and subject; keep only peaks overlapping a TSS
subsetByOverlaps(peaks, TSSs)
# GRanges object with 3 ranges and 1 metadata column:
#     seqnames    ranges strand |        type
#        <Rle> <IRanges>  <Rle> | <character>
# [1]     chr1   320-418      * |       peaks
# [2]     chr1 1221-1317      * |       peaks
# [3]     chr1 1329-1372      * |       peaks
# -------
#   seqinfo: 1 sequence from an unspecified genome; no seqlengths

# counting overlaps between query and a subject
# for each range in a query, count how many ranges in subject it overlaps with
countOverlaps(query = peaks, subject = TSSs)
# [1] 1 0 0 1 1 0 0 0

# swap query and subject to see how many peaks each TSS overlaps with
countOverlaps(query = TSSs, subject = peaks)
# [1] 1 2 0
# add above counts to original query object in 'n_peaks' column:
TSSs$n_peaks <- countOverlaps(query = TSSs, subject = peaks)
TSSs
# GRanges object with 3 ranges and 2 metadata columns:
#     seqnames    ranges strand |        type   n_peaks
#        <Rle> <IRanges>  <Rle> | <character> <integer>
# [1]     chr1   348-368      + |        TSSs         1
# [2]     chr1 1314-1334      + |        TSSs         2
# [3]     chr1 2722-2742      + |        TSSs         0
# -------
#   seqinfo: 1 sequence from an unspecified genome; no seqlengths

## Find nearest range from a subject for each range in a query
# peaks[8] and TSSs[3] are very near to each other, but not quite overlapping

nearest(peaks, TSSs)
# [1] 1 1 2 2 2 2 3 3

TSSs[nearest(peaks, TSSs)]
# GRanges object with 8 ranges and 2 metadata columns:
#     seqnames    ranges strand |        type   n_peaks
#        <Rle> <IRanges>  <Rle> | <character> <integer>
# [1]     chr1   348-368      + |        TSSs         1
# [2]     chr1   348-368      + |        TSSs         1
# [3]     chr1 1314-1334      + |        TSSs         2
# [4]     chr1 1314-1334      + |        TSSs         2
# [5]     chr1 1314-1334      + |        TSSs         2
# [6]     chr1 1314-1334      + |        TSSs         2
# [7]     chr1 2722-2742      + |        TSSs         0
# [8]     chr1 2722-2742      + |        TSSs         0
# -------
#   seqinfo: 1 sequence from an unspecified genome; no seqlengths

# calculate distance to nearest between ranges in a query and ranges in a subject
distanceToNearest(peaks, TSSs)
# Hits object with 8 hits and 1 metadata column:
#     queryHits subjectHits |  distance
#     <integer>   <integer> | <integer>
# [1]         1           1 |         0
# [2]         2           1 |       143
# [3]         3           2 |       421
# [4]         4           2 |         0
# [5]         5           2 |         0
# [6]         6           2 |       517
# [7]         7           3 |       189
# [8]         8           3 |         3 #==== worth considering an overlap
# -------
#   queryLength: 8 / subjectLength: 3

# add column to peaks GRanges object
peaks$distance_to_nearest_TSS <- mcols(distanceToNearest(peaks, TSSs))$distance
peaks
# GRanges object with 8 ranges and 2 metadata columns:
#     seqnames    ranges strand |        type distance_to_nearest_TSS
#        <Rle> <IRanges>  <Rle> | <character>               <integer>
# [1]     chr1   320-418      * |       peaks                       0
# [2]     chr1   512-567      * |       peaks                     143
# [3]     chr1   843-892      * |       peaks                     421
# [4]     chr1 1221-1317      * |       peaks                       0
# [5]     chr1 1329-1372      * |       peaks                       0
# [6]     chr1 1852-1909      * |       peaks                     517
# [7]     chr1 2489-2532      * |       peaks                     189
# [8]     chr1 2746-2790      * |       peaks                       3
# -------
#   seqinfo: 1 sequence from an unspecified genome; no seqlengths

## GInteractions Class ------------------------------------------------
library(InteractionSet)
# define two parallel GRanges objects (i.e. of same length)
gr_first <- GRanges(c(
  'chr1:1-100', 
  'chr1:1001-2000', 
  'chr1:5001-6000', 
  'chr1:8001-9000', 
  'chr1:7001-8000'  
))
gr_second <- GRanges(c(
  'chr1:1-100', 
  'chr1:3001-4000', 
  'chr1:8001-9000', 
  'chr1:7001-8000', 
  'chr2:13000-14000'  
))

# can "bind" two objects of same length (5) from above 
library(InteractionSet)
gi <- GInteractions(gr_first, gr_second)
gi
# GInteractions object with 5 interactions and 0 metadata columns:
#     seqnames1   ranges1     seqnames2     ranges2
#         <Rle> <IRanges>         <Rle>   <IRanges>
# [1]      chr1     1-100 ---      chr1       1-100
# [2]      chr1 1001-2000 ---      chr1   3001-4000
# [3]      chr1 5001-6000 ---      chr1   8001-9000
# [4]      chr1 8001-9000 ---      chr1   7001-8000
# [5]      chr1 7001-8000 ---      chr2 13000-14000
# -------
#   regions: 7 ranges and 0 metadata columns
# seqinfo: 2 sequences from an unspecified genome; no seqlengths

gi[1] # possible to have interactions joining two identical anchors
gi[4] # Not advised to have "first" end after "second" end along chromosome
gi[5] # inter-chromosomal (i.e. trans) interactions

## GInteractions specific slots include anchors and regions

# anchors of a single genomic interaction refer to the two ends of interaction
anchors(gi) # this extracts two sets of anchors ("first" and "second")
# $first
# GRanges object with 5 ranges and 0 metadata columns:
#       seqnames    ranges strand
#          <Rle> <IRanges>  <Rle>
#   [1]     chr1     1-100      *
#   [2]     chr1 1001-2000      *
#   [3]     chr1 5001-6000      *
#   [4]     chr1 8001-9000      *
#   [5]     chr1 7001-8000      *
#   -------
#   seqinfo: 2 sequences from an unspecified genome; no seqlengths
# 
# $second
# GRanges object with 5 ranges and 0 metadata columns:
#       seqnames      ranges strand
#          <Rle>   <IRanges>  <Rle>
#   [1]     chr1       1-100      *
#   [2]     chr1   3001-4000      *
#   [3]     chr1   8001-9000      *
#   [4]     chr1   7001-8000      *
#   [5]     chr2 13000-14000      *
#   -------
#   seqinfo: 2 sequences from an unspecified genome; no seqlengths

# can query for "first" or "second" set of anchors directly
anchors(gi, "first")
anchors(gi, "second")

## Regions
# regions function returns the regions associated with GInteractions object
regions(gi)
# GRanges object with 7 ranges and 0 metadata columns:
#   seqnames      ranges strand
# <Rle>   <IRanges>  <Rle>
#   [1]     chr1       1-100      *
#   [2]     chr1   1001-2000      *
#   [3]     chr1   3001-4000      *
#   [4]     chr1   5001-6000      *
#   [5]     chr1   7001-8000      *
#   [6]     chr1   8001-9000      *
#   [7]     chr2 13000-14000      *
#   -------
#   seqinfo: 2 sequences from an unspecified genome; no seqlengths

## GInteractions methods

# Metadata can be added directly to Ginteractions object
mcols(gi)
#  DataFrame with 5 rows and 0 columns
mcols(gi) <- data.frame(
  idx = seq(1, length(gi)),
  type = c("cis", "cis", "cis", "cis", "trans")
)
gi
# GInteractions object with 5 interactions and 2 metadata columns:
#     seqnames1   ranges1     seqnames2     ranges2 |       idx        type
#         <Rle> <IRanges>         <Rle>   <IRanges> | <integer> <character>
# [1]      chr1     1-100 ---      chr1       1-100 |         1         cis
# [2]      chr1 1001-2000 ---      chr1   3001-4000 |         2         cis
# [3]      chr1 5001-6000 ---      chr1   8001-9000 |         3         cis
# [4]      chr1 8001-9000 ---      chr1   7001-8000 |         4         cis
# [5]      chr1 7001-8000 ---      chr2 13000-14000 |         5       trans
# -------
#   regions: 7 ranges and 0 metadata columns
# seqinfo: 2 sequences from an unspecified genome; no seqlengths

# metadata columns can be added directly to regions of a GInteractions object***
regions(gi)
# GRanges object with 7 ranges and 0 metadata columns:
#   seqnames      ranges strand
# <Rle>   <IRanges>  <Rle>
#   [1]     chr1       1-100      *
#   [2]     chr1   1001-2000      *
#   [3]     chr1   3001-4000      *
#   [4]     chr1   5001-6000      *
#   [5]     chr1   7001-8000      *
#   [6]     chr1   8001-9000      *
#   [7]     chr2 13000-14000      *
#   -------
#   seqinfo: 2 sequences from an unspecified genome; no seqlengths

regions(gi)$binID <- seq_along(regions(gi))
regions(gi)$type <- c("P", "P", "P", "E", "E", "P", "P")
regions(gi)
# GRanges object with 7 ranges and 2 metadata columns:
#     seqnames      ranges strand |     binID        type
#        <Rle>   <IRanges>  <Rle> | <integer> <character>
# [1]     chr1       1-100      * |         1           P
# [2]     chr1   1001-2000      * |         2           P
# [3]     chr1   3001-4000      * |         3           P
# [4]     chr1   5001-6000      * |         4           E
# [5]     chr1   7001-8000      * |         5           E
# [6]     chr1   8001-9000      * |         6           P
# [7]     chr2 13000-14000      * |         7           P
# -------
#   seqinfo: 2 sequences from an unspecified genome; no seqlengths

## Sorting GInteractions
gi
# GInteractions object with 5 interactions and 2 metadata columns:
#     seqnames1   ranges1     seqnames2     ranges2 |       idx        type
#         <Rle> <IRanges>         <Rle>   <IRanges> | <integer> <character>
# [1]      chr1     1-100 ---      chr1       1-100 |         1         cis
# [2]      chr1 1001-2000 ---      chr1   3001-4000 |         2         cis
# [3]      chr1 5001-6000 ---      chr1   8001-9000 |         3         cis
# [4]      chr1 8001-9000 ---      chr1   7001-8000 |         4         cis
# [5]      chr1 7001-8000 ---      chr2 13000-14000 |         5       trans
# -------
#   regions: 7 ranges and 2 metadata columns
# seqinfo: 2 sequences from an unspecified genome; no seqlengths

sort(gi)
# GInteractions object with 5 interactions and 2 metadata columns:
#     seqnames1   ranges1     seqnames2     ranges2 |       idx        type
#         <Rle> <IRanges>         <Rle>   <IRanges> | <integer> <character>
# [1]      chr1     1-100 ---      chr1       1-100 |         1         cis
# [2]      chr1 1001-2000 ---      chr1   3001-4000 |         2         cis
# [3]      chr1 5001-6000 ---      chr1   8001-9000 |         3         cis
# [4]      chr1 7001-8000 ---      chr2 13000-14000 |         5       trans
# [5]      chr1 8001-9000 ---      chr1   7001-8000 |         4         cis
# -------
#   regions: 7 ranges and 2 metadata columns
# seqinfo: 2 sequences from an unspecified genome; no seqlengths

## Swapping GInteractions anchors
# "first" and "second" anchors can be sorted for individual interactions
# a.k.a. "pairs swapping"
swapAnchors(gi)
# GInteractions object with 5 interactions and 2 metadata columns:
#     seqnames1   ranges1     seqnames2     ranges2 |       idx        type
#         <Rle> <IRanges>         <Rle>   <IRanges> | <integer> <character>
# [1]      chr1     1-100 ---      chr1       1-100 |         1         cis
# [2]      chr1 1001-2000 ---      chr1   3001-4000 |         2         cis
# [3]      chr1 5001-6000 ---      chr1   8001-9000 |         3         cis
# [4]      chr1 7001-8000 ---      chr1   8001-9000 |         4         cis
# [5]      chr1 7001-8000 ---      chr2 13000-14000 |         5       trans
# -------
#   regions: 7 ranges and 2 metadata columns
# seqinfo: 2 sequences from an unspecified genome; no seqlengths

#=== "sorting" reorganizes all rows (interactions)
#=== "swapping" anchors reorganizes "first" and "second" for each interaction independently

## GInteractions distance method
# "Distance" for genomic interactions refers to genomic distance bewteen anchors
# of a single interaction
gi
# GInteractions object with 5 interactions and 2 metadata columns:
#     seqnames1   ranges1     seqnames2     ranges2 |       idx        type
#         <Rle> <IRanges>         <Rle>   <IRanges> | <integer> <character>
# [1]      chr1     1-100 ---      chr1       1-100 |         1         cis
# [2]      chr1 1001-2000 ---      chr1   3001-4000 |         2         cis
# [3]      chr1 5001-6000 ---      chr1   8001-9000 |         3         cis
# [4]      chr1 8001-9000 ---      chr1   7001-8000 |         4         cis
# [5]      chr1 7001-8000 ---      chr2 13000-14000 |         5       trans
# -------
#   regions: 7 ranges and 2 metadata columns
# seqinfo: 2 sequences from an unspecified genome; no seqlengths
pairdist(gi)
# [1]    0 2000 3000 1000   NA
#=== for "trans" inter-chromosomal interactions, or interactions with anchors on 
# different chromosomes, notion of genomic distance is meaningless, hence NA

?`Interaction-overlaps`
## GInteractions overlap methods
# "overlaps" for genomic interactions computed for following contexts:
# 1. overlap between any of the two anchors of an interaction with a genomic range
# 2. overlap between anchors of an interaction with anchors of another interaction
# 3. spanning of the interaction "across" a genomic range

# 1. overlap between any of the two anchors of an interaction with a genomic range
gr <- GRanges(c("chr1:7501-7600", "chr1:8501-8600"))
findOverlaps(query = gi, subject = gr)
# Hits object with 4 hits and 0 metadata columns:
#      queryHits subjectHits
#      <integer>   <integer>
# [1]         3           2
# [2]         4           1
# [3]         4           2
# [4]         5           1
# -------
#   queryLength: 5 / subjectLength: 2

countOverlaps(gi, gr)
# [1] 0 0 1 2 1

subsetByOverlaps(gi, gr)
# GInteractions object with 3 interactions and 2 metadata columns:
#     seqnames1   ranges1     seqnames2     ranges2 |       idx        type
#         <Rle> <IRanges>         <Rle>   <IRanges> | <integer> <character>
# [1]      chr1 5001-6000 ---      chr1   8001-9000 |         3         cis
# [2]      chr1 8001-9000 ---      chr1   7001-8000 |         4         cis
# [3]      chr1 7001-8000 ---      chr2 13000-14000 |         5       trans
# -------
#   regions: 7 ranges and 2 metadata columns
# seqinfo: 2 sequences from an unspecified genome; no seqlengths
#=== the order matters here
countOverlaps(gr, gi)
# [1] 2 2
#=== %over% operator can be used here
gi %over% gr
# [1] FALSE FALSE  TRUE  TRUE  TRUE
gi[gi %over% gr]#=== equivalent to `subsetByOverlaps(gi, gr)`
# GInteractions object with 3 interactions and 2 metadata columns:
#     seqnames1   ranges1     seqnames2     ranges2 |       idx        type
#         <Rle> <IRanges>         <Rle>   <IRanges> | <integer> <character>
# [1]      chr1 5001-6000 ---      chr1   8001-9000 |         3         cis
# [2]      chr1 8001-9000 ---      chr1   7001-8000 |         4         cis
# [3]      chr1 7001-8000 ---      chr2 13000-14000 |         5       trans
# -------
#   regions: 7 ranges and 2 metadata columns
# seqinfo: 2 sequences from an unspecified genome; no seqlengths

# 2. overlap between anchors of an interaction with anchors of another interaction
gi2 <- GInteractions(
  GRanges("chr1:1081-1090"), 
  GRanges("chr1:3401-3501")
)
gi %over% gi2
# [1] FALSE  TRUE FALSE FALSE FALSE
# above: one interaction in Set-1 has its two anchors overlapping anchors from an interaction in Set-2

gi3 <- GInteractions(
  GRanges("chr1:1-1000"), 
  GRanges("chr1:3401-3501")
)
gi %over% gi3
# [1] FALSE FALSE FALSE FALSE FALSE
# single interactions - anchors from query need to overlap to pair of anchors from subject

# 3. spanning of the interaction "across" a genomic range
gi <- swapAnchors(gi) #=== Make sure anchors are correctly sorted
gi <- sort(gi) #=== Make sure interactions are correctly sorted
gi <- gi[!is.na(pairdist(gi))] #=== Remove inter-chromosomal interactions
spanning_gi <- GRanges(
  seqnames = seqnames(anchors(gi)[[1]]), 
  ranges = IRanges(
    start(anchors(gi)[[1]]), 
    end(anchors(gi)[[2]])
  )
)
spanning_gi 

# GRanges object with 4 ranges and 0 metadata columns:
#        seqnames    ranges strand
#          <Rle> <IRanges>  <Rle>
#   [1]     chr1     1-100      *
#   [2]     chr1 1001-4000      *
#   [3]     chr1 5001-9000      *
#   [4]     chr1 7001-9000      *
#   -------
#   seqinfo: 2 sequences from an unspecified genome; no seqlengths
spanning_gi %over% gr
# [1] FALSE FALSE  TRUE  TRUE

## Accessing example Hi-C files
# download example contact files
library(HiContactsData)
coolf <- HiContactsData('yeast_wt', 'mcool')
# example files for other formats
hicf <- HiContactsData('yeast_wt', 'hic')
hicpromatrixf <- HiContactsData('yeast_wt', 'hicpro_matrix')
hicproregionsf <- HiContactsData('yeast_wt', 'hicpro_bed')
pairsf <- HiContactsData('yeast_wt', 'pairs.gz')

# verify content of some of the files
#=== HiC-Pro generates a tab-separated `regions.bed` file
readLines(hicproregionsf, 25)
# [1] "I\t0\t1000"      "I\t1000\t2000"   "I\t2000\t3000"   "I\t3000\t4000"   "I\t4000\t5000"  
# [6] "I\t5000\t6000"   "I\t6000\t7000"   "I\t7000\t8000"   "I\t8000\t9000"   "I\t9000\t10000" 
# [11] "I\t10000\t11000" "I\t11000\t12000" "I\t12000\t13000" "I\t13000\t14000" "I\t14000\t15000"
# [16] "I\t15000\t16000" "I\t16000\t17000" "I\t17000\t18000" "I\t18000\t19000" "I\t19000\t20000"
# [21] "I\t20000\t21000" "I\t21000\t22000" "I\t22000\t23000" "I\t23000\t24000" "I\t24000\t25000"

#=== we can see that pairs are tab-separated
readLines(pairsf, 25)
# [1] "## pairs format v1.0"                                                             
# [2] "#sorted: chr1-pos1-chr2-pos2"                                                     
# [3] "#columns: readID chr1 pos1 chr2 pos2 strand1 strand2 frag1 frag2"                 
# [4] "#chromsize: I 230218"                                                             
# [5] "#chromsize: II 813184"                                                            
# [6] "#chromsize: III 316620"                                                           
# [7] "#chromsize: IV 1531933"                                                           
# [8] "#chromsize: V 576874"                                                             
# [9] "#chromsize: VI 270161"                                                            
# [10] "#chromsize: VII 1090940"                                                          
# [11] "#chromsize: VIII 562643"                                                          
# [12] "#chromsize: IX 439888"                                                            
# [13] "#chromsize: X 745751"                                                             
# [14] "#chromsize: XI 666816"                                                            
# [15] "#chromsize: XII 1078177"                                                          
# [16] "#chromsize: XIII 924431"                                                          
# [17] "#chromsize: XIV 784333"                                                           
# [18] "#chromsize: XV 1091291"                                                           
# [19] "#chromsize: XVI 948066"                                                           
# [20] "#chromsize: Mito 85779"                                                           
# [21] "NS500150:527:HHGYNBGXF:3:21611:19085:3986\tII\t105\tII\t48548\t+\t-\t1358\t1681"  
# [22] "NS500150:527:HHGYNBGXF:4:13604:19734:2406\tII\t113\tII\t45003\t-\t+\t1358\t1658"  
# [23] "NS500150:527:HHGYNBGXF:2:11108:25178:11036\tII\t119\tII\t687251\t-\t+\t1358\t5550"
# [24] "NS500150:527:HHGYNBGXF:1:22301:8468:1586\tII\t160\tII\t26124\t+\t-\t1358\t1510"   
# [25] "NS500150:527:HHGYNBGXF:4:23606:24037:2076\tII\t169\tII\t39052\t+\t+\t1358\t1613"

## ContactFile fudnamentals
# ContactFile object establishes connection with disk-stored Hi-C file (e.g. .cool or .pairs)
# these are only connections - they do not contain actual data
library(HiCExperiment)

CoolFile(coolf) # creates connection to .(m)cool file
HicFile(hicf) # creates connection to .hic file
HicproFile(hicpromatrixf, hicproregionsf) # creates connection to output files from HiC-Pro
PairsFile(pairsf) # creates connection to pairs file

## ContactFile slots (i.e. pieces of info)
cf <- CoolFile(coolf) # acce
cf
# access individual slots
resolution(cf)
pairsFile(cf)
metadata(cf)

## ContactFile methods
availableResolutions(cf)
# resolutions(5): 1000 2000 4000 8000 16000
availableChromosomes(cf)
# Seqinfo object with 16 sequences from an unspecified genome:
#   seqnames seqlengths isCircular genome
# I            230218       <NA>   <NA>
#   II           813184       <NA>   <NA>
#   III          316620       <NA>   <NA>
#   IV          1531933       <NA>   <NA>
#   V            576874       <NA>   <NA>
#   ...             ...        ...    ...
# XII         1078177       <NA>   <NA>
#   XIII         924431       <NA>   <NA>
#   XIV          784333       <NA>   <NA>
#   XV          1091291       <NA>   <NA>
#   XVI          948066       <NA>   <NA>

## HiCExperiment Class
# GInteractions class represent genomic interactions
# ContactFile class establishes connection with disk-stored Hi-C files

# use import method to create a HiCExperiment Object from a ContactFile
cf <- CoolFile(coolf)
hic <- import(cf)
hic # gives summary of data stored in object, not actual data stored
# `HiCExperiment` object with 8,757,906 contacts over 12,079 regions 
# -------
#   fileName: "/home/scatsac/.cache/R/ExperimentHub/283eb5e62932f_7752" 
# focus: "whole genome" 
# resolutions(5): 1000 2000 4000 8000 16000
# active resolution: 1000 
# interactions: 2945692 
# scores(2): count balanced 
# topologicalFeatures: compartments(0) borders(0) loops(0) viewpoints(0) 
# pairsFile: N/A 
# metadata(0):

# import works for other types of ContactFile (HicFile, HicproFile, PairsFile)
hf <- HicFile(hicf)
hic <- import(hf)
hic
# `HiCExperiment` object with 13,681,280 contacts over 12,165 regions 
# -------
#   fileName: "/home/scatsac/.cache/R/ExperimentHub/283eb7e4ae486_7836" 
# focus: "whole genome" 
# resolutions(5): 1000 2000 4000 8000 16000
# active resolution: 1000 
# interactions: 2965693 
# scores(2): count balanced 
# topologicalFeatures: compartments(0) borders(0) loops(0) viewpoints(0) 
# pairsFile: N/A 
# metadata(0):

# PairsFile returns an object that is a representation of Hi-C "pairs" (i.e., GInteractions)
pf <- PairsFile(pairsf)
pairs <- import(pf)
pairs

# GInteractions object with 471364 interactions and 3 metadata columns:
#     seqnames1   ranges1     seqnames2   ranges2 |     frag1     frag2  distance
#         <Rle> <IRanges>         <Rle> <IRanges> | <numeric> <numeric> <integer>
# [1]        II       105 ---        II     48548 |      1358      1681     48443
# [2]        II       113 ---        II     45003 |      1358      1658     44890
# [3]        II       119 ---        II    687251 |      1358      5550    687132
# [4]        II       160 ---        II     26124 |      1358      1510     25964
# [5]        II       169 ---        II     39052 |      1358      1613     38883
# ...       ...       ... ...       ...       ... .       ...       ...       ...
# [471360]        II    808605 ---        II    809683 |      6316      6320      1078
# [471361]        II    808609 ---        II    809917 |      6316      6324      1308
# [471362]        II    808617 ---        II    809506 |      6316      6319       889
# [471363]        II    809447 ---        II    809685 |      6319      6321       238
# [471364]        II    809472 ---        II    809675 |      6319      6320       203
# -------
#   regions: 549331 ranges and 0 metadata columns
# seqinfo: 1 sequence from an unspecified genome; no seqlengths

## Customize import to parse relevant data to the study
# focus: used to parse data for specific genome location
# resolution: choose resolution to parase contact matrix at
hic <- import(cf, focus = 'II', resolution = 2000)
regions(hic)
# GRanges object with 407 ranges and 4 metadata columns:
#   seqnames        ranges strand |    bin_id    weight   chr    center
#                  <Rle>     <IRanges>  <Rle> | <numeric> <numeric> <Rle> <integer>
# II_1_2000           II        1-2000      * |       116       NaN    II      1000
# II_2001_4000        II     2001-4000      * |       117       NaN    II      3000
# II_4001_6000        II     4001-6000      * |       118       NaN    II      5000
# II_6001_8000        II     6001-8000      * |       119       NaN    II      7000
# II_8001_10000       II    8001-10000      * |       120 0.0461112    II      9000
# ...      ...           ...    ... .       ...       ...   ...       ...
# II_804001_806000       II 804001-806000      * |       518 0.0493107    II    805000
# II_806001_808000       II 806001-808000      * |       519 0.0611355    II    807000
# II_808001_810000       II 808001-810000      * |       520       NaN    II    809000
# II_810001_812000       II 810001-812000      * |       521       NaN    II    811000
# II_812001_813184       II 812001-813184      * |       522       NaN    II    812592
# -------
#   seqinfo: 16 sequences from an unspecified genome

## Interacting with HiCExperiment data
# HiCExperiment object wraps together a ContactFile and GIinteractions
yeast_hic <- contacts_yeast(full = TRUE)
yeast_hic
# `HiCExperiment` object with 8,757,906 contacts over 763 regions 
# -------
#   fileName: "/home/scatsac/.cache/R/ExperimentHub/283eb5e62932f_7752" 
# focus: "whole genome" 
# resolutions(5): 1000 2000 4000 8000 16000
# active resolution: 16000 
# interactions: 267709 
# scores(2): count balanced 
# topologicalFeatures: compartments(0) borders(0) loops(0) viewpoints(0) centromeres(16) 
# pairsFile: N/A 
# metadata(0):

# imported interactions can be directly exposed using interactions function
interactions(yeast_hic)
# GInteractions object with 267709 interactions and 4 metadata columns:
#     seqnames1       ranges1     seqnames2       ranges2 |   bin_id1   bin_id2     count  balanced
#         <Rle>     <IRanges>         <Rle>     <IRanges> | <numeric> <numeric> <numeric> <numeric>
# [1]         I       1-16000 ---         I       1-16000 |         0         0      2836 1.0943959
# [2]         I       1-16000 ---         I   16001-32000 |         0         1      2212 0.9592069
# [3]         I       1-16000 ---         I   32001-48000 |         0         2      1183 0.4385242
# [4]         I       1-16000 ---         I   48001-64000 |         0         3       831 0.2231192
# [5]         I       1-16000 ---         I   64001-80000 |         0         4       310 0.0821255
# ...       ...           ... ...       ...           ... .       ...       ...       ...       ...
# [267705]       XVI 896001-912000 ---       XVI 912001-928000 |       759       760      3565  1.236371
# [267706]       XVI 896001-912000 ---       XVI 928001-944000 |       759       761      1359  0.385016
# [267707]       XVI 912001-928000 ---       XVI 912001-928000 |       760       760      3534  2.103988
# [267708]       XVI 912001-928000 ---       XVI 928001-944000 |       760       761      3055  1.485794
# [267709]       XVI 928001-944000 ---       XVI 928001-944000 |       761       761      4308  1.711565
# -------
#   regions: 763 ranges and 4 metadata columns
# seqinfo: 16 sequences from an unspecified genome

# regions and anchors work on HiCExperiment objects similar to GInteractions
regions(yeast_hic)
# GRanges object with 763 ranges and 4 metadata columns:
#   seqnames        ranges strand |    bin_id     weight   chr    center
# <Rle>     <IRanges>  <Rle> | <numeric>  <numeric> <Rle> <integer>
#   I_1_16000        I       1-16000      * |         0  0.0196442     I      8000
# I_16001_32000        I   16001-32000      * |         1  0.0220746     I     24000
# I_32001_48000        I   32001-48000      * |         2  0.0188701     I     40000
# I_48001_64000        I   48001-64000      * |         3  0.0136679     I     56000
# I_64001_80000        I   64001-80000      * |         4  0.0134860     I     72000
# ...      ...           ...    ... .       ...        ...   ...       ...
# XVI_880001_896000      XVI 880001-896000      * |       758 0.00910873   XVI    888000
# XVI_896001_912000      XVI 896001-912000      * |       759 0.01421350   XVI    904000
# XVI_912001_928000      XVI 912001-928000      * |       760 0.02439992   XVI    920000
# XVI_928001_944000      XVI 928001-944000      * |       761 0.01993237   XVI    936000
# XVI_944001_948066      XVI 944001-948066      * |       762        NaN   XVI    946033
# -------
#   seqinfo: 16 sequences from an unspecified genome

anchors(yeast_hic)
# $first
# GRanges object with 267709 ranges and 4 metadata columns:
#   seqnames        ranges strand |    bin_id    weight   chr    center
# <Rle>     <IRanges>  <Rle> | <numeric> <numeric> <Rle> <integer>
#   [1]        I       1-16000      * |         0 0.0196442     I      8000
# [2]        I       1-16000      * |         0 0.0196442     I      8000
# [3]        I       1-16000      * |         0 0.0196442     I      8000
# [4]        I       1-16000      * |         0 0.0196442     I      8000
# [5]        I       1-16000      * |         0 0.0196442     I      8000
# ...      ...           ...    ... .       ...       ...   ...       ...
# [267705]      XVI 896001-912000      * |       759 0.0142135   XVI    904000
# [267706]      XVI 896001-912000      * |       759 0.0142135   XVI    904000
# [267707]      XVI 912001-928000      * |       760 0.0243999   XVI    920000
# [267708]      XVI 912001-928000      * |       760 0.0243999   XVI    920000
# [267709]      XVI 928001-944000      * |       761 0.0199324   XVI    936000
# -------
#   seqinfo: 16 sequences from an unspecified genome
# 
# $second
# GRanges object with 267709 ranges and 4 metadata columns:
#   seqnames        ranges strand |    bin_id    weight   chr    center
# <Rle>     <IRanges>  <Rle> | <numeric> <numeric> <Rle> <integer>
#   [1]        I       1-16000      * |         0 0.0196442     I      8000
# [2]        I   16001-32000      * |         1 0.0220746     I     24000
# [3]        I   32001-48000      * |         2 0.0188701     I     40000
# [4]        I   48001-64000      * |         3 0.0136679     I     56000
# [5]        I   64001-80000      * |         4 0.0134860     I     72000
# ...      ...           ...    ... .       ...       ...   ...       ...
# [267705]      XVI 912001-928000      * |       760 0.0243999   XVI    920000
# [267706]      XVI 928001-944000      * |       761 0.0199324   XVI    936000
# [267707]      XVI 912001-928000      * |       760 0.0243999   XVI    920000
# [267708]      XVI 928001-944000      * |       761 0.0199324   XVI    936000
# [267709]      XVI 928001-944000      * |       761 0.0199324   XVI    936000
# -------
#   seqinfo: 16 sequences from an unspecified genome

## Bins and seqinfo
seqinfo(yeast_hic)
# Seqinfo object with 16 sequences from an unspecified genome:
#   seqnames seqlengths isCircular genome
# I            230218       <NA>   <NA>
#   II           813184       <NA>   <NA>
#   III          316620       <NA>   <NA>
#   IV          1531933       <NA>   <NA>
#   V            576874       <NA>   <NA>
#   ...             ...        ...    ...
# XII         1078177       <NA>   <NA>
#   XIII         924431       <NA>   <NA>
#   XIV          784333       <NA>   <NA>
#   XV          1091291       <NA>   <NA>
#   XVI          948066       <NA>   <NA>

bins(yeast_hic)
# GRanges object with 763 ranges and 2 metadata columns:
#   seqnames        ranges strand |    bin_id     weight
# <Rle>     <IRanges>  <Rle> | <numeric>  <numeric>
#   I_1_16000        I       1-16000      * |         0  0.0196442
# I_16001_32000        I   16001-32000      * |         1  0.0220746
# I_32001_48000        I   32001-48000      * |         2  0.0188701
# I_48001_64000        I   48001-64000      * |         3  0.0136679
# I_64001_80000        I   64001-80000      * |         4  0.0134860
# ...      ...           ...    ... .       ...        ...
# XVI_880001_896000      XVI 880001-896000      * |       758 0.00910873
# XVI_896001_912000      XVI 896001-912000      * |       759 0.01421350
# XVI_912001_928000      XVI 912001-928000      * |       760 0.02439992
# XVI_928001_944000      XVI 928001-944000      * |       761 0.01993237
# XVI_944001_948066      XVI 944001-948066      * |       762        NaN
# -------
#   seqinfo: 16 sequences from an unspecified genome

#=== bins are not equivalent to regions of HiCExperiment
# bins refer to all possible regions of a HiCExperiment

## Scores - the frequency for each genomic interaction
head(scores(yeast_hic))
# List of length 2
# names(2): count balanced

head(scores(yeast_hic, "count"))
# [1] 2836 2212 1183  831  310  159

head(scores(yeast_hic, "balanced"))
# [1] 1.09439586 0.95920688 0.43852417 0.22311917 0.08212549 0.03345221

# calling interactions(hic) returns GInteractions with scores
interactions(yeast_hic)
# GInteractions object with 267709 interactions and 4 metadata columns:
#     seqnames1       ranges1     seqnames2       ranges2 |   bin_id1   bin_id2     count  balanced
#         <Rle>     <IRanges>         <Rle>     <IRanges> | <numeric> <numeric> <numeric> <numeric>
# [1]         I       1-16000 ---         I       1-16000 |         0         0      2836 1.0943959
# [2]         I       1-16000 ---         I   16001-32000 |         0         1      2212 0.9592069
# [3]         I       1-16000 ---         I   32001-48000 |         0         2      1183 0.4385242
# [4]         I       1-16000 ---         I   48001-64000 |         0         3       831 0.2231192
# [5]         I       1-16000 ---         I   64001-80000 |         0         4       310 0.0821255
# ...       ...           ... ...       ...           ... .       ...       ...       ...       ...
# [267705]       XVI 896001-912000 ---       XVI 912001-928000 |       759       760      3565  1.236371
# [267706]       XVI 896001-912000 ---       XVI 928001-944000 |       759       761      1359  0.385016
# [267707]       XVI 912001-928000 ---       XVI 912001-928000 |       760       760      3534  2.103988
# [267708]       XVI 912001-928000 ---       XVI 928001-944000 |       760       761      3055  1.485794
# [267709]       XVI 928001-944000 ---       XVI 928001-944000 |       761       761      4308  1.711565
# -------
#   regions: 763 ranges and 4 metadata columns
# seqinfo: 16 sequences from an unspecified genome

## topologicalFeatures - identified genomic structures
# by default, the four empty topologicalFeatures are:
# compartments, borders, loops, viewpoints

topologicalFeatures(yeast_hic)
# names(5): compartments borders loops viewpoints centromeres

topologicalFeatures(yeast_hic, 'centromeres')
# GRanges object with 16 ranges and 0 metadata columns:
#   seqnames        ranges strand
# <Rle>     <IRanges>  <Rle>
#   [1]        I 151583-151641      +
#   [2]       II 238361-238419      +
#   [3]      III 114322-114380      +
#   [4]       IV 449879-449937      +
#   [5]        V 152522-152580      +
#   ...      ...           ...    ...
# [12]      XII 151366-151424      +
#   [13]     XIII 268222-268280      +
#   [14]      XIV 628588-628646      +
#   [15]       XV 326897-326955      +
#   [16]      XVI 556255-556313      +
#   -------
#   seqinfo: 17 sequences (1 circular) from R64-1-1 genome

## pairsFile
# PairsFile can be created and associated with corresponding HiCExperiment object
# this allows for more accurate estimation of contact distribution
pairsFile(yeast_hic) <- pairsf
pairsFile(yeast_hic)
readLines(pairsFile(yeast_hic), 25)
# [1] "## pairs format v1.0"                                                             
# [2] "#sorted: chr1-pos1-chr2-pos2"                                                     
# [3] "#columns: readID chr1 pos1 chr2 pos2 strand1 strand2 frag1 frag2"                 
# [4] "#chromsize: I 230218"                                                             
# [5] "#chromsize: II 813184"                                                            
# [6] "#chromsize: III 316620"                                                           
# [7] "#chromsize: IV 1531933"                                                           
# [8] "#chromsize: V 576874"                                                             
# [9] "#chromsize: VI 270161"                                                            
# [10] "#chromsize: VII 1090940"                                                          
# [11] "#chromsize: VIII 562643"                                                          
# [12] "#chromsize: IX 439888"                                                            
# [13] "#chromsize: X 745751"                                                             
# [14] "#chromsize: XI 666816"                                                            
# [15] "#chromsize: XII 1078177"                                                          
# [16] "#chromsize: XIII 924431"                                                          
# [17] "#chromsize: XIV 784333"                                                           
# [18] "#chromsize: XV 1091291"                                                           
# [19] "#chromsize: XVI 948066"                                                           
# [20] "#chromsize: Mito 85779"                                                           
# [21] "NS500150:527:HHGYNBGXF:3:21611:19085:3986\tII\t105\tII\t48548\t+\t-\t1358\t1681"  
# [22] "NS500150:527:HHGYNBGXF:4:13604:19734:2406\tII\t113\tII\t45003\t-\t+\t1358\t1658"  
# [23] "NS500150:527:HHGYNBGXF:2:11108:25178:11036\tII\t119\tII\t687251\t-\t+\t1358\t5550"
# [24] "NS500150:527:HHGYNBGXF:1:22301:8468:1586\tII\t160\tII\t26124\t+\t-\t1358\t1510"   
# [25] "NS500150:527:HHGYNBGXF:4:23606:24037:2076\tII\t169\tII\t39052\t+\t+\t1358\t1613"  

## Importing a PairsFile in a GInteractions object
import(pairsFile(yeast_hic), format = 'pairs')
# GInteractions object with 471364 interactions and 3 metadata columns:                                    
#     seqnames1   ranges1     seqnames2   ranges2 |     frag1     frag2  distance
#         <Rle> <IRanges>         <Rle> <IRanges> | <numeric> <numeric> <integer>
# [1]        II       105 ---        II     48548 |      1358      1681     48443
# [2]        II       113 ---        II     45003 |      1358      1658     44890
# [3]        II       119 ---        II    687251 |      1358      5550    687132
# [4]        II       160 ---        II     26124 |      1358      1510     25964
# [5]        II       169 ---        II     39052 |      1358      1613     38883
# ...       ...       ... ...       ...       ... .       ...       ...       ...
# [471360]        II    808605 ---        II    809683 |      6316      6320      1078
# [471361]        II    808609 ---        II    809917 |      6316      6324      1308
# [471362]        II    808617 ---        II    809506 |      6316      6319       889
# [471363]        II    809447 ---        II    809685 |      6319      6321       238
# [471364]        II    809472 ---        II    809675 |      6319      6320       203
# -------
#   regions: 549331 ranges and 0 metadata columns
# seqinfo: 1 sequence from an unspecified genome; no seqlengths

