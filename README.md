Plots and Scripts for: Many purported pseudogenes in bacterial genomes
are bonafide genes
================
Nicholas P. Cooley, Department of Biomedical Informatics, University of
Pittsburgh
2024-02-06

# Pseudogenes!

This github repo contains data and scripts necessary to recreate the
plots present in the manuscript **Many purported pseudogenes in
bacterial genomes are bonafide genes**. The data present in this repo
are mostly lightweight summary tables capable of fitting within github
size restrictions. Scripts used to generate large initial data sets on
the Open Science Grid are present in the `OSG_Jobs` folder, while the
data sets generated from those jobs have been deposited on zenodo under
the following DOIs:

- [10.5281/zenodo.8360505](https://zenodo.org/record/8360505) -
  Assemblies generated Figure 3 in this repo.
- [10.5281/zenodo.8361514](https://zenodo.org/record/8361514) -
  Annotations and parsed data for Figure 3 in this repo.
- [10.5281/zenodo.8356318](https://zenodo.org/record/8356318) - Summary
  data that is too large for github, but necessary to knit this README
  and it’s associated figures, enumerated in the `.gitignore` file for
  this repo.
- [10.5281/zenodo.8366931](https://zenodo.org/record/8366931) -
  Assemblies generated for Figure 4 in this repo.
- [10.5281/zenodo.8378433](https://zenodo.org/record/8378433) -
  Annotations and parsed data for Figure 4 in this repo, as well as
  assemblies annotations and parsed data for Figure 5 in this repo.
- [10.5281/zenodo.10621232](https://zenodo.org/records/10621233) -
  Assemblies, annotations and comparison data used in the generation of
  figure 6.
- [10.5281/zenodo.10622275](https://zenodo.org/records/10622276) -
  Assemblies and raw assembly data, part 1.
- [10.5281/zenodo.10625339](https://zenodo.org/records/10625340) -
  Assemblies and raw assembly data, part 2.

Pseudogenizations due to internal stops or frameshifts can represent one
of at least three separate phenomena in assembled genomes, 1) recent
evolutionary changes that can serve as an observational marker of how
pressure is affecting functions and tools within a genome, 2) an error
introduced into an assembly via error modes inherent to the sequencing
platform or the assembly process, or 3) an inaccurate annotation of a
programmed frameshift or non-canonical amino acid inclusion in lieu of a
stop codon. Without confirmation such as Sanger sequencing, it can be
unclear which of these options any individual pseudogene actually
represents. The wide variety of platform and assembler choices available
to data submitters additionally introduces the possibility for
stochasticity in the rates at which pseudogenes are `TRUE` or `FALSE`
depending on the combination of choices made in data collection and
generation.

It would take an enormous effort to wholesale Sanger sequence even a
modest number of the pseudogenes present in RefSeq or Genbank. It is not
even clear that that type of experiment is necessary. However, some
interrogations of the diverse data present in RefSeq, GenBank, and the
SRA are possible and potentially useful. Metadata can be scraped from
the SRA and we can generate direct observations of how relative counts
of pseudogenes are related to extractable pieces of data, such as
reported assembler, platforms (sequencing technology) for available SRA
runs, submission year, reported assembly status, contig N50 over total
length, and genus. These direct observations can be coupled with causal
inference via Tetrad to predict causal links between metadata categories
and relative pseudogene counts.

We can additionally reassemble available reads under a variety of
conditions, and assemble reads simulated under varying coverages and
qualities to interrogate factors that can affect the pseudogene content
of finished assemblies.

PNGs of the manuscript figures are embedded in this document below,
while better quality PDFs are included in the the
`README_files/figure-gfm/` folder of this repo.

### Figure 1:

A reasonable *a priori* expectation is that two assemblies that
represent only very recently diverged genomes should contain very
similar numbers of pseudogenes. We can plot reported count differences
between pseudogenes by type against ANI and show that there are cases of
very closely related genomes with considerably different pseudogene
repertoires by count.

    # 
    # frameshift spline fit at 100 == 0.215859503163376 
    #  internal stop split fit at 100 == 0.0822019556556262

<div class="figure" style="text-align: center">

<img src="README_files/figure-gfm/pg-incongruency-1.png" alt="Pseudogenes often show orthologous relationships with non-pseudogenes"  />
<p class="caption">
Pseudogenes often show orthologous relationships with non-pseudogenes
</p>

</div>

We can calculate the number of very close pairs (pairs with an ANI \>=
99.9), where the rates of pseudogenization imply that one partner has
more pseudogenes than expected. The code snippet below shows how this
can be accomplished for both internal stops and frameshifts.

``` r
# from the data file: InputData/Counts_Orthos_v02.RData

# match up the total coding counts to the pair partners
mat1 <- match(x = dat3$id1,
              table = as.integer(rownames(adjusted_counts)))
mat2 <- match(x = dat3$id2,
              table = as.integer(rownames(adjusted_counts)))

dat3$features1 <- adjusted_counts$all_coding[mat1]
dat3$features2 <- adjusted_counts$all_coding[mat2]

size <- ifelse(test = dat3$is1 > dat3$is2,
               yes = dat3$features1,
               no = dat3$features2) # features in assembly with more IS
rate <- ifelse(test = dat3$is1 < dat3$is2,
               yes = dat3$is1/dat3$features1,
               no = dat3$is2/dat3$features2) # rate of IS in assembly with fewer IS
num <- ifelse(test = dat3$is1 > dat3$is2,
              yes = dat3$is1,
              no = dat3$is2) # number of IS in assembly with more IS

reject <- num > qbinom(0.99, size, rate)
# sum(dat3$ANI >= 99.9)
cat("Rejections by Internal Stops == ")
```

    # Rejections by Internal Stops ==

``` r
mean(reject[dat3$ANI >= 99.9])
```

    # [1] 0.03647652

``` r
size <- ifelse(test = dat3$fs1 > dat3$fs2,
               yes = dat3$features1,
               no = dat3$features2) # features in assembly with more FS
rate <- ifelse(test = dat3$fs1 < dat3$fs2,
               yes = dat3$fs1/dat3$features1,
               no = dat3$fs2/dat3$features2) # rate of FS in assembly with fewer FS
num <- ifelse(test = dat3$fs1 > dat3$fs2,
              yes = dat3$fs1,
              no = dat3$fs2) # number of FS in assembly with more FS

reject <- num > qbinom(0.99, size, rate)
# sum(dat3$ANI >= 99.9)
cat("Rejections by Frameshifts == ")
```

    # Rejections by Frameshifts ==

``` r
mean(reject[dat3$ANI >= 99.9])
```

    # [1] 0.112534

### Figure 2:

The incomplete nature of the metadata available in public repositories
makes direct modeling from that metadata unwise. Causal inference can be
applied to that data however, and we can plot out how distributions of
relative pseudogene counts by label groups appear to be dissimilar.

<div class="figure" style="text-align: center">

<img src="README_files/figure-gfm/inference-and-observation-1.png" alt="Causal inference and observational distributions of pseudogenes"  />
<p class="caption">
Causal inference and observational distributions of pseudogenes
</p>

</div>

### Tetrad Table for Figure 2:

| p1          | p2              | weight | type |
|:------------|:----------------|-------:|:-----|
| genus       | partials        |   1.00 | –\>  |
| genus       | stops           |   1.00 | –\>  |
| submit year | assembler       |   1.00 | o-\> |
| submit year | genus           |   1.00 | o-\> |
| technology  | partials        |   1.00 | o-\> |
| technology  | stops           |   1.00 | o-\> |
| cov         | stops           |   0.62 | o-\> |
| cov         | frameshifts     |   0.37 | o-\> |
| cov         | assembly status |   0.01 | o-\> |
| technology  | frameshifts     |   0.01 | o-\> |

### Figure 3:

Reassembly of reads on the SRA from within a single genus using
different assemblers shows unique distributions as assembler is varied.

    # [1] "242 total source reads with completed reassemblies for all chosen assemblers."

<div class="figure" style="text-align: center">

<img src="README_files/figure-gfm/influence-of-assembler-1.png" alt="Within a species, different assemblers provide unique distributions of pseudogenes"  />
<p class="caption">
Within a species, different assemblers provide unique distributions of
pseudogenes
</p>

</div>

### Figure 4:

In cases where multiple sequencing runs are deposited on the SRA for a
given biosample, those reads can be individually reassembled and their
relative pseudogene counts compared. This requires some fairly stringent
metadata filtering during data selection and after assembly which are
enumerated in the manuscript. Those steps are taken to ensure that the
two generated assemblies actually represent the same information.

<div class="figure" style="text-align: center">

<img src="README_files/figure-gfm/matched-sra-run-deviation-1.png" alt="When multiple SRA runs for the same biosample are available deviations between run can be queried"  />
<p class="caption">
When multiple SRA runs for the same biosample are available deviations
between run can be queried
</p>

</div>

### Figure 5:

Reads can be simulated under a variety of conditions and the resulting
assemblies can be used to generate models for pseudogenization rates
under those conditions. Further visualizations of this data are present
in the supporting information.

<div class="figure" style="text-align: center">

<img src="README_files/figure-gfm/predicted-behavior-1.png" alt="Assemblies generated from simulated reads provide an opportunity to model pseudogenes as an outcome for coverage (quality not shown)"  />
<p class="caption">
Assemblies generated from simulated reads provide an opportunity to
model pseudogenes as an outcome for coverage (quality not shown)
</p>

</div>

    # 
    #  0.672571234596583 percent of biosamples in genbank have associated SRA reads

    # 
    #  0.989657046482021 percent of unique SRA biosamples have associated Illumina reads

### Figure 6:

Reassembly of reads on the SRA that imply extremely deep sequencing
coverage under stepwise subsampling can be used to mimic the data
generated for simulated reads.

<div class="figure" style="text-align: center">

<img src="README_files/figure-gfm/reviewer-request01-1.png" alt="Assemblies generated from subsamples of available real reads on the SRA show increases in pseudogenization as subsamples become more sparse"  />
<p class="caption">
Assemblies generated from subsamples of available real reads on the SRA
show increases in pseudogenization as subsamples become more sparse
</p>

</div>

``` r
# ls()
sessionInfo()
```

    # R version 4.3.1 (2023-06-16)
    # Platform: x86_64-apple-darwin20 (64-bit)
    # Running under: macOS Sonoma 14.1.2
    # 
    # Matrix products: default
    # BLAS:   /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRblas.0.dylib 
    # LAPACK: /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
    # 
    # locale:
    # [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    # 
    # time zone: America/New_York
    # tzcode source: internal
    # 
    # attached base packages:
    #  [1] grid      parallel  stats4    stats     graphics  grDevices utils    
    #  [8] datasets  methods   base     
    # 
    # other attached packages:
    #  [1] plotrix_3.8-4       igraph_1.5.1        VennDiagram_1.7.3  
    #  [4] futile.logger_1.4.3 pdftools_3.4.0      magick_2.8.1       
    #  [7] SynExtend_1.14.0    DECIPHER_2.30.0     RSQLite_2.3.2      
    # [10] Biostrings_2.70.1   GenomeInfoDb_1.38.0 XVector_0.42.0     
    # [13] IRanges_2.36.0      S4Vectors_0.40.1    BiocGenerics_0.48.0
    # [16] knitr_1.45         
    # 
    # loaded via a namespace (and not attached):
    #  [1] bit_4.0.5               highr_0.10              compiler_4.3.1         
    #  [4] qpdf_1.3.2              crayon_1.5.2            Rcpp_1.0.11            
    #  [7] blob_1.2.4              bitops_1.0-7            yaml_2.3.7             
    # [10] fastmap_1.1.1           GenomeInfoDbData_1.2.11 DBI_1.1.3              
    # [13] rlang_1.1.1             cachem_1.0.8            xfun_0.40              
    # [16] bit64_4.0.5             memoise_2.0.1           cli_3.6.1              
    # [19] formatR_1.14            magrittr_2.0.3          futile.options_1.0.1   
    # [22] zlibbioc_1.48.0         digest_0.6.33           rstudioapi_0.15.0      
    # [25] askpass_1.2.0           vctrs_0.6.4             evaluate_0.22          
    # [28] lambda.r_1.2.4          RCurl_1.98-1.12         rmarkdown_2.25         
    # [31] pkgconfig_2.0.3         tools_4.3.1             htmltools_0.5.6.1
