
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mrQC

<!-- badges: start -->

<!-- badges: end -->

mrQC is an R package for checking the accuracy of meta- and summary-data
for outcome traits prior to their use in Mendelian randomisation
analyses. For example, the package provides tools for checking that the
reported effect allele and effect allele frequency columns are correct.
It also checks for possible errors in the summary data that might
introduce bias into downstream Mendelian randomisation analyses.

## Installation

To install the latest version of mrQC, perform as normal:

``` r
install.packages("devtools")
devtools::install_github("MRCIEU/mrQC")
```

## General overview

This package exploits three groups of SNPs: 1) **A MAF 1KG reference
set**. This is a set of 2297 SNPs that have the same minor allele across
the 1000 genomes super populations and that have a minor allele
frequency between 0.1 and 0.3. These SNPs are used to check for allele
frequency errors; 2) **GWAS catalog top hits**. These are SNPs in the
GWAS catalog that are strongly associated with the outcome trait of
interest; and 3) the **exposure SNPs** or genetic instruments for the
exposure of interest.

The package performs 4 types of checks:

1)  confirm the identity of the effect allele frequency column. This
    step relies on comparison of the MAF reference set between the
    outcome study and 1000 genomes reference dataset.
2)  confirm the identity of the effect allele column. This step relies
    on comparison of GWAS hits between the outcome study and the GWAS
    catalog.  
3)  assess the ancestral origins of the study. This step relies on
    comparison of the MAF reference set between the outcome study and
    1000 genomes reference dataset.  
4)  for case-control studies, identify potential errors or issues in the
    effect sizes. This step can use any of the SNP sets.

## Step 1. Extract and format the outcome summary data

The first step is to make a list of SNP rsids. Let’s assume our outcome
trait of interest is glioma. Let’s identify the top GWAS hits for glioma
in the GWAS catalog. We will search using both the reported trait name
as well as the trait experimental factor ontology (EFO). It is advisable
to search for both because not all studies in the GWAS catalog have
up-to date EFO annotations. Searching for both therefore maximises the
number of retrieved hits.

``` r
library(mrQC)
snplist<-make_snplist(efo="glioma",trait="glioma",ref1000G_superpops=FALSE) #search on EFO abd reported trait 
head(snplist)
length(snplist)
```

We see that searching on the glioma EFO and glioma trait retrieves 54
unique SNP rsids. We could also have searched on the glioma efo\_id,
which retrieves the same SNPs as searching on the “glioma” efo.

``` r
snplist<-make_snplist(efo_id="EFO_0005543",ref1000G_superpops=FALSE) 
```

In the above examples we set ref1000G\_superpops to FALSE. We now set
this to TRUE, which will result in a SNP list that includes the “MAF
reference set”. This is a set of 2297 SNPs that have the same minor
allele across all 1000 genomes super populations and a minor allele
frequency of
0.1-0.3.

``` r
snplist<-make_snplist(efo="glioma",trait="glioma",ref1000G_superpops=TRUE) 
head(snplist)
length(snplist)
```

We can also define a set of “exposure SNPs” and include this in list of
rsids. For example, let’s assume we are conducting a Mendelian
randomisation study to assess the effect of polyunsaturated fatty acid
exposure on risk of glioma. Let’s define the “exposure SNPs”, or genetic
instrument, as SNPs associated with polyunsaturated fatty acids with a P
value \<5e-8. Let’s use the ieugwasr package to extract “exposure SNPs”
for polyunsaturated fatty acids from the Open GWAS project and also add
these to the other SNPs.

``` r
instruments<-ieugwasr::tophits(id="met-d-PUFA",pval = 5e-08)
snplist<-make_snplist(efo = "glioma",trait="glioma",ref1000G_superpops=TRUE,snplist_user=instruments$rsid)
head(snplist)
length(snplist)
```

Our SNP list now contains the rsids for: 1) the GWAS catalog top hits,
2) the MAF reference set and 3) the “exposure SNPs”. Next, we extract
the summary associations statistics for these SNPs from the glioma
outcome dataset. The example below has already been restricted to the
SNPs of interest to save space. The extract\_snps function only works on
MAC or linux operating systems.

``` r
File<-system.file("extdata", "glioma_test_dat.txt", package = "mrQC")
gli<-extract_snps(snplist=snplist,path_to_target_file=File,path_to_target_file_sep="\t")
dim(gli)
head(gli)
```

In the above example, the summary data for glioma was stored locally on
our machine in a tab separated text file. Alternatively, we could have
sourced the outcome summary data from the Open GWAS project
(<https://gwas.mrcieu.ac.uk/>). For example, to extract summary data for
thyroid cancer we could
run:

``` r
snplist<-make_snplist(efo = "thyroid carcinoma",trait="thyroid carcinoma",ref1000G_superpops=TRUE,snplist_user=instruments$rsid)
head(snplist)
length(snplist)
thy <- ieugwasr::associations(id="ieu-a-1082", variants=snplist,proxies=0)  
dim(thy)
thy
```

The make\_snplist() function returns a warning that no GWAS hits were
found when searching for thyroid carcinoma as the reported trait.
However the search was able to identify hits searching on the efo for
thyroid carcinoma.

Returning to the glioma example, having extracted the summary data for
the SNP rsids of interest, we now need to format the summary data. This
is to get the data into the expected format for the QC
functions.

``` r
Dat<-format_data(dat=gli,outcome="Glioma",population="European",pmid=22886559,study="GliomaScan",ncase="cases",ncontrol="controls",UKbiobank=FALSE,rsid="Locus",effect_allele="Allele1",other_allele="Allele2",or="OR",lci="OR_95._CI_l",uci="OR_95._CI_u",eaf="eaf.controls",p="p",efo="glioma")
head(Dat)
```

Now we are ready to perform some quality checks on the summary data

## Step 2. Check that the effect allele frequency column is correct

Next we create some plots to visualise potential problems with the
effect allele frequency column. We do this by comparing allele frequency
in the outcome glioma dataset to the 1000 genomes super populations.
Let’s restrict the comparison to the European super population, since
we know that the glioma dataset was derived from a European ancestry
population.

``` r
Plot1<-make_plot_maf(ref_1000G="EUR",target_dat=Dat)
Plot1
```

SNPs with a red colour are SNPs with incompatible minor allele
frequencies, i.e. the allele frequencies are above 0.5 in the target
dataset but less than 0.5 in the 1000 genomes dataset. In this example,
all SNPs are flagged as problematic and there is a strong inverse
correlation in minor allele frequency between the datasets. This
indicates that the reported effect allele frequency column, in the
outcome dataset, does not correspond to the reported effect allele. The
strong inverse correlation implies that the effect allele column
actually refers to the non-effect allele.

Next we construct a similar plot but this time comparing allele
frequency with all 1000 genomes super
populations.

``` r
Plot2<-make_plot_maf(ref_1000G=c("AFR","AMR","EAS","EUR","SAS","ALL"),target_dat=Dat)
Plot2
```

All SNPs across all super populations are flagged as problematic. This
illustrates that the function can identify problematic SNPs regardless
of the ancestry of the outcome dataset. We also see a strong inverse
correlation in the comparison with the European 1000 genomes super
population. This is not surprising, since we know that the the glioma
GWAS results were generated in a European ancestry population. This
illustrates that the strength of the correlation in MAF between the
datasets can also support inferences about the ancestral background of
the outcome dataset, although the function was not designed with this
objective in mind. A more efficient approach would be to select SNPs
with a much wider range of variation in minor allele frequency that
chosen here.

## Step 3. Check that the effect allele column is correct

We next compare effect alleles between the glioma outcome dataset and
the GWAS catalog, in order identify potential mis-specification of the
effect allele
column.

``` r
Plot3<-make_plot_gwas_catalog(dat=Dat,efo=unique(Dat$efo),trait="glioma")
Plot3
```

We see that there are three groups of SNPs: those flagged as showing no
effect size conflict; those with moderate effect size conflict; and
those with high effect size conflict. The distinction between moderate
and high effect size conflict is arbitrary but is specified to make
allowance for chance deviations between datasets. If the effect sizes
are in opposite directions, the effect size conflict flag is set to
moderate or high. If the two-sided P value in both datasets is ≤0.0001
then the flag is upgraded to high. In addition, if the summary
association statistics in the target dataset and GWAS catalog are
derived from the same publication, and effect sizes are in opposite
directions, the conflict flag is set to high regardless of P value.

Overall, the plot shows that the majority of SNPs show effect size
conflicts. In other words, alleles associated with higher risk of glioma
in the GWAS catalog tend to be associated with lower risk of glioma in
the outcome dataset. This indicates that the reported effect allele
column corresponds to the non-effect allele. Taken together with the
previous allele frequency plots, showing an inverse correlation in MAF
between the outcome dataset and 1000 genomes European superpopulation,
this strongly indicates that the effect and non-effect allele columns
have been incorrectly specified. The reported non-effect allele column
is very likely the effect allele column.

We can also make a plot comparing effect allele frequency between the
outcome dataset and the GWAS catalog, which we show in the next example.

``` r

Plot4<-make_plot_gwas_catalog(dat=Dat,plot_type="plot_eaf",efo=unique(Dat$efo),trait=unique(Dat$outcome))
Plot4
```

We see an inverse correlation in effect allele frequency (EAF) between
the outcome glioma dataset and the GWAS catalog in European ancestry
studies. This is the expected pattern when the effect allele column has
been incorrectly specified. In general, the EAF conflict flag is set to
moderate or high if EAF is not consistent between the outcome dataset
and the GWAS catalog (e.g. is \<0.5 in the outcome dataset but is \>0.5
in the GWAS catalog or vice versa). For conflicting SNPs, the flag is
further upgraded to high if effect allele frequency is \>0.6 or \<0.4.
This makes allowance for chance deviations in allele frequency for SNPs
with minor allele frequency close to
0.5.

## Step 4. Check whether the reported effect size corresponds to log odds ratios

We next compare the predicted log odds ratio to the reported effect
size, in order to identify other potential errors or issues. This step
is applicable to summary data that has been derived from a case-control
genome-wide association study. Since the function to derive the
predicted log odds ratio can be a bit slow, we restrict the glioma
example to just the first 20 SNPs.

``` r
Dat1<-Dat[1:5,]
Pred<-predict_lnor_sh(dat=Dat1)
Plot5<-make_plot_predlnor(dat=Pred)
Plot5
```

The plot shows a strong positive correlation between the predicted log
odds ratios and the reported effect size, an intercept close to zero and
a slope that is \>0.8. When the predicted log odds ratios and reported
effect sizes are identical, the intercept should be 0 and the slope
should be 1.

We can also plot the relative bias, i.e. the percentage deviation of the
predicted log odds ratio from the reported effect size.

``` r
Plot6<-make_plot_predlnor(dat=Pred,bias=TRUE)
Plot6
```

Overall the relative bias seems small and mostly varies from -10.9% to
-13.5%. Since genetic effect sizes tend to be small, a relative bias of
10% is very small in absolute terms.

In the next example we show a dataset that returns a slope and intercept
very different from
expectation

``` r
snplist<-make_snplist(efo = "kidney cancer",snplist_user=instruments$rsid)
ukb <- data.frame(ieugwasr::associations(id="ukb-b-1316", variants=snplist,proxies=0),stringsAsFactors=FALSE)
Ukb<-format_data(dat=ukb,outcome="Kidney cancer",population="European",pmid="ukb-b-1316",ncase=1114,ncontrol=461896,study="UKB",UKbiobank=TRUE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=145,efo = "kidney cancer")
Pred<-predict_lnor_sh(dat=Ukb[1:5,])
Plot5_2<-make_plot_predlnor(dat=Pred)
Plot5_2
```

Although there is a strong positive correlation, the slope is 400, when
we expect a slope of 1. In fact, further investigation revealed that
Open GWAS dataset ukb-b-1316 was generated using a linear model. In
other words, the results were derived from a model where kidney cancer
case-control status (controls coded 1 and cases coded 2) was regressed
on SNP genotype (additively coded). The effect size from this model can
be interpreted as the absolute risk of kidney cancer per copy of the
effect allele. We can transform this into a log odds ratio scale using
the transform\_betas() function.

``` r
Ukb2<-transform_betas(dat=Ukb,effect="lnor",effect.se="se")
Pred2<-predict_lnor_sh(dat=Ukb2[1:5,])
Plot5_3<-make_plot_predlnor(dat=Pred2)
Plot5_3
```

After transforming the reported effect size to a log odds ratio scale,
we now see a slope close to 1 for the relationship with the predicted
log odds ratio. More generally, we suggest that any dataset with an
unsual intercept or with a slope very different from 1 (e.g. \<0.8 or
\>1.2) should be investigated by the user for potential problems. In the
glioma example, the slope and intercept look reasonably close to what
we’d expect. In the kidney cancer example, the slope was very
different from 1 because the reported effect sizes had not been
generated in a logistic regression model. Other factors that could cause
the slope to differ from 1 include: 1) the impact of covariate
adjustment in the original GWAS, 2) deviations from Hardy Weinberg
equilibrium, 3) mismatches between reported and actual allele frequency
for some or all SNPs or 4) mismatches between the reported and effective
SNP-level sample sizes for some or all
SNPs.

## Step 5. Check whether the reported P values correspond to the reported effect sizes

Let’s now return back to the glioma dataset. Let’s generate some ZZ
plots, in order to flag SNPs with P values that don’t coincide with
their reported effect sizes. The zz\_plot() function compares Zp scores
(inferred from the reported P values) to Zb scores (inferred from the
reported effect size and standard error).

``` r
Plot7<-zz_plot(dat=Dat)
Plot7
```

The correlation between the Zp and Zb scores is 1, indicating very
strong concordance between the reported P values and reported effect
sizes in the glioma dataset for the selected SNPs. In the next example,
we highlight a dataset where there is discordance between the reported P
values and effect sizes.

``` r
instruments<-ieugwasr::tophits(id="met-d-PUFA")
snplist<-make_snplist(efo = "lung carcinoma",snplist_user=instruments$rsid)
luc <- ieugwasr::associations(id="ieu-a-966", variants=snplist,proxies=0)  
Dat<-format_data(dat=data.frame(luc,stringsAsFactors=F),outcome="Lung cancer",population="European",pmid=24880342,ncase=11348,ncontrol=15861,study="ILCCO",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",se="se",eaf="eaf",p="p",efo = "lung carcinoma")
Plot7_2<-zz_plot(dat=Dat)
Plot7_2
```

The correlation between the Zp and Zb scores is less than 1, suggesting
discordance between the reported P values and reported effect sizes,
which is also clear from visual inspection of the plot. In this example,
we have included three SNP sets (the GWAS catalog hits, the MAF
reference set and the “exposure SNPs”. Let’s restrict the comparison to
only the “exposure SNPs” and the GWAS catalog top
hits.

``` r
snplist<-make_snplist(efo = "lung carcinoma",snplist_user=instruments$rsid,ref1000G_superpops=FALSE)
luc <- ieugwasr::associations(id="ieu-a-966",variants=snplist,proxies=0)  
Dat<-format_data(dat=data.frame(luc,stringsAsFactors=F),outcome="Lung cancer",population="European",pmid=24880342,ncase=11348,ncontrol=15861,study="ILCCO",rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",se="se",eaf="eaf",p="p",efo = "lung carcinoma")
Plot7_3<-zz_plot(dat=Dat)
Plot7_3
```

The correlation between the Zp and Zb scores is still less than 1,
suggesting some discordance between the reported P values and reported
effect sizes, which is also clear from visual inspection of the plot.
For example, 1 SNP has a Z score close to 15 when estimated from the
reported P value but a Z score close to 5 when estimated from the
reported effect size and standard error.

## Step 6. Combine all plots into a single report

Next we combine all the plots into a single report.

``` r
Plot_list2<-ls()[grep("Plot[0-9]",ls())] 
Plot_list<-lapply(1:length(Plot_list2),FUN=function(x) eval(parse(text=Plot_list2[x])))
combine_plots(Plot_list=Plot_list,out_file="~/qc_report.png")
```

![title of image](/man/figures/README-qc_report.png)
