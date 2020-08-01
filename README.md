# non-competitive Allele Specific Expression analysis (ncASE)

ncASE is a workflow initially developed by Maria Gutin, and further developed by Ana Pinharanda, Andrew Taverner, Cheyenne Payne, Molly Schumer, and Patrick Reilly.

The core of ncASE is `ncASE.R`, which performs the filtering of sites and summarization/normalization of counts from an ASE experiment. The central idea of the method is, rather than competitively mapping reads to one reference or the other (which may introduce reference bias in the process), instead you map to both parental references (hence "noncompetitive"), and filter sites that appear subject to reference bias.

## Dependencies
1. R
1. tidyverse (R package)
1. awk (GNU awk preferred, POSIX awk untested)

## Dependencies for input preparation
### Parallelization if not using a job engine:
1. GNU Parallel
### Transcriptome mapping (BWA-MEM) pipeline
1. BWA (can be obtained via `git clone https://github.com/lh3/bwa --recursive`)
1. Samtools (mainly used with versions 1.0+, unknown compatibility below that)
1. HTSlib (a dependency of Samtools, so usually the version should match that of Samtools)
### Genome+Annotation mapping (STAR TranscriptomeBAM) pipeline
1. STAR (can be obtained via `git clone https://github.com/alexdobin/STAR --recursive`)
### Post-mapping processing
1. Picard (can be obtained via `git clone https://github.com/broadinstitute/picard --recursive`)
### Variant calling (for ascertainment of allele depths)
1. BCFtools (for variant calling, tested with various versions 1.3+)
1. GATK (mainly tested with 3.4, and IR and HC tasks should work with newer)

## Algorithm details

Filters for sites:



Weighting of allele depths to get normalized counts:



## Input preparation recipes

One possible approach to generating non-competitive allele depths uses BWA-MEM to map to parental transcriptomes, and BCFtools for variant calling. This approach has been used for an ASE analysis of *D. yakuba* inversion arrangements.

An easy way to do this is using the [Pseudoreference Pipeline](https://github.com/YourePrettyGood/PseudoreferencePipeline):

```bash
#Prepare your parental transcriptomes in identical coordinate spaces


#Index the parental transcriptomes:
PIPELINEDIR="[path to Pseudoreference Pipeline installation directory]"
${PIPELINEDIR}/indexDictFai.sh [parent1 transcriptome] BWA
${PIPELINEDIR}/indexDictFai.sh [parent2 transcriptome] BWA

#Prepare a metadata file with information about your samples and
# the references to map them to.
#The format is tab-separated, with columns:
#1) Prefix for output files (may contain a relative path which must exist)
#2) Reference to map to
#3) Path to read 1 FASTQ (may be relative or absolute, can be gzipped)
#4) Path to read 2 FASTQ (optional, use if you have PE data)
#For here, we assume your metadata file is called ncASE_mapping_metdata.tsv
# and that you have 9 libraries to map to 2 transcriptomes, so 18 lines

#Now map all your samples to each transcriptome:
#We use GNU Parallel here, but you could use SBATCH arrays or whatnot
NUMTHREADS=1 #The number of threads to use per sample
parallel -j18 --eta "${PIPELINEDIR}/localArrayCall_v2.sh {1} MAP ncASE_mapping_metadata.tsv ${NUMTHREADS} only_bwa" ::: {1..18}

#Now call variants with BCFtools:
parallel -j18 --eta "${PIPELINEDIR}/localArrayCall_v2.sh {1} MPILEUP ncASE_mapping_metadata.tsv ${NUMTHREADS} no_markdup,no_IR" ::: {1..18}

#Finally, generate the .snps files from the bgzipped VCFs:
#We do something slightly special here with GNU parallel, which is
# colsep with ::::. Four colons means to read lines from the file,
# and colsep splits each line by its value, meaning that {1} will
# get substituted with the value in the first column of the file.
#You can take a different approach if you want to for cluster job
# engines.
parallel -j18 --eta --colsep="\t" "gzip -dc {1}_nomarkdup_mpileupcall.vcf.gz | [ncASE directory]/vcf_to_ncASE.awk > {1}.snps" :::: ncASE_mapping_metadata.tsv
```
`vcf_to_ncASE.awk` is currently written for BCFtools VCFs, but can be modified to use the `FORMAT/AD` tag of GATK VCFs produced by GenotypeGVCFs. Adding this flexibility is on the to-do list.

## ncASE usage

Once the .snps files are generated, provide them to `ncASE.R` in a particular order: sample1_par1 sample1_par2 sample2_par1 sample2_par2 ... threshold readlength outputfilename allow0DP

Yes, the command will look rather ugly if your .snps filenames are long. We haven't recoded it to use a file of filenames instead, yet. That's on the to-do list.

## Rescuing homozygous ref sites and including extreme SNPs

In the original ncASE code (and the awk code for preparation of the .snps files), sites called as homozygous ref were omitted from the SNPs files. These sites do have the potential to be the strongest indicators of ASE, but are also the most likely cases for reference bias. If you do wish to include these sites in the calculations of per-transcript per-parent read count estimates output by `ncASE.R`, you will need to perform some extra steps before running `ncASE.R`:

1. `union_AIMs.awk` must be used with all .snps files as input. This will generate a TSV of all sites with non-homozygous-ref genotypes in at least one sample. We use this as a pre-filter during the rescuing process, rather than generating .snps files including all 0/0 sites, and then needing to filter out sites where all samples are 0/0, as these sites should be quite abundant. Hence, we save memory, space, and time by taking this union.

2. `rescue_hom_ref.awk` takes the output of `union_AIMs.awk` as its first argument, and the uncompressed BCFtools VCF as its second argument. The output is a .snps file that includes sites with all the biallelic genotypes (but only those sites in the union), and with the alt allele inferred from the union (since 0/0 calls have . in the ALT column in a standard VCF).

You can then use these "rescued" .snps files as input to `ncASE.R`, with the last argument being set to 1 if you want to allow sites with 0 allele depth.

If you would like to reproduce results of the original code, skip the rescue process and set `allow0DP` to 0. If you would like to use such extreme sites, go through the rescue process and set `allow0DP` to 1.

## Example rescue usage on Dyak inversion ASE data:

```bash
#We continue the process after running vcf_to_ncASE.awk on each VCF:
#Take the union of segregating sites:
/usr/bin/time -v [path to ncASE]/union_AIMs.awk Dyak_NY73PB_backbone/Dyak_*/*.snps 2> union_AIMs_Dyak_NY73PB_backbone.stderr > Dyak_NY73PB_backbone_union.aims
#Rescue segregating sites with homref genotypes for each sample:
parallel -j18 --eta '/usr/bin/time -v [path to ncASE]/rescue_hom_ref.awk Dyak_NY73PB_backbone_union.aims <(gzip -dc Dyak_NY73PB_backbone/Dyak_cross{3}_rep{2}_{1}/Dyak_{4}_cross{3}_rep{2}_{1}_nomarkdup_realigned_mpileupcall.vcf.gz) 2> logs/rescue_Dyak_{4}_cross{3}_rep{2}_{1}.stderr > Dyak_{4}_cross{3}_rep{2}_{1}_rescued.snps' ::: carcass head ::: {1..3} ::: {1..3} ::: NY73PB Tai18E2
#Now run ncASE.R while allowing sites with 0 allele depth:
#Command is quite long
/usr/bin/time -v [path to ncASE]/ncASE.R Dyak_NY73PB_cross1_rep1_carcass_rescued.snps Dyak_Tai18E2_cross1_rep1_carcass_rescued.snps Dyak_NY73PB_cross1_rep2_carcass_rescued.snps Dyak_Tai18E2_cross1_rep2_carcass_rescued.snps Dyak_NY73PB_cross1_rep3_carcass_rescued.snps Dyak_Tai18E2_cross1_rep3_carcass_rescued.snps Dyak_NY73PB_cross2_rep1_carcass_rescued.snps Dyak_Tai18E2_cross2_rep1_carcass_rescued.snps Dyak_NY73PB_cross2_rep2_carcass_rescued.snps Dyak_Tai18E2_cross2_rep2_carcass_rescued.snps Dyak_NY73PB_cross2_rep3_carcass_rescued.snps Dyak_Tai18E2_cross2_rep3_carcass_rescued.snps Dyak_NY73PB_cross3_rep1_carcass_rescued.snps Dyak_Tai18E2_cross3_rep1_carcass_rescued.snps Dyak_NY73PB_cross3_rep2_carcass_rescued.snps Dyak_Tai18E2_cross3_rep2_carcass_rescued.snps Dyak_NY73PB_cross3_rep3_carcass_rescued.snps Dyak_Tai18E2_cross3_rep3_carcass_rescued.snps 0.10 150 Dyak_allCrosses_carcass_on_Dyak_NY73PB_f0.10_L150_allow0DP_ncASE.tsv 1 2> logs/ncASE_allCrosses_Dyak_NY73PB_backbone_carcass_f0.10_L150_allow0DP.stderr > logs/ncASE_allCrosses_Dyak_NY73PB_backbone_carcass_f0.10_L150_allow0DP.stdout
```

### Notes on rescue:

The rescue process has only been minimally tested. Due to the distinction between the genotype filter in `vcf_to_ncASE.awk` and the allele depth filter in `ncASE.R`, different combinations of rescuing and setting `allow0DP` will give different results.

Rescue recovers not only sites with 0 depth for one of the two alleles, but also other extreme imbalances (e.g. ~100 for ref, ~1 for alt). Thus the result of skipping rescue is likely *not* the same as the result of performing rescue but setting `allow0DP` to 0.

The `allow0DP` toggle switches the within-parent allele depth filter from a logical AND (implemented by multiplying the ref and alt depths together, testing for > 0) to a logical OR (implemented as a vectorized OR in R, testing for > 0 for the two alleles). Thus under `allow0DP=1`, the only way this filter eliminates a site is if the allele depths are both 0. 
