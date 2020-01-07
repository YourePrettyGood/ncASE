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

Once the .snps files are generated, provide them to `ncASE.R` in a particular order: sample1_par1 sample1_par2 sample2_par1 sample2_par2 ... threshold readlength

Yes, the command will look rather ugly if your .snps filenames are long. We haven't recoded it to use a file of filenames instead, yet. That's on the to-do list.
