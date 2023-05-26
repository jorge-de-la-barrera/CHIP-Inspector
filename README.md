# Overview
The accumulation of somatic mutations in hematopoietic stem cells can lead to a clonal expansion of these mutant cells, a phenomenon referred to as clonal hematopoiesis of indeterminate potential or CHIP (Steensma et al., 2015). Recent studies describe CHIP as a novel cardiovascular risk factor independent of traditional risk factors (Fuster and Walsh, 2018; Jaiswal et al., 2017; Páramo Fernández, 2018). Identifying individuals with CHIP can have significant implications for the management of cardiovascular diseases. To facilitate the study of CHIP's influence on these conditions, a bioinformatics tool called CHIP Inspector has been created.

# CHIP-Inspector

CHIP Inspector is a bioinformatics tool for identification of somatic mutations potentially causative of clonal hematopoiesis with high sensitivity from deep-sequencing NGS data. It can be used in cross-sectional and longitudinal studies to characterize clonal expansion and the influence of CHIP driver mutations on pathological processes. It orchestrates the entire bioinformatics workflow (coded as wdl scripts) , spanning from the processing of raw sequencing data to the generation of candidate variants ready to be manually curated. CHIP-Inspector uses the following third-party software:

## At the sample level:
- fgbio (https://fulcrumgenomics.github.io/fgbio/): Used for grouping reads by Unique Molecular Identifiers (UMIs), calling and filtering consensus reads, and clipping them if they overlap.
- Picard (https://broadinstitute.github.io/picard/: Performs tasks such as marking duplicates, collecting hybrid-selection (HS) and library quality metrics, and other sequence manipulations like converting SAM to FASTQ, sorting BAM files, and merging BAM files.
- Samtools(http://www.htslib.org/): Utilized for tasks such as retrieving barcodes and performing BAM flagstat (to get alignment statistics), converting BAM to CRAM (compressed BAM format), and generating indexes for BAM and CRAM files.
- BWA (https://bio-bwa.sourceforge.net/) : Used for mapping both raw reads and consensus reads to a reference genome.
- GATK (https://gatk.broadinstitute.org/): Employed for calling genomic variants using tools like Mutect2 (for somatic variant calling) and HaplotypeCaller (for germline variant calling).

## At the cohort level:
- Bcftools (https://samtools.github.io/bcftools) : Used for merging somatic calls from individual single-sample VCF (Variant Call Format) files into a joint VCF file representing the cohort. It can also normalize and index the resulting VCF file.
- VEP (https://www.ensembl.org/info/docs/tools/vep/index.html): Used for annotating variants by providing functional and genomic information about the variants based on various databases and prediction algorithms.
