# Neo Epitope Prediction
The provided code is inteded to be used as a part of the [Neo Epitope Prediction Workflow](https://github.com/csam5596/EpitopePrediction)
but can also be used as a standalone analytical tool to derive neo epitopes from mutational data.

Neo Epitopes are generated from somatic Single Nucleotide Variations, Indels, and Frameshift Mutations
while also considering variant phasing using [Phaser](https://github.com/secastel/phaser) as well as taking germline context into account.
All epitopes are evaluated in terms of their immunogenicity using a Random Forest model based on BLOMAP encoding.

For more information see: https://github.com/csam5596/MasterThesis


## Requirements
As a part of the analysis the following tools are required to be preinstalled:
  - grep
  - xargs
  - Python (required to run [Phaser](https://github.com/secastel/phaser))
  - samtools
  - bgzip
  - tabix
  - [vep](https://www.ensembl.org/info/docs/tools/vep/index.html)
  - [NetMHCpan](http://www.cbs.dtu.dk/services/NetMHCpan/)

## Arguments

### Required

| Argument     | Type   | Description |
|--------------|--------|-------------|
|--dna         | String | Tumor DNA - aligned, sorted and indexed (**.bam**)  |
|--rna         | String | Tumor RNA - aligned, sorted and indexed (**.bam**) |
|--somatic     | List   | Space separated list of variant caller-files containing somatic variants (**.vcf** \| **.vcf.gz**) |
|--vc          | List   | Space separated list of names of variant callers used for each vcf file given with --somatic. Each can be any of **Mutect2**, **Strelka**, **Varscan** |
|--germline    | String | File containing germline variants called by GATK HaplotypeCaller (**.vcf** \| **.vcf.gz**) |
|--proteins    | String | Ensembl reference protein file (**.fa**) </br>(see https://www.ensembl.org/info/data/ftp/index.html)|
|--gtf         | String | File containing genomic positions (**.gtf**) </br>(see https://www.ensembl.org/info/data/ftp/index.html)|
|--ref         | String | Reference Genome (**.fa**) </br>(see https://www.ensembl.org/info/data/ftp/index.html)|
|--hla         | String | Output file derived from Optitype |
|--tpm         | String | Output file derived from Kallisto (**.tsv**) |
|--phaser      | String | Path to **phaser.py** (see [Phaser](https://github.com/secastel/phaser))|
|--vepCacheDir | String | Path to cache for Variant Effect Predictor |

For both dna and rna sequence files an index file is required with the same name followed by the **.bai** extension.

### Optional Filter Parameters

| Argument      | Type    | Default | Description |
|---------------|---------|--------:|-------------|
|--minVcCount   | Integer | 2       | Minimum number of variant callers required to confirm a mutation |
|--maxNormalVaf | Float   | 0.05    | Maximum VAF for variant in normal sample |
|--minTumorVaf  | Float   | 0.05    | Minimum VAF for variant in tumor sample |
|--minRcAbs     | Integer | 1       | Minimum absolute number of reads required to confirm a variant |
|--minRcRel     | Float   | 0.2     | Minimum fraction of reads spanning a mutated position required to confirm a variant |
|--minTpmCount  | Float   | 1.0     | Minimum TPM for a gene to be considered expressed |
|--dnaMapQ      | Integer | 1       | Minimum DNA mapping quality for phaser (Depends on aligner used) |
|--rnaMapQ      | Integer | 255     | Minimum RNA mapping quality for phaser (Depends on aligner used) |



### Flags

| Argument   | Description |
|------------|-------------|
| --help, -h | Produce help message |
| --dirty    | Keep temporary files |

## Output

As output a single file: **epitopes.tsv** is generated containing the following columns:

| Parameter             | Type    | Description |
|-----------------------|---------|-------------|
| Peptide               | String  | Sequence of identified neo Epitope |
| Variants              | String  | Comma separated list of mutations that are contained in the sequence. Each mutation is given in the form: </br> **chr\<Id\>:\<Position\>:\<Reference\>:\<Mutated\>** |
| VariantTypes          | String  | Comma separated list of variant types. Each can be one of: </br> **SNV**, **Insertion**, **Deletion**, **Frameshift** |
| VariantCallers        | String  | Comma separated list of number of variant callers that confirmed the respective mutation |
| DNA VAF Normal        | String  | Comma separated list of VAFs in normal DNA sample (Each averaged over the reporting variant callers) |
| DNA VAF Tumor         | String  | Comma separated list of variant VAFs in tumor DNA sample (Each averaged over the reporting variant callers) |
| RNA VAF               | String  | Comma separated list of variant VAFs in tumor RNA |
| Read Count DNA Normal | String  | Comma separated list of variant read counts in normal DNA sample (Each averaged over the reporting variant callers) |
| Read Count DNA Tumor  | String  | Comma separated list of variant read counts in tumor DNA sample (Each averaged over the reporting variant callers) |
| Read Count RNA Tumor  | String  | Comma separated list of variant read counts in tumor RNA sample |
| gene                  | String  | Name of the gene the variants occurred in |
| Protein               | String  | Ensembl-ID of the mutated protein |
| HLA                   | String  | HLA type for which binding affinity is reported |
| BindingCore           | String  | 9-mer directly bound to MHC complex |
| ICore                 | String  | Interaction core |
| RawPredictionScore    | Float   | Raw binding affinity prediction score |
| Affinity (nM)         | Float   | Binding affinity reported in nano molar |
| %Rank                 | Float   | Rank of binding affinity |
| Exp                   | Integer | NetMHCpan associated value |
| Strength              | String  | **WB** if %Rank \< 2, **SB** if %Rank < 0.5, empty otherwise
| TPM                   | Float   | Transcripts per million of the respective gene |
| expressed             | String  | **yes** if TPM \> minTpmCount - **no** otherwise |
| ImmunogenicityScore   | Float   | Raw predicted immunogenicity score ranging from **0** to **1**
| Immunogenic           | String  | **yes** if ImmunogenicityScore \> 0.5 - **no** otherwise 
