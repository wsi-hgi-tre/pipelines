# Solita TRE Pipelines

## Workflows

### SAIGE

**NOTE** `saige.wdl` is the WDL 1.0 workflow, which is easier to
understand; `saige-d2.wdl` is a WDL `draft-2`-compliant port, without
subworkflows and hacks, for the TRE.

#### Containers

* `eu.gcr.io/fg-qmul-testing-master/saige:0.44.5`

#### Inputs

```json
{
  "SAIGE.phenotypes":    ["Array", "of", "binary", "and", "continuous", "traits"],

  "SAIGE.plinkPrefix":   "PLINK file prefix for creating the GRM",
  "SAIGE.phenoFile":     "Phenotype file",
  "SAIGE.covariants":    ["Array", "of", "covariant", "columns"],

  "SAIGE.SPATestFOFN":   "File of chromosome-SPA filename pairs, tab-delimited",
  "SAIGE.sampleFile":    "File of IDs of samples in the dosage file"
}
```

##### `phenotypes` Array

Each phenotype should be prefixed with `b_` or `c_` for binary or
continuous traits, respectively.

##### `SPATestFOFN`

The SPA Test FOFN should be a tab-delimited file with fields:

1. Chromosome
2. Path to 8-bit BGEN (for autosomal chromosomes) or VCF (for allosomal
   chromosomes) input file

**NOTE** On the Solita TRE, the paths must be fully qualified Google
Cloud Storage bucket addresses.

Example:

```
1	gs://my-bucket/path/to/my-chr1.bgen
X	gs://my-bucket/path/to/my-chrX.vcf
```
