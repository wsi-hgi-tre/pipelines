# Solita TRE Pipelines

## Workflows

### SAIGE

**Note** `saige-d2.wdl` is the WDL `draft-2`-compliant workflow, without
any subworkflows and including the hacks needed to work around the TRE's
assumptions. Additional commentary has been included to help clarify
things.

#### Containers

* `eu.gcr.io/fg-qmul-containers/saige:0.44.5`

#### Inputs

```json
{
  "SAIGE.phenotypeFile": "File of traits (one per line)",

  "SAIGE.plinkPrefix":   "PLINK file prefix for creating the GRM",
  "SAIGE.phenoFile":     "Phenotype file",
  "SAIGE.covarFile":     "File of covariants (one per line)",

  "SAIGE.SPATestFOFN":   "File of chromosome-SPA filename pairs, tab-delimited",
  "SAIGE.sampleFile":    "File of IDs of samples in the dosage file"
}
```

##### `phenotypeFile` File

Each phenotype should be prefixed with `b_` or `c_` for binary or
continuous traits, respectively.

##### `SPATestFOFN` File

The SPA Test FOFN should be a tab-delimited file with fields:

1. Chromosome;
2. Path to 8-bit BGEN (for autosomal chromosomes) or VCF (for allosomal
   chromosomes) input file.

**Note** On the Solita TRE, the paths must be fully qualified Google
Cloud Storage bucket addresses.

Example:

```
1	gs://my-bucket/path/to/my-chr1.bgen
X	gs://my-bucket/path/to/my-chrX.vcf
```

**Note** For allosomal chromosomes, a sibling TABIX index, suffixed with
`.tbi`, is expected alongside the VCF file. These are not explicitly
enumerated as part of the input.

#### Outputs

For each phenotype (i.e., the contents of the `phenotypeFile` file), a
gzip'd TSV file is produced, named `step3-PHENOTYPE.gwas.gz`.

**Note** The allosomal GWAS output excludes the rsID field, so a dummy
entry (copied from the previous field) is added in its place.
