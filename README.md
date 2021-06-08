# Solita TRE Pipelines

## Workflows

### SAIGE

**NOTE** `saige.wdl` is the WDL 1.0 workflow; `saige-d2.wdl` is a WDL
`draft-2`-compliant port, without subworkflows, for the TRE.

#### Inputs

```json
{
  "SAIGE.plinkPrefix":   "PLINK file prefix for creating the GRM",
  "SAIGE.phenoFile":     "Phenotype file",
  "SAIGE.autosomeBGENs": "File of bgen filenames, per autosomal chromosome",
  "SAIGE.sampleFile":    "File of IDs of samples in the dosage file"
}
```
