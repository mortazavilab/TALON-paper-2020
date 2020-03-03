# GENCODE_v29
Human GENCODE annotation

## Combine with SIRVs and ERCCs
```
cat gencode.v29.annotation.gtf ../SIRV_ERCC/ERCC-spikein_reformatted.gtf \
    ../SIRV_ERCC/SIRV-spikein_gene-sep_reformatted.gtf > gencode.v29.SIRV.ERCC.annotation.gtf  
```
