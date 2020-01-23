# Reference files for TranscriptClean

## Create human and mouse splice junction files based on GENCODE
```
./make_gencode_v29_splice_reference.sh
```

## Create a combined human + SIRV splice junction file
```
cat gencode_v29_SJs.tsv ../SIRV_ERCC/SIRV_SJs.tsv  > gencode_v29_SIRV_SJs.tsv
```
