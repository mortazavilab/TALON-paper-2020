# CAGE
To externally validate our long read transcript starts, we compare them to 5' end cap sites derived from the CAGE assay.

## FANTOM5
For FANTOM, we have one file of robust peaks measured across a variety of human samples rather than individual files per cell line. The peaks are mapped to hg19, so we need to use the UCSC genome browser liftover tool to convert them to hg38. I already downloaded this tool and the necessary chain file when working on the RNA-PET data.
### Download the CAGE data
```
./download_FANTOM5_CAGE.sh
```
### LiftOver to hg38
#### Download Liftover
```
# Program:
Linux: wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
Mac OSX: wget http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/liftOver

# Chain file:
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
gunzip hg19ToHg38.over.chain.gz
```
Run Liftover on the CAGE peaks
```
./liftover_FANTOM.sh
```
