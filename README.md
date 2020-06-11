# feox_timeseries

## Files for making the phyloseq object

Rscript: FeOx_phyloseq.R

OTU table: seqtab.rds

Taxa table: taxa.rds

Metadata: FeOx_demux_metadata.csv

## Files for BEEM-static analysis (network graphs)

https://github.com/lch14forever/BEEM-static

### Shows how to go from a phyloseq object to the necessary input for BEEM-static package

Rscript: FeOx_beemstatic.R

Phyloseq object: feox-phyloseq.rds

## Converting phyloseq object into format to use with Web-gLV

Tool: https://web.rniapps.net/webglv/

Rscript: phyloseq_to_web_gLV.R

Phyloseq object: feox-phyloseq.rds

Output TSV table: gLV_otu_table.tsv

## Make network out of phyloseq object using tidyverse and ggplot2

Rscript: FeOx_networks.R
Rscript for faceting by time: FeOx_TimeNetworks.R
