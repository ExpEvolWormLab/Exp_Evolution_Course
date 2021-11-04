

## Downloaded the main genomic data file from:
#https://datadryad.org/stash/downloads/file_stream/14282

#Awk command to filter on CHR value / kept 2L only
#awk '($1=="2L")'  Dsim_F0-F60_Q20_polymorphic_CMH_FET_blockID.sync > Dsim_F0-F60_Q20_polymorphic_CMH_FET_blockID_CHR_2L.sync

# Load the file in R

Dsim_F0-F60_Q20_polymorphic_CMH_FET_blockID_CHR_2L.sync

load()
