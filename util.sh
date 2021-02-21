#!/bin/bash

###############################################################################
# salmon index transcriptome
salmon_index_dir_name=salmon_index_mm10
path_to_salmon_index=~/Downloads/${salmon_index_dir_name}
mkdir -p $path_to_salmon_index
echo "Downloading salmon transcriptome index file"
wget -O ~/Downloads/${salmon_index_dir_name}.tgz http://refgenomes.databio.org/v2/asset/mm10_cdna/salmon_index/archive?tag=default
tar -xvzf ~/Downloads/${salmon_index_dir_name}.tgz -C $path_to_salmon_index
salmon_index=${path_to_salmon_index}/salmon_index

# salmon_index=/home/prime/Downloads/salmon_index_mouse/salmon_index
# http://refgenomes.databio.org/v2/asset/mm10_cdna/salmon_index/archive?tag=default
# download other indexes here http://refgenomes.databio.org/
###############################################################################
