#!/usr/bin/env bash

# Create conda environment
conda env create -f bin/netResponseEnv.yml

# Download and unzip large files
cd ./database
curl -O http://iid.ophid.utoronto.ca/static/download/mouse_annotated_PPIs.txt.gz ; gunzip mouse_annotated_PPIs.txt.gz
curl -O http://iid.ophid.utoronto.ca/static/download/human_annotated_PPIs.txt.gz ; gunzip human_annotated_PPIs.txt.gz
curl -O https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv
curl -O https://www.grnpedia.org/trrust/data/trrust_rawdata.mouse.tsv
curl -O https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt
curl -O https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/MOUSE_10090_idmapping.dat.gz ; gunzip MOUSE_10090_idmapping.dat.gz
curl -O https://www.informatics.jax.org/downloads/reports/MRK_List2.rpt
curl -O https://s3.amazonaws.com/data.clue.io/repurposing/downloads/repurposing_drugs_20200324.txt
curl -OL https://github.com/EwaldLab/2019_Prkd1/raw/master/Georgess_Prkd1_2019_Data_Final.xlsx
curl -O https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3941052/bin/supp_jcb.201306088_JCB_201306088_TableS1.xlsx
cd ..

#activate conda environment
source activate netResponseEnv

#create Twist1 biological knowledge
python ./bin/createTwist1BiologicalKnowledge.py
