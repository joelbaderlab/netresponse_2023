#!/usr/bin/env bash

#activate conda environment
source activate netResponseEnv

#network directory
network_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

#run NetResponse analysis
python ./bin/netResponse.py analysis -s mouse $network_dir Twist1

#compare NetResponse to biological knowledge using dissemination assay results from Georgess et al. (2020)
python ./bin/netResponse.py comparison -s mouse $network_dir Twist1
