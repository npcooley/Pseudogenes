#!/bin/bash

# running inside a singularity container, ENV commands in the dockefile aren't run, so need to export the PATH
# This path needs to match where the executable was installed in the dockerfile
# export TMPDIR=$(pwd)
# export NCBI_SETTINGS="$(pwd)/user-settings.mkfg"
# vdb-config -o n NCBI_SETTINGS
# echo "Aexyos" | vdb-config -i
# export PATH=$PATH:/ANIcalculator_v1

# ani calc is hypothetically correctly setup in the new docker container
# tar -zxvf ANIcalculator_v1.tgz
# export PATH=/srv/ANIcalculator_v1:$PATH
# chmod +x ANIcalculator_v1/*

Rscript JobScript.R ${1} ${2} ${3} ${4}

if [ -e Result*.RData ]
then
  exit 0
else
  exit 1
fi
