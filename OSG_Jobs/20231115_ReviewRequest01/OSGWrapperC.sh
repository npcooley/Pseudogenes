#!/bin/bash

# hypothetically we're landing in the srv directory of the container,
# but the directory structure above that should still exist and be accessible?

pwd
ls -lh

# new ID for jobmap B
echo ${1}
# input file
echo ${2}
# rename for output file
echo ${3}

# unpack PGAP reference data
tar xzf /srv/input-2022-04-14.build6021.tgz

ls -lh

# rename assembly
mv "${2}" assembly.fna.gz

# initial fasta is gzipped, unzip it
gunzip -v assembly.fna.gz

# rename the submol file
mv "${6}" submol.yaml

# cwltool -h
cwltool --timestamps --disable-color --preserve-entire-environment --outdir /srv/output /pgap/pgap/pgap.cwl /srv/controller.yaml

# rm "${6}"
rm assembly.fna submol.yaml

ANNOT=$(find . -name annot.gff)

echo ${ANNOT}
echo ${7}

mv "${ANNOT}" "${7}"
gzip "${7}"

ls -lh
