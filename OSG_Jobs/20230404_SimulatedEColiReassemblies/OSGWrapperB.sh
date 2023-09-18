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
mv "${2}" Assembly.fna.gz

# initial fasta is gzipped, unzip it
gzip -dv Assembly.fna.gz

# cwltool -h
cwltool --timestamps --disable-color --preserve-entire-environment --outdir /srv/output /pgap/pgap/pgap.cwl /srv/Controller.yaml

# rm "${6}"
rm Assembly.fna

ANNOT=$(find . -name annot.gff)

echo ${ANNOT}
echo ${3}

mv "${ANNOT}" "${3}"

# mv /svr/output/annot.gff "${4}"

ls -lh
