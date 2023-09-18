#!/bin/bash

# hypothetically we're landing in the srv directory of the container,
# but the directory structure above that should still exist and be accessible?

pwd
ls -lh

echo ${1}
echo ${2}
echo ${3}
echo ${4}
echo ${5}
echo ${6}

# unpack PGAP reference data
tar xzf /srv/input-2022-04-14.build6021.tgz

ls -lh

gzip -dv "${5}"

# cwltool -h
cwltool --timestamps --disable-color --preserve-entire-environment --outdir /srv/output /pgap/pgap/pgap.cwl "/srv/${2}"

rm "${6}"

ANNOT=$(find . -name annot.gff)

echo ${ANNOT}
echo ${4}

mv "${ANNOT}" "${4}"

# mv /svr/output/annot.gff "${4}"

ls -lh
