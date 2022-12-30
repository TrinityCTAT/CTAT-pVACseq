#!/bin/bash

set -e

#VERSION=`cat VERSION.txt`

#docker build --no-cache -t mbrown/pvactools:devel .
# docker build -f Dockerfile -t mbrown/pvactools:devel .
docker build --build-arg CACHEBUST=$(date +%s) -f Dockerfile -t brownmp/pvactools:devel .
