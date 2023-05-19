#!/bin/bash
GTF="ref/chr22/chr22.gtf"
ref="ref/chr22/chr22.fa"

mkdir -p test/ref/rsem

rsem-prepare-reference \
        --gtf $GTF \
        -p 1 \
        $ref \
        test/ref/rsem/chr22
