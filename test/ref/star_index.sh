#!/bin/bash

genome="ref/chr22/chr22.fa"
annotation="ref/chr22/chr22.gtf"

STAR \
        --runThreadN 1 \
        --runMode genomeGenerate \
        --genomeDir test/ref/star/chr22 \
        --genomeFastaFiles $genome \
        --sjdbGTFfile $annotation \
	--genomeSAindexNbases 11
