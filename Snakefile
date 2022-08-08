configfile: "config/config.yaml"
configfile: "config/params.yaml"

rule all:
	input:
		expand("output/{name}/final/{file}", 
			name = ["sim6A", "sim6B", "sim6C"],
			file = ["pseudoalignments.bam", "phased.vcf", "expressionMatrix.tsv"])


# ===========================================
# Generate genomes of both alleles
# ===========================================

rule generatingAltGenome:
	input:
		config["ref_genome"]
	output: 
		genome = "output/{prefix}/genome/allele{num}.fa", 
		vcf = "output/{prefix}/genome/allele{num}.refseq2simseq.SNP.vcf"
	params: 
		snp_count=45000
	log: "output/{prefix}/log/generatingAltGenome_{num}.log"
	shell:
		"""
		mkdir -p output/{wildcards.prefix}/genome
		
		perl scripts/simuG/simuG.pl \
			-refseq {input} \
			-snp_count {params.snp_count} \
			-prefix output/{wildcards.prefix}/genome/allele{wildcards.num} &> {log}
		#-gene_gff /projects/koneill_scratch/pog_tmp/ref/Homo_sapiens.GRCh38.100.chr.gff3.gz 
		#-coding_partition_for_snp_simulation coding &> {log} 
		sed 's/^>/>chr/' output/{wildcards.prefix}/genome/allele{wildcards.num}.simseq.genome.fa > {output.genome}
		awk -v OFS="\t" '/^#/ {{print $0; next}} {{$1 = "chr"$1 ; print $0}}' {output.vcf} > {output.vcf}.tmp && mv {output.vcf}.tmp {output.vcf}
		"""	


# ============================================
# Generating allelic RNA read
# ============================================
#rule sortBedFile:
#	input: config["annotation_bedfile"]
#	output: "output/{prefix}/simRNA/annotation.sorted.bed"
#	log: "output/{prefix}/log/sortBedFile.log"
#	shell:
#		"""
#		mkdir -p output/{wildcards.prefix}/simRNA
#		sort -k 1,1 -k 2,2n {input} > {output} 2> {log}
#		"""
#
#rule generateExpressionProfile:
#	input: "output/{prefix}/simRNA/annotation.sorted.bed"
#	output: "output/{prefix}/simRNA/explvprofile.txt"
#	log: "output/{prefix}/log/generateExpressionProfile.log"
#	shell:
#		"""
#		scripts/RNASeqReadSimulator/genexplvprofile.py -e 4,4 {input} > {output} 2> {log}
#		"""
#
#rule allelicExpressionProfile: 
#	input: "output/{prefix}/simRNA/explvprofile.txt"
#	output: 
#		"output/{prefix}/simRNA/allele1.explvprofile.txt",
#		"output/{prefix}/simRNA/allele2.explvprofile.txt",
#		"output/{prefix}/final/ASEGene.txt"
#	params:
#		ASEgene = 1000
#	log: "output/{prefix}/log/allelicExpressionProfile.log"
#	shell:
#		"""
#		mkdir -p output/{wildcards.prefix}/final
#		scripts/generateAllelicExpressionProfile.R -e {input} -n {params.ASEgene} -o output/{wildcards.prefix}/simRNA &> {log}
#		mv output/{wildcards.prefix}/simRNA/knownASE.txt output/{wildcards.prefix}/final/ASEGene.txt
#		"""
#
#rule generateAlleleBed: 
#	input: 
#		bedfile="output/{prefix}/simRNA/annotation.sorted.bed",
#		expressionProfile="output/{prefix}/simRNA/allele{num}.explvprofile.txt"
#	output: "output/{prefix}/simRNA/allele{num}.bed"
#	params:
#		numRead=50000000//2,
#		readLength=200
#	log: "output/{prefix}/log/generateAlleleBed_{num}.log"
#	shell:
#		"""
#		mkdir -p output/{wildcards.prefix}/simRNA/reads/
#		scripts/RNASeqReadSimulator/gensimreads.py \
#			-e {input.expressionProfile} \
#			-n {params.numRead} \
#			-l {params.readLength} \
#			-p 200,20 \
#			{input.bedfile} > {output} 2> {log}
#		"""
#
#rule generateAlleleReadFasta:  
#	input: 
#		genome = "output/{prefix}/genome/allele{num}.fa",
#		bed = "output/{prefix}/simRNA/allele{num}.bed"
#	output: "output/{prefix}/simRNA/reads/allele{num}.fa"
#	params:
#		readLength=100
#	log: "output/{prefix}/log/generateAlleleReadFasta_{num}.log"
#	shell:
#		"""
#		scripts/RNASeqReadSimulator/getseqfrombed.py \
#			-f A -r 0.1 \
#			-l {params.readLength} \
#			{input.bed} \
#			{input.genome} > {output} 2> {log}
#		"""
#
#rule splitFasta:
#        input: "output/{prefix}/simRNA/reads/allele{num}.fa"
#        output:
#                "output/{prefix}/simRNA/reads/allele{num}_1.fa",
#                "output/{prefix}/simRNA/reads/allele{num}_2.fa"
#        log: "output/{prefix}/log/splitFasta_{num}.log"
#        shell:
#                """
#                scripts/RNASeqReadSimulator/splitfasta.py -o output/{wildcards.prefix}/simRNA/reads/allele{wildcards.num} < {input}
#                """
#
#rule catReads:
#        input:
#                R1 = expand("output/{{prefix}}/simRNA/reads/allele{num}_1.fa", num=[1,2]),
#                R2 = expand("output/{{prefix}}/simRNA/reads/allele{num}_2.fa", num=[1,2])
#        output:
#                R1 = "output/{prefix}/final/rnaReads_R1.fa",
#                R2 = "output/{prefix}/final/rnaReads_R2.fa"
#        log: "output/{prefix}/log/catReads.log"
#        shell:
#                """
#                cat {input.R1} > {output.R1} 2> {log}
#                cat {input.R2} > {output.R2} 2> {log}
#                """

#########################
# Generate allelic read 2
#########################

#rule getTranscriptome:
#	input: 
#		genome = "output/{prefix}/genome/allele{num}.fa",
#		transcriptomeBed = config["annotation_bedfile"]
#	output: "output/{prefix}/transcriptome/allele{num}.cDNA.fa"
#	singularity: "docker://quay.io/biocontainers/bedtools:2.23.0--h5b5514e_6"
#	log: "output/{prefix}/log/getTranscriptome_{num}.log"
#	shell:
#		"""
#		mkdir -p output/{wildcards.prefix}/transcriptome
#		bedtools getfasta -fi {input.genome} -bed {input.transcriptomeBed} -name -fo {output} &> {log}
#		"""	

rule getTranscriptome:
	input:
		genome = "output/{prefix}/genome/allele{num}.fa",
		gtf = config["annotation_gtf"]
	output: "output/{prefix}/transcriptome/allele{num}.cDNA.fa"
	singularity: "docker://quay.io/biocontainers/gffread:0.12.7--hd03093a_1"
	log: "output/{prefix}/log/getTranscriptome_{num}.log"
	shell:
		"""
		mkdir -p output/{wildcards.prefix}/transcriptome
		gffread \
			-w {output} \
			-g {input.genome} {input.gtf}
		"""

rule cutNTranscriptome:
	input: "output/{prefix}/transcriptome/allele{num}.cDNA.fa"
	output: "output/{prefix}/transcriptome/allele{num}.cDNA.cutN.fa"
	singularity: "docker://quay.io/biocontainers/seqtk:1.3--h7132678_4"
	log: "output/{prefix}/log/cutNTranscriptome_{num}.log"
	shell:
		"seqtk cutN {input} > {output} 2> {log}"

rule generateExpressionProfile:
	input: "output/{prefix}/transcriptome/allele1.cDNA.cutN.fa"
	output: 
		"output/{prefix}/simRNA/allele1/expressionProfile.rds",
		"output/{prefix}/simRNA/allele2/expressionProfile.rds",
		"output/{prefix}/simRNA/ASEGene.txt"
	params:
		numASE = 100,
		foldChange = 6,
		percentExpressed = 0.8,
		depth = 10
	log: "output/{prefix}/log/generateExpressionProfile.log"
	shell:
		"""
		mkdir -p output/{wildcards.prefix}/simRNA/allele1
		mkdir -p output/{wildcards.prefix}/simRNA/allele2

		scripts/polyesterGenerateProfile.R \
			-t {input} \
			-n {params.numASE} \
			-f {params.foldChange} \
			-e {params.percentExpressed} \
			-d {params.depth} \
			-o output/{wildcards.prefix}/simRNA &> {log}
		"""

rule simulateReadPolyester: 
	input: 
		expressionProfile = "output/{prefix}/simRNA/allele{num}/expressionProfile.rds",
		transcriptome = "output/{prefix}/transcriptome/allele{num}.cDNA.cutN.fa"
	output: 
		temp("output/{prefix}/simRNA/allele{num}/sample_01_1.fasta"),
		temp("output/{prefix}/simRNA/allele{num}/sample_01_2.fasta")
	log: "output/{prefix}/log/simulateReadPolyester_{num}.log"
	shell:
		"""
		scripts/polyester.R \
			-t {input.transcriptome} \
			-e {input.expressionProfile} \
			-s $(($RANDOM * {wildcards.num})) \
			-o output/{wildcards.prefix}/simRNA/allele{wildcards.num} &> {log}
		"""

rule catReads:
	input: 
		R1 = expand("output/{{prefix}}/simRNA/allele{num}/sample_01_1.fasta", num=[1,2]),
		R2 = expand("output/{{prefix}}/simRNA/allele{num}/sample_01_2.fasta", num=[1,2])
	output: 
		R1 = "output/{prefix}/final/rnaReads_R1.fa",
		R2 = "output/{prefix}/final/rnaReads_R2.fa"
	log: "output/{prefix}/log/catReads.log"
	shell:
		"""
		cat {input.R1} > {output.R1} 2> {log}
		cat {input.R2} > {output.R2} 2> {log}
		"""


# ===========================================
# RNA-seq analysis 
# ===========================================
rule fasta2fastq:
	input: 
		R1 = "output/{prefix}/final/rnaReads_R1.fa",
		R2 = "output/{prefix}/final/rnaReads_R2.fa"
	output: 
		R1 = temp("output/{prefix}/rnaSeq/rnaReads_R1.fq"),
		R2 = temp("output/{prefix}/rnaSeq/rnaReads_R2.fq")
	singularity: "docker://quay.io/biocontainers/seqtk:1.3--h7132678_4"
	log: "output/{prefix}/log/fasta2fastq.log"
	shell:
		"""
		mkdir -p output/{wildcards.prefix}/rnaSeq	
		seqtk seq -F 'J' {input.R1} > {output.R1} 2> {log}
		seqtk seq -F 'J' {input.R2} > {output.R2} 2> {log}
		"""
	
rule kalliso:
	input:
		index = config["kalliso_index"],
		gtf_annotation = config["gtf_annotation"],
		chrom_length = config["chrom_length"],
		reads = expand("output/{{prefix}}/rnaSeq/rnaReads_R{num}.fq", num=[1,2])
	output: 
		"output/{prefix}/final/pseudoalignments.bam",
		"output/{prefix}/rnaSeq/abundance.h5"
	singularity: "docker://quay.io/biocontainers/kallisto:0.48.0--h0d531b0_1"
	log: "output/{prefix}/log/kalliso.log"
	threads: 72
	shell:
		"""
		kallisto quant \
			-i {input.index} \
			-g {input.gtf_annotation} \
			-c {input.chrom_length} \
        		--genomebam \
        		-o output/{wildcards.prefix}/rnaSeq \
			-t {threads} \
        		{input.reads} &> {log}
		mv output/{wildcards.prefix}/rnaSeq/pseudoalignments.bam output/{wildcards.prefix}/final/pseudoalignments.bam 
		mv output/{wildcards.prefix}/rnaSeq/pseudoalignments.bam.bai output/{wildcards.prefix}/final/pseudoalignments.bam.bai
		"""

rule expressionMatrix:
	input: "output/{prefix}/rnaSeq/abundance.h5"
	output: "output/{prefix}/final/expressionMatrix.tsv"
	log: "output/{prefix}/log/expressionMatrix.log"
	shell:
		"""
		scripts/expressionMatrix.R -p {wildcards.prefix} -a {input} -o output/{wildcards.prefix}/rnaSeq &> {log}
		mv output/{wildcards.prefix}/rnaSeq/expressionMatrix.tsv {output}
		"""

# =======================================
# Phasing
# =======================================

rule convertToPhasedVCF:
	input:
		alleleOne = "output/{prefix}/genome/allele1.refseq2simseq.SNP.vcf",
		alleleTwo = "output/{prefix}/genome/allele2.refseq2simseq.SNP.vcf"
	output:
		alleleOne = temp("output/{prefix}/genome/allele1.phased.vcf"),
		alleleTwo = temp("output/{prefix}/genome/allele2.phased.vcf")
	log: "output/{prefix}/log/convertToPhasedVCF.log"
	shell:	
		"""
		awk -v OFS="\t" '/^##/ {{print $0;next}} /^#/{{print $0, "FORMAT", "SAMPLE"; next}} \
			{{print $1, $2, $3, $4, $5, $6, "PASS", $8, "GT", "0|1"}}' {input.alleleOne} > {output.alleleOne} 2> {log}
		awk -v OFS="\t" '/^##/ {{print $0;next}}  /^#/{{print $0, "FORMAT", "SAMPLE"; next}} \
			{{print $1, $2, $3, $4, $5, $6, "PASS", $8, "GT", "1|0"}}' {input.alleleTwo} > {output.alleleTwo} 2> {log}
		"""	

rule mergePhasedVCF:
	input:
		alleleOne = "output/{prefix}/genome/allele1.phased.vcf",
		alleleTwo = "output/{prefix}/genome/allele2.phased.vcf"
	output:
		"output/{prefix}/final/phased.vcf"
	params:
		index = config["picardIndex"] 
	singularity: "docker://quay.io/biocontainers/picard:2.27.1--hdfd78af_0"
	log: "output/{prefix}/log/mergePhasedVCF.log"
	shell:
		"""
		picard MergeVcfs \
        		I={input.alleleOne} \
        		I={input.alleleTwo} \
        		O={output} \
			D={params.index} 2> {log}
		"""
