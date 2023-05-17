configfile: "config/refPaths.yaml"
configfile: "config/params.yaml"
rule all:
	input:
		expand("output/{name}/final/{file}", 
			name = config["sample"],
			file = ["rnaReadAlignment.bam", "phased.vcf.gz", "expression_matrix.tsv"]),
		expand("output/{name}/params.txt", 
			name = config["sample"])

# ===========================================
# Save parameter
# ===========================================
rule saveParams:
	output: "output/{prefix}/params.txt"
	params:
		normalSNP = config["snps"]["normalSNPcount"],
		tumorSNP = config["snps"]["tumorSNPcount"],
		numASE = config["numASE"],
		tumorContent = config["tumorContent"],
		rnaReadLength = config["rnaReads"]["readLength"],
		rnaPercentExpressed = config["rnaReads"]["percentExpressed"],
		rnaDepth = config["rnaReads"]["depth"],
		rnaASEfoldChange = config["rnaReads"]["ASEfoldchange"],
		simulatePhasing = config["simulatePhasing"],
		wasp = config["WASPfilter"],
		longRead = config["longRead"]["numReads"],
		refGenome = config["ref_genome"],
		nanosim_model = config["nanosim_model"]
	log: "output/{prefix}/log/params.txt"
	shell:
		"""
		date >> {output}
		printf "{wildcards.prefix}\n------------------------\n\n" >> {output} 2> {log}
		echo "Reference Genome: {params.refGenome}" >> {output}
		echo "" >> {output}
		echo "Number of SNP in normal: {params.normalSNP}" >> {output}
		echo "Number of SNP in tumor: {params.tumorSNP}" >> {output}
		echo "Tumor Content: {params.tumorContent}" >> {output}
		echo "" >> {output}
		echo "Number of ASE genes: {params.numASE}" >> {output}
		echo "ASE fold change: {params.rnaASEfoldChange}" >> {output}
		echo "" >> {output}
		echo "RNA Read Parameter:" >> {output}
		echo "    Read Length: {params.rnaReadLength}" >> {output}
		echo "    Percentage gene expressed: {params.rnaPercentExpressed}" >> {output}
		echo "    RNA depth: {params.rnaDepth}" >> {output}
		echo "Long Read parameter:" >> {output}
		echo "    Number of reads: {params.longRead}" >> {output}
		echo "    Nanosim model: {params.nanosim_model}" >> {output}
		echo "" >> {output}
		echo "Simulate phasing: {params.simulatePhasing}" >> {output}
		echo "Use WASP filter: {params.wasp}" >> {output}
		"""



# ===========================================
# Generate genomes of both alleles
# ===========================================
rule generatingRefGenome:
	input:
		config["ref_genome"]
	output: 
		genome = "output/{prefix}/genome/ref_allele{num}.fa",
		vcf = "output/{prefix}/genome/ref_allele{num}.refseq2simseq.SNP.vcf"
	params: 
		snp_count=config["snps"]["normalSNPcount"]
	log: "output/{prefix}/log/generatingRefGenome_{num}.log"	
	singularity: "docker://quay.io/biocontainers/simug:1.0.0--hdfd78af_1"
	shell:
		"""
		mkdir -p output/{wildcards.prefix}/genome
	
		simuG \
			-refseq {input} \
			-snp_count {params.snp_count} \
			-prefix output/{wildcards.prefix}/genome/ref_allele{wildcards.num} &> {log}
		mv output/{wildcards.prefix}/genome/ref_allele{wildcards.num}.simseq.genome.fa {output.genome}
		"""	


rule generatingAltGenome:
	input:
		"output/{prefix}/genome/ref_allele{num}.fa"
	output: 
		genome = "output/{prefix}/genome/alt_allele{num}.fa", 
		vcf = "output/{prefix}/genome/alt_allele{num}.refseq2simseq.SNP.vcf"
	params: 
		snp_count=config["snps"]["tumorSNPcount"]
	log: "output/{prefix}/log/generatingAltGenome_{num}.log"
	singularity: "docker://quay.io/biocontainers/simug:1.0.0--hdfd78af_1"
	shell:
		"""
		mkdir -p output/{wildcards.prefix}/genome
	
		simuG \
			-refseq {input} \
			-snp_count {params.snp_count} \
			-prefix output/{wildcards.prefix}/genome/alt_allele{wildcards.num} &> {log}
	
		mv output/{wildcards.prefix}/genome/alt_allele{wildcards.num}.simseq.genome.fa {output.genome}
		"""	

rule snpRegion:
	input: 
		expand("output/{{prefix}}/genome/alt_allele{num}.refseq2simseq.SNP.vcf", num = [1, 2])
	output: 
		"output/{prefix}/genome/SNP.txt"
	params:
		geneBed = config["gene_annotation"]
	singularity: "docker://quay.io/biocontainers/bedtools:2.23.0--h5b5514e_6"
	log: "output/{prefix}/log/snpRegion.log"
	shell:
		"""
		bedtools intersect -u -a {params.geneBed} -b {input} > {output} 2> {log}
		"""

#########################
# Generate allelic read 
#########################

rule getTranscriptome:
	input:
		genome = "output/{prefix}/genome/{genomeType}_allele{num}.fa",
		gtf = config["annotation_gtf"]
	output: temp("output/{prefix}/transcriptome/{genomeType}_allele{num}.cDNA.fa")
	singularity: "docker://quay.io/biocontainers/gffread:0.12.7--hd03093a_1"
	log: "output/{prefix}/log/getTranscriptome_{genomeType}_{num}.log"
	shell:
		"""
		mkdir -p output/{wildcards.prefix}/transcriptome
		gffread \
			--no-pseudo \
			-w {output} \
			-g {input.genome} {input.gtf}
		"""

rule cutNTranscriptome:
	input: "output/{prefix}/transcriptome/{genomeType}_allele{num}.cDNA.fa"
	output: temp("output/{prefix}/transcriptome/{genomeType}_allele{num}.cDNA.cutN.fa")
	log: "output/{prefix}/log/cutNTranscriptome_{genomeType}_{num}.log"
	shell:
		"sed '/^>/! s/[^AGCT]//g' {input} > {output} 2> {log}"

rule generateExpressionProfile:
	input: 
		genome = "output/{prefix}/transcriptome/ref_allele1.cDNA.cutN.fa",
		Snp = "output/{prefix}/genome/SNP.txt"
	output: 
		"output/{prefix}/simRNA/allele1/alt_expressionProfile.rds",
		"output/{prefix}/simRNA/allele2/alt_expressionProfile.rds",
		"output/{prefix}/simRNA/allele1/ref_expressionProfile.rds",
		"output/{prefix}/simRNA/allele2/ref_expressionProfile.rds",
		"output/{prefix}/simRNA/expressedGene.tsv"
	params:
		numASE = config["numASE"],
		foldChange = config["rnaReads"]["ASEfoldchange"],
		percentExpressed = config["rnaReads"]["percentExpressed"],
		depth = config["rnaReads"]["depth"],
		tumorContent = config["tumorContent"],
		readLength = config["rnaReads"]["readLength"],
		gene2transcript = "ref/ensembl100_transcript2gene.tsv"
	log: "output/{prefix}/log/generateExpressionProfile.log"
	singularity: "docker://glenn032787/asesim_rcontainer:1.0"
	shell:
		"""
		mkdir -p output/{wildcards.prefix}/simRNA/allele1
		mkdir -p output/{wildcards.prefix}/simRNA/allele2
		mkdir -p output/{wildcards.prefix}/final

		Rscript scripts/polyesterGenerateProfile.R \
			-t {input.genome} \
			-n {params.numASE} \
			-f {params.foldChange} \
			-e {params.percentExpressed} \
			-d {params.depth} \
			-c {params.gene2transcript} \
			-s {input.Snp} \
			-u {params.tumorContent} \
			-l {params.readLength} \
			-o output/{wildcards.prefix}/simRNA &> {log}
		
		cp output/{wildcards.prefix}/simRNA/expressedGene.tsv output/{wildcards.prefix}/final/expressedGene.tsv
		"""

rule simulateReadPolyester: 
	input: 
		expressionProfile = "output/{prefix}/simRNA/allele{num}/{genomeType}_expressionProfile.rds",
		transcriptome = "output/{prefix}/transcriptome/{genomeType}_allele{num}.cDNA.cutN.fa"
	output: 
		temp("output/{prefix}/simRNA/allele{num}/{genomeType}/sample_01_1.fasta"),
		temp("output/{prefix}/simRNA/allele{num}/{genomeType}/sample_01_2.fasta")
	log: "output/{prefix}/log/simulateReadPolyester_{genomeType}_{num}.log"
	params:
		readLength = config["rnaReads"]["readLength"]
	singularity: "docker://glenn032787/asesim_rcontainer:1.0"
	shell:
		"""
		mkdir -p output/{wildcards.prefix}/simRNA/allele{wildcards.num}/{wildcards.genomeType}
		Rscript scripts/polyester.R \
			-t {input.transcriptome} \
			-e {input.expressionProfile} \
			-r {params.readLength} \
			-o output/{wildcards.prefix}/simRNA/allele{wildcards.num}/{wildcards.genomeType} &> {log}
		"""

if config["tumorContent"] == 1.0:
	tc = ["alt"]
else:
	tc = ["ref", "alt"]

rule catReads:
	input: 
		R1 = expand("output/{{prefix}}/simRNA/allele{num}/{type}/sample_01_1.fasta", num=[1,2], type=tc),
		R2 = expand("output/{{prefix}}/simRNA/allele{num}/{type}/sample_01_2.fasta", num=[1,2], type=tc)
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
	
if not config["WASPfilter"]:
	rule kalliso:
		input:
			index = config["kalliso_index"],
			gtf_annotation = config["annotation_gtf"],
			chrom_length = config["chrom_length"],
			reads = expand("output/{{prefix}}/rnaSeq/rnaReads_R{num}.fq", num=[1,2])
		output: 
			"output/{prefix}/final/rnaReadAlignment.bam",
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
			mv output/{wildcards.prefix}/rnaSeq/pseudoalignments.bam output/{wildcards.prefix}/final/rnaReadAlignment.bam
			mv output/{wildcards.prefix}/rnaSeq/pseudoalignments.bam.bai output/{wildcards.prefix}/final/pseudoalignments.bam.bai
			"""

	rule expressionMatrix:
		input: "output/{prefix}/rnaSeq/abundance.h5"
		output: "output/{prefix}/final/expression_matrix.tsv"
		log: "output/{prefix}/log/expressionMatrix.log"
		singularity: "docker://glenn032787/asesim_rcontainer:1.0"
		shell:
			"""
			Rscript scripts/expressionMatrix.R -p {wildcards.prefix} -a {input} -o output/{wildcards.prefix}/rnaSeq &> {log}
			mv output/{wildcards.prefix}/rnaSeq/expressionMatrix.tsv {output}
			"""
# =======================================
# STAR alignment
# =======================================
if config["WASPfilter"]:
	rule star:
		input:
			r1="output/{prefix}/rnaSeq/rnaReads_R1.fq",
			r2="output/{prefix}/rnaSeq/rnaReads_R2.fq",
			phased_vcf="output/{prefix}/final/phased.vcf.gz"
		output:
			"output/{prefix}/star/starAligned.sortedByCoord.out.bam",
			temp("output/{prefix}/star/starAligned.toTranscriptome.out.bam")
		singularity: "docker://quay.io/biocontainers/star:2.7.10a--h9ee0642_0"
		params:
			ref="/projects/glchang_prj/ref/star/hg38"
		log: "output/{prefix}/log/star.log"
		threads: 20
		shell:
			"""
			mkdir -p output/{wildcards.prefix}/star
			zcat {input.phased_vcf} > output/{wildcards.prefix}/star/phased.vcf
			STAR \
				--genomeDir {params.ref} \
				--runThreadN {threads} \
				--readFilesIn  {input.r1} {input.r2} \
				--outFileNamePrefix output/{wildcards.prefix}/star/star \
				--outSAMtype BAM SortedByCoordinate \
				--outSAMunmapped Within \
				--outSAMattributes Standard \
				--waspOutputMode SAMtag \
				--varVCFfile output/{wildcards.prefix}/star/phased.vcf \
				--quantMode TranscriptomeSAM \
				--twopassMode Basic \
				--twopass1readsN -1 &> {log}
			"""


	rule waspFilter:
		input: "output/{prefix}/star/starAligned.sortedByCoord.out.bam"
		output: "output/{prefix}/final/rnaReadAlignment.bam"
		singularity: "docker://quay.io/biocontainers/samtools:1.15.1--h1170115_0"
		log: "output/{prefix}/log/waspFilter.log"
		shell:
			"""
			samtools view -h {input} | grep -e '^@' -e 'vW:i:1' | samtools view -b -S > {output} 2> {log}
			samtools index {output} &> {log}
			"""

	rule rsem:
		input: "output/{prefix}/star/starAligned.toTranscriptome.out.bam"
		output: "output/{prefix}/rsem/star.genes.results"
		params:
			ref="/projects/glchang_prj/ref/RSEM/hg38_no_alt"
		singularity: "docker://quay.io/biocontainers/rsem:1.3.3--pl5321hecb563c_4"
		threads: 72
		log: "output/{prefix}/log/rsem.log"
		shell:
			"""
			mkdir -p output/{wildcards.prefix}/rsem
			rsem-calculate-expression \
				--alignments \
				-p {threads} \
				--paired-end \
				{input} \
				{params.ref} \
				output/{wildcards.prefix}/rsem/star &> {log}
			"""

	rule rsemExpressionMatrix:
		input: "output/{prefix}/rsem/star.genes.results"
		output: "output/{prefix}/final/expression_matrix.tsv"
		log: "output/{prefix}/log/rsemExpressionMatrix.log"
		singularity: "docker://glenn032787/asesim_rcontainer:1.0"
		shell:
			"Rscript scripts/generateExpressionMatrix.R -i {input} -o output/{wildcards.prefix}/final -s {wildcards.prefix} &> {log}"


# =======================================
# Phasing
# =======================================
if not config["simulatePhasing"]:
	rule convertToPhasedVCF:
		input:
			alleleOne = "output/{prefix}/genome/ref_allele1.refseq2simseq.SNP.vcf",
			alleleTwo = "output/{prefix}/genome/ref_allele2.refseq2simseq.SNP.vcf"
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
			"output/{prefix}/final/phased.vcf.gz"
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

# ====================================
# Nanopore phasing
# ==================================== 
if config["simulatePhasing"]:
	rule nanosim:
		input: "output/{prefix}/genome/ref_allele{num}.fa"
		output: 
			temp("output/{prefix}/longRead/allele{num}_aligned_reads.fasta"),
			temp("output/{prefix}/longRead/allele{num}_unaligned_reads.fasta")
		threads: 35
		params:
			numReads = config["longRead"]["numReads"],
			model = config["nanosim_model"]
		singularity: "docker://glenn032787/nanosim:1.5"
		log: "output/{prefix}/log/nanosim_{num}.log"
		shell:
			"""
			mkdir -p output/{wildcards.prefix}/longRead

			simulator.py genome \
				-rg {input} \
				-o output/{wildcards.prefix}/longRead/allele{wildcards.num} \
				-n {params.numReads} \
				-t {threads} \
				-b guppy \
				-c {params.model} &> {log}
				
			"""

	rule catLongRead:
		input: expand("output/{{prefix}}/longRead/allele{num}_{type}_reads.fasta", num = [1,2], type=["aligned", "unaligned"])
		output: temp("output/{prefix}/longRead/nanoporeReads.fa")
		shell:
			"cat {input} > {output}"


	rule longReadAlignment:
		input: "output/{prefix}/longRead/nanoporeReads.fa"
		output: temp("output/{prefix}/longRead/nanoporeAlignment.sam")
		singularity: "docker://quay.io/biocontainers/minimap2:2.24--h7132678_1"
		params:
			ref = "/gsc/resources/Homo_sapiens_genomes/hg38_no_alt/genome/fasta/hg38_no_alt.fa"
		log: "output/{prefix}/log/longReadAlignment.log"
		threads: 3
		shell:
			"""
			minimap2 \
				-ax map-ont \
				-t {threads} \
				{params.ref} \
				{input} > {output} 2> {log}
			"""

	rule longReadSamtoBam:
		input: "output/{prefix}/longRead/nanoporeAlignment.sam"
		output: 
			sorted_bam = "output/{prefix}/longRead/nanoporeAlignment.sorted.bam",
			bam = temp("output/{prefix}/longRead/nanoporeAlignment.bam")
		singularity: "docker://quay.io/biocontainers/sambamba:0.8.2--h98b6b92_2" 
		threads: 72
		log: "output/{prefix}/log/longReadSamToBam.log"
		shell: 
			"""	
			sambamba view -S -t {threads} -f bam {input} >  {output.bam} 2> {log}
			sambamba sort -t {threads} {output.bam} &> {log}
			"""

	rule clair3:
		input: "output/{prefix}/longRead/nanoporeAlignment.sorted.bam"
		output: 
			"output/{prefix}/final/phased.vcf.gz"
		threads: 72
		singularity: "docker://hkubal/clair3:latest"
		params:
			ref = config["ref_genome"]
		log: "output/{prefix}/log/clair3.log"
		shell:
			"""
			mkdir -p output/{wildcards.prefix}/clair3
			output=$(realpath output/{wildcards.prefix}/clair3)
			input=$(realpath {input})

			run_clair3.sh \
				--bam_fn=${{input}} \
				--ref_fn={params.ref} \
				--threads={threads} \
				--platform="ont" \
				--model_path=/opt/models/r941_prom_sup_g5014 \
				--output=${{output}} \
				--enable_phasing &> {log}
			
			cp output/{wildcards.prefix}/clair3/phased_merge_output.vcf.gz output/{wildcards.prefix}/final/phased.vcf.gz
			"""











