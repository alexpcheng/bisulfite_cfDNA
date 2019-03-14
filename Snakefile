## CONFIG FILES ##
configfile: 'config.yaml'

## PARAMETERS AND FILE PATHS ##
ADAPTOR_SEQUENCES = config['ADAPTOR_SEQUENCES']
METHREF = config['METHREF']
GENOMEREF = config['GENOMEREF']
BWAMETH = config['BWAMETH']
POM = config['POM']
TRIM = config['TRIM']
AUTOSOMALCHROMO = config['AUTOSOMALCHROMO']
METHPIPETOMR = config['METHPIPETOMR']
METHPIPEBSRATE = config['METHPIPEBSRATE']

## RULES ##
# all: list all final output files for downstream analysis
rule all:
	input:
	#from fqc
		expand('V1/fastqc/{sample}{read}_fastqc.html', sample=config["SAMPLES"], read=config["READS"]),
	#from trim
		# only creates temp files
	#from alignment
		expand('V1/aligned/{sample}{bwaext}', sample=config["SAMPLES"], bwaext=[".bam", "_mapped_all_chr.bam", "_mapped_autosomal.bam"]),
		expand('V1/aligned/stats/{sample}{ext}', sample=config["SAMPLES"], ext=[".flagstat", ".numreads"]),
		expand('V2/Lengths/{sample}_aligned.lengths', sample=config["SAMPLES"]),
	#from BSconversion
		expand('V2/conversion_rates/{sample}.bsrate.txt', sample=config["SAMPLES"]),
	#from tidyUP
		"V2/conversion_rates/conversion_summary.txt",
		"V2/depth/coverage_summary.txt",
	#from HQ
		expand('V1/aligned/HQ/{sample}_mapped_autosomal.hq{ext}', sample=config["SAMPLES"], ext=[".bam", "_CpG.bedGraph"]),
		expand('V1/depth/{sample}.depth', sample=config["SAMPLES"]),
		expand('V1/binned_samples/{sample}', sample=config["SAMPLES"])

# fqc: fastQC report on R1 and R2 reads for each sample
rule fqc:
	input:
		r1 = 'Data/{sample}_R1.fastq.gz',
		r2 = 'Data/{sample}_R2.fastq.gz'
	output:
		"V1/fastqc/{sample}_R1_fastqc.html",
		"V1/fastqc/{sample}_R2_fastqc.html"
	threads: 2
	shell:
		"""
		fastqc {input.r1} {input.r2} -t {threads} --outdir V1/fastqc
		"""

# trim: using TRIMMOMATIC software to remove adaptor sequences
rule trim:
	input:
	       	r1 = 'Data/{sample}_R1.fastq.gz',
        	r2 = 'Data/{sample}_R2.fastq.gz'
	output:
        	r1p = temp("V1/{sample}_R1_trim.fastq"),
        	r2p = temp("V1/{sample}_R2_trim.fastq"),
		r1u = temp("V1/{sample}_R1_unpaired.fastq"),
		r2u = temp("V1/{sample}_R2_unpaired.fastq")
	threads: 3
	log: 'logs/trim/{sample}.trim.log'
	shell:
		"""
		java -jar {TRIM} PE -phred33 -threads {threads} \
                	{input.r1} {input.r2} \
			{output.r1p} {output.r1u} {output.r2p} {output.r2u} \
                	{ADAPTOR_SEQUENCES} &>{log}
		"""

rule alignment:
	input:
                r1p = "V1/{sample}_R1_trim.fastq",
                r2p = "V1/{sample}_R2_trim.fastq"
	output:
                bam = "V1/aligned/{sample}.bam",
		mapped_all_chr="V1/aligned/{sample}_mapped_all_chr.bam",
		mapped_autosomal= "V1/aligned/{sample}_mapped_autosomal.bam",
		flagstat="V1/aligned/stats/{sample}.flagstat",
                numreads="V1/aligned/stats/{sample}.numreads",
		lengths=temp("V2/Lengths/{sample}_aligned.lengths")
	threads: 5
	log:'logs/alignment/{sample}.alignment.log'
	shell:
		"""
		({BWAMETH} --reference {METHREF} -t {threads} {input.r1p} {input.r2p} |
		samtools view -@ {threads} -b -F 1024,512 - > {output.bam}) &>{log}

		samtools view -@ {threads} -b -f 3 -F 256 {output.bam} |
			samtools sort -n -@ {threads} -o - | samtools fixmate -@ {threads} -m - - |
			samtools sort -@ {threads} -o - | samtools markdup -r - {output.mapped_all_chr}

		samtools index {output.mapped_all_chr}

		samtools view -@ {threads} -b {output.mapped_all_chr} {AUTOSOMALCHROMO} | samtools sort -@ {threads} - -o {output.mapped_autosomal}
		samtools index {output.mapped_autosomal}

		samtools view {output.mapped_autosomal} | awk '{{if ($9>0 && $9<1000) print $9}}' > {output.lengths}
                Rscript Bin/HistogramFragmentLengths.R {output.lengths} V2/Lengths/{wildcards.sample}_FragsHistogram.pdf
                samtools flagstat {output.mapped_autosomal} >  {output.flagstat}
		NUM_MAPPED=`cat {output.flagstat} | grep "mapped (" | cut -f1 -d' '`
                echo $NUM_MAPPED > {output.numreads}

		mkdir -p V1/aligned/downsampled
		"""

rule BSconversion:
        input:
               #'V1/aligned/HQ/{sample}_mapped_autosomal.hq.bam'
              	"V1/aligned/{sample}_mapped_autosomal.bam"
        output:
               	mr=temp("V1/{sample}.mr"),
                bsrate="V2/conversion_rates/{sample}.bsrate.txt"
        threads: 1
        shell:
              	"""
                {METHPIPETOMR} -m general -o {output.mr} -L 500 {input[0]}
                {METHPIPEBSRATE} -c {GENOMEREF} -o {output.bsrate} {output.mr}
                X=$(head -1 {output.bsrate})
                echo -e "{wildcards.sample}\t$X" >> V2/conversion_rates/summary.txt
                """
rule tidyUP:
	input:
		expand("V2/conversion_rates/{sample}.bsrate.txt", sample=config["SAMPLES"]),
		expand("V1/depth/{sample}.depth", sample=config["SAMPLES"])
	output:
		conversion="V2/conversion_rates/conversion_summary.txt",
		coverage="V2/depth/coverage_summary.txt"
	threads: 1
	shell:
		"""
		sort V2/conversion_rates/summary.txt > {output.conversion}
		cat V1/depth/*.depth | sort > {output.coverage}
		"""

rule HQ:
	input:
		allQ='V1/aligned/{sample}_mapped_autosomal.bam',
		refFile='Methylation_References_snaked/MethylMatrix/all_chr/all_groups_merged_overlap.metilene'
	output:
		HQ='V1/aligned/HQ/{sample}_mapped_autosomal.hq.bam',
		CpGHQ_bg='V1/aligned/HQ/{sample}_mapped_autosomal.hq_CpG.bedGraph',
		CpGHQ_bg_tmp=temp('V1/aligned/HQ/{sample}_mapped_autosomal.hq_CpG.bedGraph.tmp'),
		binned_CpG_DMR='V1/binned_samples/{sample}',
		depth='V1/depth/{sample}.depth'
	threads: 1
	shell:
		"""
		file={output.HQ}
		out=${{file::-4}}
		samtools view -b -q 40 {input.allQ} > {output.HQ}
                {POM} extract {METHREF} {output.HQ} -o $out

		tail -n +2 {output.CpGHQ_bg} | awk -F'\t' '{{print $1,$2,$3,NR,$4}}' | sort-bed - > {output.CpGHQ_bg_tmp}

		bedmap --echo --fraction-map 1 --mean --delim '\t' \
			{input.refFile} \
			{output.CpGHQ_bg_tmp} > {output.binned_CpG_DMR}

		depth=$(samtools depth {output.HQ} | awk '{{sum+=$3}} END {{print sum/2684573069}}')
		echo -e "{wildcards.sample}\t$depth" > {output.depth}
		"""
