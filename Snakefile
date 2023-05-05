ext=['html', 'zip']
SAMPLES=['R1_001', 'R2_001']

rule all:
    input:
        "fastqc_raw/multiqc_report.html"
        "fastqc_raw_trimmed/multiqc_report.html"
        "counts.txt"

rule fastqc_raw:
    input: "data/{sample}.fastq"
    output: 
        "fastqc_raw/{sample}_fastqc.html",
        "fastqc_raw/{sample}_fastqc.zip"
    shell:'''
    FastQC/fastqc -o fastqc_raw {input}
    '''
rule trim_barcode:
    input:
        r1 = "data/R1_001.fastq",
        r2 = "data/R2_001.fastq",
    output:
        r1 = "trimmed/R1_001.fastq",
        r2 = "trimmed/R2_001.fastq",
    shell:
        "bbduk.sh in1={input.r1} in2={input.r2} out1={output.r1} out2={output.r2} ref=data/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=10"
 
rule multiqc_raw:
    input: expand("fastqc_raw/{sample}_fastqc.html", sample = SAMPLES)
    output: "fastqc_raw/multiqc_report.html"
    shell:'''
    multiqc -o fastqc_raw fastqc_raw
    '''

rule fastqc_raw_trim:
    input: "trimmed/{sample}.fastq"
    output:
        "fastqc_raw_trimmed/{sample}_fastqc.html",
        "fastqc_raw_trimmed/{sample}_fastqc.zip"
    shell:'''
    FastQC/fastqc -o fastqc_raw_trimmed {input}
    '''

rule multiqc_raw_trim:
    input: expand("fastqc_raw_trimmed/{sample}_fastqc.html", sample = SAMPLES)
    output: "fastqc_raw_trimmed/multiqc_report.html"
    shell:'''
    multiqc -o fastqc_raw_trimmed fastqc_raw_trimmed
    '''

rule star_index:
    input:
        "data/chr19_20Mb.fa",
        "data/chr19_20Mb.gtf",
        "trimmed/R1_001.fastq",
        "trimmed/R2_001.fastq"
    output:
        bam = directory("genome")
    shell:
        "mkdir -p {output} && "
        "STAR --runThreadN 1 "
        "--runMode genomeGenerate "
        "--genomeDir {output} "
        "--genomeFastaFiles data/chr19_20Mb.fa "
        "--sjdbGTFfile data/chr19_20Mb.gtf "
        "--genomeSAindexNbases 11 "
        "--readFilesIn trimmed/R1_001.fastq trimmed/R2_001.fastq "

rule star_align:
    input:
        "genome",
        "trimmed/R1_001.fastq",
        "trimmed/R2_001.fastq"
    output:
        directory("aligned")
    shell:
        "mkdir -p aligned && "
        "STAR --runThreadN 1 "
        "--genomeDir {input[0]} "
        "--readFilesIn {input[1]} {input[2]} "
        "--outSAMtype BAM SortedByCoordinate "
        "--outFileNamePrefix aligned/ "
        "&& samtools index aligned/Aligned.sortedByCoord.out.bam"

rule featurecounts:
    input:
        "aligned/Aligned.sortedByCoord.out.bam",
        "data/chr19_20Mb.gtf"
    output:
        "counts.txt"
    params:
        strand = "1"
    shell:
        "featureCounts -p -t exon -g gene_id -a {input[1]} -o {output} {input[0]} -s {params.strand}"
        