#configfile: "singleconfig.yaml"
workdir: "/home/myan/scratch/private/practice_data/popgenomics/rubi.gbs"

SRR, = glob_wildcards("00.raw.fq.single/" + "{srr}.fastq.gz")

print(SRR)

rule all:
    input:
        #expand("01.fastp.filtered.single/" + "{srr}_clean.fastq.gz",srr=SRR),
        #"reference/genome_index/Pr.1.bt2",
        #expand("02.sam/" + "{srr}.sam",srr=SRR),
        #expand("02.sorted.bam/" + "{srr}.sorted.bam.bai",srr=SRR),
        "03.bcf/" + "raw.vcf"


rule a_runfastp:
    input:
        read01 = "00.raw.fq.single/" + "{srr}.fastq.gz"
    output:
        read01 = "01.fastp.filtered.single/" + "{srr}_clean.fastq.gz",
        json = "01.fastp.report.single/" + "{srr}.json",
        html = "01.fastp.report.single/" + "{srr}.html"
    threads:
        8
    params:
        "-f 10"
    shell:
        """
        fastp -i {input.read01} -o {output.read01} {params} \
        -w {threads} -j {output.json} -h {output.html}
        """

rule b_bowtie2index:
    input:
        ref = "reference/genome/Pr.fna"
    output:
        index = "reference/genome_index/Pr.1.bt2"
    params:
        "reference/genome_index/Pr"
    threads:
        8
    shell:
        """
        bowtie2-build {input.ref} {params}
        """

rule c_bowtie2align:
    input:
        read01 = rules.a_runfastp.output.read01
    output:
        sam = "02.sam/" + "{srr}.sam"
    threads:
        8
    params:
        index = "reference/genome_index/Pr",
        others = "-q --very-sensitive --no-unal --local --rg-id {srr} --rg SM:{srr}"
    shell:
        """
        bowtie2 -x {params.index} -U {input.read01} -S {output.sam} {params.others} -p {threads}
        """

rule d_samtoolsview:
    input:
        sam = rules.c_bowtie2align.output.sam
    output:
        bam = "02.bam/" + "{srr}.bam"
    threads:
        2
    shell:
        """
        samtools view -@ {threads} -bS -o {output.bam} {input.sam}
        """

rule e_samtoolssort:
    input:
        bam = rules.d_samtoolsview.output.bam
    output:
        sorted = "02.sorted.bam/" + "{srr}.sorted.bam"
    threads:
        2
    shell:
        """
        samtools sort -@ {threads} -O bam {input.bam} -o {output.sorted}
        """

rule f_samtoolsindex:
    input:
        sorted = rules.e_samtoolssort.output.sorted
    output:
        bai = "02.sorted.bam/" + "{srr}.sorted.bam.bai"
    threads:
        2
    shell:
        """
        samtools index {input.sorted}
        """

rule g_bcftools:
    input:
        bam = expand("02.sorted.bam/" + "{srr}.sorted.bam",srr=SRR),
        ref = 'reference/genome/Pr.fna'
    output:
        vcf = "03.bcf/" + "raw.vcf"
    threads:
        8
    shell:
        """
        bcftools mpileup -O b -f {input.ref} --threads 8 -q 20 \
        -Q 30 {input.bam} | bcftools call --ploidy 2 -m -v -o {output.vcf} --threads 8
        """
    