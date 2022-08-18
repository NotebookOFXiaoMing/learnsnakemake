#configfile: "singleconfig.yaml"
workdir: "/home/myan/scratch/private/practice_data/popgenomics/rubi.gbs"

SRR,FRR = glob_wildcards("00.raw.fq.paired/" + "{srr}_{frr}.fastq.gz")

print(SRR)

rule all:
    input:
        #expand("01.fastp.filtered.single/" + "{srr}_clean.fastq.gz",srr=SRR),
        #"reference/genome_index/Pr.1.bt2",
        #expand("02.sam/" + "{srr}.sam",srr=SRR),
        expand("02.sorted.bam/" + "{srr}.sorted.bam.bai",srr=SRR)


rule a_runfastp:
    input:
        read01 = "00.raw.fq.paired/" + "{srr}_1.fastq.gz",
        read02 = "00.raw.fq.paired/" + "{srr}_2.fastq.gz"
    output:
        read01 = "01.fastp.filtered.paired/" + "{srr}_clean_1.fastq.gz",
        read02 = "01.fastp.filtered.paired/" + "{srr}_clean_2.fastq.gz",
        json = "01.fastp.report.paired/" + "{srr}.json",
        html = "01.fastp.report.paired/" + "{srr}.html"
    threads:
        8
    params:
        "-f 10"
    shell:
        """
        fastp -i {input.read01} -I {input.read02} -o {output.read01} -O {output.read02} \
        {params} -w {threads} -j {output.json} -h {output.html}
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
        read01 = rules.a_runfastp.output.read01,
        read02 = rules.a_runfastp.output.read02
    output:
        sam = "02.sam/" + "{srr}.sam"
    threads:
        8
    params:
        index = "reference/genome_index/Pr",
        others = "-q --very-sensitive --no-unal --local --rg-id {srr} --rg SM:{srr}"
    shell:
        """
        bowtie2 -x {params.index} -1 {input.read01} -2 {input.read02} -S {output.sam} {params.others} -p {threads}
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