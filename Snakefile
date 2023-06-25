import os
import pandas as pd
from collections import Counter
SAMPLES = []
DIRECTION=["1","2"]
GROUPS=set()
sample_sheet=pd.read_csv("samples.tsv", sep="\t")

count={}
syn=[]
#for a in list(sample_sheet["Type"]):
#    if a not in count:
#        count[a]=1
#        syn.append(a+str(count[a]))
#    elif a in count:
#        count[a]=count[a]+1
#        syn.append(a+str(count[a]))
#sample_sheet["working"]=syn

def check_symlink(file1, file2):
    try:
        os.symlink(file1, file2)
    except FileExistsError:
        print("Link for " + file1 + " is already present in 01_raw")

SAMPLES=list(sample_sheet["ID"])
#WORKING_SAMPLES=list(sample_sheet["working"])
#sample_sheet.to_csv("sample_working.tsv",index=False, sep="\t")


for a in sample_sheet['ID']:
    os.makedirs("../01_raw/" +a, exist_ok=True)
    check_symlink(sample_sheet[sample_sheet["ID"]==a]["forward reads"].values[0], "../01_raw/"+a +"/"+ sample_sheet[sample_sheet["ID"]==a]["ID"].values[0] +"_1P.fastq.gz")
    check_symlink(sample_sheet[sample_sheet["ID"]==a]["reverse reads"].values[0], "../01_raw/"+a +"/"+ sample_sheet[sample_sheet["ID"]==a]["ID"].values[0] +"_2P.fastq.gz")
    os.makedirs("../02_trimmed/" + a +"/fastqc", exist_ok=True)
    os.makedirs("../03_star/"+a, exist_ok=True)
    #os.makedirs("../04_Bedtools_count/"+b, exist_ok=True)
    #os.makedirs("../05_featureCount/sub_count/"+"b", exist_ok=True)
    #os.makedirs("../06_Overview_Mapping", exist_ok=True)

#fastqc_output1=[]
#fastqc_output2=[]
#for a in SAMPLES:
#    fastqc_output1.append("../01_raw/fastqc/"+sample_sheet[sample_sheet["ID"]==a]["forward reads"].values[0].replace(".fastq.gz","_fastqc.html"))
#    fastqc_output1.append("../01_raw/fastqc/"+sample_sheet[sample_sheet["ID"]==a]["reverse reads"].values[0].replace(".fastq.gz","_fastqc.html"))
#    fastqc_output2.append("../02_trimmed/fastqc/" +a+ "_trimmed_1P_fastqc.html")
#    fastqc_output2.append("../02_trimmed/fastqc/" +a+ "_trimmed_2P_fastqc.html")

rule all:
    input:
        #expand("../03_mapped/{sample}/{sample}Aligned.sortedByCoord.out.mkdup.unique.bam",sample=WORKING_SAMPLES),
        #expand("../04_Bedtools_count/{sample}/{sample}_intronCount.txt", sample=WORKING_SAMPLES),
        expand('../03_star/{sample}/{sample}Aligned.sortedByCoord.out.bam', sample=SAMPLES),
        #expand('02_trimmed/{sample}_trimmed_1P.fastq.gz', sample=WORKING_SAMPLES),
        expand('../01_raw/{sample}/fastqc/{sample}_{direction}P_fastqc.html',direction=DIRECTION, sample=SAMPLES),
        expand('../02_trimmed/{sample}/fastqc/{sample}_trimmed_{direction}P_fastqc.html', direction=DIRECTION, sample=SAMPLES),
        #expand("../05_featureCount/sub_count/{sample}/FeatureCountTable", sample=WORKING_SAMPLES),
        #expand("../03_mapped/{sample}/{sample}_sub{sub}.bam",sub=SUBSAMPLE, sample=WORKING_SAMPLES),
        #fastqc_output2,
        #["../05_featureCount/FeatureCountTable"],
        #["../06_Overview_Mapping/Mapping_overview.html"]

rule fastqc1:
    input:
        r = '../01_raw/{sample}/{sample}_{direction}P.fastq.gz',
    threads: 2
    priority: 50
    output:
        '../01_raw/{sample}/fastqc/{sample}_{direction}P_fastqc.html'
    conda:
        "envs/fastqc.yaml"
    shell:
        'fastqc -o ../01_raw/{wildcards.sample}/fastqc -t {threads} --extract {input.r}'

rule trimming:
    input:
        r1= '../01_raw/{sample}/{sample}_1P.fastq.gz',
        r2= '../01_raw/{sample}/{sample}_2P.fastq.gz',
    output:
        p1="../02_trimmed/{sample}/{sample}_trimmed_1P.fastq.gz",
        u1="../02_trimmed/{sample}/{sample}_trimmed_1U.fastq.gz",
        p2="../02_trimmed/{sample}/{sample}_trimmed_2P.fastq.gz",
        u2="../02_trimmed/{sample}/{sample}_trimmed_2U.fastq.gz",
    conda:
        "envs/trimmomatic.yaml"
    threads: 4
    priority: 50
    shell:
        'trimmomatic PE -threads {threads} \
        {input.r1} {input.r2} \
        {output.p1} {output.u1} \
        {output.p2} {output.u2} \
        ILLUMINACLIP:data/TruSeq3-PE.fa:2:30:10:2:true \
        LEADING:15 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:36 CROP:75'


rule fastqc2:
    input:
        r1 = "../02_trimmed/{sample}/{sample}_trimmed_1P.fastq.gz",
        r2 = "../02_trimmed/{sample}/{sample}_trimmed_2P.fastq.gz"
    threads: 2
    output:
        r1_trim='../02_trimmed/{sample}/fastqc/{sample}_trimmed_1P_fastqc.html',
        r2_trim='../02_trimmed/{sample}/fastqc/{sample}_trimmed_2P_fastqc.html'
    conda:
        "envs/fastqc.yaml"
    priority: 50
    shell:
        'fastqc -o ../02_trimmed/{wildcards.sample}/fastqc -t {threads} --extract {input.r1}; \
         fastqc -o ../02_trimmed/{wildcards.sample}/fastqc -t {threads} --extract {input.r2}'


rule StarMapping:
    input:
        r1 = "../02_trimmed/{sample}/{sample}_trimmed_1P.fastq.gz",
        r2 = "../02_trimmed/{sample}/{sample}_trimmed_2P.fastq.gz",
        Genom="/mnt/ceph/ressources/01_refGenomes/Hsapiens/GRCh37/RefSeq/fasta/star_2.7.9",
        #sampleID="{sample}",
    output:
        "../03_star/{sample}/{sample}Aligned.sortedByCoord.out.bam",
    conda:
        "envs/star.yaml"
    threads:
        8
    priority: 50
    shell:
        'STAR \
            --genomeDir {input.Genom} \
            --readFilesIn {input.r1} {input.r2} \
            --runThreadN {threads} \
            --outFileNamePrefix ../03_star/{wildcards.sample}/{wildcards.sample} \
            --readFilesCommand zcat \
            --alignIntronMax 500000 \
            --alignMatesGapMax 500000 \
            --outBAMcompression 0 \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMprimaryFlag OneBestScore \
            --outFilterMultimapNmax 100 \
            --outFilterMismatchNoverLmax 0.05 \
            --chimSegmentMin 15 \
            --chimOutType WithinBAM \
            --chimScoreMin 1 \
            --chimScoreJunctionNonGTAG 0 \
            --chimJunctionOverhangMin 15 \
            --chimSegmentReadGapMax 3 \
            --alignSJstitchMismatchNmax 5 -1 5 5'
