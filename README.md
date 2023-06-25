## Module 1 of the Snakemake RNA Pipeline for H. sapiens

* needs a working conda installation and conda in path

* Snakemake needs to be installed in the working enviroment

* sample_working.tsv provides a example of the needed sampleseheet

* use the module like this:

```
$ cd rna_module_1
```

the current working dirtectory needs to be the rna_module_1 folder 

```
snakemake --use-conda --cores 8 --conda-prefix "folder where conda should store the envs data (is optional)"

```


* the pipeline creates a folder structure and writes all the output in the parent directory of the current working directory 

#Folder Structure
```
.
├── 01_raw
├── 02_trimmed
├── 03_mapped
└── rna_module_1
    ├── data
    │   ├── NexteraPE-PE.fa
    │   ├── TruSeq2-PE.fa
    │   ├── TruSeq2-SE.fa
    │   ├── TruSeq3-PE-2.fa
    │   ├── TruSeq3-PE.fa
    │   └── TruSeq3-SE.fa
    ├── envs
    │   ├── bedtools.yaml
    │   ├── fastqc.yaml
    │   ├── mapping_viz.yaml
    │   ├── sambamba.yaml
    │   ├── samtools.yaml
    │   ├── samtools.yaml.save
    │   ├── star.yaml
    │   ├── subread.yaml
    │   └── trimmomatic.yaml
    ├── samples.tsv
    ├── sample_working.tsv
    ├── scripts
    │   ├── mapping.py
    │   └── sub_sampleVIZ.py
    └── Snakefile
