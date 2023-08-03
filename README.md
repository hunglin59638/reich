# Reich &middot; [![GitHub license](https://img.shields.io/badge/license-GPL3-brightgreen.svg)](https://github.com/chanzuckerberg/czid-web/blob/master/LICENSE) ![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg)  
<img align="left" width="150" height="60" src="pics/cover.png">  
  

## Introduction
Reich is a pathogen detection pipeline, inspired by [CZID](https://czid.org/).  
It has backend and fronted modules. The backend module is based on [Nextflow](https://www.nextflow.io/), and the frontend module is based on Rshiny.

Currently, only the backend module is available, and the frontend module is still under development.

**This project is intended for learning and educational purposes only. It is not suitable for clinical or research applications.**  


## Installation
For local installation, follow these steps:

```bash
git clone https://github.com/hunglin59638/reich.git``
cd reich
```

you can install the pipeline by using conda:
```bash
conda env create -f backend/environment.yml
```
Alternatively, you can install from the lock file using conda:
```bash
conda env create -f backend/conda.lock
```

## Usage

First, you need to build the database for pathogen detection. Using the following command to build the database:

```bash
backend/script/build_db.py --out_dir /path/to/reich_db
```
It will download NT blastdb and convert to fasta. Then, use human nt fasta to build bowtie2 index.   
The fasta excluded human nt will be used to non-human pathogen detection.
Hisat2 index for human removal will be download from [hisat2](https://daehwankimlab.github.io/hisat2/download/).

It will take about 4-6 hours to download and build the database.  
The database size is about 180G. Please make sure you have enough disk space.

Expected output:
```
/path/to/reich_db
├── human_nt.fna
├── non_human_nt.fna
├── blastdb
├── human
│   ├── bowtie
│   └── grch38_tran_hisat
└── taxdump
```

Then, you can run the pipeline by using the following command:

```bash
nextflow run backend/main.nf --reads '/path/to/reads/*.fastq.gz' --db_dir /path/to/reich_db --out_dir /path/to/output --threads 12
```

The memory usage is about 70-80G. Please make sure you have enough memory.



## Output

Expected output:
```
/path/to/output
├── report/${sample_name}.report.tsv
```

The output is a tsv file, which contains the following columns:  

| Column | Description |
| --- | --- |
|Taxon|Taxon name|
|Score|Aggregate Score|
|Z score| Z score|
|rPM| Reads per million|
|r|reads mapped to the taxon|
|%id| average identity|
|L| average alignment length|

The detail of the columns can be found [here](https://chanzuckerberg.zendesk.com/hc/en-us/articles/360034790574-Single-Sample-Report-Table#score).