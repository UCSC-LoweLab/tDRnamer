## tDRnamer

### About
Transfer RNA, well-studied as the key adapter molecule in protein translation, is now 
recognized as having a plethora of additional regulatory functions when further processed into 
shorter RNAs. Many groups have contributed to this rapidly advancing field, inevitably leading 
to a constellation of diverse yet often overlapping names attributed to tRNA-derived fragments 
(e.g., tRFs, tsRNAs, tiRNAs, SHOT-RNAs, tRNA halves, etc.). In order to promote the growth of 
the research community and aid the recognition of structural/processing/functional 
relationships among the many forms of tRNA-derived RNAs (tDRs), the tRNA community proposes 
a consistent, uniform naming system integrating three main parts: (i) the prefix "tDR" (tRNA-Derived 
RNA), (ii) the portion of the mature tRNA from which the fragment is derived (using 
standardized Sprinzl tRNA base numbering), and (iii) the Genomic tRNA Database Identifier 
(GtRNAdb ID) of the tRNA gene(s) or transcript(s). 

tDRnamer provides a consistent, stable name with annotations for submitted tDR sequences based 
on the described naming system. It consists of a standalone version that can be run on a Linux/Unix 
platform with a large data set including pre-processed small RNA-seq reads (FASTQ file format) 
and tDR sequences (FASTA file format). Researchers who have a smaller set of tDR sequences or 
prefer to use a point-and-click environment can utilize the web-based version at 
[http://trna.ucsc.edu/tDRnamer/](http://trna.ucsc.edu/tDRnamer/) that does not require any software 
installation. Moreover, tDRnamer can take tDR names that are formatted according to our naming 
system and obtain the corresponding sequences with related annotations.

### System requirements
The standalone version of tDRnamer requires to be run on a Linux/Unix system with at least 8 cores and 
16 GB memory. If small RNA sequencing data is used as input, we do not recommend using tDRnamer on a 
regular desktop or laptop. The following dependencies have been tested, use newer versions at 
your own risk.

#### Dependencies
* Python 2.7
* pysam 0.15.3
  * Older versions have a memory leak, make sure you have an updated version
* bowtie2 2.3.5
* NCBI BLAST+ 2.3 or higher
* EMBOSS 6.6
* samtools 1.9
* Infernal 1.1 or higher


### Using Docker Image
To eliminate the need of installing dependencies, you can download the Docker image from our 
[DockerHub repository](https://hub.docker.com/r/ucsclowelab/tDRnamer) using the command
```
docker pull ucsclowelab/tDRnamer
```

### Using Conda Enviroment
In addition to Docker you can alternatively use a Conda environment using the command
```
conda env create -f tDRnamer_env.yaml
```

### Details and Reference Manual
For more details on the tDR naming system and how to use tDRnamer, please check out the 
[user manual](http://trna.ucsc.edu/tDRnamer/docs/).
