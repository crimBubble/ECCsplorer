# ECCsplorer pipeline #
-------------------------------------------------------------------------------

![fig1_graph_abstract](https://user-images.githubusercontent.com/45664228/112352923-80acf380-8ccb-11eb-94a1-788a917bd89a.png)


## Authors
Ludwig Mann, Kathrin M. Seibt, Beatrice Weber, Tony Heitkam

Institute of Botany, Technische Universität Dresden, 01069 Dresden, Germany

## What is ECCsplorer

The ECCsplorer is a bioinformatics pipeline for the automated detection of extrachromosomal circular DNA from paired-end read data of amplified circular DNA.

![fig2_workflow](https://user-images.githubusercontent.com/45664228/112352695-46dbed00-8ccb-11eb-9a01-e06a44b0203f.png)

## How to use ECCsplorer

ECCsplorer is started over the linux command line interface.

Minimum data requirements

1. Paired-end read data (FASTA/Q) from amplified (phi29 polymerase) circular DNA from the organism of interest and a reference (genome) sequence (FASTA).
or

2. Paired-end read data (FASTA/Q) from amplified (phi29 polymerase) circular DNA from the organism of interest and paired-end read data (FASTA/Q) of a control (either non-amplified, or amplified from different treatment or organism).

Recommended data

3. Paired-end read data (FASTA/Q) from amplified (phi29 polymerase) circular DNA from the organism of interest, paired-end read data (FASTA/Q) of a control (either non-amplified, or amplified from different treatment or organism) and a reference (genome) sequence (FASTA).

### ECCsplorer command line options

```bash
python ECCsplorer.py <read_files> [optional arguments]
```
Short read paired-end sequencing data might be in FASTA or (compressed) FASTQ format in separated files. For usage specify either second/control read data or a reference (genome) sequence (for details see manual):

```bash
python ECCsplorer.py readsA1.fa/q readsA2.fa/q readsB1.fa/q readsB2.fa/q
```
or 
```bash
python ECCsplorer.py readsR1.fa/q readsR2.fa/q -ref sequence.fa
```

The first command will run the clustering module the second command will run the mapping module. Best results are achieved by inclusion of the comparative module by providing control read data and reference sequence (for details see manual), using:

```bash
python ECCsplorer.py readsA1.fa/q readsA2.fa/q readsB1.fa/q readsB2.fa/q -ref sequence.fa [optional arguments] 
```
For an overview on all optional arguments use:

```bash
 python ECCsplorer.py –h/--help.
```
## How to install ECCsplorer

### Installation using conda

```bash
git clone https://github.com/crimBubble/ECCsplorer

cd ECCsplorer
```

Installation of dependencies can be done using a [conda environment](https://docs.conda.io/projects/conda/en/latest/). The required environment can be prepared using the following command (includes dependencies for RepeatExplorer2, solving the environment may take a while):

```bash
conda env create -f environment.yml
```

Activate the environment:

```bash
conda activate eccsplorer
```

Install RepeatExplorer2 following the [manual](https://bitbucket.org/petrnovak/repex_tarean/src/devel/) without(!) creating a new conda environment.

Edit the "Locations of 3rd party tools (PATH)" for Trimmomatic and RepeatExploer2 in *ECCsplorer/lib/config.py* to match your installation. 

### Manual installation

Once you have downloaded the pipeline and installed python packages, R libraries and 3rd party tools, the ECCsplorer.py script can be run in place. Optionally, the pipeline directory can be added to the ```$PATH``` environment variable.

Python (3.4 or higher) is required with the additional packages:

- numpy
- biopython
- scipy
- pyRserve

Install the additional python packages using the pip:
```bash
$cd /.../ECCsplorer # pipeline directory 
$pip install -r requirements.txt
$pip3 install -r requirements.txt 
```

R (3.4.3 or higher) is required with the additional libraries:

- ggplot2
- ggrepel
- grid & gridExtra
- dplyr

Install the additional libraries using: 
```r
install.packages(c("ggplot2", "ggrepel", "grid", "gridExtra", "dplyr"))
```

### Third party tool requiered

Install 3rd party tools following their installation instructions:

- [NCBI blast+ (v2.2.28 or higher)](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) (optional)
    - needed for using -trm/--trim_reads option and using compressed files as input.
- [segemehl (including haarz)](https://www.bioinf.uni-leipzig.de/Software/segemehl/)
- [samtools (1.9 or higher)](https://github.com/samtools/samtools) 
- [bedtools (v2.28.0 or higher)](https://bedtools.readthedocs.io/en/latest/content/installation.html)
- [RepeatExplorer2](https://bitbucket.org/petrnovak/repex_tarean/src/devel/)
- [seqtk](https://github.com/lh3/seqtk) (optional)
    - improves performance for file convertion, but not essential.

After installing 3rd party tools you might need to add them to the ```$PATH``` environment variable or specify the location of their executables in the lib/config.py file (for details see manual). 
The ECCsplorer pipeline will check if all 3rd party tools are available before starting its modules.

## How to cite

Mann L, Seibt KM, Weber B, Heitkam T. (2021). ECCsplorer: a pipeline to detect extrachromosomal circular DNA (eccDNA) from next-generation sequencing data. bioRxiv 2021.06.08.447410 (preprint). https://doi.org/10.1101/2021.06.08.447410 
