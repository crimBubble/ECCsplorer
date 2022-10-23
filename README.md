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

Find detailed instructions on how to use the ECCsplorer pipeline here: [Mini-workshop (usage of ECCsplorer)](https://github.com/crimBubble/ECCsplorer/blob/master/tutorials/Mini-workshop.md)

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

see [Installation instructions](https://github.com/crimBubble/ECCsplorer/blob/master/tutorials/Installation_instructions.md).

## How to cite

Mann, L., Seibt, K.M., Weber, B. , Heitkam, T. ECCsplorer: a pipeline to detect extrachromosomal circular DNA (eccDNA) from next-generation sequencing data. BMC Bioinformatics 23, 40 (2022). https://doi.org/10.1186/s12859-021-04545-2
