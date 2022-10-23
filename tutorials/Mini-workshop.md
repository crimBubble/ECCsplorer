# Mini-workshop #
-------------------------------------------------------------------------------

**Follow these instructions to learn how to use the ECCsplorer pipeline.**

The example commands shown in the following assume that the detailed instructions for [(2.) installation using conda envs](https://github.com/crimBubble/ECCsplorer/blob/master/tutorials/Installation_instructions.md) were used (*see* Note 0). Note that the instructions are meant to be used step-by-step (some commands require previous ones to be run before, e.g. directory changes).

Open a command line window (terminal/powershell/etc.) and when working on a server (or VM) login with ```ssh```:

```bash
ssh <user>@<server-address>
# <user>@<server-address>'s password:
# <user>@<machine>:~$
```

1. Installation
2. Input data
3. Command options (--help)
4. Running the pipeline
5. Output data


### 1. Installation

Find detailed installation instructions [here](https://github.com/crimBubble/ECCsplorer/blob/master/tutorials/Installation_instructions.md).

Activate the ECCsplorers conda environment (to make sure all dependencies and 3rd party tools are available):

```bash
conda activate eccsplorer
```

### 2. Input data

0. Circle-Seq (circSeq) or mobilome-seq data (required, *see* Note 1)

You find example mobilome-seq testdata files in the ```/.../ECCsplorer/testdata/``` directory called ```aDNA_R<1/2>.fastq```. Find a description of the test data [here](https://github.com/crimBubble/ECCsplorer/blob/master/testdata/README.md). Check if the two files are available with:

```bash
cd ~/miniconda3/envs/eccsplorer/bin/ECCsplorer/testdata
ls
```
Out: ```README.md  aDNA_R1.fastq  aDNA_R2.fastq  gDNA_R1.fastq  gDNA_R2.fastq```

1. Control sequencing data (option 1, *see* Note 2)

You find example control testdata files in the ```/.../ECCsplorer/testdata/``` directory called ```gDNA_R<1/2>.fastq``` similar to the mobilome-seq testdata files.

2. Reference genome assembly (option 2, *see* Note 3)

Download the example reference assembly using [E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK179288/):

```bash
efetch -db nucleotide -id CM009438.1 -seq_start 8971216 -seq_stop 9147030 -format fasta > chr1.fa
efetch -db nucleotide -id CM009440.1 -seq_start 44915437 -seq_stop 45118630 -format fasta > chr2.fa
efetch -db nucleotide -id CM009444.1 -seq_start 23920431 -seq_stop 24120625 -format fasta > chr3.fa

cat chr1.fa chr2.fa chr3.fa | awk '/^>/{print ">chr" ++i; next}{print}' > RefGenomeSeq.fa

rm chr1.fa chr2.fa chr3.fa
```

3. Reference sequence database of potential eccDNA candidates (optional, only used for annotation during clustering when input option 1 available)

Download the example reference database:

```bash
efetch -db nucleotide -id JX455085.1 -format fasta > RefSeq_DB.fa
```

Check if all input files are available in the ```/.../ECCsplorer/testdata/``` directory:

```bash
ls
```
Out: ```README.md  RefGenomeSeq.fa  RefSeq_DB.fa  aDNA_R1.fastq  aDNA_R2.fastq  gDNA_R1.fastq  gDNA_R2.fastq```

Copy the testdata to the ```home``` directory for easier usage.

```bash
cd ~
mkdir ecc_WS
cp -r miniconda3/envs/eccsplorer/bin/ECCsplorer/testdata/ ecc_WS/
```

### 3. Command options

Find a detailed description of all options in the [quick-start manual](https://github.com/crimBubble/ECCsplorer/blob/master/Quickstart_Manual.pdf).

```bash
eccsplorer -h
```
Out:
```
usage: eccsplorer [-h] [-ref <file>] [-out <directory>] [-trm <option>]
                  [-img <option>] [-dsa <txt>] [-dsb <txt>] [-rgs <int>]
                  [-cnt <int>] [-win <int>] [-tax <tax>] [-log]
                  [-d <DB> [<DB> ...]] [-cpu <int>] [-m <option>]
                  <file1A> <file2A> [<file1B>] [<file2B>]

ECCsplorer v0.9b: detecting extrachromosomal circular DNAs (eccDNA) from short read sequencing data.

positional arguments:
  <file1A>              Paired-end reads file1 of data set A (R1). Required.
  <file2A>              Paired-end reads file2 of data set A (R2). Required.
  <file1B>              Paired-end reads file1 of data set B (R1). Recommended.
  <file2B>              Paired-end reads file2 of data set B (R2). Recommended.

optional arguments:
  -h, --help            show this help message and exit
  -ref <file>, --reference_genome <file>
                        Reference genome sequence in FASTA format.
                        With single chromosomes named as chr1, chr2, ...chrN
  -out <directory>, --output_dir <directory>
                        Name your project output directory.
                        Default: eccpipe (old content will be partially overwritten !)
  -trm <option>, --trim_reads <option>
                        Read trimming with trimmomatic v0.38.
                        Strongly recommended, for usage specify adapter option:
                        - nex (Nextera),
                        - tru2 (TruSeq2), tru3 (TruSeq3), tru3-2 (TruSeq3-2)
                        - custom (see trimmomatic manual, name UserAdapter.fa)
  -img <option>, --image_format <option>
                        Choose your desired image format.
                        Options: png (default), jpeg, bmp, tiff, pdf
  -dsa <txt>, --preA <txt>
                        Set readID prefix for data set A. Max. 10 characters. Default = TR
                        Used for comparative analysis.
  -dsb <txt>, --preB <txt>
                        Set readID prefix for data set B. Equal length as -dsa. Default = CO
                        Used for comparative analysis.
  -rgs <int>, --genome_size <int>
                        Set the genome size of your organism in base pairs [bp].
                        Only needed if -cnt set to "auto". Default = None
  -cnt <int>, --read_count <int>
                        Number of reads to use for clustering, if not set max. available reads are used.
                        Set "auto" to use 0.1x genome coverage (only with -ref or -rgs set)
                        Note: for mapping analysis max. available reads are used.
  -win <int>, --window_size <int>
                        Window size for mapping analysis.
                        Used for peak detection and visualization.
                        Smaller window size increases memory usage. Default = 100
  -tax <tax>, --taxon <tax>
                        Use this option to specify taxon using:
                        vir for Viridiplantea (default) or met for Metazoa.
  -log                  Use this option to print logging to file.
                        If not set logging is only printed to stdout.
  -d <DB> [<DB> ...], --database <DB> [<DB> ...]
                        Fasta file for custom BLASTn (annotation database).
                        Existing database might be used.
                        Usage of multiple databases separated by space possible.
  -cpu <int>, --max_threads <int>
                        Specify max. threads to use.
                        Default = max. available cpu threads are used.
  -m <option>, --mode <option>
                        Choose mode to run.
                        Options:
                        all (default, run all modules)
                        map (run only mapping module)
                        clu (run only clustering module)
                        PRExer (only run preparation module)

Thanks for using ECCsplorer.
```

Additional parameters can be changed in the ```/.../ECCsplorer/lib/config.py```. Find a detailed description and recommendations in the [quick-start manual](https://github.com/crimBubble/ECCsplorer/blob/master/Quickstart_Manual.pdf).

### 4. Running the pipeline

Running the ECCsplorer pipeline using the example testdata (described in 2. Input data) and the following setup:
```
Output directory:          testrun                      # in the current working directory
Prefix data set A:         TR                           # default setting
File 1 data set A (f1a):   testdata/aDNA_R1.fastq       # example mobilome-seq testdata (forward reads)
File 2 data set A (f2a):   testdata/aDNA_R2.fastq       # example mobilome-seq testdata (reverse reads)
Prefix data set B:         CO                           # default setting
File 1 data set B (f1b):   testdata/gDNA_R1.fastq       # example control testdata (forward reads)
File 2 data set B (f2b):   testdata/gDNA_R2.fastq       # example control testdata (reverse reads)
Reference genome sequence: testdata/RefGenomeSeq.fa     # example reference assembly
Custom BLAST+ database:    ['testdata/RefSeq_DB.fa']    # example reference database
Taxon:                     vir                          # default setting
Read trimming option:      tru3                         # Trimming using the truSeq3 adapter
Mapping window size:       100                          # default setting
Genome size:               ---                          # default setting
User read count:           1000                         # use 1k reads from each input read data file
Run mode:                  all                          # default setting
Image format:              png                          # default setting
Max threads used:          8                            # default setting (max. available)
Logging to file:           Yes                          # log the command line output
```

To run the pipeline with the setup described above execute the following command (run time < 5 min):

```bash
cd ~/ecc_WS

eccsplorer --output_dir testrun testdata/aDNA_R1.fastq testdata/aDNA_R2.fastq testdata/gDNA_R1.fastq testdata/gDNA_R2.fastq --reference_genome testdata/RefGenomeSeq.fa --database testdata/RefSeq_DB.fa --trim_reads tru3 --read_count 1000 -log

```

Out (successful run, last 5 lines):
```
[1] "Starting comparative scatter plot."
[1] "Finished comparative scatter plot."
2022-10-22 14:21:11,603 - [exit_msg] INFO: Exiting...
2022-10-22 14:21:11,804 - [r_shutdown] INFO: Shutting down Rserve.
2022-10-22 14:21:11,804 - [exit_msg] INFO: Thanks for using ECCsplorer!
```

Find detailed information on the analysis performed within the pipeline in the [ECCsplorer article's method section](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04545-2#Sec2) or in the [quick-start manual](https://github.com/crimBubble/ECCsplorer/blob/master/Quickstart_Manual.pdf).

### 5. Output data

Find a detailed description of all output files in the [quick-start manual](https://github.com/crimBubble/ECCsplorer/blob/master/Quickstart_Manual.pdf).

Some of the output data are HTML and image files. When working on a server the output needs to be transferred to a computer with a graphical user interface (GUI). This can be done using e.g. [```sftp```](https://en.wikipedia.org/wiki/SSH_File_Transfer_Protocol), [```scp```](https://en.wikipedia.org/wiki/Secure_copy_protocol), or [```rsync```](https://en.wikipedia.org/wiki/Rsync). In the example we will use ```sftp```. Additionally, the ECCsplorer output includes many small files. To speed up the transfer process the output directory can be archived before (*see* Note 4).

Archive files:
```bash
zip -r ecc_output.zip testrun/
```

Disconnect from a server (```exit```) and navigate (in the terminal/shell/etc.) to the location where you want to place the output files. 

Login to server using ```sftp```, transfer files and disconnect:
```bash
sftp <user>@<server-address>
# <user>@<server-address>'s password:
# sftp>
cd ecc_WS
ls
# <date-time>.log  ecc_output.zip  formatdb.log  testdata  testrun
get ecc_output.zip
exit
```

Find a detailed description of all output files, image layouts and special file formats in the [quick-start manual](https://github.com/crimBubble/ECCsplorer/blob/master/Quickstart_Manual.pdf). Find additional explanations to the ```clustering_results``` on the [RepeatExplorer2's website](http://repeatexplorer.org/).

**You completed the mini-workshop. Thanks for trying the ECCsplorer pipeline.**


### Notes

0. If you used a different way of installation make sure to adjust the commands to match you file locations.

1. Circle-Seq/mobilome-seq data is paired-end read sequencing data in FASTA/Q format from amplified (phi29 polymerase) circular DNA with the linear DNA fraction removed (exonuclease) before amplification. For more detailed information on the method see [Mann *et al.* (2022)](https://doi.org/10.1186/s12859-021-04545-2) [Supp. data 2](https://static-content.springer.com/esm/art%3A10.1186%2Fs12859-021-04545-2/MediaObjects/12859_2021_4545_MOESM2_ESM.docx) or [Lanciano *et al.* (2021)](https://doi.org/10.1007/978-1-0716-1134-0_7)

2. Control data is is paired-end read sequencing data in FASTA/Q format to be compared with the Circle-Seq/mobilome-seq data. It can be WGS data or from a different treatment (depending on the research question).

3. The reference assembly should be a chromosome-level and high-quality assembly.

4. The ECCsplorer output includes various files which are used as input for the different modules. The actual result files are located in the ```eccpipe_results``` directory.