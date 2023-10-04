# Installation instructions #
-------------------------------------------------------------------------------

Here you find detailed installation instructions for:

1. Installation using conda envs (quick-guide)
2. Installation using conda envs (detailed step-by-setp instructions)
3. Manual installation (old)


### 1. Installation using conda envs (quick-guide)

```bash
git clone https://github.com/crimBubble/ECCsplorer

cd ECCsplorer
```

Installation of dependencies can be done using a [conda environment](https://docs.conda.io/projects/conda/en/latest/). The required environment can be prepared using the following command (includes dependencies for RepeatExplorer2, solving the environment may take a while):

```bash
conda env create -f environment.yml
```

If the creation of the environment fails due to dependency errors, try:

```bash
CONDA_CHANNEL_PRIORITY=flexible conda env create -f environment.yml
```

Activate the environment:

```bash
conda activate eccsplorer
```

Install RepeatExplorer2 following the [manual](https://bitbucket.org/petrnovak/repex_tarean/src/devel/) without(!) creating a new conda environment.

Edit the "Locations of 3rd party tools (PATH)" for Trimmomatic and RepeatExploer2 in *ECCsplorer/lib/config.py* to match your installation.


### 2. Installation using conda envs (detailed step-by-setp instructions)

The following instructions are a step-by-step tutorial with no previous knowledge required. The installation is based on the use of conda environments. If you have already installed conda on your machine continue from step 2.

These installation instructions are tested for Linux Ubuntu 20.04 LTS but should work on any similar system with minimal hardware requirements = CPU: 6 cores or higher; RAM: 16 GB; ROM: 50 GB. Depending on your data the minimal requirements to successfully run the ECCsplorer pipeline might be substantially higher. 

1. Install the environment manager conda. (Note: Here we use miniconda3, but anaconda3 will work as well. You only will have to adjust the commands accordingly.)
The first command line downloads the latest version of miniconda3. The next line makes the installation file runnable. The third line runs the installation file. You will have to confirm several installation steps by either pressing enter or typing ‘yes’ and pressing enter. After the installation finished you apply the changes to the terminal using the last command.
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh -u
source ~/.bashrc
```
After the installation your command line should look similar to the following (Note the base-statement in the beginning): ```(base) user@machine:~/Downloads$```

2. Install Git a version control software. It is here used to download the most recent version of the ECCsplorer pipeline.
```bash
conda install -c conda-forge git -y
```
3. Install Mamba an environment solving tool. Solving environments is also possible with conda itself but mamba is a lot faster.
```bash
conda install -c conda-forge mamba -y
```
4. Download the ECCsplorer pipeline from GitHub and place the files in the appropriate directory (Note: adjust the commands if you use anaconda3 instead of miniconda3. This includes to adjust the ```TRIM_ADAPTER_PATH``` variable in the ```lib/config.py``` file). The ECCsplorer files are downloaded with the git command and the environment is created with mamba. The last steps place the ECCsplorer.py file in the ```$PATH``` and make it runnable.
```bash
cd ~/Documents
git clone https://github.com/crimBubble/ECCsplorer
cd ECCsplorer
mamba env create -f environment.yml
```
If the creation of the environment fails due to dependency errors, try:

```bash
CONDA_CHANNEL_PRIORITY=flexible mamba env create -f environment.yml
```

```bash
cd ~/miniconda3/envs/eccsplorer/bin
mv ~/Documents/ECCsplorer/ ECCsplorer
ln -s ECCsplorer/ECCsplorer.py eccsplorer
chmod +x eccsplorer
```
5. Download and install the RepeatExplorer2 pipeline (Note: You need administrator “sudo” rights to install the necessary dependencies). You again will have to confirm several installation steps by either pressing enter or typing ‘y’ and pressing enter.
```bash
cd ~/miniconda3/envs/eccsplorer/bin
git clone https://bitbucket.org/petrnovak/repex_tarean.git
cd repex_tarean
conda activate eccsplorer
make
sudo dpkg --add-architecture i386
sudo apt-get update
sudo apt-get install libc6:i386 libncurses5:i386 libstdc++6:i386
cd .. 
ln -s repex_tarean/seqclust seqclust
chmod +x seqclust
conda deactivate
```
6. Download the test dataset to confirm your ECCsplorer installation.
```bash
conda activate eccsplorer
cd ~/miniconda3/envs/eccsplorer/bin/ECCsplorer/testdata
efetch -db nucleotide -id CM009438.1 -seq_start 8971216 -seq_stop 9147030 -format fasta > chr1.fa
efetch -db nucleotide -id CM009440.1 -seq_start 44915437 -seq_stop 45118630 -format fasta > chr2.fa
efetch -db nucleotide -id CM009444.1 -seq_start 23920431 -seq_stop 24120625 -format fasta > chr3.fa
cat chr1.fa chr2.fa chr3.fa | awk '/^>/{print ">chr" ++i; next}{print}' > RefGenomeSeq.fa
rm chr1.fa chr2.fa chr3.fa
efetch -db nucleotide -id JX455085.1 -format fasta > RefSeq_DB.fa
conda deactivate
```
7. Confirm the successful installation by running the ECCsplorer pipeline with the test dataset.
```bash
conda activate eccsplorer
cd ~/miniconda3/envs/eccsplorer/bin/ECCsplorer/
eccsplorer testdata/aDNA_R1.fastq testdata/aDNA_R2.fastq testdata/gDNA_R1.fastq testdata/gDNA_R2.fastq --reference_genome testdata/RefGenomeSeq.fa --database testdata/RefSeq_DB.fa --output_dir testrun --trim_reads tru3 --read_count 1000
conda deactivate
```

If the installation was successful, the pipline will finish with the statement: ```INFO: Thanks for using ECCsplorer!```. You can find the exemplary output files in the ```test_run``` directory.

8. (Optional) Install seqtk a fast and lightweight tool for processing sequences. This speeds up the running time for some steps in the pipeline.
```bash
conda install -c bioconda seqtk -y
```


### 3. Manual installation

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
