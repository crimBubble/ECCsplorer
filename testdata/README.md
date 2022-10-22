# ECCsplorer test-data and test-run

## Test-data

The test-data is semi-artificial data. It contains reference sequences from real world data and simulated reads. 

Three random parts of the sugar beet (Beta vulgaris) reference genome [EL10](https://www.ncbi.nlm.nih.gov/assembly/GCA_002917755.1) ([Funk *et al.*, 2018](https://doi.org/10.1111/tpj.13977)) serve as reference sequence. 
The sugar beet LTR-retrotransposon [*Beetle7*](https://www.ncbi.nlm.nih.gov/nuccore/408362947) ([Weber *et al.*, 2013](https://dx.doi.org/10.1186%2F1759-8753-4-8)) is used as simulated circular DNA as well as annotation reference. 
The test-data reads have been simulated from both references (genome and annotation reference¹) using [dwgsim](https://github.com/nh13/DWGSIM) with standard parameters.

Download test-data references using [E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK179288/) from NCBI:

```{bash}
cd /.../ECCsplorer/testdata

efetch -db nucleotide -id CM009438.1 -seq_start 8971216 -seq_stop 9147030 -format fasta > chr1.fa
efetch -db nucleotide -id CM009440.1 -seq_start 44915437 -seq_stop 45118630 -format fasta > chr2.fa
efetch -db nucleotide -id CM009444.1 -seq_start 23920431 -seq_stop 24120625 -format fasta > chr3.fa

cat chr1.fa chr2.fa chr3.fa | awk '/^>/{print ">chr" ++i; next}{print}' > RefGenomeSeq.fa

rm chr1.fa chr2.fa chr3.fa

efetch -db nucleotide -id JX455085.1 -format fasta > RefSeq_DB.fa
```

or download the files manually from NCBI and rename and concatenate the chromosomes.

  chr1... CM009438.1:8971216-9147030
  chr2... CM009440.1:44915437-45118630
  chr3... CM009444.1:23920431-24120625
  
  annotation_reference... JX455085.1
  
## Test-run

For testing the ECCsplorer pipeline installation run:

```{bash}
python3 ECCsplorer.py testdata/aDNA_R1.fastq testdata/aDNA_R2.fastq testdata/gDNA_R1.fastq testdata/gDNA_R2.fastq --reference_genome testdata/RefGenomeSeq.fa --database testdata/RefSeq_DB.fa --output_dir testdata --trim_reads tru3 --read_count 1000
```

Note 1: As the annotation reference is a LTR-retrotransposon the reads have been simulated considering LTR-LTR-junctions and solo-LTR-junctions for the circularization.
