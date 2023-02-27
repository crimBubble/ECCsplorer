#!/usr/bin/python3

"""
Configuration file for the ECCsplorer pipeline:

    - Locations (PATH) of 3rd party bioinformatic tools  (modify if needed)
    - Pipeline parameters                                (modify if needed)
    - 3rd party tools parameters                         (can be modified)
    - Pipeline settings                                  (might be modified, not recommended)
    - 3rd party tools settings                           (do not modify)
    - Commands                                           (do not modify)
    - Pipeline values, dictionaries, and lists           (do not modify)
    - Miscellaneous                                      (do not modify)
"""

import os

#######################################
# Locations of 3rd party tools (PATH) #
#######################################

# TODO: PATH values can be modified if necessary

# Trimmomatic
TRIM_PATH = 'trimmomatic'
TRIM_ADAPTER_PATH = '~/miniconda3/envs/eccsplorer/share/trimmomatic/adapters/'

# seqtk
SEQTK_PATH = 'seqtk'

# BLAST+
BLASTn_PATH = 'blastn'
BLASTx_PATH = 'blastx'
BLASTMAKEDB_PATH = 'makeblastdb'
BLAST_ALIAS_PATH = 'blastdb_aliastool'

# segemehl (and haarz)
SEGEMEHL_PATH = 'segemehl.x'
HAARZ_PATH = 'haarz.x'

# SAMtools
SAMTOOLS_PATH = 'samtools'
git config --global --add safe.directory /mnt/Storage_3TB/Ludwig_Storage/Dauerdaten/Projects/ECCsplorer/working
# BEDtools
BEDTOOLS_PATH = 'bedtools'

# RepeatExplorer2
REPEATEXPLORER_PATH = 'seqclust'

################################
# Pipeline parameters #
################################

# TODO: Additional settings can be changed

# Image output (plotted with R ggplot2) #

# image resolution in dpi
IMAGE_RES = 300
# image dimensions in px
IMAGE_WIDTH = 1625
IMAGE_HEIGHT = 925
# Font size
IMAGE_POINTS = 10

# Rserve #

# Rserve port
PORT = 6311

#############################
# 3rd party tool parameters #
#############################

# TODO: Parameters might be changed
#  be sure to read manual of each 3rd party tool before changing parameters

# Trimmomatic
TRIM_ILLUMINACLIP = ':2:30:10'
TRIM_SET = 'LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:35'

# BLASTn
BLASTn_HSP = '-max_hsps 25'

# segemehl (and haarz)
SEGEMEHL_ACCURACY = '--accuracy 95'
HAARZ_QUAL = '--minsplit 5 --minqual 1'

#####################
# Pipeline settings #
#####################

# TODO: Settings might be changed
#  Note that changing settings might drastically impact results and comparability

# Mapping module #

# Potential eccDNA candidate parameters
MAX_eccDNA_LENGTH = 35000
MIN_eccDNA_LENGTH = 100
# Merge potential eccDNA candidate regions within x bp
MERGE_CLOSE_REGIONS = 1000
# Default window size for rough coverage calculation
WINDOW_SIZE = 100
# DR: Signal to noise ratio (min. coverage = % of max. coverage detected)
BACKGROUND_PERC = 0.05
# Peak finder parameters
PEAK_THRESHOLD = None
PEAK_DISTANCE = 20
PEAK_PROMINENCE = 1
# Enrichment threshold
ENRICH_THRESHOLD = 2.0
# Maximal amount of candidates to be analyzed in detail (decreasing enrichment score)
MAX_CAND_CNT = 500

# Clustering module #

# Minimal probe percentage for cluster to be candidate
REPEX_ECC_PROPORTION = 0.8

#########################################
# 3rd party tools settings #
#########################################

# TODO: Changing 3rd party tools settings can corrupt pipeline

# Trimmomatic #

# Trimmomatic adapter dictionary
TRIM_ADAPTER = {"nex": "NexteraPE-PE.fa",
                "tru2": "TruSeq2-PE.fa",
                "tru3": "TruSeq3-PE.fa",
                "tru3-2": "TruSeq3-PE-2.fa",
                "custom": "UserAdapter.fa"}

# Blast+ #

# BLAST database
BLASTDB = []  # list with path to fasta files
BLASTDB_is_nucl = True
# BLAST table format
BLASTn_OUTFMT = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend',
                 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen']
BLASTn_ADD_PARAMETERS = '-outfmt \"6 ' + ' '.join(BLASTn_OUTFMT) + '\"' + ' ' + BLASTn_HSP

# SAMtools #

# SAMtool flags: reads not mapped in proper paired
SAMFLAGS_DR_not = [2]
# SAMtool flags: reads mapped in unusual orientation (reverse-forward)
SAMFLAGS_DR = [83, 163]

# RepeatExplorer2 #

# Conda environment name
CONDA_ENV = 'repeatexplorer'
# Taxon dictionary
REPEX_TAX = {"vir": "VIRIDIPLANTAE3.0",
             "met": "METAZOA3.0"}
# Approximate genome size in bp
REPEX_GENOME_SIZE = None  # standard_value = None
# Genome coverage to use for RepeatExplorer2 runs
REPEX_FOLD_COV = 0.1  # standard_value = 0.1


############
# Commands #
############

# TODO: Changing commands can corrupt pipeline

# R (with Rserve)
RSERVE_CMD = 'R CMD Rserve --RS-port {port_} -q --no-save'

# Trimmomatic
TRIM_CMD = '{path_} PE -threads {cpu_} {R1_} {R2_} ' \
           '{out_}/{pre_}-R1_trim.fq {out_}/{pre_}-R1_utrim.fq ' \
           '{out_}/{pre_}-R2_trim.fq {out_}/{pre_}-R2_utrim.fq ' \
           'ILLUMINACLIP:' + TRIM_ADAPTER_PATH + '/{adapter_}' + TRIM_ILLUMINACLIP + ' ' \
           + TRIM_SET

# seqtk
SEQTK_CMD = '{path_} seq -A {fq_}'

# BLAST+
BLAST_CMD = '{path_} -query {fasta_} -out {m6_} -db {database_} {add_}'

# segemehl (and haarz)
SEGEMEHL_CMD = '{path_} --index {refidx_} --database {ref_} ' \
               '--query {r1_} --mate {r2_} ' \
               '--threads {cpu_} ' \
               '--outfile {out_} --nomatchfilename {unout_} ' \
               '--splits --briefcigar --MEOP ' + SEGEMEHL_ACCURACY
HAARZ_CMD = '{path_} split --files {files_} ' + HAARZ_QUAL + ' > {out_}'

# BEDtools
BEDTOOLS_MEAN = '{bedtools_} coverage -mean -a {win_} -b {map_}'

# RepeatExplorer2
REPEX_CMD = '{path_} ' \
            '--paired ' \
            '--prefix_length {preflen_} ' \
            '--output_dir {out_} ' \
            '--taxon {tax_} ' \
            '--cpu {cpu_} ' \
            '{in_} ' \
            '--cleanup --keep_names --options ILLUMINA'

############################################
# Pipeline values, dictionaries, and lists #
############################################

# TODO: Changing pipeline values, dictionaries, and lists can corrupt pipeline

# Version
VERSION = 'v0.9b'
PRJ_NAME = 'eccpipe'  # absolute path for R os.path.abspath(os.path.dirname(PRJ_NAME))

# Pipeline running modes
PIPELINE_MODE = 'all'
PIPELINE_MODES = ['all', 'map', 'clu', 'PRExer']

# Multithreading
CPU = os.cpu_count()

# Seed for pseudo random selections
PIPELINE_SEED = 12

# Image types
IMAGE_TYPE = 'png'  # default (more options)
IMAGE_EXT = {'png': '.png',
             'jpeg': '.jpeg',
             'bmp': '.bmp',
             'tiff': '.tiff',
             'pdf': '.pdf'}

# Directories
PIPE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir))
CURRENT_DIR = os.getcwd()
DIR_TREE = {
    'reference_data': 'reference_data',
    'read_data': 'read_data',
    'results': 'eccpipe_results',
    'mapping_results': 'eccpipe_results/mapping_results',
    'mapping_candidates': 'eccpipe_results/mapping_results/candidates',
    'clustering_results': 'eccpipe_results/clustering_results',
    'comparative_results': 'eccpipe_results/comparative_results',
    'comparative_candidates': 'eccpipe_results/comparative_results/candidates'
}

# Reference and FASTA fileformat
REFERENCE_PATH = None
FASTA_EXT = ['.f', '.fa', '.fas', '.fasta', '.fna']
FASTA_ALPHABET = ['A', 'T', 'G', 'C', 'N', 'a', 't', 'g', 'c', 'n']

#################
# Miscellaneous #
#################

# TODO: Changing miscellaneous can corrupt pipeline

# Rserve
CONNECTOR = None
