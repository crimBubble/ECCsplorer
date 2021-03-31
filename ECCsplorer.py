#!/usr/bin/env python3

"""
ECCsplorer: a tool for detection of extrachromosomal circular DNA from NGS data
Copyright (C) 2020 Ludwig Mann

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import argparse
import atexit
import logging
import os
import re
import shlex
import shutil
import subprocess
import sys
import textwrap

from datetime import datetime
from lib import config
from lib import eccDNA_Rcodes, eccPrepare
from lib import eccMapper, eccClusterer, eccComparer

LOGGER = logging.getLogger(__name__)


class DirectoryType(object):
    """
    this class is similar to argparse.FileType
    for mode 'w' creates and check the access to the directory
    for mode 'r' check the presence of the dictory and accesibility
    """

    def __init__(self, mode='r'):
        self._mode = mode

    def __call__(self, string):
        if self._mode == 'w':
            try:
                os.makedirs(string, exist_ok=True)
            except FileExistsError:
                raise argparse.ArgumentTypeError(
                    ("Cannot create directory, '{}' is a file".format(string)))
            if os.access(string, os.W_OK):
                return string
            else:
                raise argparse.ArgumentTypeError(
                    "Directory '{}' is not writable".format(string))
        if self._mode == 'r':
            if not os.path.isdir(string):
                raise argparse.ArgumentTypeError(
                    "'{}' is not a directory".format(string))
            if os.access(string, os.R_OK):
                return string
            else:
                raise argparse.ArgumentTypeError(
                    "Directory '{}' is not readable".format(string))


def exit_err():
    LOGGER.error('Sorry, something went wrong.')


def exit_msg():
    LOGGER.info('Exiting...')

    eccDNA_Rcodes.r_shutdown(port=config.PORT, logger=LOGGER)

    LOGGER.info('Thanks for using ECCsplorer!')

    os._exit(0)


def get_argparser():
    """
    handle user command line input options, get help for usage with -h, --help
    """

    # TODO: check if DB title can be used instead of fasta file in -d option

    parser = argparse.ArgumentParser(description='ECCsplorer {0}: detecting extrachromosomal circular DNAs (eccDNA)'
                                                 ' from short read sequencing data.'.format(config.VERSION),
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog='Thanks for using ECCsplorer.')

    # input files
    parser.add_argument('f1a',
                        help='Paired-end reads file1 of data set A (R1). Required.',
                        metavar="<file1A>",
                        nargs=1,
                        type=argparse.FileType())

    parser.add_argument('f2a',
                        help='Paired-end reads file2 of data set A (R2). Required.',
                        metavar="<file2A>",
                        nargs=1,
                        type=argparse.FileType())

    parser.add_argument('f1b',
                        help='Paired-end reads file1 of data set B (R1). Recommended.',
                        metavar="<file1B>",
                        nargs='?',
                        default=None,
                        type=argparse.FileType())

    parser.add_argument('f2b',
                        help='Paired-end reads file2 of data set B (R2). Recommended.',
                        metavar="<file2B>",
                        nargs='?',
                        default=None,
                        type=argparse.FileType())

    parser.add_argument('-ref', '--reference_genome',
                        help=('Reference genome sequence in FASTA format.\n'
                              'With single chromosomes named as chr1, chr2, ...chrN'),
                        default=None,
                        metavar="<file>",
                        type=argparse.FileType())

    # input options
    parser.add_argument('-out', '--output_dir',
                        help=('Name your project output directory.\n'
                              'Default: {} (old content will be partially overwritten !)'
                              .format(config.PRJ_NAME)),
                        required=False, metavar="<directory>",
                        type=DirectoryType('w'),
                        default='{}'.format(config.PRJ_NAME))

    parser.add_argument('-trm', '--trim_reads',
                        help=('Read trimming with trimmomatic v0.38. \n'
                              'Strongly recommended, for usage specify adapter option: \n'
                              '- nex (Nextera), \n'
                              '- tru2 (TruSeq2), tru3 (TruSeq3), tru3-2 (TruSeq3-2) \n'
                              '- custom (see trimmomatic manual, name UserAdapter.fa)'),
                        required=False,
                        metavar="<option>")

    parser.add_argument('-img', '--image_format',
                        help=('Choose your desired image format. \n'
                              'Options: png (default), jpeg, bmp, tiff, pdf'),
                        required=False,
                        metavar="<option>",
                        default='{}'.format(config.IMAGE_TYPE))

    parser.add_argument('-dsa', '--preA',
                        help=('Set readID prefix for data set A. Max. 10 characters. Default = TR \n'
                              'Used for comparative analysis.'),
                        default='TR',
                        metavar="<txt>")

    parser.add_argument('-dsb', '--preB',
                        help=('Set readID prefix for data set B. Equal length as -dsa. Default = CO \n'
                              'Used for comparative analysis.'),
                        default='CO',
                        metavar="<txt>")

    parser.add_argument('-rgs', '--genome_size',
                        help='Set the genome size of your organism in base pairs [bp].\n'
                             'Only needed if -cnt set to "auto". Default = {}'.format(config.REPEX_GENOME_SIZE),
                        required=False,
                        type=int,
                        metavar="<int>",
                        default=config.REPEX_GENOME_SIZE)

    parser.add_argument('-cnt', '--read_count',
                        help='Number of reads to use for clustering, if not set max. available reads are used.\n'
                             'Set "auto" to use 0.1x genome coverage (only with -ref or -rgs set)\n'
                             'Note: for mapping analysis max. available reads are used.',
                        required=False,
                        metavar="<int>",
                        default=None)

    parser.add_argument('-win', '--window_size',
                        help=('Window size for mapping analysis. \n'
                              'Used for peak detection and visualization.\n'
                              'Smaller window size increases memory usage. Default = {}'.format(config.WINDOW_SIZE)),
                        default=config.WINDOW_SIZE,
                        type=int,
                        metavar="<int>")

    parser.add_argument('-tax', '--taxon',
                        help=('Use this option to specify taxon using:\n'
                              'vir for Viridiplantea (default) or met for Metazoa.'),
                        metavar="<tax>",
                        default='vir')

    parser.add_argument('-log',
                        help=('Use this option to print logging to file.\n'
                              'If not set logging is only printed to stdout.'),
                        action='store_true')

    parser.add_argument('-d', '--database',
                        help=('Fasta file for custom BLASTn (annotation database).\n'
                              'Existing database might be used.\n'
                              'Usage of multiple databases separated by space possible.'),
                        metavar='<DB>',
                        nargs='+',
                        type=str,
                        default=None)

    parser.add_argument('-cpu', '--max_threads',
                        help=('Specify max. threads to use. \n'
                              'Default = max. available cpu threads are used.'),
                        default=config.CPU,
                        type=int,
                        metavar="<int>")

    parser.add_argument('-m', '--mode',
                        help=('Choose mode to run. \n'
                              'Options: \n'
                              'all (default, run all modules)\n'''
                              'map (run only mapping module)\n'
                              'clu (run only clustering module)\n' 
                              'PRExer (only run preparation module)'),
                        required=False,
                        metavar="<option>",
                        default='{}'.format(config.PIPELINE_MODE))

    return parser


def setup_logging(args):
    """
    start logging on logging level info
    and write logging to file if set by user
    :param args: user input
    """

    if args.log:
        log2console = logging.StreamHandler()
        log2console.setLevel(logging.INFO)
        log2console.setFormatter(logging.Formatter('%(asctime)s - [%(funcName)-s] %(levelname)-s: %(message)s'))

        LOGGER.addHandler(log2console)

        logfile = '{:%Y%m%d%H%M}.log'.format(datetime.today())

        logging.basicConfig(filename=logfile,
                            format='%(asctime)s - %(name)s - %(funcName)s - %(levelname)s -\n%(message)s\n',
                            filemode='w',
                            level=logging.INFO)

    else:
        logging.basicConfig(format='%(asctime)s - [%(funcName)-s] %(levelname)-s: %(message)s',
                            level=logging.INFO)

    # starting string of pipeline and logging (all input options)
    starting_str = 'Starting ECCsplorer pipeline with...\n' \
                   'Output directory:          {output_dir_}\n' \
                   'Prefix data set A:         {preA_}\n' \
                   'File 1 data set A (f1a):   {f1a_}\n' \
                   'File 2 data set A (f2a):   {f2a_}\n' \
                   'Prefix data set B:         {preB_}\n' \
                   'File 1 data set B (f1b):   {f1b_}\n' \
                   'File 2 data set B (f2b):   {f2b_}\n' \
                   'Reference genome sequence: {reference_genome_}\n' \
                   'Custom BLAST+ database:    {database_}\n' \
                   'Taxon:                     {taxon_}\n' \
                   'Read trimming option:      {trim_reads_}\n' \
                   'Mapping window size:       {window_size_}\n' \
                   'Genome size:               {genome_size_}\n' \
                   'User read count:           {user_cnt_}\n' \
                   'Run mode:                  {mode_}\n' \
                   'Image format:              {img_}\n' \
                   'Max threads used:          {cpu_}\n' \
                   'Logging to file:           {log_}' \
        .format(output_dir_=args.output_dir,
                preA_=args.preA,
                f1a_=args.f1a[0].name,
                f2a_=args.f2a[0].name,
                preB_=args.preB,
                f1b_=(args.f1b.name if args.f1b is not None else '---'),
                f2b_=(args.f2b.name if args.f2b is not None else '---'),
                reference_genome_=(args.reference_genome.name if args.reference_genome is not None else '---'),
                database_=args.database,
                taxon_=args.taxon,
                trim_reads_=args.trim_reads,
                window_size_=args.window_size,
                genome_size_=args.genome_size if args.genome_size is not None else '---',
                user_cnt_=args.read_count if args.read_count is not None else 'not set, using max. available reads',
                mode_=args.mode,
                img_=args.image_format,
                cpu_=args.max_threads,
                log_='Yes' if args.log else 'No')
    LOGGER.info(starting_str.strip())


def basic_checkups(args):
    """
    check availability of bioinformatic tools, and if necessary their versions.
    check for segemehl/RepeatExplorer depending on run mode (--mode).
    start and check pyRserve connection.
    :return connector
    """
    LOGGER.info('Performing basic checkups.')

    # python
    python_version = (3, 4)
    if sys.version_info < python_version:
        error_msg = '\npython 3.4 or higher is required!\n'
        LOGGER.error(error_msg)
        raise Exception(error_msg)

    # R_version = '3.6.0'
    if not shutil.which('R'):
        error_msg = textwrap.dedent('\nR not found!\n\n'
                                    'Please make sure R is installed and in PATH.\n')
        LOGGER.error(error_msg)
        raise Exception(error_msg)

    # segemehl_version = '0.3.4'
    if args.mode in ['all', 'map']:
        if not shutil.which(config.SEGEMEHL_PATH):
            error_msg = textwrap.dedent('\nsegemehl not found!\n\n'
                                        'Please make sure segemehl is installed.\n'
                                        'If not in PATH modify location: '
                                        'config > SEGEMEHL_PATH.\n')
            LOGGER.error(error_msg)
            raise Exception(error_msg)

        # haarz_version = '0.3.4'
        if not shutil.which(config.HAARZ_PATH):
            error_msg = textwrap.dedent('\nhaarz not found!\n\n'
                                        'Please make sure haarz is installed.\n'
                                        'If not in PATH modify location: '
                                        'config > HAARZ_PATH.\n')
            LOGGER.error(error_msg)
            raise Exception(error_msg)

    # repeatexplorer
    # TODO: activate conda env within ECCsplorer (possible?)

    if args.mode in ['all', 'clu']:

        if not shutil.which(config.REPEATEXPLORER_PATH):
            error_msg = textwrap.dedent('\nRepeatExplorer2 not found!\n\n'
                                        'Please make sure RepeatExplorer2 is installed.\n'
                                        'If you installed RepeatExplorer2 using conda, '
                                        'please activate the environment\n'
                                        'before starting the ECCsplorer pipeline.\n'
                                        'If not in PATH modify location: '
                                        'config > REPEATEXPLORER_PATH.\n')
            LOGGER.error(error_msg)
            raise Exception(error_msg)

    # blast+
    if not shutil.which(config.BLASTn_PATH):
        error_msg = textwrap.dedent('\nBLAST+ not found!\n\n'
                                    'Please make sure BLAST+ is installed.\n'
                                    'If not in PATH modify location: '
                                    'config > BLAST(n/x)_PATH; BLASTMAKEDB.\n')
        LOGGER.error(error_msg)
        raise Exception(error_msg)

    if args.mode in ['all', 'map']:

        # SAMtools
        samtools_version = "1.9"
        if shutil.which(config.SAMTOOLS_PATH):
            is_version = subprocess.check_output((config.SAMTOOLS_PATH, '--version'), universal_newlines=True)
            is_version = is_version.split('\n')[0]
            is_version = is_version.replace('v', '').split(' ')[-1]
            for idx in range(samtools_version.count(".") + 1):
                if int(is_version.split(".")[idx]) < int(samtools_version.split(".")[idx]):
                    error_msg = '\nsamtools {} or higher is required!\n'.format(samtools_version)
                    LOGGER.error(error_msg)
                    raise Exception(error_msg)
        else:
            error_msg = textwrap.dedent('\nsamtools not found!\n\n'
                                        'Please make sure samtools is installed.\n'
                                        'If not in PATH modify location: '
                                        'config > SAMTOOLS_PATH.\n')
            LOGGER.error(error_msg)
            raise Exception(error_msg)

        # BEDtools
        bedtools_version = 2.28
        if shutil.which(config.BEDTOOLS_PATH):
            is_version = subprocess.check_output((config.BEDTOOLS_PATH, '--version'), universal_newlines=True)
            is_version = is_version.split('\n')[0]
            is_version = is_version.replace('v', '').split(' ')[-1]
            is_version = is_version.split('.')
            is_version = float(str(is_version[0] + '.' + is_version[1]))
            if is_version < bedtools_version:
                error_msg = '\nbedtools {} or higher is required!\n'.format(bedtools_version)
                LOGGER.error(error_msg)
                raise Exception(error_msg)
        else:
            error_msg = textwrap.dedent('\nbedtools not found!\n\n'
                                        'Please make sure bedtools is installed.\n'
                                        'If not in PATH modify location: '
                                        'config > BEDTOOLS_PATH.\n')
            LOGGER.error(error_msg)
            raise Exception(error_msg)

    # Trimmomatic
    if args.trim_reads is not None:
        if not shutil.which(config.TRIM_PATH):
            error_msg = textwrap.dedent('\nTrimmomatic not found!\n\n'
                                        'Please make sure Trimmomatic is installed.\n'
                                        'If not in PATH modify location: '
                                        'config > TRIM_PATH.\n')
            LOGGER.error(error_msg)
            raise Exception(error_msg)

    conn = eccDNA_Rcodes.r_connection(LOGGER)

    LOGGER.info('Basic checkups passed successfully.')

    return conn


def input_checkups(args, current_read_files):
    """
    check user input options
    """

    # check sufficient mode (--mode)
    if args.mode is not None:
        if args.mode not in config.PIPELINE_MODES:
            error_msg = textwrap.dedent('Please specify a valid run mode.\n'
                                        'Options: {}'.format(' '.join(config.PIPELINE_MODES)))
            LOGGER.error(error_msg)
            raise Exception(error_msg)

    # minimal data requirements
    if args.mode in ['all', 'map']:

        if args.reference_genome is None:
            error_msg = textwrap.dedent('Please specify a reference genome sequence (-ref/--reference_genome).\n')
            LOGGER.error(error_msg)
            raise Exception(error_msg)

    if args.mode in ['all', 'clu']:

        if args.f1b is None or args.f2b is None:
            error_msg = textwrap.dedent('Please specify a control data set (f1b, f2b).\n')
            LOGGER.error(error_msg)
            raise Exception(error_msg)

    # image output format option (--img_format)
    if args.image_format is not None:
        if args.image_format not in config.IMAGE_EXT:
            error_msg = textwrap.dedent('Please specify a valid image format.\n'
                                        'Options: {}'.format(' '.join(config.IMAGE_EXT.keys())))
            LOGGER.error(error_msg)
            raise Exception(error_msg)

    # reference sequence (--reference_genome)
    if args.reference_genome is not None:
        with open(args.reference_genome.name) as file:
            first_char = []
            for line in file:
                first_char.append(line[0])
            for i in range(0, len(first_char)):
                if first_char[i] == '>':
                    if first_char[i + 1] in config.FASTA_ALPHABET:
                        file_is_fasta = True
                    else:
                        file_is_fasta = False
                        LOGGER.info('Line: {line_} = {letter_}'
                                    .format(line_=i + 1, letter_=first_char[i + 1]))
                        break

        if not file_is_fasta:
            error_msg = 'Reference genome sequence file is not in FASTA format.'
            LOGGER.error(error_msg)
            raise Exception(error_msg)

    # trimming (--trim_reads)
    trm_opts = config.TRIM_ADAPTER.keys()

    if args.trim_reads is not None:
        if args.trim_reads not in trm_opts:
            error_msg = 'Possible trimming options: {}'.format(trm_opts)
            LOGGER.error(error_msg)
            raise Exception(error_msg)

    # read count (--read_count)
    if args.read_count is not None:
        if args.read_count != 'auto':
            try:
                args.read_count = int(args.read_count)
            except ValueError:
                error_msg = 'Read count option must be an integer or "auto".'
                LOGGER.error(error_msg)
                raise Exception(error_msg)

    # prefix (--preA/--preB)
    prefix_msg = ''
    if len(args.preA) > 10:
        prefix_msg = 'Prefix A is too long.'
    if len(args.preB) > 10:
        prefix_msg = 'Prefix B is too long.'
    if len(args.preA) != len(args.preB):
        prefix_msg = 'Prefixes are unequal in length.'
    if args.preA == args.preB:
        prefix_msg = 'Prefixes can not be the same.'

    if prefix_msg != '':
        LOGGER.error(prefix_msg)
        raise Exception(prefix_msg)

    # taxon (--taxon)
    if args.mode in ['all', 'clu']:
        tax_opts = config.REPEX_TAX.keys()

        if args.taxon not in tax_opts:
            error_msg = 'Possible taxon options: {}'.format(tax_opts)
            LOGGER.error(error_msg)
            raise Exception(error_msg)

    # check if read files are already fasta (quick check)
    files_are_fasta = False

    for file in current_read_files:
        if file is not None:
            file_ext = str.lower(os.path.splitext(os.path.basename(file))[-1])
            if file_ext in config.FASTA_EXT:
                with open(file) as fasta_file:
                    head = [next(fasta_file) for x in range(10)]
                    for i in range(10):
                        if head[i][0] == '>':
                            for j in range(len(head[i + 1])):
                                if head[i + 1][j] in config.FASTA_ALPHABET:
                                    files_are_fasta = True

    return files_are_fasta


def basic_setup(args):
    """
    setup directories, reference files and databases
    :return: config.REFERENCE_PATH, config.BLASTDB, config.PRJ_NAME
    """
    # make directories and sub directories
    LOGGER.info('Creating directory structure.')

    # TODO: only set up necessary folders depending on run mode

    config.PRJ_NAME = args.output_dir

    def create_dir(dir_keys):
        for dir_key in dir_keys:
            dir_path = os.path.join(config.CURRENT_DIR, config.PRJ_NAME, config.DIR_TREE[dir_key])
            os.makedirs(dir_path, exist_ok=True)
            LOGGER.debug('Created directory: {}'.format(config.DIR_TREE[dir_key]))

    create_dir(['read_data'])

    if args.mode != 'PRExer':
        create_dir(['reference_data', 'results'])

    if args.mode in ['all', 'map']:
        create_dir(['mapping_results', 'mapping_candidates'])

    if args.mode in ['all', 'clu']:
        create_dir(['clustering_results'])

    if args.mode in ['all']:
        create_dir(['comparative_results', 'comparative_candidates'])

    # setup reference sequence
    ref_out_dir = os.path.join(config.CURRENT_DIR, config.PRJ_NAME, config.DIR_TREE['reference_data'])

    if args.reference_genome is not None:
        try:
            config.REFERENCE_PATH = shutil.copy(args.reference_genome.name, ref_out_dir)
        except shutil.SameFileError:
            LOGGER.info('Reference genome file already inplace.')
            config.REFERENCE_PATH = os.path.join(config.CURRENT_DIR, args.reference_genome.name.strip('./'))

        # check for existing chromosome size file or create it
        refseq_size_file = os.path.splitext(config.REFERENCE_PATH)[0] + '_chrSize.txt'

        if os.path.isfile(refseq_size_file) and os.path.getsize(refseq_size_file) > 0:
            LOGGER.info('Chromosome size file found.')
        else:
            LOGGER.info('Creating chromosome size file.')
            from Bio import SeqIO
            refseq_chr_name = []
            refseq_chr_size = []
            with open(config.REFERENCE_PATH) as refseq:
                for record in SeqIO.parse(refseq, "fasta"):
                    refseq_chr_name.append(record.id)
                    refseq_chr_size.append(len(record))
            with open(refseq_size_file, 'wt') as size_file:
                for i in range(len(refseq_chr_name)):
                    line = '{0}\t{1}\n'.format(refseq_chr_name[i], refseq_chr_size[i])
                    size_file.write(line)

    # setup blast database
    default_db = None
    # search for RepeatExplorer blastDB (also for map mode)
    try:
        default_db = shutil.which(config.REPEATEXPLORER_PATH)
        default_db = os.path.realpath(default_db)
        default_db = default_db.replace('seqclust', 'databases/dna_database_masked.fasta')
        LOGGER.info('RepeatExplorer database location: {}'.format(default_db))
    except TypeError:
        LOGGER.info('RepeatExplorer database not found')
    blast_db_list = config.BLASTDB
    if default_db is not None:
        blast_db_list.append(default_db)

    if args.database is not None:
        for db in args.database:
            blast_db_list.append(db)

    blast_db_list_updated = []

    for db in blast_db_list:
        db_does_not_exist = True
        if os.path.isfile(db) and os.path.getsize(db) > 0:
            if os.path.isfile(db + '.nhr') and os.path.getsize(db + '.nhr') > 0:
                if os.path.isfile(db + '.nin') and os.path.getsize(db + '.nin') > 0:
                    if os.path.isfile(db + '.nsq') and os.path.getsize(db + '.nsq') > 0:
                        LOGGER.info('BLAST+ database found {}.'.format(db))
                        db_does_not_exist = False

        if db_does_not_exist:
            LOGGER.info('Creating BLAST+ database: {}.'.format(db))

            try:
                db = shutil.copy(db, ref_out_dir)
            except shutil.SameFileError:
                LOGGER.info('Database fasta file already inplace.')
                db = os.path.join(config.CURRENT_DIR, db.strip('./'))

            cmd = '{path_} -in {blastdb_} -dbtype {type_}'.format(path_=config.BLASTMAKEDB_PATH,
                                                                  blastdb_=db,
                                                                  type_='nucl' if config.BLASTDB_is_nucl else 'prot')

            makeblastdb = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, universal_newlines=True)
            out = str.format(makeblastdb.communicate()[0]).lstrip('\n').rstrip('\n')
            LOGGER.info('{}'.format(out))

        blast_db_list_updated.append(db)

    cmd = '{path_} -dblist "{dblist_}" -dbtype {type_} ' \
          '-out {out_}/combinedDB.fas -title combinedDB'.format(path_=config.BLAST_ALIAS_PATH,
                                                                dblist_=' '.join(blast_db_list_updated),
                                                                type_='nucl' if config.BLASTDB_is_nucl else 'prot',
                                                                out_=ref_out_dir)

    blast_alias = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = str(blast_alias.communicate()[1], encoding='ascii').lstrip('\n').rstrip('\n')
    LOGGER.info('{}'.format(out))

    config.BLASTDB = '{out_}/combinedDB.fas'.format(out_=ref_out_dir)
    LOGGER.info('Using combined BLAST+ databases (combinedDB.fas) containing: \n'
                '{}'.format('\n'.join(blast_db_list_updated)))

    # image output format
    config.IMAGE_TYPE = args.image_format

    # max cpus
    config.CPU = args.max_threads

    return config.REFERENCE_PATH, config.BLASTDB, config.PRJ_NAME, config.IMAGE_TYPE, config.CPU


def main():
    """
    Run and coordinate pipeline
    """

    # get user input handle
    parser = get_argparser()
    args = parser.parse_args()

    setup_logging(args=args)

    atexit.register(exit_err)

    # basic checkup (are all tools installed)
    if args.mode != 'PRExer':
        config.CONNECTOR = basic_checkups(args=args)

    # startup parameters
    four_read_files = False
    analysis_errors = False

    # file paths
    repex_prepared_file = ''
    sum_mapper_candidate_fas = ''
    sum_mapper_win_coverage = ''

    current_read_files = [args.f1a[0].name,
                          args.f2a[0].name,
                          (args.f1b.name if args.f1b is not None else None),
                          (args.f2b.name if args.f2b is not None else None)]

    read_file_dict = {args.preA: [current_read_files[0], current_read_files[1]]}

    if args.f1b is not None and args.f2b is not None:
        read_file_dict[args.preB] = [current_read_files[2], current_read_files[3]]
        four_read_files = True

    # basic checkup (is user input appropriate)
    files_are_fasta = input_checkups(args=args, current_read_files=current_read_files)

    # basic setup
    config.REFERENCE_PATH, config.BLASTDB, config.PRJ_NAME, config.IMAGE_TYPE, config.CPU = basic_setup(args=args)

    LOGGER.info('Starting pipeline modules.')

    current_key_ext = ''  # key extension changes when read files are processed

    # prepare files for mapping

    LOGGER.info('Starting trimming. This might take a while.')

    obj_prepare = eccPrepare.eccTrimConvert(read_file_dict, current_read_files, args.trim_reads, LOGGER)

    if args.trim_reads is not None and not files_are_fasta:
        read_file_dict, current_read_files = obj_prepare.read_trimming()
        current_key_ext = '-trim'
    elif args.trim_reads is None:
        LOGGER.info('Trimming skipped.')
    elif files_are_fasta:
        LOGGER.info('Trimming skipped: Files are in FASTA format.')

    if not files_are_fasta:
        read_file_dict, current_read_files = obj_prepare.converter()
        current_key_ext = '-conv'
    else:
        LOGGER.info('Not converting: Files are in FASTA format.')

    if four_read_files:
        sorted_fasta_files = read_file_dict[args.preA + current_key_ext] + read_file_dict[args.preB + current_key_ext]
    else:
        sorted_fasta_files = read_file_dict[args.preA + current_key_ext]

    # prepare files for clustering
    if four_read_files or args.mode == 'PRExer':
        obj_prexer = eccPrepare.eccSampleConcat(sorted_fasta_files=sorted_fasta_files,
                                                pre_a=args.preA, pre_b=args.preB,
                                                log=LOGGER,
                                                genome_size=args.genome_size,
                                                read_count=args.read_count,
                                                four_read_files=four_read_files)
        repex_prepared_file = obj_prexer.prepare_for_clustering()

    if args.mode == 'PRExer':
        LOGGER.info('Read file preparation for manual RepeatExplorer run or mapping finished successfully.')
        exit_msg()

    LOGGER.info('Creating Rscript for visualization purpose. Might be edited to fit individual needs')
    with open(os.path.join(config.PRJ_NAME, config.DIR_TREE['results'], 'viz.R'), 'w+') as viz_file:
        viz_file.write(('#!/usr/bin/env Rscript\n\n'
                        '#Summarized Rscripts for visualization\n{libs_}\n'
                        'options(warn = -1)').format(libs_=eccDNA_Rcodes.load_r_libs))

    if args.reference_genome is not None and args.mode in ['all', 'map']:
        obj_mapper = eccMapper.eccMapper(seq_files=sorted_fasta_files,
                                         ref_seq_file=config.REFERENCE_PATH,
                                         prefix=[args.preA, args.preB],
                                         win_size=args.window_size,
                                         log=LOGGER,
                                         conn=config.CONNECTOR)
        sum_mapper_win_coverage, sum_mapper_candidate_fas, analysis_errors = obj_mapper.mapper_coordinator()

    # exit(1)  # test mapping module

    if four_read_files and args.mode in ['all', 'clu']:
        obj_clusterer = eccClusterer.eccClusterer(seq_file=repex_prepared_file,
                                                  prefix=[args.preA, args.preB],
                                                  taxon=args.taxon,
                                                  conn=config.CONNECTOR,
                                                  log=LOGGER)
        obj_clusterer.cluster_coordinator()

    if os.path.isfile(sum_mapper_candidate_fas) and os.path.getsize(sum_mapper_candidate_fas) != 0:

        if four_read_files:
            obj_comp = eccComparer.eccComparer(map_sum_cov_file=sum_mapper_win_coverage,
                                               mapper_fasta_files=sum_mapper_candidate_fas,
                                               prefix=[args.preA, args.preB],
                                               ref_seq_file=config.REFERENCE_PATH,
                                               window_size=args.window_size,
                                               conn=config.CONNECTOR,
                                               log=LOGGER)
            obj_comp.compare_coordinator()

    viz = True

    if analysis_errors:
        viz = False

    if viz:
        LOGGER.info('Visualizing results. This may take a while!')
        os.system('Rscript {rscript_} --no-save'.format(
            rscript_=os.path.join(config.PRJ_NAME, config.DIR_TREE['results'], 'viz.R')))
    else:
        LOGGER.info('No visualization of results!')

    exit_msg()


if __name__ == '__main__':
    main()
