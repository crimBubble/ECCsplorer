#!usr/bin/env python3

import math
import multiprocessing
import numpy as np
import os
import re
import shlex
import subprocess
import time
import statistics

from Bio import SeqIO
from functools import partial
from scipy import optimize

from lib import config


def converting(in_file, out):
    """
    Static, for multiprocessing as LOGGER raises errors w/ multiprocessing
    :param in_file: read file to convert from FASTA to FASTQ
    :param out: directory of read_data
    :return file: converted file (path)
    """

    if in_file is not None:

        file_name = os.path.splitext(os.path.basename(in_file))[0]

        fasta_out = '{0}/{1}.fa'.format(out, file_name)

        if os.path.isfile(fasta_out) and os.path.getsize(fasta_out) > 0:

            print('Converted {0}.FASTA file found. Skipping converting. (ProcessID {1})'.format(file_name, os.getpid()))

        else:

            try:
                print('Trying to convert FASTQ to FASTA using "seqtk": {0} (ProcessID {1})'.format(file_name, os.getpid()))

                cmd = config.SEQTK_CMD.format(path_=config.SEQTK_PATH,
                                              fq_=in_file)

                with open(fasta_out, 'w+') as out_file:

                    subprocess.call(shlex.split(cmd), stdout=out_file)

                print('Converted {0}.FASTQ to {0}.FASTA using "seqtk" (ProcessID {1})'.format(file_name, os.getpid()))

            except FileNotFoundError:

                print('Install "seqtk" to improve converting speed!')

                SeqIO.convert(in_file, 'fastq', fasta_out, 'fasta-2line')

                print('Converted {0}.FASTQ to {0}.FASTA (ProcessID {1})'.format(file_name, os.getpid()))

        return fasta_out


def get_best_read_length(in_file):
    """
    Static, for multiprocessing as LOGGER raises errors w/ multiprocessing
    :param in_file: FASTA files (text file)
    :return min_loss: number of reads and optimal read length
    """

    if in_file is not None:
        file_name = os.path.splitext(os.path.basename(in_file))[0]

        print('Checking reads in {0}.FASTA (ProcessID {1})'.format(file_name, os.getpid()))

        # TODO: user input option for read length to use (80 - 250 bp, warning if not in this range: 'continue? y/n)

        with open(in_file, 'r') as file:

            read_length_list = []
            # longest_read = 0
            # shortest_read = 1000

            for record in SeqIO.parse(file, "fasta"):
                read_length_list.append(len(record.seq))
                # if len(record.seq) > longest_read:
                #    longest_read = len(record.seq)
                # if len(record.seq) < shortest_read:
                #    shortest_read = len(record.seq)

            def cumulative_bases(set_read_length):
                """Function to calculate data loss @ specific read length"""
                lost_bases = 0

                set_read_length = math.floor(set_read_length)

                print('{0}.FASTA current checked read lenght: {1}[bp]'.format(file_name, set_read_length))

                for read_len_val in read_length_list:
                    if read_len_val < set_read_length:
                        lost_bases = lost_bases + read_len_val
                    elif read_len_val > set_read_length:
                        lost_bases = lost_bases + (read_len_val - set_read_length)
                return lost_bases

            min_loss = math.floor(optimize.fmin(cumulative_bases,
                                                x0=np.mean(np.array(read_length_list)),
                                                maxiter=100, ftol=10000, xtol=10000))

            max_bases = sum(read_length_list)
            perc_loss = round(cumulative_bases(min_loss) / max_bases * 100, 2)

            # left_read_count = 0
            # for read_len_val in read_length_list:
            #    if read_len_val >= min_loss:
            #        left_read_count += 1

            print('Optimal read length for {}.FASTA: '.format(file_name), min_loss,
                  '(Cumulative bases: {}bp; Loss: {}%)'.format(max_bases, perc_loss))

            return min_loss


def get_max_read_count(in_file, best_read_length):
    """
    Static, for multiprocessing as LOGGER raises errors w/ multiprocessing
    :param in_file: paired-end FASTA files (text file)
    :param best_read_length: calculated optimal read length for minimal data loss
    :return left_read_pair_count: number of reads and optimal read length
    """

    if in_file is not None:
        file_name = in_file[2]

        print('Counting paired-end reads >{2}bp in {0}.FASTAs (ProcessID {1})'
              .format(file_name, os.getpid(), best_read_length))

        left_read_pair_count = 0

        with open(in_file[0], 'r') as f1, open(in_file[1], 'r') as f2:
            for record1, record2 in zip(SeqIO.parse(f1, "fasta"), SeqIO.parse(f2, "fasta")):
                if len(record1.seq) >= best_read_length and len(record2.seq) >= best_read_length:
                    left_read_pair_count += 1

            return [file_name, left_read_pair_count]


def prexing_reads(read_file_set, best_read_length, read_count_set, read_count_counted, out):
    """
    Sub sample reads, cutoff read sequence to get unified read length, add prefix and suffix to read id
    :param out: temporary read file (to be concatenated later)
    :param read_count_counted: real read count
    :param read_file_set: paired-end read data
    :param best_read_length: optimal read length with minimal data loss (get_best_read_length)
    :param read_count_set: maximal available reads @ optimal length for paired-end data (get_max_read_count)
    :return interlaced_reads: interlaced reads to be concatenated
    """

    np.random.seed(config.PIPELINE_SEED)

    interlaced_reads = []
    temp_fasta = os.path.join(out, 'REPEATEXPLORER_' + read_file_set[2] + '.fa.tmp')

    reads_in_data_set = read_count_counted[read_file_set[2]]  # ~ left_read_pair_count

    random_indexes = np.random.choice(range(reads_in_data_set), read_count_set, replace=False)
    random_indexes = set(random_indexes)

    i = 0
    percentag_points = np.around(np.linspace(0, read_count_set, 34), 0)
    count = 0

    print('{}: Adding prefix and suffix to read identifiers.'.format(read_file_set[2]))

    with open(temp_fasta, 'w+') as fasta_file:
        print('{}: Creating temporary read file.'.format(read_file_set[2]))
        fasta_file.write('')

    with open(read_file_set[0]) as f1, open(read_file_set[1]) as f2:
        for record1, record2 in zip(SeqIO.parse(f1, "fasta"), SeqIO.parse(f2, "fasta")):
            if len(record1.seq) >= best_read_length and len(record2.seq) >= best_read_length:
                if i in random_indexes:
                    record1.__init__(seq=record1.seq[0:best_read_length],
                                     id=read_file_set[2] + '_' + record1.id + '_#0/1',
                                     description='')
                    record2.__init__(seq=record2.seq[0:best_read_length],
                                     id=read_file_set[2] + '_' + record2.id + '_#0/2',
                                     description='')

                    interlaced_reads.append(record1)
                    interlaced_reads.append(record2)

                    if count in percentag_points:
                        print('{pre_}: {perc_}%'.format(pre_=read_file_set[2], perc_=int(count / read_count_set * 100)))
                        with open(temp_fasta, 'a') as fasta_file:
                            SeqIO.write(interlaced_reads, fasta_file, 'fasta-2line')
                        interlaced_reads.clear()

                    count += 1
                    random_indexes.remove(i)
                i += 1

    print('{pre_}: 100%'.format(pre_=read_file_set[2]))
    with open(temp_fasta, 'a') as fasta_file:
        SeqIO.write(interlaced_reads, fasta_file, 'fasta-2line')
    interlaced_reads.clear()

    return temp_fasta


class eccTrimConvert:
    """
    Trim and convert reads.
    From user input (FASTQ) to segemehl input (FASTA).
    """

    def read_trimming(self):

        """
        Read trimming using trimmomatic.
        :return read_file_dict, current_read_files:
        """

        for key in self.read_file_dict:
            trimmed_file = '{out_}/{pre_}-R1_trim.fq'.format(out_=self.out, pre_=key)

            if os.path.isfile(trimmed_file) and os.path.getsize(trimmed_file) > 0:
                self.LOGGER.info('Trimmed file found. Skipping trimming.')

            else:
                if self.read_file_dict[key][0] is not None and self.read_file_dict[key][1] is not None:

                    self.LOGGER.info('Trimming using: {f1_} and {f2_} as read files.'
                                     .format(f1_=self.read_file_dict[key][0],
                                             f2_=self.read_file_dict[key][1]))

                    cmd = config.TRIM_CMD.format(path_=config.TRIM_PATH,
                                                 cpu_=config.CPU,
                                                 pre_=key,
                                                 R1_=self.read_file_dict[key][0],
                                                 R2_=self.read_file_dict[key][1],
                                                 out_=self.out,
                                                 pipe_dir_=config.PIPE_DIR,
                                                 adapter_=config.TRIM_ADAPTER[self.user_adapter])

                    trim = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    out = str(trim.communicate()[1], encoding='ascii').lstrip('\n').rstrip('\n')
                    self.LOGGER.info('{}\nFinished!'.format(out))

        for prefix in self.prefix:
            if self.read_file_dict[prefix][0] is not None and self.read_file_dict[prefix][1] is not None:
                file = self.out + '/' + prefix + '-R{}_trim.fq'
                self.read_file_dict[prefix + '-trim'] = [file.format('1'), file.format('2')]
                self.current_read_files.append(file.format('1'))
                self.current_read_files.append(file.format('2'))
            else:
                self.current_read_files.append(None)
                self.current_read_files.append(None)

        return self.read_file_dict, self.current_read_files

    def converter(self):
        """
        Converting FASTQ files to FASTA using SeqIO.convert and multiprocess
        :return read_file_dict, current_read_files:
        """

        start = time.time()
        pool = multiprocessing.Pool(config.CPU)
        func = partial(converting, out=self.out)
        in_files = self.current_read_files[-4:]
        converted_files = pool.map(func, in_files)
        self.LOGGER.info('Converted files:\n' + '\n'.join(filter(None, converted_files)))
        pool.close()
        pool.join()
        end = time.time()
        self.LOGGER.info('Converting took {:0.2f}s'.format(end - start))

        read_file_dict_update = {}

        for file in filter(None, converted_files):
            file_name = os.path.splitext(os.path.basename(file))[0]
            for key, value in self.read_file_dict.items():
                for element in value:
                    if re.search(file_name, element):
                        new_key = key.split('-')[0] + '-conv'
                        if new_key not in read_file_dict_update.keys():
                            read_file_dict_update[new_key] = [file]
                        else:
                            read_file_dict_update[new_key].append(file)

        self.current_read_files = self.current_read_files + converted_files
        self.read_file_dict.update(read_file_dict_update)

        return self.read_file_dict, self.current_read_files

    def __init__(self, read_file_dict, current_read_files, adapter, log):
        self.read_file_dict = read_file_dict
        self.current_read_files = current_read_files
        self.prefix = list(self.read_file_dict.keys())
        self.user_adapter = adapter
        self.out = config.PRJ_NAME + '/' + config.DIR_TREE['read_data']
        self.LOGGER = log


class eccSampleConcat:
    """
    Unify read length, sub sample and concatenate reads.
    From processed user input (FASTA) to RepeatExplorer input (concatenated FASTA).
    """

    def prepare_for_clustering(self):

        prexer_fasta = os.path.join(self.out, 'REPEATEXPLORER_READY.fa')

        if os.path.isfile(prexer_fasta) and os.path.getsize(prexer_fasta) > 0:

            self.LOGGER.info('REPEATEXPLORER_READY.fa file found. Skipping read preparation for clustering.')

            return prexer_fasta

        self.LOGGER.info('Starting read preparation for cluster analysis.')

        self.LOGGER.info('Checking reads.')
        # get the optimal read length for each read file with minimal bases loss, use minimum overall
        start = time.time()
        pool = multiprocessing.Pool(config.CPU)
        best_read_length = pool.map(get_best_read_length, self.sorted_read_files)
        best_read_length = min(i for i in best_read_length if i is not None)
        pool.close()
        pool.join()
        end = time.time()
        self.LOGGER.info('Checking reads I/II took {:0.2f}s'.format(end - start))
        self.LOGGER.info('Optimal usable read length: {0}bp'.format(best_read_length))

        # read sets for calculation of max available reads
        extended_read_set_a = [self.sorted_read_files[0], self.sorted_read_files[1], self.prefix[0]]
        read_sets = [extended_read_set_a]

        if self.four_read_files:
            extended_read_set_b = [self.sorted_read_files[2], self.sorted_read_files[3], self.prefix[1]]
            read_sets.append(extended_read_set_b)

        # get the max available read count for each read file @ optimal read length, use minimum overall
        start = time.time()
        pool = multiprocessing.Pool(config.CPU)
        func = partial(get_max_read_count, best_read_length=best_read_length)
        reads = pool.map(func, read_sets)
        read_count_counted = {}
        for item in reads:
            read_count_counted[item[0]] = item[1]
        max_read_count = min(i for i in read_count_counted.values() if i is not None)
        pool.close()
        pool.join()
        end = time.time()
        self.LOGGER.info('Checking reads II/II took {:0.2f}s'.format(end - start))
        self.LOGGER.info('Maximal usable read count: {} @{}bp'.format(max_read_count, best_read_length))

        if max_read_count == 0:
            read_count_error = 'Read count is 0. You my need to edit the ESTIMATED_READ_LENGTH value in config file!'
            self.LOGGER.error(read_count_error)
            raise Exception(read_count_error)

        # sub sample, add prefix and suffix, concatenate
        if self.read_count_set is None:

            self.read_count_set = max_read_count

        elif self.read_count_set == 'auto':
            if config.REFERENCE_PATH is None and self.genome_size is None:

                self.read_count_set = max_read_count

                self.LOGGER.info('No genome size found. Automatically using max. available reads.')

            else:

                if self.genome_size is None:

                    refseq_size_file = os.path.splitext(config.REFERENCE_PATH)[0] + '_chrSize.txt'

                    genome_size = 0

                    with open(refseq_size_file) as size_file:
                        for line in size_file:
                            line = line.strip().split('\t')
                            genome_size += int(line[1])

                    self.LOGGER.info('Estimate genome size using given reference genome file.')

                else:
                    genome_size = self.genome_size

                self.read_count_set = math.floor(genome_size * config.REPEX_FOLD_COV / best_read_length)
                self.LOGGER.info('Using {} reads equal to {}x genome coverage.'.
                                 format(self.read_count_set, config.REPEX_FOLD_COV))

        if self.read_count_set > max_read_count:
            self.read_count_set = max_read_count
            self.LOGGER.info('User set read count is greater than available read count, '
                             'using max. available reads.')

        # prexing read file sets with set read count @ optimal read length
        self.LOGGER.info('Sub sampling read files using {} paired-end reads.\n'.format(self.read_count_set))
        start = time.time()

        pool = multiprocessing.Pool(config.CPU)

        func = partial(prexing_reads, best_read_length=best_read_length,
                       read_count_set=self.read_count_set,
                       read_count_counted=read_count_counted,
                       out=self.out)
        interlaced_fastas = pool.map(func, read_sets)

        pool.close()
        pool.join()

        try:
            self.LOGGER.info('Concatenating files.')
            # interlaced_reads = interlaced_fastas[0] + interlaced_fastas[1]
            os.system('cat {file1_} {file2_} > {out_}'.format(file1_=interlaced_fastas[0],
                                                              file2_=interlaced_fastas[1],
                                                              out_=prexer_fasta))
            try:
                self.LOGGER.info('Removing temporary fasta files.')
                os.remove(interlaced_fastas[0])
                os.remove(interlaced_fastas[1])
            except OSError:
                self.LOGGER.error('Temporary files could not be found. Delete manually!')

        except IndexError:
            self.LOGGER.info('Removing temporary fasta files.')
            os.system('mv -f {src_} {trg_}'
                      .format(src_=interlaced_fastas[0],
                              trg_=prexer_fasta))

        end = time.time()
        self.LOGGER.info('Prexing reads took {:0.2f}s'.format(end - start))

        return prexer_fasta

    def __init__(self, pre_a, pre_b, sorted_fasta_files, log, genome_size, read_count, four_read_files):
        self.sorted_read_files = sorted_fasta_files
        self.prefix = [pre_a, pre_b]
        self.out = config.PRJ_NAME + '/' + config.DIR_TREE['read_data']
        self.genome_size = genome_size
        self.read_count_set = read_count
        self.four_read_files = four_read_files
        self.LOGGER = log
