#!usr/bin/env python3

import numpy as np
import os
import subprocess
import shlex
import time

from Bio import SeqIO
from collections import OrderedDict
from functools import partial
from io import StringIO
from multiprocessing import Pool
from re import search
from scipy.signal import find_peaks, peak_widths
from operator import itemgetter

from lib import config
from lib import eccHTML_templates
from lib import eccDNA_Rcodes


def subp_bedtools_to_arr(window_file, mapping_file):
    """
     run 'bedtools coverage -mean' command and return as numpy array
    """

    cmd = config.BEDTOOLS_MEAN.format(bedtools_=config.BEDTOOLS_PATH,
                                      win_=window_file,
                                      map_=mapping_file)

    proc = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE)

    win_cov_bed = proc.communicate()[0].decode('utf-8')  # binary output to string

    win_cov_bed = np.array(list(win_cov_bed.replace('\t', '\n').split('\n')))  # string to array
    win_cov_bed = np.delete(win_cov_bed, -1)  # remove empty cell at the end from last '\n'

    win_cov_bed = win_cov_bed.reshape((win_cov_bed.shape[0] // 4), 4)  # list to table

    return win_cov_bed


def get_rough_coverage(region, mapping, mapping_dict, files_name_peak_dict, res_out):
    """
     create coverage files (mean depth over region = a or b)
     a) using window size from user input (default = 100, might also be changed in config)
     b) using calculated potential regions of eccDNA origin
     executed in parallel
    """

    if mapping is not None:
        multi_id = os.path.basename(region.split('win')[0].replace('x', ''))

        print('Calculating mean coverage of {identifier_}.MAPPING (ProcessID {id_}; MultiID {mid_})'
              .format(identifier_=mapping_dict[mapping], id_=os.getpid(), mid_=multi_id))

        win_cov_bed = subp_bedtools_to_arr(region, mapping)

        # add data set identifier at beginning of array
        head = np.array(['chr', 'start', 'end', mapping_dict[mapping]])
        head = head.reshape((head.shape[0] // 4), 4)
        win_cov_bed = np.append(head, win_cov_bed, axis=0)

        # Start peak calculation

        if res_out is not None:

            """
            peak identification using scipys peak_finder
            parameter might be changed in config (threshold and distance)
            convert to .bed-like format and save directly
            """

            # windowed coverage bed file = self.win_coverage[:, sum_cov_col]

            print('Calculating peaks in {identifier_}.COVERAGE (ProcessID {id_})'
                  .format(identifier_=mapping_dict[mapping], id_=os.getpid()))

            # calculate peaks in coverage
            win_peak, _ = find_peaks(win_cov_bed[:, -1][1:].astype('float'),
                                     threshold=config.PEAK_THRESHOLD, distance=config.PEAK_DISTANCE)

            # calculate peak width coordinates ~ candidate regions
            win_peak_width = peak_widths(win_cov_bed[:, -1][1:].astype('float'),
                                         win_peak, rel_height=0.9)

            # initialize region.bed with first peak boundaries
            try:
                regions = np.array([np.append(win_cov_bed[int(win_peak_width[2][0]), :2],
                                              win_cov_bed[int(win_peak_width[3][0]), 2:3])])
                # fill region.bed with all other peaks boundaries
                for i in range(1, win_peak.size):
                    regions = np.vstack((regions, np.append(win_cov_bed[int(win_peak_width[2][i]), :2],
                                                            win_cov_bed[int(win_peak_width[3][i]), 2:3])))

                # save and sort peak region file
                region_name = os.path.splitext(os.path.basename(region))[0]

                temp_file = res_out + '/temp/' + str(mapping_dict[mapping]) + '_' + str(region_name) + '.tmp'

                np.savetxt(temp_file, regions, delimiter='\t', fmt='%s')

            except IndexError:
                print('No peaks in {identifier_}.COVERAGE (ProcessID {id_})'
                      .format(identifier_=mapping_dict[mapping], id_=os.getpid()))

        # End peak calculation

        # remove position info and store coverage only
        win_cov_bed = win_cov_bed[:, -1][1:].astype('float')
        # win_cov_bed = np.append([mapping_dict[mapping]], win_cov_bed, axis=0)

        # return coverage per window and region as tuple with identifier
        return multi_id, win_cov_bed


def get_overall_coverage(mapping, region, mapping_dict):
    """
     create coverage files (mean depth over region = a or b)
     a) using window size from user input (default = 100, might also be changed in config)
     b) using calculated potential regions of eccDNA origin
     executed in parallel
    """

    if mapping is not None:
        print('Calculating mean coverage of {identifier_}.MAPPING (ProcessID {id_})'
              .format(identifier_=mapping_dict[mapping], id_=os.getpid()))

        # calculate coverage over each a) window b) region
        win_cov_bed = subp_bedtools_to_arr(region, mapping)

        # add data set identifier at beginning of array
        head = np.array(['chr', 'start', 'end', mapping_dict[mapping]])
        head = head.reshape((head.shape[0] // 4), 4)
        win_cov_bed = np.append(head, win_cov_bed, axis=0)

        # remove position info and store coverage only
        win_cov_bed = win_cov_bed[:, -1][1:].astype('float')
        win_cov_bed = np.append([mapping_dict[mapping]], win_cov_bed, axis=0)

        # return coverage per window or region
        return win_cov_bed


def analyze_candidate_region(candidate_key, region_dict, candidate_file_names_dict, all_mappings,
                             mapping_dict):
    """
    save sequence of candidate as FASTA
    get coverage depth of candidate per base for each mapping (summary file)
    """

    print('Analyzing {identifier_} (ProcessID {id_})'
          .format(identifier_=candidate_key, id_=os.getpid()))

    dir_abs = os.path.abspath(os.path.join(config.CURRENT_DIR, config.PRJ_NAME))
    dir_out = os.path.join(dir_abs, config.DIR_TREE['mapping_candidates'], candidate_key)
    os.makedirs(dir_out, exist_ok=True)

    region_in = region_dict[candidate_key][0]
    region_in = str(region_in).replace(':', '\\t').replace('-', '\\t')

    reg_cov_bed = []
    out_file_fasta = candidate_file_names_dict[candidate_key][0]
    out_file_sum = candidate_file_names_dict[candidate_key][1]

    if os.path.isfile(out_file_fasta) and os.path.getsize(out_file_fasta) > 0:
        print('{identifier_}: Fasta file found. Skipping blast.'
              .format(identifier_=candidate_key))

    else:
        with open(out_file_fasta, 'wt') as out:
            out.write(str('>' + candidate_key + '|' + region_dict[candidate_key][0]) + '\n')
            out.write(str(region_dict[candidate_key][1]) + '\n')

        os.system(config.BLAST_CMD.format(path_=config.BLASTn_PATH,
                                          fasta_=out_file_fasta,
                                          m6_=candidate_file_names_dict[candidate_key][4],
                                          database_=config.BLASTDB,
                                          add_=config.BLASTn_ADD_PARAMETERS))

    if os.path.isfile(out_file_sum) and os.path.getsize(out_file_sum) > 0:
        print('{identifier_}: Summary file found. Skipping analyzing.'
              .format(identifier_=candidate_key))

    else:

        for i in range(len(all_mappings)):
            j = 0
            while j <= 5:
                try:
                    curr_bed = subprocess.check_output('echo \'{region_}\' | bedtools coverage -d -a stdin -b {in_}'
                                                       .format(region_=region_in, in_=all_mappings[i]),
                                                       shell=True, universal_newlines=True)
                    j = 12
                except:
                    print('{identifier_} ({map_}): Bedtools error, trying again.'
                          .format(identifier_=candidate_key, map_=mapping_dict[all_mappings[i]]))

                    time.sleep(30)

                    j += 1

            if j != 12:
                curr_bed = 'chr\t0\t0\t0\t0\n'
                print('{identifier_}: Bedtools error, writing command to missed_analysis.sh.'
                      .format(identifier_=candidate_key))
                for file in os.listdir(os.path.join(dir_abs, config.DIR_TREE['mapping_results'])):
                    if file.endswith('missed_analysis.sh'):
                        print(file)
                        with open(os.path.join(dir_abs, config.DIR_TREE['mapping_results'], file), 'a') as missed_file:
                            missed_file.write('echo \'{region_}\' | bedtools coverage -d -a stdin -b {in_} | '
                                              'sed -e \'1 i\\chr\\tstart\\tend\\tpos\\t{map_}\' > {out_}\n'
                                              .format(region_=region_in,
                                                      in_=all_mappings[i],
                                                      out_=dir_out + '/' + mapping_dict[all_mappings[i]] + '_rerun.bed',
                                                      map_=mapping_dict[all_mappings[i]]))

            # bed to np.array and delete last empty line from last \n
            curr_bed = np.array(list(curr_bed.replace('\t', '\n').split('\n')))
            curr_bed = np.delete(curr_bed, -1)

            # reshape array from list-like to table-like format
            columns = 5
            curr_bed = curr_bed.reshape((curr_bed.shape[0] // columns), columns)

            if i == 0:
                reg_cov_bed = np.array(['chr', 'start', 'end', 'pos', mapping_dict[all_mappings[i]]])
                reg_cov_bed = reg_cov_bed.reshape(1, 5)
                reg_cov_bed = np.append(reg_cov_bed, curr_bed, axis=0)

            else:
                header = np.array(['chr', 'start', 'end', 'pos', mapping_dict[all_mappings[i]]])
                header = header.reshape(1, 5)
                curr_bed = np.append(header, curr_bed, axis=0)
                curr_bed = curr_bed[:, -1]
                curr_bed = curr_bed.reshape((curr_bed.shape[0]), 1)
                reg_cov_bed = np.append(reg_cov_bed, curr_bed, axis=1)

        np.savetxt(out_file_sum, reg_cov_bed, delimiter='\t', fmt='%s')


class eccMapper:
    """
    Align reads against reference sequence using the segemehl tool.
    Detect read peaks, split reads and discordant mapping reads.
    Extract and annotate (blast) regions of potential eccDNA origin.
    Create html report.
    """

    def run_mapping(self):
        # map reads from each paired-end read set against reference genome sequence

        i = 0
        j = 0
        while i < len(self.fastas_reads):
            if self.fastas_reads[i] is not None and self.fastas_reads[i + 1] is not None:

                if not os.path.isfile(self.fname_base[j] + '.sam') or os.path.getsize(self.fname_base[j] + '.sam') == 0:
                    # run segemehl tool generating output alignment
                    cmd = config.SEGEMEHL_CMD.format(path_=config.SEGEMEHL_PATH,
                                                     refidx_=self.idx_refseq,
                                                     ref_=self.fasta_refseq,
                                                     r1_=self.fastas_reads[i],
                                                     r2_=self.fastas_reads[i + 1],
                                                     cpu_=config.CPU,
                                                     out_=self.fname_base[j] + '.sam',
                                                     unout_=self.fname_base[j] + '_unmapped.sam')

                    run_segemehl = subprocess.Popen(shlex.split(cmd), stderr=subprocess.PIPE, universal_newlines=True)

                    out = str.format(run_segemehl.communicate()[1]).lstrip('\n').rstrip('\n')

                    self.LOGGER.info('{}'.format(out))

                else:
                    self.LOGGER.info('Mapping file found. Using existing file.')

                # rewrite .sngl.bed file eliminating uncommon errors (start > end) and rename to match convention
                self.LOGGER.info('Checking split read file.')
                os.system('awk \'$2 > $3 {{ temp = $3; $3 = $2; $2 = temp }} 1\' OFS=\'\\t\' {sngl_bed_} > tmp '
                          '&& mv -f tmp {bed_}'
                          .format(sngl_bed_=self.fname_base[j] + '.sngl.bed',
                                  bed_=self.fname_base[j] + '_aligned-SR.bed'))

                # convert .sam alignments into .bed alignments
                self.LOGGER.info('Converting alignment from SAM to BED.')
                os.system('samtools view -S {in_} -@ {cpu_} -u | '  # .bam
                          'samtools sort -@ {cpu_} -O BAM - | '  # sorted .bam
                          'bedtools bamtobed -i stdin > {out_}'  # .bed
                          .format(in_=self.fname_base[j] + '.sam',
                                  out_=self.fname_base[j] + '.bed',
                                  cpu_=config.CPU))

                # get discordant mapping reads from alignment using SAM-flags
                self.LOGGER.info('Gathering discordant mapping reads.')
                single_flag_dr_read_files = []

                for flag in config.SAMFLAGS_DR_not:
                    self.LOGGER.info('Searching for reads not mapped in proper pair.')

                    curr_flag_file = self.fname_base[j] + '_aligned-DR-nF' + str(flag) + '.bed'

                    os.system('samtools view -S -G {flags_} {in_} -@ {cpu_} -u | '
                              'samtools sort -@ {cpu_} -O BAM - | '
                              'bedtools bamtobed -i - > {out_}'
                              .format(flags_=str(flag),
                                      in_=self.fname_base[j] + '.sam',
                                      out_=curr_flag_file,
                                      cpu_=config.CPU))

                    single_flag_dr_read_files.append(curr_flag_file)

                for flag in config.SAMFLAGS_DR:
                    self.LOGGER.info('Searching for reads mapped in unusual orientation (rev-for, {}).'.format(flag))

                    curr_flag_file = self.fname_base[j] + '_aligned-DR-F' + str(flag) + '.bed'

                    os.system('samtools view -S -f {flags_} {in_} -@ {cpu_} -u | '
                              'samtools sort -@ {cpu_} -O BAM - | '
                              'bedtools bamtobed -i - > {out_}'
                              .format(flags_=str(flag),
                                      in_=self.fname_base[j] + '.sam',
                                      out_=curr_flag_file,
                                      cpu_=config.CPU))

                    single_flag_dr_read_files.append(curr_flag_file)

                self.LOGGER.info('Removing duplicates from discordant mapping reads.')

                os.system('cat {files_} |'
                          ' bedtools groupby -g 1,2,3,4,5,6 -c 4 -o count_distinct |'
                          ' sort -k1,1 -k2,2n > {out_}'
                          .format(files_=' '.join(single_flag_dr_read_files),
                                  out_=self.fname_base[j] + '_aligned-DR' + '.bed'))

                self.LOGGER.info('Finished gathering discordant mapping reads.')

                # create statistics file from .sam alignments
                self.LOGGER.info('Collecting statistics from alignment file (SAM).')

                os.system('samtools stats -s -@ {cpu_} --reference {ref_} {in_} |'
                          'grep ^SN | cut -f 2- > {out_}'
                          .format(ref_=self.fasta_refseq,
                                  in_=self.fname_base[j] + '.sam',
                                  out_=self.fname_base[j] + '_alignment-stats.txt',
                                  cpu_=config.CPU))

                # import relevant information from .sam statistics file
                with open(self.fname_base[j] + '_alignment-stats.txt') as stats:
                    for line in stats:
                        if search('reads mapped:', line):
                            line = line.split('\t')[-1]
                            self.maps_readcnt[j] = int(line.replace('\n', ''))
                        if search('bases mapped \(cigar\):', line):
                            line = line.split('\t')[-2]
                            self.maps_basecnt[j] = int(line.replace('\n', ''))
                        if search('average length:', line):
                            line = line.split('\t')[-1]
                            self.maps_avgreadlength[j] = int(line.replace('\n', ''))

                # create list with .bed alignment file names and associated identifiers (headers in summary)
                self.LOGGER.info('Summarizing mapping.')

                self.beds_maps.extend([self.fname_base[j] + '.bed',
                                       self.fname_base[j] + '_aligned-SR.bed',
                                       self.fname_base[j] + '_aligned-DR.bed'])
                self.id_maps.extend([self.prefix[j] + '_map-all',
                                     self.prefix[j] + '_map-SR',
                                     self.prefix[j] + '_map-DR'])

                self.LOGGER.info('Finished mapping: {}'.format(self.prefix[j]))

                i += 2
                j += 1
            else:
                i += 2
                j += 1

    def run_splitread_detect(self):

        # summarize split read mappings using segemehls haarz function

        self.LOGGER.info('Calculating regions from split reads.')

        fname_haarz_out = self.fname_base[0] + '_haarz-SR.bed'

        os.system(config.HAARZ_CMD.format(path_=config.HAARZ_PATH,
                                          files_=self.fname_base[0] + '_aligned-SR.bed',
                                          out_=fname_haarz_out))

        # create .bed file containing the potential eccDNA regions according to split reads
        # remove too large and small regions (might be changed in config)
        # merge regions in close surroundings

        self.LOGGER.info('Merging and cleaning up regions.')

        os.system('sed \'1d\' {in_} | '
                  'bedtools sort -i stdin | '
                  'awk \'($3-$2)<={max_} && ($4)>0\' | '
                  'bedtools merge -d {merge_} -i stdin |'
                  'awk \'($3-$2)<={max_} && ($3-$2)>={min_}\' '
                  '> {out_}'
                  .format(in_=fname_haarz_out,
                          out_=self.fname_base[0] + '_regions-SR.bed',
                          max_=config.MAX_eccDNA_LENGTH,
                          min_=config.MIN_eccDNA_LENGTH,
                          merge_=config.MERGE_CLOSE_REGIONS))

        # extend regions by its own size to cover upstream and downstream regions, not in use
        # os.system('awk \'($4=($3-$2))\' {in_} | '
        #          'awk \'($2=$2-$4)\' | '
        #          'awk \'($3=$3+$4)\' | '
        #          'awk \'BEGIN{{OFS="\\t"}}{{print $0}}\' > {out_}'
        #          .format(in_=self.files_name_base[0] + '_regions-SR.bed',
        #                  out_=self.files_name_base[0] + '_extended-regions-SR.bed'))

    def run_discordantread_detect(self):

        # summarize discordant mapped reads detected with SAM-flags

        self.LOGGER.info('Calculating genome coverage from discordant mapping reads.')

        bed_braph = self.fname_base[0] + '_aligned-DR_graph.bed'

        os.system('bedtools genomecov -bga -i {map_DR_} -g {ref_size_} > {graph_DR_}'
                  .format(map_DR_=self.fname_base[0] + '_aligned-DR.bed',
                          ref_size_=self.txt_refseq_size,
                          graph_DR_=bed_braph))

        self.LOGGER.info('Merging and cleaning up regions.')

        cmd_getmaxdepth = 'awk -v max=0 \'{{if($4>max){{max=$4}}}}END{{print max}}\' {}'.format(bed_braph)

        max_coverage = subprocess.check_output(shlex.split(cmd_getmaxdepth))

        min_coverage_allowed = int(float(max_coverage) * config.BACKGROUND_PERC)

        os.system('awk \'$4 > {min_cov_}\' {in_} | '
                  'bedtools merge -d {merge_} -i stdin | '
                  'awk \'($3-$2)<={max_} && ($3-$2)>={min_}\' '
                  '> {out_}'
                  .format(in_=bed_braph,
                          min_cov_=min_coverage_allowed,
                          merge_=config.MERGE_CLOSE_REGIONS,
                          max_=config.MAX_eccDNA_LENGTH,
                          min_=config.MIN_eccDNA_LENGTH,
                          out_=self.fname_base[0] + '_regions-DR.bed'))

    def extract_candidate_regions(self, all_ident, DR_ident):

        # extract eccDNA candidate regions by comparing SR regions with peak regions (all, DR)

        # check .bed files
        bed_files = [self.fname_base[0] + '_regions-SR.bed', self.dict_idmaps_peakregion[all_ident],
                     self.dict_idmaps_peakregion[DR_ident]]
        for file in bed_files:
            os.system('awk \'($2<$3)\' {bed_} > tmp && mv tmp {bed_}'
                      .format(bed_=file))

        # high confident regions containing split reads, discordant reads, and high overall read coverage (3/3)
        os.system('bedtools intersect -u -a {SR_} -b {all_} | bedtools intersect -u -a stdin -b {DR_} > {out_}'
                  .format(SR_=self.fname_base[0] + '_regions-SR.bed',
                          all_=self.dict_idmaps_peakregion[all_ident],
                          DR_=self.dict_idmaps_peakregion[DR_ident],
                          out_=self.fname_base[0] + '_hiconf-ECC-REGIONS.bed'))

        # confident regions containing split reads, and discordant reads, or high overall read coverage (2/3)

        out_file_SRall = self.fname_base[0] + '_lowconf-ECC-regions_SR-all.bed'
        out_file_SRDR = self.fname_base[0] + '_lowconf-ECC-regions_SR-DR.bed'
        out_file_DRall = self.fname_base[0] + '_lowconf-ECC-regions_DR-all.bed'

        # only SR and high overall read coverage, not DR
        os.system('bedtools intersect -u -a {SR_} -b {all_} | bedtools intersect -v -a stdin -b {DR_} > {out_}'
                  .format(SR_=self.fname_base[0] + '_regions-SR.bed',
                          all_=self.dict_idmaps_peakregion[all_ident],
                          DR_=self.dict_idmaps_peakregion[DR_ident],
                          out_=out_file_SRall))

        # only SR and DR, not all
        os.system('bedtools intersect -u -a {SR_} -b {DR_} | bedtools intersect -v -a stdin -b {all_} > {out_}'
                  .format(SR_=self.fname_base[0] + '_regions-SR.bed',
                          all_=self.dict_idmaps_peakregion[all_ident],
                          DR_=self.dict_idmaps_peakregion[DR_ident],
                          out_=out_file_SRDR))

        # only all and DR, not SR
        os.system('bedtools intersect -u -a {all_} -b {DR_} | bedtools intersect -v -a stdin -b {SR_} > {out_}'
                  .format(SR_=self.fname_base[0] + '_regions-SR.bed',
                          all_=self.dict_idmaps_peakregion[all_ident],
                          DR_=self.dict_idmaps_peakregion[DR_ident],
                          out_=out_file_DRall))

        out_file = self.fname_base[0] + '_lowconf-ECC-regions.bed'

        os.system('cat {noDR_} {noALL_} {noSR_} | sort  -k1,1 -k2,2n > {out_}'
                  .format(noDR_=out_file_SRall,
                          noALL_=out_file_SRDR,
                          noSR_=out_file_DRall,
                          out_=out_file))

        bed_ecc_region = self.fname_base[0] + '_hiconf-ECC-REGIONS.bed'

        if os.path.getsize(bed_ecc_region) == 0:

            bed_ecc_region = self.fname_base[0] + '_lowconf-ECC-regions.bed'

        os.system('bedtools intersect -c -a {ref_win_} -b {regions_} | '
                  'sed -e \'1 i\\chr\\tstart\\tend\\thiconf\' > {out_}'
                  .format(ref_win_=self.bed_refseq_win,
                          regions_=bed_ecc_region,
                          out_=self.dir_result + '/temp/report_hiconf_win.tmp'))

    def html_summarize(self):

        # create overview summary

        manhattan_rel = os.path.join('mapping_results', self.prefix[0] + '{}'
                                     .format('-' + self.prefix[1] if len(self.fname_base) == 2 else '') +
                                     '_chr_manhattan-plot' + config.IMAGE_EXT[config.IMAGE_TYPE])

        with open(os.path.join(config.PRJ_NAME, config.DIR_TREE['results'], 'eccMap_summary.html'), 'a') as html_file:
            html_file.write(eccHTML_templates.eccMap_HTML_sum.format(manhattan_=manhattan_rel,
                                                                     pre1_=self.prefix[0],
                                                                     reads1_=self.maps_readcnt[0],
                                                                     bases1_=self.maps_basecnt[0],
                                                                     pre2_=self.prefix[1],
                                                                     reads2_=self.maps_readcnt[1],
                                                                     bases2_=self.maps_basecnt[1], ) + '\n')
            html_file.write(eccHTML_templates.eccMap_HTML_th + '\n')

        for key in self.candidate_list:

            blast_res = ''
            enr_score = ''

            with open(self.candidate_file_names_dict[key][4]) as blast_file:
                score = 0
                for line in blast_file:
                    line = line.strip().split('\t')
                    new_score = int(float(line[11]))
                    if new_score > score:
                        score = new_score
                        blast_file_rel = os.path.join('mapping_results', 'candidates', key, key + '_blast.m6')
                        blast_res = '<a href="{}">{}</a>'.format(blast_file_rel, str(line[1]))

            with open(self.csv_sumreg_nrm) as enr_file:
                for line in enr_file:
                    if search(key, line):
                        line = line.replace('\n', '').split('\t')
                        enr_score = line[9]

            with open(os.path.join(config.PRJ_NAME, config.DIR_TREE['results'], 'eccMap_summary.html'),
                      'a') as html_file:
                html_file.write(eccHTML_templates.eccMap_HTML_td.format(
                    cand_=key,
                    plot_=self.candidate_file_names_dict[key][5],
                    seq_=self.region_dict[key][1],
                    enr_=enr_score,
                    len_=self.region_dict[key][2],
                    pos_=self.region_dict[key][0],
                    blast_=blast_res) + '\n')

    def mapper_coordinator(self):
        """
        input: trimmed paired-end read files from probe and control
        process: map reads against reference genome sequence, calculate eccDNA candidate regions
        output: overview for whole genome and detailed information about eccDNA candidate regions
        """

        # check for existing index data
        if os.path.isfile(self.idx_refseq) and os.path.getsize(self.idx_refseq) > 0:
            self.LOGGER.info('Index file for segemehl mapping found.')
        else:
            self.LOGGER.info('Creating index file for segemehl mapping.')

            cmd = '{0} -x {1} -d {2}'.format(config.SEGEMEHL_PATH, self.idx_refseq, self.fasta_refseq)

            create_idx = subprocess.Popen(shlex.split(cmd), stderr=subprocess.PIPE, universal_newlines=True)

            out = str.format(create_idx.communicate()[1]).lstrip('\n').rstrip('\n')

            self.LOGGER.info('{}'.format(out))

        # check for existing chromosome size file and read in, or create one
        if os.path.isfile(self.txt_refseq_size) and os.path.getsize(self.txt_refseq_size) > 0:
            with open(self.txt_refseq_size, 'r') as chr_size_file:
                for line in chr_size_file:
                    line = line.replace('\n', '').split('\t')
                    self.refseq_chrs.append(str(line[0]))
                    self.refseq_chrsizes.append(int(line[1]))
            self.LOGGER.info('Chromosome size file found.')
        else:
            self.LOGGER.info('Creating chromosome size file.')
            with open(self.fasta_refseq) as reference_sequences:
                for record in SeqIO.parse(reference_sequences, "fasta"):
                    self.refseq_chrs.append(record.id)
                    self.refseq_chrsizes.append(len(record))
            with open(self.txt_refseq_size, 'wt') as chr_size_file:
                for i in range(len(self.refseq_chrsizes)):
                    line = '{0}\t{1}\n'.format(self.refseq_chrs[i], self.refseq_chrsizes[i])
                    chr_size_file.write(line)

        # cytoBand file not yet needed (no visualization w/ circos)
        # check for existing cytoBand file for circlize visualization
        # if os.path.isfile(self.ref_seq_cyt) and os.path.getsize(self.ref_seq_cyt) > 0:
        #     self.LOGGER.info('CytoBand file found.')
        # else:
        #     self.LOGGER.info('Creating cytoBand file.')
        #     with open(self.ref_seq_cyt, 'wt') as chr_cyt_file:
        #         for i in range(len(self.ref_seq_size_list)):
        #             line = '{0}\t0\t{1}\t{0}\tgeng\n'.format(self.ref_seq_chr_list[i], self.ref_seq_size_list[i])
        #             chr_cyt_file.write(line)

        # check for existing reference genome sequence window file
        if os.path.isfile(self.bed_refseq_win) and os.path.getsize(self.bed_refseq_win) > 0:
            self.LOGGER.info('Reference genome sequence window file found.')
        else:
            self.LOGGER.info('Creating Reference genome sequence window file.')
            os.system('bedtools makewindows -g {in_} -w {win_} > {out_}'
                      .format(in_=self.txt_refseq_size,
                              win_=self.int_winsize,
                              out_=self.bed_refseq_win))

        # run mapping using the segemehl tool
        self.LOGGER.info('Start: Map reads against reference genome sequence.')

        self.run_mapping()

        self.LOGGER.info('Finished: Mapped reads against reference genome sequence!')

        # get candidate regions from split reads (run haarz)
        self.LOGGER.info('Start: Summarize mapped split reads (SR) and calculate SR regions.')

        self.run_splitread_detect()

        self.maps_readcnt = np.array(self.maps_readcnt)
        self.maps_basecnt = np.array(self.maps_basecnt)
        self.maps_avgreadlength = np.array(self.maps_avgreadlength)

        self.LOGGER.info('Finished: Summarized mapped SR and calculated SR regions.')

        # get candidate regions from discordant mapping reads
        self.LOGGER.info('Start: Summarize discordant mapped reads (DR) and calculate DR regions.')

        self.run_discordantread_detect()

        self.LOGGER.info('Start: Summarized mapped DR and calculated DR regions.')

        # calculating coverage files and peaks
        for i in range(len(self.beds_maps)):  # beds_maps are 3 or 6
            self.dict_maps_id[self.beds_maps[i]] = self.id_maps[i]

        for i in range(len(self.id_maps)):
            self.dict_idmaps_peakregion[self.id_maps[i]] = '{}'.format(self.fname_base[0]
                                                                       if search(self.prefix[0],
                                                                                 self.id_maps[i])
                                                                       else self.fname_base[1]) \
                                                           + '_regions-' \
                                                           + self.id_maps[i].split('_map-')[-1] \
                                                           + '.bed'

        # check if coverage and peak calculations has been done in run before
        if not os.path.isfile(self.csv_sum_raw) or os.path.getsize(self.csv_sum_raw) == 0:

            with open(self.bed_refseq_win) as window_file:
                for line in window_file:
                    line = line.replace('\n', '').split('\t')
                    self.win_coverage.append(line)

            self.win_coverage = np.array(self.win_coverage)

            self.LOGGER.info('Calculating rough coverage and find peaks.')
            start = time.time()

            os.system('cd {dir_ref_} && split -n l/{cpu_} -a 2 -d --additional-suffix win{win_}-part.bed {refseq_}'
                      .format(dir_ref_=self.dir_refseq,
                              cpu_=config.CPU,
                              win_=self.int_winsize,
                              refseq_=self.bed_refseq_win))

            window_file_parts = []

            for file in os.listdir(self.dir_refseq):
                if file.endswith('-part.bed'):
                    window_file_parts.append(os.path.join(self.dir_refseq, file))

            func = partial(get_rough_coverage,
                           mapping=self.beds_maps[0],
                           mapping_dict=self.dict_maps_id,
                           res_out=self.dir_result,
                           files_name_peak_dict=self.dict_idmaps_peakregion)

            pool = Pool(config.CPU)

            tuples_coverages = pool.map(func, window_file_parts)  # part info, array with coverage info

            # print(list_of_tuples_coverages[2][1][0:10])  # indices: list, tuple, array
            # tuples with part info to sorted summarized array
            self.sum_coverages_tr = sorted(tuples_coverages, key=itemgetter(0))[0][1][:]

            for i in range(1, len(tuples_coverages)):
                self.sum_coverages_tr = np.append(self.sum_coverages_tr,
                                                  sorted(tuples_coverages, key=itemgetter(0))[i][1][:])

            # reshape from row to column
            self.sum_coverages_tr = np.reshape(self.sum_coverages_tr, (len(self.sum_coverages_tr), 1),)

            # combine mapping identifier with coverage value
            self.sum_coverages_tr = np.append(np.array([[self.dict_maps_id[self.beds_maps[0]]]]),
                                              self.sum_coverages_tr, axis=0)

            pool.close()
            pool.join()

            if len(self.beds_maps) == 6:

                func = partial(get_rough_coverage,
                               mapping=self.beds_maps[3],
                               mapping_dict=self.dict_maps_id,
                               res_out=None,
                               files_name_peak_dict=self.dict_idmaps_peakregion)

                pool = Pool(config.CPU)

                tuples_coverages = pool.map(func, window_file_parts)  # part info, array with coverage info

                # print(list_of_tuples_coverages[2][1][0:10])  # indices: list, tuple, array
                # tuples with part info to sorted summarized array
                self.sum_coverages_co = sorted(tuples_coverages, key=itemgetter(0))[0][1][:]

                for i in range(1, len(tuples_coverages)):
                    self.sum_coverages_co = np.append(self.sum_coverages_co,
                                                      sorted(tuples_coverages, key=itemgetter(0))[i][1][:])

                # reshape from row to column
                self.sum_coverages_co = np.reshape(self.sum_coverages_co, (len(self.sum_coverages_co), 1), )

                # combine mapping identifier with coverage value
                self.sum_coverages_co = np.append(np.array([[self.dict_maps_id[self.beds_maps[3]]]]),
                                                  self.sum_coverages_co, axis=0)

                pool.close()
                pool.join()

                self.win_coverage = np.concatenate((self.win_coverage, self.sum_coverages_tr, self.sum_coverages_co),
                                                   axis=1)
            else:
                self.win_coverage = np.concatenate((self.win_coverage, self.sum_coverages_tr), axis=1)

            # save coverage summary table
            np.savetxt(self.csv_sum_raw, self.win_coverage, delimiter='\t', fmt='%s')

            window_file_parts = []

            for prefix in self.prefix:
                for file in os.listdir(os.path.join(self.dir_result, 'temp')):
                    if file.startswith(prefix):
                        if file.endswith('-part.tmp'):
                            window_file_parts.append(os.path.join(self.dir_result, 'temp', file))

                os.system('cat {ins_} | sort -k1,1 -k2,2n > {out_}'
                          .format(ins_=' '.join(x for x in window_file_parts),
                                  out_=os.path.join(self.dir_result, prefix + '_regions-all.bed')))

            end = time.time()
            self.LOGGER.info('Coverage calculation and peak finding took {:0.2f}s.'.format(end - start))

        else:
            self.LOGGER.info('Coverage and peak files found.')

        # Extract high confident and confident regions for eccDNA candidates
        self.LOGGER.info('Extracting eccDNA candidate regions.')
        start = time.time()
        self.extract_candidate_regions(self.id_maps[0], self.id_maps[2])
        end = time.time()
        self.LOGGER.info('Extracting eccDNA candidate regions took {:0.2f}s.'.format(end - start))

        self.LOGGER.info('Normalizing coverage data and editing Rscript for visualization.')
        start = time.time()

        # Rearrange data and convert to Reads per million mapped reads (RPM)
        self.conn.voidEval(eccDNA_Rcodes.Rconvert_genome)
        self.conn.r.convgenome(self.csv_sum_raw, self.maps_basecnt, self.csv_sum_nrm)

        # Create plotting script: mean coverage over genome
        png_path = self.fname_base[0] + '{}' \
            .format('-' + self.prefix[1] if len(self.fname_base) == 2 else '') \
                   + '_chr_manhattan-plot' + config.IMAGE_EXT[config.IMAGE_TYPE]

        with open(os.path.join(config.PRJ_NAME, config.DIR_TREE['results'], 'viz.R'), 'a') as viz_file:
            viz_file.write('''\n\n#Rscript for visualization: mapping Manhattan plot
                        \nfile.in <- "{file_in_}"\nfile.out <- "{png_}"\nwin.size <- {win_}\nfile.add <- "{file_temp_}"
                        \n{cmd_manhattan_}\nmanhattan(file.in, file.out, win.size, file.add)
                        '''.format(file_in_=self.csv_sum_nrm,
                                   png_=png_path,
                                   win_=self.int_winsize,
                                   file_temp_=self.dir_result + '/temp/report_hiconf_win.tmp',
                                   cmd_manhattan_=eccDNA_Rcodes.rmanhattan_plot()))

        end = time.time()
        self.LOGGER.info('Normalizing took {:0.2f}s.'.format(end - start))

        # Extract sequences from (high) confident regions
        # Blast sequences
        # Get and plot (with annotations) mappings from regions and regions + surroundings

        bed_ecc_regions = self.fname_base[0] + '_hiconf-ECC-REGIONS.bed'

        if os.path.getsize(bed_ecc_regions) == 0:
            self.LOGGER.info('No high confidants eccDNA candidate regions found.\n'
                             'Continuing with low confidants eccDNA candidate regions.')

            bed_ecc_regions = self.fname_base[0] + "_lowconf-ECC-regions.bed"

        if os.path.getsize(bed_ecc_regions) == 0:

            self.LOGGER.info('No potential eccDNA candidate regions found. Sorry.')

            self.analysis_errors = True

            with open(self.fasta_sumreg, 'wt') as out_file:
                out_file.write('>no candidate found by mapping\n')
                out_file.write('NNNNNN\n')

            return self.csv_sum_nrm, self.fasta_sumreg, self.analysis_errors

        # create base columns for region coverage file (chr, start, end) and make it an numpy array
        with open(bed_ecc_regions) as region_file:
            for line in region_file:
                line = line.replace('\n', '').split('\t')
                self.region_coverage.append(line)
        self.region_coverage = np.array(self.region_coverage)

        # Calculate the mean coverage for each high confident region
        # -> position of region and mean coverage of all, SR and DR per data set
        self.LOGGER.info('Calculating detailed mean coverage over eccDNA candidate regions.')

        if not os.path.isfile(self.csv_sumreg_raw) or os.path.getsize(self.csv_sumreg_raw) == 0:

            start = time.time()

            # reuse get_rough_coverage w/ high confident region file
            region = bed_ecc_regions
            func = partial(get_overall_coverage, mapping_dict=self.dict_maps_id, region=region)

            pool = Pool(config.CPU)
            self.sum_region_coverages = np.array(pool.map(func, self.beds_maps))
            self.sum_region_coverages = np.swapaxes(self.sum_region_coverages, 0, 1)
            self.sum_region_coverages = np.concatenate((self.region_coverage, self.sum_region_coverages), axis=1)
            pool.close()
            pool.join()
            end = time.time()
            self.LOGGER.info('Coverage calculation took {:0.2f}s.'.format(end - start))

            np.savetxt(self.csv_sumreg_raw, self.sum_region_coverages, delimiter='\t', fmt='%s')

        else:
            self.LOGGER.info('Coverage and peak files found.')

        # Normalize coverage data and calculate fold enrichment for each high confident region w/ R

        if self.conn is not None:
            self.conn.voidEval(eccDNA_Rcodes.Rconvert_region)
            self.conn.r.convregion(self.csv_sumreg_raw, self.maps_basecnt, self.csv_sumreg_nrm)

            self.conn.voidEval(eccDNA_Rcodes.Renrichment)
            self.conn.r.enrichment(self.csv_sumreg_nrm, self.csv_sum_nrm)

        eccCand_fasta = subprocess.check_output(['bedtools', 'getfasta',
                                                 '-fi', self.fasta_refseq,
                                                 '-bed', bed_ecc_regions],
                                                universal_newlines=True)

        with open(self.csv_sumreg_nrm) as sum_reg_file:
            for line in sum_reg_file:
                line = line.replace('\n', '').split('\t')
                self.region_dict[line[-1]] = '{chr_}:{start_}-{end_}'.format(chr_=line[0],
                                                                             start_=line[1],
                                                                             end_=line[2])
        if 'id' in self.region_dict:
            del self.region_dict['id']

        self.region_dict = OrderedDict(sorted(self.region_dict.items()))

        for record in SeqIO.parse(StringIO(eccCand_fasta), 'fasta'):
            for cand, pos in self.region_dict.items():
                if record.id == pos:
                    record.description = record.id
                    record.id = cand
                    self.region_dict[cand] = (record.description, record.seq, len(record.seq))

        with open(self.fasta_sumreg, 'wt') as out_file:
            for key in self.region_dict.keys():
                out_file.write(str('>' + key + '|' + self.region_dict[key][0]) + '\n')
                out_file.write(str(self.region_dict[key][1]) + '\n')
                self.candidate_list.append(key)

        self.LOGGER.info('Analyzing eccDNA candidate regions.')
        start = time.time()

        with open(self.dir_result + '/missed_analysis.sh', 'w+') as missed_file:
            missed_file.write('# Sometimes bedtools fails due to an unknown reason.\n'
                              '# The missed runs will be collected here. Rerun manually if files are missing.\n')

        def sortby(x):
            try:
                return int(x.split('_')[1])
            except ValueError:
                return float('inf')

        self.candidate_list.sort(key=sortby)
        self.candidate_list = self.candidate_list[:config.MAX_CAND_CNT]

        for cand in self.candidate_list:
            base_name = os.path.join(self.dir_canddt, cand, cand)
            out_file_fasta = base_name + '.fasta'
            out_file_sum = base_name + '_coverage-summary_RAW.csv'
            out_file_sum_converted = base_name + '_coverage-summary_NORM.csv'
            out_file_png = base_name + '_collage-plot' + config.IMAGE_EXT[config.IMAGE_TYPE]
            rel_file_png = os.path.join('mapping_results/candidates', cand,
                                        cand + '_collage-plot' + config.IMAGE_EXT[config.IMAGE_TYPE])
            out_file_m6 = base_name + '_blast.m6'
            self.candidate_file_names_dict[cand] = (out_file_fasta, out_file_sum, out_file_sum_converted, out_file_png,
                                                    out_file_m6, rel_file_png)

        func = partial(analyze_candidate_region,
                       region_dict=self.region_dict,
                       candidate_file_names_dict=self.candidate_file_names_dict,
                       all_mappings=self.beds_maps,
                       mapping_dict=self.dict_maps_id)

        pool = Pool(config.CPU)
        pool.map(func, self.candidate_list)
        pool.close()
        pool.join()

        count_lines = 0

        for file in os.listdir(self.dir_result):
            if file.endswith('missed_analysis.sh'):
                with open(os.path.join(self.dir_result, file)) as missed_file:
                    for line in missed_file:
                        count_lines += 1

        if count_lines >= 3:
            self.LOGGER.info('Try to rerun missed. Missed eccDNA candidate analysis will result in '
                             'visualization errors.\nManual editing of either the viz.R file or the '
                             'eccCand_0_coverage-summary_NORM.csv needed.\n'
                             'Missed bedtools commands are summed up in the missed_analysis.sh.')
            self.analysis_errors = True
            try:
                os.system('{}'.format(self.dir_result + '/missed_analysis.sh'))
            except:
                print('')

        end = time.time()
        self.LOGGER.info('Analyzing took {:0.2f}s.'.format(end - start))

        # Normalize coverage data and calculate fold enrichment for each high confident region w/ R
        if self.conn is not None:

            self.LOGGER.info('Starting normalization of eccDNA candidates.')
            start = time.time()

            self.conn.voidEval(eccDNA_Rcodes.Rconvert_cand)
            for candidate in self.candidate_list:
                self.conn.r.convcand(self.candidate_file_names_dict[candidate][1], self.maps_basecnt,
                                     self.candidate_file_names_dict[candidate][2])

            if len(self.fastas_reads) == 4:
                with open(os.path.join(config.PRJ_NAME, config.DIR_TREE['results'], 'viz.R'), 'a') as viz_file:
                    viz_file.write('\n\n#Rscript for visualization: detailed mapping plot of eccCandidates\n\n'
                                   '{cmd_r_}\n\n'
                                   'candidate.list <- c({candidates_})\n'
                                   'ext.in <- "{ext_csv_}"\n'
                                   'ext.out <-"{ext_png_}"\n'
                                   'dir.candidate <- "{dir_canddt_}"\n\n'
                                   'for (candidate in candidate.list)\n'
                                   '{{\n'
                                   'file.in <- paste0(dir.candidate, "/", candidate, "/", candidate, ext.in)\n'
                                   'file.out <- paste0(dir.candidate, "/", candidate, "/", candidate, ext.out)\n'
                                   'linemulti(candidate, file.in, file.out)\n}}'
                                   .format(cmd_r_=eccDNA_Rcodes.rline_multiplot(),
                                           candidates_=', '.join('"{}"'.format(str(x)) for x in self.candidate_list),
                                           ext_csv_='_coverage-summary_NORM.csv',
                                           ext_png_='_collage-plot' + config.IMAGE_EXT[config.IMAGE_TYPE],
                                           dir_canddt_=self.dir_canddt))

            else:
                with open(os.path.join(config.PRJ_NAME, config.DIR_TREE['results'], 'viz.R'), 'a') as viz_file:
                    viz_file.write('\n\n#Rscript for visualization: detailed mapping plot of eccCandidates\n\n'
                                   '{cmd_r_}\n\n'
                                   'candidate.list <- c({candidates_})\n'
                                   'ext.in <- "{ext_csv_}"\n'
                                   'ext.out <-"{ext_png_}"\n'
                                   'dir.candidate <- "{dir_canddt_}"\n\n'
                                   'for (candidate in candidate.list)\n'
                                   '{{\n'
                                   'file.in <- paste0(dir.candidate, "/", candidate, "/", candidate, ext.in)\n'
                                   'file.out <- paste0(dir.candidate, "/", candidate, "/", candidate, ext.out)\n'
                                   'linemulti_noco(candidate, file.in, file.out)\n}}'
                                   .format(cmd_r_=eccDNA_Rcodes.rline_multiplot_noco(),
                                           candidates_=', '.join('"{}"'.format(str(x)) for x in self.candidate_list),
                                           ext_csv_='_coverage-summary_NORM.csv',
                                           ext_png_='_collage-plot' + config.IMAGE_EXT[config.IMAGE_TYPE],
                                           dir_canddt_=self.dir_canddt))

            end = time.time()
            self.LOGGER.info('Normalization took {:0.2f}s.'.format(end - start))

        self.LOGGER.info('Creating eccMap summary HTML.')

        with open(os.path.join(config.PRJ_NAME, config.DIR_TREE['results'], 'eccMap_summary.html'), 'w+') as html_file:
            html_file.write(eccHTML_templates.eccHTML_style + '\n')

        self.html_summarize()

        return self.csv_sum_nrm, self.fasta_sumreg, self.analysis_errors

    def __init__(self, seq_files, ref_seq_file, prefix, win_size, log, conn):
        # internal input
        self.fastas_reads = seq_files
        self.fasta_refseq = ref_seq_file
        self.prefix = prefix
        self.int_winsize = win_size
        self.fname_refseq = os.path.splitext(os.path.basename(self.fasta_refseq))[0]
        self.LOGGER = log
        self.conn = conn

        # internal temporary variables
        self.refseq_chrs = []
        self.refseq_chrsizes = []

        self.fnames_reads_fas = []
        for file in self.fastas_reads:
            if file is not None:
                self.fnames_reads_fas.append(os.path.splitext(os.path.basename(file))[0])

        self.beds_maps = []  # List of bed file names of mappings [all, SR, DR]
        self.id_maps = []  # List of IDs: prefix + _map- + [all, SR, DR]
        self.dict_maps_id = {}
        self.maps_readcnt = [0, 0]
        self.maps_basecnt = [0, 0]
        self.maps_avgreadlength = [0, 0]
        self.maps_avgreaddepth = []  # average depth per base to calculate enrichment w/o control data

        self.win_coverage = [['chr', 'start', 'end']]
        self.region_coverage = [['chr', 'start', 'end']]
        self.win_peak_dict = {}
        self.win_peak_width_dict = {}

        self.sum_coverages_tr = []
        self.sum_coverages_co = []
        self.sum_region_coverages = []

        self.dict_idmaps_peakregion = {}

        self.region_dict = {}
        self.candidate_list = []
        self.candidate_file_names_dict = {}

        # output directories
        self.dir_temp = os.path.join(config.PRJ_NAME, config.DIR_TREE['mapping_results'], 'temp')
        os.makedirs(self.dir_temp, exist_ok=True)
        self.dir_abs = os.path.abspath(os.path.join(config.CURRENT_DIR, config.PRJ_NAME))
        self.dir_refseq = os.path.join(self.dir_abs, config.DIR_TREE['reference_data'])
        self.dir_result = os.path.join(self.dir_abs, config.DIR_TREE['mapping_results'])
        self.dir_canddt = os.path.join(self.dir_abs, config.DIR_TREE['mapping_candidates'])

        # output files: reference related files
        self.idx_refseq = self.dir_refseq + '/' + self.fname_refseq + '.idx'
        self.txt_refseq_size = self.dir_refseq + '/' + self.fname_refseq + '_chrSize' + '.txt'
        self.txt_refseq_cyto = self.dir_refseq + '/' + self.fname_refseq + '_cytoBand' + '.txt'
        self.bed_refseq_win = self.dir_refseq + '/' + self.fname_refseq + '_win{}'.format(self.int_winsize) + '.bed'

        # output files: file name base
        self.fname_base = []
        if len(self.fastas_reads) == 4:
            for pre in self.prefix:
                self.fname_base.append(self.dir_result + '/' + pre)
        else:
            self.fname_base.append(self.dir_result + '/' + self.prefix[0])

        # output files: summary files
        self.fname_sum = self.fname_base[0] + '{}'.format('-' + self.prefix[1] if len(self.fname_base) == 2 else '')

        self.csv_sum_raw = self.fname_sum + '_summary_coverages.csv'
        self.csv_sum_nrm = self.fname_sum + '_summary_coverages_normalized.csv'
        self.csv_sumreg_raw = self.fname_sum + '_summary_region-coverages.csv'
        self.csv_sumreg_nrm = self.fname_sum + '_summary_region-coverages_normalized.csv'
        self.fasta_sumreg = self.fname_sum + '_ECC-SEQUENCES.fasta'

        # internal outpu
        self.analysis_errors = False
