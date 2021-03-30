#!/usr/bin/env python3

import numpy as np
import os
import shutil
import subprocess
import time

from collections import OrderedDict
from functools import partial
from multiprocessing import Pool
from re import search

from lib import config
from lib import eccHTML_templates
from lib import eccDNA_Rcodes


def blast_cand(cl, repex_dir, comp_dir, cand_blastdb):
    # blastn all cluster candidates against a database of map candidates
    print('Blasting candidate cluster: CL{}'.format(cl))

    cl_dir = 'dir_CL' + str(cl).zfill(4)
    fasta = os.path.join(repex_dir, 'seqclust/clustering/clusters',
                         cl_dir, 'reads.fas.CL{}.contigs'.format(cl))
    m6 = os.path.join(comp_dir, 'blast_CL{}.m6'.format(cl))

    if not os.path.isfile(fasta):
        # follow new RepeatExplorer file naming convention
        fasta = os.path.join(repex_dir, 'seqclust/clustering/clusters',
                             cl_dir, 'contigs.fasta'.format(cl))

    if os.path.isfile(fasta) and os.path.getsize(fasta) > 0:
        os.system(config.BLAST_CMD.format(path_=config.BLASTn_PATH,
                                          fasta_=fasta,
                                          m6_=m6,
                                          database_=cand_blastdb,
                                          add_=config.BLASTn_ADD_PARAMETERS))

    if os.path.isfile(m6) and os.path.getsize(m6) == 0:
        print('CL{} had no hit.'.format(cl))
        os.remove(m6)


class eccComparer:
    """
    Compare the found candidates from eccMapper with candidates from eccClusterer and find matches
    """

    def compare_coordinator(self):

        # create blastDB from eccMapper candidates
        self.LOGGER.info('Creating BLAST+ database from eccMapper candidates.')

        shutil.copy(self.map_fas_files, self.blastdb_fas)

        os.system('{path_} -in {blastdb_} -dbtype nucl'
                  .format(path_=config.BLASTMAKEDB_PATH,
                          blastdb_=self.blastdb_fas))

        self.LOGGER.info('Blast comparison of candidates from eccMapper and eccClusterer.')
        start = time.time()
        pool = Pool(config.CPU)
        func = partial(blast_cand, repex_dir=self.dir_clu, comp_dir=self.dir_res, cand_blastdb=self.blastdb_fas)
        pool.map(func, self.cl_list)
        pool.close()
        pool.join()
        end = time.time()
        self.LOGGER.info('Blast calculation took {:0.2f}s.'.format(end - start))

        # create a dictionary with mapper candidates as keys and from blast corresponding clusters as values
        self.LOGGER.info('Reporting candidates from mapping and clustering.')
        for file in os.listdir(self.dir_res):
            if file.endswith('.m6'):
                with open(os.path.join(self.dir_res, file), newline='\n') as f:
                    cluster = file.replace('blast_CL', '').replace('.m6', '')
                    self.cand_cl_list.append(cluster)
                    for line in f:
                        line = line.strip().replace('|', '\t').split('\t')
                        map_cand = line[1]
                        map_pos = line[2].replace('-', ':').split(':')
                        if map_cand not in self.map_cand_dict.keys():
                            self.map_cand_dict[map_cand] = []
                            self.map_cand_dict[map_cand].append(int(cluster))
                            self.map_fas_pos_bed += '{chr_}\t{start_}\t{end_}\n' \
                                .format(chr_=map_pos[0], start_=map_pos[1], end_=map_pos[2])
                            self.map_cand_pos_array = np.vstack((self.map_cand_pos_array, np.array(
                                [map_cand, map_pos[0], int(map_pos[1]), int(map_pos[2])])))
                        else:
                            if int(cluster) not in self.map_cand_dict[map_cand]:
                                self.map_cand_dict[map_cand].append(int(cluster))
                                self.map_cand_dict[map_cand] = sorted(self.map_cand_dict[map_cand])

        map_cand_pos_file = self.dir_map + '/temp/map_cand_pos_array.tmp'

        np.savetxt(map_cand_pos_file, self.map_cand_pos_array, delimiter='\t', fmt='%s')

        with open(os.path.join(config.PRJ_NAME, config.DIR_TREE['results'], 'viz.R'), 'a') as viz_file:
            viz_file.write('\n\n#Rscript for visualization: detailed comparative plot of eccCandidates\n\n'
                           '{cmd_r_}\n\n'
                           'cmp.candidate.list <- c({candidates_})\n'
                           'ext.in <- "{ext_csv_}"\n'
                           'ext.out <-"{ext_png_}"\n'
                           'dir.candidate.map <- "{dir_mapcanddt_}"\n'
                           'dir.candidate.cmp <- "{dir_cmpcanddt_}"\n'
                           'path.base.clu.blast <- "{path_clu_blast_}"\n'
                           'blast.header <- c({blast_})\n'
                           'dict.candmap.clu <- vector(mode="list", length=length(cmp.candidate.list))\n'
                           'names(dict.candmap.clu) <- cmp.candidate.list\n'
                           .format(cmd_r_=eccDNA_Rcodes.rcomparative_plot(),
                                   candidates_=', '.join('"{}"'.format(str(x)) for x in self.map_cand_dict.keys()),
                                   dir_mapcanddt_=os.path.join(self.dir_map, 'candidates'),
                                   dir_cmpcanddt_=self.dir_cand,
                                   path_clu_blast_=os.path.join(self.dir_res, 'blast_CL'),
                                   ext_csv_='_coverage-summary_NORM.csv',
                                   ext_png_='_comparative-plot' + config.IMAGE_EXT[config.IMAGE_TYPE],
                                   blast_=', '.join('"{}"'.format(str(x)) for x in config.BLASTn_OUTFMT),
                                   count_cands_=len(self.map_cand_dict.keys())))

        for cand in self.map_cand_dict.keys():
            os.makedirs(os.path.join(self.dir_cand, cand), exist_ok=True)
            self.LOGGER.info('Creating comparative plot scripts: {}'.format(cand))

            # png_path = os.path.join(self.cand_out, cand,
            #                        cand + '_comparative-plot' + config.IMAGE_EXT[config.IMAGE_TYPE])
            cluster_list = self.map_cand_dict[cand]
            # cluster_blast_path = os.path.join(self.res_out, 'blast_CL')

            with open(os.path.join(config.PRJ_NAME, config.DIR_TREE['results'], 'viz.R'), 'a') as viz_file:
                viz_file.write('dict.candmap.clu${candidate_} <- c({clusters_})\n'
                               .format(candidate_=cand,
                                       clusters_=', '.join(map(str, cluster_list))))
                '''{cmd_r_}\n\n'
                                'candidate.list <- c({candidates_})\n'
                                'ext.in <- "{ext_csv_}"\n'
                                'ext.out <-"{ext_png_}"\n'
                                'dir.candidate.map <- "{dir_mapcanddt_}"\n'
                                
                                'dict.candmap.clu <- vector(mode="list", length=length(candidate.list)\n'
                                'names(dict.candmap.clu) <- candidate.list\n'


                                'dir.candidate.clu <- "{dir_clucanddt_}"\n'
                                'dir.candidate.cmp <- "{dir_cmpcanddt_}"\n'
                                'for (candidate in candidate.list)\n'
                                '{{\n'
                                'file.in.map <- paste0(dir.candidate.map, "/", candidate, "/", candidate, ext.in)\n'
                                'file.in.clu <- paste0(dir.candidate.map, "/", candidate, "/", candidate, ext.in)\n'

                                'file.out <- paste0(dir.candidate, "/", candidate, "/", candidate, ext.out)\n'
                                'comp_plot(candidate, file.in, file.out)\n}}'

                                'file.in <- "{file_in_}"\n'
                                'file.out <- "{file_out_}"\n'
                                'file.add1 <- c({file_add1_})\n'
                                'file.add2 <- "{file_add2_}"\n'
                                'file.add3 <- c({file_add3_})\n'
                                'cand.key <- "{key_}"\n'
                                '\n{cmd_r_}\n'
                                'comp_plot(file.in, file.out, file.add1, file.add2, file.add3, cand.key)') \
                               .format(cmd_r_=eccDNA_Rcodes.rcomparative_plot(),
                                       candidates_=', '.join('"{}"'.format(str(x)) for x in self.map_cand_dict.keys()),
                                       dir_mapcanddt_=os.path.join(self.map_in, 'candidates'),
                                       ext_csv_='_coverage-summary_NORM.csv',
                                       ext_png_='_comparative-plot' + config.IMAGE_EXT[config.IMAGE_TYPE],
                                       count_cands_=len(self.map_cand_dict.keys()),
                                       file_out_=png_path,
                                       file_add1_=', '.join(map(str, cluster_list)),
                                       file_add2_=cluster_blast_path,
                                       file_add3_=', '.join('"{}"'.format(str(x)) for x in config.BLASTn_OUTFMT),
                                       key_=cand, ))'''

            # save combined blast table in mapper candidate folder
            cl_blast_tab = np.array(config.BLASTn_OUTFMT)

            self.LOGGER.info('Creating comparative blast results: {}'.format(cand))
            for i in range(len(cluster_list)):
                cl_blast_tab_add = np.loadtxt(os.path.join(self.dir_res, 'blast_CL{}.m6'.format(cluster_list[i])),
                                              dtype=np.str, delimiter='\t', skiprows=0)
                cl_blast_tab = np.vstack((cl_blast_tab, cl_blast_tab_add))

            # remove hits from other mapping candidates (e.g. one cluster hits multiple mapping candidates if they are
            # from repetitive DNAs)
            indices_other_cl = []
            for i in range(len(cl_blast_tab[:, 0])):
                if not search(cand, cl_blast_tab[i, 1]):
                    indices_other_cl.append(i)

            cl_blast_tab = np.delete(cl_blast_tab, indices_other_cl, 0)

            cand_blast_path = os.path.join(self.dir_cand, cand, cand + '_comparative-blast.m6')

            np.savetxt(cand_blast_path, cl_blast_tab, delimiter='\t', fmt='%s')

        self.LOGGER.info('Reporting data for comparative summary plots.')
        os.system('echo "{bed_}" | '
                  'bedtools intersect -c -a {ref_win_} -b stdin | '
                  'sed -e \'1 i\\chr\\tstart\\tend\\tcomp_cand\' > {out_}'
                  .format(ref_win_=self.ref_seq_win,
                          bed_=self.map_fas_pos_bed,
                          out_=self.dir_map + '/temp/report_comp_cand_win.tmp'))

        hicnfdc_path = self.dir_map + '/temp/report_hiconf_win.tmp'
        comp_cand_path = self.dir_map + '/temp/report_comp_cand_win.tmp'
        png_path = os.path.join(self.dir_res, 'comparative-summary_manhattan' + config.IMAGE_EXT[config.IMAGE_TYPE])

        self.LOGGER.info('Creating Rscripts for visualization of comparative summary plots.')

        with open(os.path.join(config.PRJ_NAME, config.DIR_TREE['results'], 'viz.R'), 'a') as viz_file:
            viz_file.write('\nfor (candidate in cmp.candidate.list)\n'
                           '{{\n'
                           'file.in <- paste0(dir.candidate.map, "/", candidate, "/", candidate, ext.in)\n'
                           'file.out <- paste0(dir.candidate.cmp, "/", candidate, "/", candidate, ext.out)\n'
                           'comp_plot(file.in, file.out, dict.candmap.clu[[candidate]], path.base.clu.blast, '
                           'blast.header, candidate)\n}}'
                           .format())

        with open(os.path.join(config.PRJ_NAME, config.DIR_TREE['results'], 'viz.R'), 'a') as viz_file:
            viz_file.write(('\n\n#Rscript for visualization: comparative Manhattan plot\n\n'
                            'file.in <- "{file_in_}"\n'
                            'file.out <- "{file_out_}"\n'
                            'file.add1 <- "{file_add1_}"\n'
                            'file.add2 <- "{file_add2_}"\n'
                            'file.add3 <- "{file_add3_}"\n'
                            'win.size <- {win_}\n'
                            '\n{cmd_r_}\n'
                            'comp_manhattan(file.in, file.out, file.add1, file.add2, file.add3, win.size)')
                           .format(file_in_=self.map_sum_cov_file,
                                   file_out_=png_path,
                                   file_add1_=hicnfdc_path,
                                   file_add2_=comp_cand_path,
                                   file_add3_=map_cand_pos_file,
                                   win_=self.win_size,
                                   cmd_r_=eccDNA_Rcodes.rmanhattan_comp_plot()))

        png_path = os.path.join(self.dir_res, 'comparative-summary_scatter' + config.IMAGE_EXT[config.IMAGE_TYPE])

        with open(os.path.join(config.PRJ_NAME, config.DIR_TREE['results'], 'viz.R'), 'a') as viz_file:
            viz_file.write(('\n\n#Rscript for visualization: comparative scatter plot\n\n'
                            'file.in <- "{file_in_}"\n'
                            'file.out <- "{file_out_}"\n'
                            'file.add1 <- "{file_add1_}"\n'
                            'file.add2 <- "{file_add2_}"\n'
                            'file.add3 <- "{file_add3_}"\n'
                            'cl.list <- c({cl_})\n'
                            '\n{cmd_r_}\n'
                            'comp_scatter(file.in, file.out, file.add1, file.add2, file.add3, cl.list)') \
                           .format(file_in_=self.cl_sum_file,
                                   file_out_=png_path,
                                   file_add1_=config.REPEX_ECC_PROPORTION,
                                   file_add2_=self.prefix[0],
                                   file_add3_=self.prefix[1],
                                   cl_=', '.join(self.cand_cl_list),
                                   cmd_r_=eccDNA_Rcodes.rscatter_comp_plot()))

        self.LOGGER.info('Comparative scatter plot done.')

        # create html summary
        self.LOGGER.info('Creating eccComp summary HTML.')
        with open(self.dir_html + '/eccComp_summary.html', 'wt') as html_file:
            html_file.write(eccHTML_templates.eccHTML_style + '\n')

        cl_cand_tab = np.loadtxt(os.path.join(self.dir_clu, 'COMPARATIVE_CLUSTER_TABLE_eccCANDIDATES.csv'),
                                 dtype=np.str, delimiter='\t', skiprows=0)

        manhattan_rel = os.path.join('comparative_results',
                                     'comparative-summary_manhattan' + config.IMAGE_EXT[config.IMAGE_TYPE])
        scatter_rel = os.path.join('comparative_results',
                                   'comparative-summary_scatter' + config.IMAGE_EXT[config.IMAGE_TYPE])

        cnt_map_reads = [0, 0]

        for file in os.listdir(self.dir_map):
            if file.endswith('_alignment-stats.txt'):
                if search(self.prefix[0], file):
                    with open(os.path.join(self.dir_map, file))as stats:
                        for line in stats:
                            if search("reads mapped:", line):
                                cnt_map_reads[0] = int(line.strip().split('\t')[-1])
                if search(self.prefix[1], file):
                    with open(os.path.join(self.dir_map, file))as stats:
                        for line in stats:
                            if search("reads mapped:", line):
                                cnt_map_reads[1] = int(line.strip().split('\t')[-1])

        cnt_cl_reads = [sum(cl_cand_tab[1:, 7].astype(int)), sum(cl_cand_tab[1:, 8].astype(int))]

        with open(self.dir_html + '/eccComp_summary.html', 'a') as html_file:
            html_file.write(eccHTML_templates.eccComp_HTML_sum.format(manhattan_=manhattan_rel,
                                                                      scatter_=scatter_rel,
                                                                      pre1_=self.prefix[0],
                                                                      readsMP1_=cnt_map_reads[0],
                                                                      readsCL1_=cnt_cl_reads[0],
                                                                      pre2_=self.prefix[1],
                                                                      readsMP2_=cnt_map_reads[1],
                                                                      readsCL2_=cnt_cl_reads[1]) + '\n')
            html_file.write(eccHTML_templates.eccComp_HTML_th.format(pre1_=self.prefix[0]) + '\n')

        self.map_cand_dict = OrderedDict(sorted(self.map_cand_dict.items()))

        for cand in self.map_cand_dict.keys():

            for file in os.listdir(self.dir_map):
                if file.endswith('_summary_region-coverages_normalized.csv'):
                    with open(os.path.join(self.dir_map, file)) as tab:
                        for line in tab:
                            if search(cand, line):
                                line = line.strip().split('\t')
                                length = int(line[2]) - int(line[1])
                                enrich = line[9]

            cand_fasta = os.path.join(self.dir_map, 'candidates', cand, cand + '.fasta')
            blast_map = ''

            with open(os.path.join(self.dir_map, 'candidates', cand, cand + '_blast.m6')) as blast_file:
                score = 0
                for line in blast_file:
                    line = line.strip().split('\t')
                    new_score = int(float(line[11]))
                    if new_score > score:
                        score = new_score
                        blast_file_rel = os.path.join(self.dir_map, 'candidates', cand, cand + '_blast.m6')
                        blast_map = '<a href="{}">{}</a>'.format(blast_file_rel, str(line[1]))

            comp_plot_rel = os.path.join('comparative_results', 'candidates', cand,
                                         cand + '_comparative-plot' + config.IMAGE_EXT[config.IMAGE_TYPE])
            cl_block = ''''''
            prop_block = ''''''
            blast_cl = ''''''
            for cluster in self.map_cand_dict[cand]:
                cl_dir = 'dir_CL' + str(cluster).zfill(4)
                cl_indexhtml = os.path.join('clustering_results', 'seqclust/clustering/clusters', cl_dir, 'index.html')
                cl_block += '<a href="{sclidx_}">{scl_}</a><br>'.format(scl_=cluster,
                                                                        sclidx_=cl_indexhtml)
                with open(os.path.join(self.dir_clu, 'COMPARATIVE_CLUSTER_TABLE_eccCANDIDATES.csv')) as tab:
                    for line in tab:
                        line = line.strip().split('\t')
                        if str(line[0]) == str(cluster):
                            prop_block += str(line[9]) + '<br>'
                            blast_cl += str(line[4]).replace('"', '') + '<br>'

            with open(self.dir_html + '/eccComp_summary.html', 'a') as html_file:
                html_file.write(eccHTML_templates.eccComp_HTML_td.format(cand_=cand,
                                                                         candidx_=cand_fasta,
                                                                         comppl_=comp_plot_rel,
                                                                         enr_=enrich,
                                                                         cl_block_=cl_block,
                                                                         prop_block_=prop_block,
                                                                         blast_map_=blast_map,
                                                                         blast_cl_=blast_cl) + '\n')

    def __init__(self, map_sum_cov_file, mapper_fasta_files, prefix, ref_seq_file, window_size, log, conn):
        self.map_fas_files = mapper_fasta_files
        self.map_fas_name = os.path.basename(self.map_fas_files)
        self.prefix = prefix
        self.ref_seq_name = os.path.splitext(os.path.basename(ref_seq_file))[0]
        self.map_sum_cov_file = map_sum_cov_file
        self.win_size = window_size
        self.LOGGER = log
        self.conn = conn

        # general directory
        self.dir_abs = os.path.abspath(os.path.join(config.CURRENT_DIR, config.PRJ_NAME))

        # input directories
        self.dir_map = os.path.join(self.dir_abs, config.DIR_TREE['mapping_results'])
        self.dir_clu = os.path.join(self.dir_abs, config.DIR_TREE['clustering_results'])
        self.dir_ref = os.path.join(self.dir_abs, config.DIR_TREE['reference_data'])

        # output directory
        self.dir_res = os.path.join(self.dir_abs, config.DIR_TREE['comparative_results'])
        self.dir_cand = os.path.join(self.dir_abs, config.DIR_TREE['comparative_candidates'])
        self.dir_html = os.path.join(self.dir_abs, config.DIR_TREE['results'])

        # internal variables
        self.blastdb_fas = os.path.join(self.dir_res, 'eccMapper_cand.fa')

        self.cl_sum_file = os.path.join(self.dir_clu, 'comparative_cluster_table.csv')
        self.cl_cand_list = os.path.join(self.dir_clu, 'comp_cl_tab_eccCandidates_list.csv')
        self.cl_list = np.loadtxt(self.cl_cand_list, dtype=np.str, delimiter='\t', skiprows=1)
        self.cand_cl_list = []

        self.map_cand_dict = {}
        self.map_cand_pos_array = np.array(['cand', 'chr', 'start', 'end'])
        self.map_fas_pos_bed = ''''''

        # reference data
        self.ref_seq_win = self.dir_ref + '/' + self.ref_seq_name + '_win{}'.format(window_size) + '.bed'
