#!/usr/bin/env python3

import numpy as np
import os
import subprocess
import shlex

from Bio import SeqIO
from lib import config
from lib import eccHTML_templates
from lib import eccDNA_Rcodes


class eccClusterer:
    """
    Cluster reads to detect highly abundant features with repetitive character.
    Tandemly arranged repetitive structures occur due to circle amplification
    """

    def cluster_coordinator(self):

        # run RepeatExplorer
        if not os.path.isfile(self.dir_res + '/index.html') or os.path.getsize(self.dir_res + '/index.html') == 0:
            # TODO: catch output in log
            self.LOGGER.info('Starting read clustering.')

            cmd = config.REPEX_CMD.format(path_=config.REPEATEXPLORER_PATH,
                                          preflen_=len(self.prefix[0]),
                                          out_=self.dir_res,
                                          tax_=config.REPEX_TAX[self.taxon],
                                          cpu_=config.CPU,
                                          in_=self.seq_file)

            create_idx = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, universal_newlines=True)

            out = create_idx.communicate()[0]  # .lstrip('\n').rstrip('\n')

            self.LOGGER.info('{}'.format(out))

        else:
            self.LOGGER.info('Clustering results found. Skipping clustering.')

        self.LOGGER.info('Summarizing clustering results.')

        # reformat cells with multiple lines from cluster output (single line per row)
        unwrap_cmd = '''sed 's/""*$/ENDE/' {in_} | sed -e ':a' -e 'N' -e '$!ba' -e 's/\\n/ /g' | 
                        sed 's/\\tENDE/\\t""\\n/g' | sed 's/ ENDE/"\\n/g' | sed 's/ENDE/"\\n/g' | 
                        sed -e 's/^[ \\t]*//' > {out_}'''

        os.system(unwrap_cmd.format(in_=os.path.join(self.dir_res, 'CLUSTER_TABLE.csv'),
                                    out_=os.path.join(self.dir_res, 'cluster_table_sngl-lined.csv')))

        # read reformatted cluster output (w/o header) and get number of analyzed clusters
        cl_tab = np.loadtxt(os.path.join(self.dir_res, 'cluster_table_sngl-lined.csv'),
                            dtype=np.str, delimiter='\t', skiprows=1)
        cl_cnt = max(cl_tab[:, 0].astype(int))

        # read comparative cluster output
        comp_cl_tab = np.loadtxt(os.path.join(self.dir_res, 'COMPARATIVE_ANALYSIS_COUNTS.csv'),
                                 dtype=np.str, delimiter='\t', skiprows=0)

        # add comparative cluster output columns to reformatted cluster output in defined order
        def add_cnt_col(cl_tab, i):
            col_cnts = comp_cl_tab[1:cl_cnt + 1, i]
            cl_tab = np.hstack((cl_tab, col_cnts.reshape(np.size(col_cnts), 1)))
            return cl_tab

        for i in range(len(comp_cl_tab[0, :])):
            if comp_cl_tab[0, i].replace('"', '') == self.prefix[0]:
                cl_tab = add_cnt_col(cl_tab, i)
        for i in range(len(comp_cl_tab[0, :])):
            if comp_cl_tab[0, i].replace('"', '') == self.prefix[1]:
                cl_tab = add_cnt_col(cl_tab, i)

        for row in cl_tab:
            reads_TR = row[7].astype(float)
            reads_CO = row[8].astype(float)
            cl_score = (reads_TR / (reads_TR + reads_CO)).round(3)
            row = np.hstack((row, cl_score))
            self.cl_sum_tab = np.vstack((self.cl_sum_tab, row))
            if cl_score > config.REPEX_ECC_PROPORTION:
                self.cl_cand_tab = np.vstack((self.cl_cand_tab, row))

        np.savetxt(self.cl_cand_file, self.cl_cand_tab, delimiter='\t', fmt='%s')
        np.savetxt(self.cl_cand_list_file, self.cl_cand_tab[:, 0], delimiter='\t', fmt='%s')
        np.savetxt(self.cl_sum_file, self.cl_sum_tab, delimiter='\t', fmt='%s')

        png_path = os.path.join(self.dir_res, 'summary_scatter' + config.IMAGE_EXT[config.IMAGE_TYPE])

        with open(os.path.join(config.PRJ_NAME, config.DIR_TREE['results'], 'viz.R'), 'a') as viz_file:
            viz_file.write('''\n\n#Rscript for visualization: clustering scatter plot
                \nfile.in <- "{file_in_}"\nfile.out <- "{png_}"\necc.prop <- {prop_}\npre.tr <- "{pre1_}"\npre.co <- "{pre2_}"
                \n{cmd_r_}\ncl_scatter(file.in, file.out, ecc.prop, pre.tr, pre.co)
                '''.format(file_in_=self.cl_sum_file,
                           png_=png_path,
                           prop_=config.REPEX_ECC_PROPORTION,
                           pre1_=self.prefix[0],
                           pre2_=self.prefix[1],
                           cmd_r_=eccDNA_Rcodes.rscatter_plot()))

        # create html summary
        print('Creating eccCL summary HTML.')
        with open(self.dir_html + '/eccCL_summary.html', 'wt') as html_file:
            html_file.write(eccHTML_templates.eccHTML_style + '\n')

        # non tab values
        scatter_rel = os.path.join('clustering_results', 'summary_scatter' + config.IMAGE_EXT[config.IMAGE_TYPE])
        cnt_cl_reads = [sum(self.cl_cand_tab[1:, 7].astype(int)), sum(self.cl_cand_tab[1:, 8].astype(int))]

        with open(self.dir_html + '/eccCL_summary.html', 'a') as html_file:
            html_file.write(eccHTML_templates.eccCL_HTML_sum.format(scatter_=scatter_rel,
                                                                    pre1_=self.prefix[0],
                                                                    reads1_=cnt_cl_reads[0],
                                                                    pre2_=self.prefix[1],
                                                                    reads2_=cnt_cl_reads[1]) + '\n')
            html_file.write(eccHTML_templates.eccCL_HTML_th.format(pre1_=self.prefix[0]) + '\n')

        for i in range(1, np.size(self.cl_cand_tab, 0)):

            cl_dir = 'dir_CL' + str(self.cl_cand_tab[i, 0]).zfill(4)
            cl_indexhtml = os.path.join('clustering_results', 'seqclust/clustering/clusters', cl_dir, 'index.html')
            scl_dir = 'dir_SC' + str(self.cl_cand_tab[i, 1]).zfill(4)
            scl_indexhtml = os.path.join('clustering_results', 'seqclust/clustering/superclusters', scl_dir,
                                         'index.html')
            cl_graph = 'CL' + str(self.cl_cand_tab[i, 0] + '.png')
            cl_graph_rel = os.path.join('clustering_results', 'seqclust/clustering/clusters', cl_dir, cl_graph)

            if os.path.isfile(os.path.join(self.dir_html, cl_graph_rel)):
                pass

            else:
                cl_graph_rel = os.path.join('clustering_results', 'seqclust/clustering/clusters', cl_dir,
                                            'graph_layout.png')

            comp_graph_rel = os.path.join('clustering_results', 'seqclust/clustering/clusters', cl_dir, 'html_files',
                                          'graph_comparative.png')
            contig_file = 'reads.fas.CL' + str(self.cl_cand_tab[i, 0]) + '.contigs'
            contig_path = os.path.join(self.dir_res, 'seqclust/clustering/clusters', cl_dir, contig_file)
            contigs = ''''''

            if os.path.isfile(contig_path):
                with open(contig_path) as f:
                    for records in SeqIO.parse(f, 'fasta'):
                        contigs += '>' + records.id + '<br>'
                        contigs += records.seq + '<br>'

            with open(self.dir_html + '/eccCL_summary.html', 'a') as html_file:
                html_file.write(eccHTML_templates.eccCL_HTML_td.format(cl_=self.cl_cand_tab[i, 0],
                                                                       clidx_=cl_indexhtml,
                                                                       scl_=self.cl_cand_tab[i, 1],
                                                                       sclidx_=scl_indexhtml,
                                                                       clgr_=cl_graph_rel,
                                                                       compgr_=comp_graph_rel,
                                                                       size_=self.cl_cand_tab[i, 3],
                                                                       prop_=self.cl_cand_tab[i, 9],
                                                                       blast_=self.cl_cand_tab[i, 4],
                                                                       seq_=contigs) + '\n')

    def __init__(self, seq_file, prefix, taxon, conn, log):

        self.seq_file = seq_file
        self.taxon = taxon
        self.prefix = prefix
        self.conn = conn
        self.LOGGER = log

        # internal temporary variables
        self.cl_cand_tab = np.array([['cluster', 'super_cluster', 'size', 'size_adjusted',
                                      'super_CL_best-hit', 'TAREAN_classification', 'CL_similarity_hits',
                                      'size_{}'.format(self.prefix[0]), 'size_{}'.format(self.prefix[1]),
                                      'proportion_{}'.format(self.prefix[0])]])
        self.cl_sum_tab = self.cl_cand_tab

        # output directory
        self.dir_abs = os.path.abspath(os.path.join(config.CURRENT_DIR, config.PRJ_NAME))
        self.dir_res = os.path.join(self.dir_abs, config.DIR_TREE['clustering_results'])
        self.dir_html = os.path.join(self.dir_abs, config.DIR_TREE['results'])

        # output files
        self.cl_cand_file = os.path.join(self.dir_res, 'COMPARATIVE_CLUSTER_TABLE_eccCANDIDATES.csv')
        self.cl_sum_file = os.path.join(self.dir_res, 'comparative_cluster_table.csv')
        self.cl_cand_list_file = os.path.join(self.dir_res, 'comp_cl_tab_eccCandidates_list.csv')

        os.system('export PYTHONHASHSEED=12')  # set seed for pseudo random clustering (repeatability)
