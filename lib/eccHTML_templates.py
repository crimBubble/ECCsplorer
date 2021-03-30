#!/usr/bin/python3

import os
from lib import config

eccHTML_style = ''''''

with open(os.path.join(config.PIPE_DIR, 'lib/html_style.css')) as style:
    for line in style:
        eccHTML_style += line

eccMap_HTML_sum = '''<h2> eccMap - summary </h2>
    <p align= left >
    <a href="{manhattan_}">
    <img src="{manhattan_}" alt="summary manhattan plot" height=450px>
    </a>
    <br>
    <br>
    <table border=0 class=dataframe>
    <tr class= firstline > 
        <th>data set</th>
        <th>mapped reads</th>
        <th>mapped bases</th>
    </tr>
    <tr>
    <td>{pre1_}</td>
    <td>{reads1_}</td>
    <td>{bases1_}</td>
    </tr>
    <tr>
    <td>{pre2_}</td>
    <td>{reads2_}</td>
    <td>{bases2_}</td>
    </tr>
    </table><br>
    </p>'''

eccMap_HTML_th = '''<h3> eccMap - candidates overview </h3>
    
    <p align= left >
    <table border=0 class=dataframe>
    <tbody> 
    <tr class= firstline > 
        <th>candidate  </th>
        <th>mapping plots  </th>
        <th>enrichment score</th>
        <th>length</th>
        <th>position</th>
        <th>blast</th>
        <th>estimated sequence</th>
    </tr>'''

eccMap_HTML_td = '''<tr>
    <td>{cand_}</td>
    <td><a href="{plot_}">
    <img src="{plot_}" alt="{cand_}" height=300px></a></td>
    <td>{enr_}</td>
    <td>{len_}</td>
    <td>{pos_}</td>
    <td>{blast_}</td>
    <td class=seqcell><div style="height:300px; overflow:auto;">
    {seq_}
    </div></td>
    </tr>'''

eccCL_HTML_sum = '''<h2> eccCL - summary </h2>
    <p align= left >
    <a href="{scatter_}">
    <img src="{scatter_}" alt="summary scatter plot" height=450px>
    </a>
    <br>
    <br>
    <table border=0 class=dataframe>
    <tr class= firstline > 
        <th>data set</th>
        <th>clustered reads</th>
    </tr>
    <tr>
    <td>{pre1_}</td>
    <td>{reads1_}</td>
    </tr>
    <tr>
    <td>{pre2_}</td>
    <td>{reads2_}</td>
    </tr>
    </table><br>
    </p>'''

eccCL_HTML_th = '''<h3> eccCL - cluster candidates overview </h3>
    <p align= left >
    <table border=0 class=dataframe>
    <tbody> 
    <tr class= firstline > 
        <th>cluster</th>
        <th>super cluster</th>
        <th>cluster graph</th>
        <th>comparative graph</th>
        <th>size</th>
        <th>{pre1_} proportion</th>
        <th>best hit annotation</th>
        <th>contigs</th>
    </tr>'''

eccCL_HTML_td = '''<tr>
    <td><a href="{clidx_}">{cl_}</a></td>
    <td><a href="{sclidx_}">{scl_}</a></td>
    <td><a href="{clgr_}"><img src="{clgr_}" alt="{cl_}" height=100px></a></td>
    <td><a href="{compgr_}"><img src="{compgr_}" alt="{cl_}" height=100px></a></td>
    <td>{size_}</td>
    <td>{prop_}</td>
    <td>{blast_}</td>
    <td class=seqcell><div style="height:100px; overflow:auto;">{seq_}</div></td>
    </tr>'''

eccComp_HTML_sum = '''<h2> eccComp - summary </h2>
    <p align= left >
    <a href="{manhattan_}"><img src="{manhattan_}" alt="summary manhattan plot" height=450px></a>
    <a href="{scatter_}"><img src="{scatter_}" alt="summary scatter plot" height=450px></a>
    <br>
    <br>
    <table border=0 class=dataframe>
    <tr class= firstline > 
        <th>data set</th>
        <th>mapped reads</th>
        <th>clustered reads</th>
    </tr>
    <tr>
    <td>{pre1_}</td>
    <td>{readsMP1_}</td>
    <td>{readsCL1_}</td>
    </tr>
    <tr>
    <td>{pre2_}</td>
    <td>{readsMP2_}</td>
    <td>{readsCL2_}</td>
    </tr>
    </table><br>
    </p>'''

eccComp_HTML_th = '''<h3> eccCL - cluster candidates overview </h3>
    <p align= left >
    <table border=0 class=dataframe>
    <tbody> 
    <tr class= firstline > 
        <th>comparative plot</th>
        <th>candidate</th>
        <th>enrichment score</th>
        <th>best hit (mapping)</th>
        <th>associated clusters</th>
        <th>{pre1_} proportion</th>
        <th>best hit (clusters)</th>
    </tr>'''

eccComp_HTML_td = '''<tr>
    <td><a href="{comppl_}"><img src="{comppl_}" alt="{comppl_}" height=300px></a></td>
    <td><a href="{candidx_}">{cand_}</a></td>
    <td>{enr_}</td>
    <td>{blast_map_}</td>
    <td>{cl_block_}</td>
    <td>{prop_block_}</td>
    <td>{blast_cl_}</td>
    </tr>'''