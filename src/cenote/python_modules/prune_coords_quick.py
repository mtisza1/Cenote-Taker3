#!/usr/bin/env python

import os
import bisect
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def prune_chunks(name, group, out_dir1, hallmark_arg):
    id_list = ['common_virus','hypothetical_protein','nonviral_gene','intergenic']
    score_list = [10,5,-3,0]
    domain_dictionary = dict(zip(id_list, score_list))

    window = 5000
    slide = 50

    contig_length1 = group['contig_length'].agg(pd.Series.mode)
    total_len = int(contig_length1.iloc[0])

    vscore_list = [0] * total_len

    group2 = group.copy()
    group2['gene_start'] = group2['gene_start'].astype(int)
    group2['gene_stop'] = group2['gene_stop'].astype(int)
    group2 = group2.sort_values('gene_start')

    for _, row in group2.iterrows():
        s = domain_dictionary[row['vscore_category']]
        a = int(row['gene_start'])
        b = int(row['gene_stop'])
        if a < 0: a = 0
        if b > total_len: b = total_len
        if a < b:
            vscore_list[a:b] = [s] * (b - a)

    n_blocks = int(((len(vscore_list) - window) / slide) + 1)
    n_windows = max(1, n_blocks + 2)

    dat_list = []
    count = 0
    count_start = -slide
    for i in range(0, n_windows * slide, slide):
        score_result = int(np.sum(vscore_list[i:i+window]))
        new_let_list = vscore_list[i:i+window]
        pf = "pass" if score_result >= 0 else "fail"
        vg = new_let_list.count(10)
        hg = new_let_list.count(5)
        bg = new_let_list.count(-3)
        ig = new_let_list.count(0)
        count += 1
        count_start += slide
        count_stop = count_start + window
        dat_list.append([count, max(0, count_start), min(total_len, count_stop), pf, score_result, vg, hg, bg, ig])

    df_0 = pd.DataFrame(dat_list, columns=[
        'Window','Position start','Position stop','Pass/Fail','Score',
        'VirusGene','HypotheticalGene','BacterialGene','Intergenic'
    ])

    x = df_0['Window']
    y = np.array(df_0['Score'])
    l = len(y)

    def smooth(arr, box_pts):
        box = np.ones(box_pts)/box_pts
        return np.convolve(arr, box, mode='same')

    smooth_val = 100
    if l <= smooth_val:
        smooth_val = max(1, int(0.5 * l))
    else:
        smooth_val = int(smooth_val)

    smoothy = smooth(y, smooth_val) if l > 0 else np.array([])

    df_0_file = os.path.join(out_dir1, name + ".window_info.tsv")
    df_0.to_csv(df_0_file, sep='\t', index=False)

    df_0['Window midpoint'] = ((df_0['Position start'] + df_0['Position stop'])/2).astype(int)

    pos_mask = (smoothy > 0) if l > 0 else np.array([])

    boundaries = []
    if l > 0:
        for i in range(l-1):
            if pos_mask[i] != pos_mask[i+1]:
                boundaries.append(i)

    merged_records = []
    for i in boundaries:
        w = int(df_0.iloc[i]['Window'])
        ps = int(df_0.iloc[i]['Position start'])
        pe = int(df_0.iloc[i]['Position stop'])
        mid = int(df_0.iloc[i]['Window midpoint'])
        sign_right = '+' if (i+1 < l and pos_mask[i+1]) else '-'
        merged_records.append([w, sign_right, ps, pe, 'none', mid])

    if l > 0:
        w_last = int(df_0.iloc[-1]['Window'])
        ps_last = int(df_0.iloc[-1]['Position start'])
        mid_last = int(df_0.iloc[-1]['Window midpoint'])
        sign_right_last = '+' if (pos_mask[-1] if l > 0 else False) else '-'
        merged_records.append([w_last, sign_right_last, ps_last, total_len + 1, 'none', mid_last])

    merged_df = pd.DataFrame(merged_records, columns=['Window','+/- to the right','Position start','Position stop','Chunk_end','Window midpoint'])

    merged_file = os.path.join(out_dir1, name + ".merged_windows.tsv")
    merged_df.to_csv(merged_file, sep='\t', index=False)

    gstart = sorted(list(group2['gene_start']))
    gstop = sorted(list(group2['gene_stop']))

    def right_cutoff(position):
        ce = int(position)
        if len(gstop) == 0:
            return ce
        rpos = bisect.bisect(gstop, ce)
        if rpos == 0:
            return gstop[0]
        return gstop[rpos-1]

    def left_cutoff(position):
        cs = int(position)
        if len(gstart) == 0:
            return cs
        lpos = bisect.bisect(gstart, cs)
        if lpos >= len(gstart):
            lpos = len(gstart) - 1
        return gstart[lpos]

    chunks = []
    i = 0
    while i < l:
        if pos_mask[i]:
            s = i
            while i + 1 < l and pos_mask[i+1]:
                i += 1
            e = i
            left_mid = int(df_0.iloc[s]['Window midpoint'])
            if e + 1 < l:
                right_mid = int(df_0.iloc[e+1]['Window midpoint'])
                leco = left_cutoff(left_mid)
                rico = right_cutoff(right_mid)
                if rico > leco:
                    left = leco
                    right = rico
                else:
                    left = int(df_0.iloc[s]['Position start'])
                    right = int(df_0.iloc[e]['Position stop'])
            else:
                leco = left_cutoff(left_mid)
                left = leco
                right = total_len + 1
            chunks.append((left, right))
        i += 1

    ddf_list = []
    for idx, (lc, rc) in enumerate(chunks):
        ddf_list.append([f"C{idx}", int(lc), int(rc)])

    chunk_df = pd.DataFrame(ddf_list, columns=["chunk_number", "left_cutoff", "right_cutoff"])
    chunk_df['contig'] = name
    chunk_df = chunk_df[['contig', 'left_cutoff', 'right_cutoff', 'chunk_number']]

    chunk_sum_file = os.path.join(out_dir1, name + ".chunks.tsv")
    chunk_df.to_csv(chunk_sum_file, sep='\t', index=False)

    if "virion" in str(hallmark_arg):
        virion_str = "Evidence_source == 'hallmark_hmm'"
    else:
        virion_str = ""
    if "rdrp" in str(hallmark_arg):
        rdrp_str = "Evidence_source == 'rdrp_hall_hmm'"
    else:
        rdrp_str = ""
    if "dnarep" in str(hallmark_arg):
        rep_str = "Evidence_source == 'rep_hall_hmm'"
    else:
        rep_str = ""

    query_str = ' | '.join(filter(None, [virion_str, rdrp_str, rep_str]))

    if query_str:
        vir_bait_table = group2[['gene_start', 'gene_stop', 'Evidence_source']].query(str(query_str)).copy()
    else:
        vir_bait_table = group2[['gene_start', 'gene_stop', 'Evidence_source']].iloc[0:0].copy()

    vir_bait_table['gene_start'] = vir_bait_table['gene_start'].astype(int)
    vir_bait_table['gene_stop'] = vir_bait_table['gene_stop'].astype(int)
    vir_bait_table['mean'] = ((vir_bait_table['gene_start'] + vir_bait_table['gene_stop'])/2).round()
    vir_bait_table_med_list = list(vir_bait_table['mean'].astype(int))

    points_list = []
    for item in vir_bait_table_med_list:
        eq = round(((item - 2500) + 50) / 50)
        if eq >= len(x):
            plot_point = (len(x) - 1)
        else:
            plot_point = eq
        points_list.append(plot_point)

    new_points_list = [1 if i <= 0 else i for i in points_list]

    df_0['smoothy'] = smooth(df_0['Score'].values, 100) if l > 0 else np.array([])

    pdf_outname = os.path.join(out_dir1, name + ".figures.pdf")
    figures = PdfPages(pdf_outname)

    plt.plot(x, y, 'o', ms=0.6)
    plt.axhline(0, 0, l)
    plt.plot(x, df_0['smoothy'], 'c', lw=2)
    plt.plot(x, df_0['smoothy'], 'y', markevery=(new_points_list), ms=11.0, marker='*')
    plt.title("Viral region calls")
    plt.xlabel('Window')
    plt.ylabel('Score')
    plt.rc('axes', titlesize=6.8)
    plt.rc('xtick', labelsize=5)
    plt.rc('ytick', labelsize=5)
    plt.rc('legend', fontsize=5)
    plt.rc('figure', titlesize=8)
    plt.grid(True)
    if l > 0:
        zero = np.zeros(l)
        idx = np.argwhere(np.diff(np.sign(zero - smooth(y,100)))).flatten()
        plt.plot(x.iloc[idx], zero[idx], 'ro', ms=5.0)
    plt.savefig(figures, format='pdf')
    plt.close()

    mycol = (["#e7ba52", "#637939", "#7b4173", "#d6616b"])
    df_0[['VirusGene','HypotheticalGene','BacterialGene','Intergenic']].plot(color=mycol)
    plt.grid(True)
    plt.xlabel('Window')
    plt.ylabel('Count')
    plt.title('Character counts')
    plt.rc('axes', titlesize=6.8)
    plt.rc('xtick', labelsize=5)
    plt.rc('ytick', labelsize=5)
    plt.rc('legend', fontsize=5)
    plt.rc('figure', titlesize=8)
    plt.savefig(figures, format='pdf')
    plt.close()

    figures.close()
