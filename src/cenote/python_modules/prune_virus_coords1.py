#!/usr/bin/env python

import os
import pandas as pd
import numpy as np
import itertools
from itertools import tee
import statistics
from statistics import mean
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import math
import collections
import bisect

def prune_chunks(name, group, out_dir1):

    #print(name + "in prune")

    ## function takes the name and data from the grouped pandas dataframe
    ## grouped_df
    ## and the output directory

    ####define scoring scheme and important variables
    id_list = ['common_virus','hypothetical_protein','nonviral_gene','intergenic']
    score_list = [10,5,-3,0]

    keys = id_list
    values = score_list
    domain_dictionary = dict(zip(keys,values))

    threshold = 0
    window = 5000
    sliiiiiide_to_the_right = 50

    count=0
    count_start=-sliiiiiide_to_the_right

    contig_length1 = group['contig_length'].agg(pd.Series.mode)
    ####

    ## make the score list based on annotations
    vscore_list = list([0] * int(contig_length1.iloc[0]))

    total_len = int(len(vscore_list))

    for index, row in group.iterrows():
        S = domain_dictionary[row['vscore_category']]
        vscore_list[row['gene_start']:row['gene_stop']] = [S] * (row['gene_stop'] - row['gene_start'])

    blocks = int(((len(vscore_list) - window) / sliiiiiide_to_the_right) + 1)
    blocks_2 = blocks + 2 #need this otherwise the last (incomplete/little) block will be cut off!

    dat_list = []
    for i in range(0, blocks_2 * sliiiiiide_to_the_right, sliiiiiide_to_the_right):
        score_result = sum(vscore_list[i:i+window])
        new_let_list = vscore_list[i:i+window]

        if score_result >= 0 :
            PF_result = "pass"
        else:
            PF_result = "fail"

        #counts for later
        VirusGene = new_let_list.count(10)
        HypotheticalGene = new_let_list.count(5)
        BacterialGene = new_let_list.count(-3)
        Intergenic = new_let_list.count(0)

        #vars for count columns
        count = count +1
        count_start += sliiiiiide_to_the_right #same as c_s = c_s + siiii...
        count_stop = count_start+window

        dat_list.append([count, count_start, count_stop, PF_result, score_result, VirusGene,
                          HypotheticalGene, BacterialGene, Intergenic])


    df_0 = pd.DataFrame(dat_list, columns=['Window', 'Position start', 'Position stop',
                                             'Pass/Fail', 'Score', 'VirusGene', 'HypotheticalGene', 
                                             'Intergenic', 'BacterialGene'])

    #Now let's make the smoothed plot

    x = df_0['Window']
    y = np.array(df_0['Score'])
    l = df_0['Window'].count()
    df_empty = pd.DataFrame(index=range(l),columns=range(1))
    for col in df_empty.columns:
        df_empty[col].values[:] = 0

    zero=df_empty[0]

    def smooth(y, box_pts):
        box = np.ones(box_pts)/box_pts
        y_smooth = np.convolve(y, box, mode='same')
        return y_smooth
    smooth_val = 100 #####we can change this if we want!


    #statement for handling short sequences (error called if len(y) < smoothing value)
    if len(y) <= smooth_val:
        smooth_val = int(0.5 * len(y))
    else:
        smooth_val = int(smooth_val)

    smoth = smooth(y,smooth_val)
    idx = np.argwhere(np.diff(np.sign(zero - smoth))).flatten()
    df_val = pd.DataFrame(zero[idx])
    df_val = df_val.reset_index()
    #we will save to figures, but first we need to do the validation steps

    #This is for validating if region is + or -
    df_val.loc[-1] = 1  # adding a row for first position
    df_val.index = df_val.index + 1  # shifting index
    df_val = df_val.sort_index()

    #last position as last row
    df_val.sort_values(by=['index']) #need to sort first otherwise +1 belwo will break things
    new_list = pd.DataFrame(df_val['index'] + 1) #df['index'][:-1] + 1 #add +1 to all for next position is +/-, except for last position, will throw erre - so it deletes it, we'll add it in later

    new_list_2 = new_list['index']

    new_y_val = list(smoth[new_list])   #find position y on smooth line

    #assigning pos / neg for that +1 position
    pos_neg_results = []
    for i in new_y_val:
        if i > 0:
            result = '+'
        else:
            result = '-'
        pos_neg_results.append(result)

    #creating dataframe for next steps
    df_val.drop(df_val.columns[len(df_val.columns)-1], axis=1, inplace=True) #to delete last column, unnamed so tricky to get rid of (?) this does it tho
    df_val['+/- to the right'] = pos_neg_results

    #append +/- and start stop coords from original table
    df_val.rename(columns={'index': 'Window'}, inplace=True)

    df_val['Window']=df_val['Window'].astype(int)
    df_0['Window']=df_0['Window'].astype(int)

    ## merging
    merged_df = df_val.merge(df_0, how = 'inner', on = ['Window'])

    merged_df = merged_df.drop(['Pass/Fail','Score','VirusGene','HypotheticalGene','Intergenic','BacterialGene'], axis = 1)
    merged_df['Chunk_end'] = 'none'
    merged_df['Window midpoint'] = merged_df.iloc[:,[2,3]].median(axis=1)
    merged_df['Window midpoint'] = merged_df['Window midpoint'].astype(int)

    #df edits to accomodate this:
    #we are duplicating the last row of the df to handle a trailing + chunk (w/ no y=0 intercept to close the chunk)
    #merged_df = merged_df.append(merged_df[-1:])
    merged_df = pd.concat([merged_df, merged_df[-1:]])
    #now need to make it read actual last stop position (this os not rounded per window like the other coords)
    merged_df = merged_df.replace(merged_df.iloc[-1][3],(total_len+1))

    #now let's get the coordinates for the > 0 'chunks'
    #iterate over for true hit testing
    def pairwise(iterable):
        "s -> (s0,s1), (s1,s2), (s2, s3), ..."
        a, b = tee(iterable)
        next(b, None)
        return zip(a, b)

    #this is to define the chunks, accounting for all the ways the graph can look
    #note: leading and trailing here mean a chunk at the start or end of the graph that

    ## use nearby ORF boundaries to refine cutoffs:

    def right_cutoff(position, groupframe):
        calced_end =  int(position)
        rbisect_pos = bisect.bisect(list(groupframe['gene_stop']), calced_end)
        rchunk_cutoff = list(groupframe['gene_stop'])[rbisect_pos-1]
        return rchunk_cutoff

    def left_cutoff(position, groupframe):
        calced_start =  int(position)
        lbisect_pos = bisect.bisect(list(groupframe['gene_start']), calced_start)
        lchunk_cutoff = list(groupframe['gene_start'])[lbisect_pos]
        return lchunk_cutoff

    ddf_list = []

    for (i1, row1), (i2, row2) in pairwise(merged_df.iterrows()):
        #for a leading chunk
        if row1['+/- to the right'] == '+' and \
            row1["Position start"] == 0 and \
            row1["Position stop"] != (total_len + 1):

                ddf = ["Chunk_" + str(i1), 
                       row1["Position start"], 
                       right_cutoff(row2["Window midpoint"], group)]
                ddf_list.append(ddf)
        #for a contained chunk
        if row1['+/- to the right'] == '+' and \
            row1["Position start"] != 0 and \
            row1["Position stop"] != (total_len + 1):
                ##making sure that the adjusted coordinates actually make sense
                leco = left_cutoff(row1["Window midpoint"], group)

                rico = right_cutoff(row2["Window midpoint"], group)

                if rico > leco:
                    ddf = ["Chunk_" + str(i1), 
                        left_cutoff(row1["Window midpoint"], group), 
                        right_cutoff(row2["Window midpoint"], group)]
                else:
                     ddf = ["Chunk_" + str(i1), row1["Window midpoint"], row2["Window midpoint"]]

                ddf_list.append(ddf)

        #3. for a trailing chunk
        if row1['+/- to the right'] == '+' and \
            row1["Position start"] != 0 and \
            row1["Position stop"] == (total_len + 1):
                
                ddf = ["Chunk_" + str(i1), 
                       left_cutoff(row1["Window midpoint"], group), 
                       row2["Position stop"]]
                ddf_list.append(ddf)

        #4. for graphs with no leading and no trailing chunk (for graphs with no y = 0 intercept -> this is is
        #a differently-defined statemnt below b/c the empty file gets appended w/ stuff above from older files when
        #it's in the loop, ALSO the criterion gets fulfilled by contained cunks which means duplicate csv rows for chunks (defined diffrently to specifiy the rules)
        if merged_df.iloc[0,1] == '+' and \
            merged_df.iloc[0,2] == 0 and \
            merged_df.iloc[0,3] == (total_len + 1): #if first column last(2nd row) == last -1 then its one chunk
                rep_list = [('Chunk_0', '0', (total_len+1))]
                ddf_list = rep_list
        else:
                ddf_list = ddf_list

    ## make chunk csv
    chunk_df = pd.DataFrame(ddf_list, columns=["chunk_number", "left_cutoff", "right_cutoff"])

    chunk_df['contig'] = name

    chunk_df = chunk_df[['contig', 'left_cutoff', 'right_cutoff', 'chunk_number']]

    chunk_sum_file = os.path.join(out_dir1, name + ".chunks.tsv")

    chunk_df.to_csv(chunk_sum_file, sep = "\t", index = False)

    ###Find optimal location on plot to place hallmark marker
    vir_bait_table = group[['gene_start', 'gene_stop', 'Evidence_source']]\
        .query("Evidence_source == 'hallmark_hmm'")
    vir_bait_table['mean'] = round(vir_bait_table[['gene_start', 'gene_stop']].mean(axis=1))
    vir_bait_table_med_list = list(vir_bait_table['mean'])

    points_list = []
    for item in vir_bait_table_med_list:
        eq = round(((item - 2500) + 50) / 50)
        if eq >= len(x):
            plot_point = (len(x) - 1) #1 because it can't = len, has to be less
        else:
            plot_point = eq

        points_list.append(plot_point)

    new_points_list = [1 if i <=0 else i for i in points_list]

    df_0['smoothy'] = smooth(df_0['Score'],100)

    #FIGURES
    pdf_outname = os.path.join(out_dir1, name + ".figures.pdf")

    figures = PdfPages(pdf_outname)

    ## first figure, shows smoothed average

    plt.plot(x, y, 'o', ms=0.6)
    plt.axhline(0, 0, l)
    plt.plot(x, df_0['smoothy'], 'c', lw=2)
    ## add halmark gene markers
    plt.plot(x, df_0['smoothy'], 'y', markevery = (new_points_list), ms=11.0, marker = '*')
    plt.title("Viral region calls")
    plt.xlabel('Window')
    plt.ylabel('Score')
    plt.rc('axes', titlesize=6.8)     # fontsize of the axes title
    plt.rc('xtick', labelsize=5)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=5)    # fontsize of the tick labels
    plt.rc('legend', fontsize=5)    # legend fontsize
    plt.rc('figure', titlesize=8)  # fontsize of the figure title
    plt.grid(True)
    idx = np.argwhere(np.diff(np.sign(zero - smooth(y,100)))).flatten()
    plt.plot(x[idx], zero[idx],  'ro', ms=5.0)


    plt.plot()
    plt.savefig(figures, format='pdf')
    plt.close()

    ## second figure shows gene type distribution
    mycol = (["#e7ba52", "#637939", "#7b4173", "#d6616b"])
    df_0[['VirusGene','HypotheticalGene','BacterialGene','Intergenic']].plot(color = mycol)
    plt.grid(True)
    plt.xlabel('Window')
    plt.ylabel('Count')
    plt.title('Character counts')
    plt.rc('axes', titlesize=6.8)     # fontsize of the axes title
    plt.rc('xtick', labelsize=5)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=5)    # fontsize of the tick labels
    plt.rc('legend', fontsize=5)    # legend fontsize
    plt.rc('figure', titlesize=8)  # fontsize of the figure title
    plt.savefig(figures, format='pdf')
    plt.close()

    figures.close()
