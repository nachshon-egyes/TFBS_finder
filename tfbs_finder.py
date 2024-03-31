import math
import os
import sys
import re
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import Bio
from Bio import motifs


# create dictionary with relevant information on each of the TFBSs
def motif_stats(tfbs_dir):
    stats_list = ['jaspar_id', 'TF', 'length', 'consensus', 'degen_consensus', 'matrix', 'pwm', 'pssm', 'max_scores', 'mean', 'dev']
    stats = {stat : [] for stat in stats_list}
    for file in os.listdir(tfbs_dir):
        with open(tfbs_dir + "/" + file) as handle:
            for m in motifs.parse(handle, "jaspar"):
                temp_max_scores = []       
                stats['jaspar_id'].append(m.matrix_id)
                stats['TF'].append(m.name)
                stats['length'].append(len(m))
                stats['consensus'].append(m.consensus)
                stats['degen_consensus'].append(m.degenerate_consensus)
                stats['matrix'].append(m.counts)
                stats['pwm'].append(m.pwm)
                stats['pssm'].append(m.pssm)
                for i in range(len(m)):
                    temp_max_scores.append(max(m.pssm[:,i].values()))
                stats['max_scores'].append(temp_max_scores)
                stats['mean'].append(np.array(temp_max_scores).mean())
            grand_mean = np.array(stats['mean']).mean()
            stats['dev'] = np.array(stats['mean']) - grand_mean
    return stats

# create dictionary with 'id:sequence' items for all test sequences provided by user
def prepare_test_seqs(seq_table, motif_stats, compare = False):

    if isinstance(seq_table, pd.DataFrame):
        if compare == True:
            seq1_list = seq_table.loc[:, 'sequence1']
            seq2_list = seq_table.loc[:, 'sequence2']
            test_seqs = list(zip(seq1_list, seq2_list))
        else:
            test_seqs = list(seq_table['sequence'])
            test_seqs = [i for i in test_seqs]  
        test_ids = list(seq_table['id'])
    elif seq_table[-5:] == 'fasta' and compare == False:
        from Bio import SeqIO
        test_seqs_raw = [seq for seq in SeqIO.parse(seq_table, 'fasta')]
        test_seqs = [str(test_seqs_raw[i].seq) for i in range(len(test_seqs_raw))]
        test_ids = [test_seqs_raw[id].id for id in range(len(test_seqs_raw))]
    else:
        raise TypeError('file type given is not valid, or a fasta file was given with compare mode')
    test_seqs_dict = dict(zip(test_ids, test_seqs))
    return test_seqs_dict

# create table with all events where TFBSs were found within the test sequences
def find_matches(stats, test_seqs_dict):
    if test_seqs_dict == {} or stats == {}:
        print('The test sequence file or JASPAR files were not provided properly.')
        sys.exit('Please try again')

    scoring_list = list()
    legend = {'A':0, 'C':1, 'G':2, 'T':3}
    for key in test_seqs_dict:
        seq = test_seqs_dict[key].upper()
        if bool(re.search("[^ACGT]", seq)):
            print('One of the test sequences has an unreadable nucleotide. Use only A C G T.')
            continue

        for motif in range(len(stats['TF'])):
            ln = stats['length'][motif]
            
            for nt in range(len(seq) - ln + 1):
                subseq = seq[nt:nt+ln]
                subseq_ix = [legend[ix] for ix in list(subseq)]
                temp_score = np.array(list(stats['pssm'][motif].values()))[subseq_ix, list(np.arange(0,ln,1))]
                    
                if temp_score.mean() > stats['mean'][motif] - 0.5:
                    norm_score = temp_score.mean() - stats['dev'][motif]
                    scoring_list.append([key, stats['jaspar_id'][motif], nt, nt+ln, subseq, stats['consensus'][motif],
                                         subseq == stats['consensus'][motif], stats['TF'][motif], round(temp_score.mean(),3),
                                         round(norm_score,3), round(stats['mean'][motif],3), round(stats['mean'][motif]-stats['dev'][motif],3)])
    scoring_df = pd.DataFrame(scoring_list,
                              columns = ['sequence_ID', 'jaspar_ID', 'start', 'end', 'subseq', 'consensus',
                                         'is_consensus', 'TF', 'score', 'norm_score', 'TFBS_score', 'norm_TFBS_score'])
    return scoring_df

# same as previous function, adjusted for compare mode
# The function includes many lines of code to enable robustness and less need for requirement compliance by the user 
def compare_matches(stats, test_seqs_dict):
    if test_seqs_dict == {} or stats == {}:
        print('The test sequence file or JASPAR files were not provided properly.')
        sys.exit('Please try again')    
    scoring_list = list()
    legend = {'A':0, 'C':1, 'G':2, 'T':3}

    for key in test_seqs_dict:
        loc = None
        seq1, seq2 = test_seqs_dict[key][0].upper(), test_seqs_dict[key][1].upper()
        if bool(re.search("[^ACGT]", seq1)) or bool(re.search("[^ACGT]", seq2)):
            print(f'Non A C T G characters found in at least one of the sequence pairs. Skipping {key} sequence pair.')
            continue

        if len(seq1) <= len(seq2): #this is in case the 2 sequences are not the same length
            seq_len = len(seq1)
        else:
            seq_len = len(seq2)
        
        for mutation_ix in range(seq_len): #find the location of the difference between both sequences
            if seq1[mutation_ix] != seq2[mutation_ix]:
                loc = mutation_ix + 1
                break
        if loc == None:
            print(f'No base substitution was found between {key} sequence pair')
            continue

        for motif in range(len(stats['TF'])):
            counter = 0
            motif_len = stats['length'][motif]

            if motif_len + loc < seq_len: #find the relevant indices for which to search the sequences for TFBSs 
                if loc > motif_len:
                    srch_rng = range(loc-motif_len, loc)
                else:
                    srch_rng = range(loc)
            else:
                if loc > motif_len:
                    srch_rng = range(loc-motif_len, seq_len-motif_len)
                else:
                    srch_rng = range(seq_len-motif_len+1)

            for nt in srch_rng:
                subseq1, subseq2 = seq1[nt:nt+motif_len], seq2[nt:nt+motif_len]
                subseq1_ix, subseq2_ix = [legend[ix] for ix in list(subseq1)], [legend[ix] for ix in list(subseq2)]
                seq1_snp, seq2_snp = legend[seq1[loc-1]], legend[seq2[loc-1]]
                temp_score1 = np.array(list(stats['pssm'][motif].values()))[subseq1_ix, list(np.arange(0,motif_len,1))]
                temp_score2 = np.array(list(stats['pssm'][motif].values()))[subseq2_ix, list(np.arange(0,motif_len,1))]
                seq1_snp_score = np.array(list(stats['pssm'][motif].values()))[seq1_snp, stats['length'][motif]-1-counter]
                seq2_snp_score = np.array(list(stats['pssm'][motif].values()))[seq2_snp, stats['length'][motif]-1-counter]
                counter += 1

                if (temp_score1.mean() > stats['mean'][motif] - 0.5) or (temp_score2.mean() > stats['mean'][motif] - 0.5):
                    norm_score1, norm_score2 = temp_score1.mean() - stats['dev'][motif], temp_score2.mean() - stats['dev'][motif]
                    scoring_list.append([key, stats['jaspar_id'][motif], subseq1, subseq2, stats['consensus'][motif],
                                         subseq1 == stats['consensus'][motif], subseq2 == stats['consensus'][motif],
                                         stats['TF'][motif], round(temp_score1.mean(),3), round(temp_score2.mean(),3), round(norm_score1,3),
                                         round(norm_score2,3), round(norm_score1 - norm_score2,3), round(seq1_snp_score - seq2_snp_score,3),
                                         round(stats['mean'][motif],3)])
    scoring_df = pd.DataFrame(scoring_list, columns = ['sequence_ID', 'jaspar_ID', 'subseq1', 'subseq2' , 'consensus',
                                         'is_consensus1', 'is_consensus2', 'TF', 'score1', 'score2',
                                         'norm_score1', 'norm_score2', 'diff_in_scores', 'diff_base_score', 'TFBS_score'])
    return scoring_df

# create a table with summary statistics on the run   
def summary_table(scoring_df):
    if scoring_df.shape[1] == 12:
        gb = scoring_df.groupby('TF')
        summary = (gb.size().to_frame(name='counts')
         .join(gb.agg({'norm_score':'sum'}).rename(columns={'norm_score':'sum_scores'}))
         .join(gb.agg({'norm_score':'mean'}).rename(columns={'norm_score':'mean_scores'}))
         .sort_values('counts', ascending = False)
         .reset_index())
    else:
        scdf = scoring_df
        scdf.loc[scoring_df['diff_in_scores'] == math.inf, 'diff_in_scores'] = 2
        scdf.loc[abs(scoring_df['diff_in_scores']) == math.inf, 'diff_in_scores'] = -2
        gb = scdf.groupby('TF')
        summary = (gb.size().to_frame(name='counts')
         .join(gb.agg({'diff_in_scores':'mean'}).rename(columns={'diff_in_scores':'mean_difference'}))
         .assign(abs_mean_diff = lambda x: abs(x.mean_difference))
         .sort_values('abs_mean_diff', ascending = False)
         .reset_index())
    return summary

# plot the number of times the top TFBSs were found / the mean change in binding affinity for top TFs    
def enrichment_plot(summary):
    for_plot = summary.iloc[0:10]
    fig, ax = plt.subplots()
    if summary.columns[2] == 'sum_scores':
        ax.bar(for_plot['TF'], for_plot['counts'])
        ax.set_ylabel('# of times TFBS found in sequences')
    else:
        for_plot = for_plot.sort_values('mean_difference')
        ax.bar(for_plot['TF'], for_plot['mean_difference'])
        ax.set_ylabel('TF binding affinity change')
    ax.tick_params(axis='x', rotation=90)
    ax.set_xlabel('Transcription factor')
    return fig

# prepare arguments to pass to main function
def sys_argvs():
    if len(sys.argv) == 3:
        argvs = [sys.argv[1], sys.argv[2]]

    elif len(sys.argv) == 2:
            argvs = [sys.argv[1], "tfbs"]

    else:
        print("parameters given are invalid. Type code file name, then file name containing sequence data,"
                       "and then directory containing JASPAR files")
        sys.exit('Please try again')
    
    return argvs

# call the functions 
def main():
    argvs = sys_argvs()
    seq_table, tfbs_dir = argvs[0], argvs[1]
    if seq_table[-3:] == 'csv':
        seq_table = pd.read_csv(seq_table)
        if 'sequence1' in list(seq_table.columns) and 'sequence2' in list(seq_table.columns):
            compare = True
    else:
        compare = False

    if os.path.exists('output'):
        import shutil
        shutil.rmtree('output')
    os.mkdir('output')
    stats = motif_stats(tfbs_dir = tfbs_dir)
    pd.DataFrame(stats).to_csv('output/TFBS_data.csv')
    test_seqs_dict = prepare_test_seqs(seq_table = seq_table, motif_stats = stats, compare = compare)
    if compare == False:
        scoring_df = find_matches(stats = stats, test_seqs_dict = test_seqs_dict)
        scoring_df.to_csv('output/matches.csv')
    else:
        scoring_df = compare_matches(stats = stats, test_seqs_dict = test_seqs_dict)
        scoring_df.to_csv('output/matches_comp.csv')
    summary = summary_table(scoring_df = scoring_df)
    summary.to_csv('output/summary.csv')
    fig = enrichment_plot(summary)
    fig.savefig('output/TF_enrichment.png', bbox_inches = 'tight', dpi = 100)

if len(sys.argv) > 1:
    main()

