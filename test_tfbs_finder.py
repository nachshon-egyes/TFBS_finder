import pytest
import math
import os
import sys
import re
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import Bio
from Bio import motifs
import tfbs_finder as tff

def test_sys_argvs():
    with pytest.raises(SystemExit):
        tff.sys_argvs()

def test_find_matches():
    test_stats = {}
    test_seqs = {}
    with pytest.raises(SystemExit):
        tff.find_matches(test_stats, test_seqs)

    test_stats = {'jaspar_id': ['MA0002.3'],'TF': ['Runx1'],'length': [9],'consensus': ['CTGTGGTTT'],
    'pssm': [{'A':[-2.02, -3.13,  -np.Inf, -2.52,  -np.Inf, -4.88, -5.64, -1.93,  0.0],
    'C':[ 1.1 ,  -np.Inf, -2.74, -1.98,  -np.Inf, -3.57, -0.32, -0.11, -1.66],
    'G':[-1.75, -6.16,  1.9 , -2.84,  1.99,  1.89, -0.99, -2.63, -0.79],
    'T':[ 0.39,  1.95, -3.24,  1.78, -5.27, -2.43,  1.42,  1.41,  1.07]}],
    'mean': [1.61],'dev': np.array([0.0])}

    test_seqs = {'id_ex':'CTGTGGTTTTTGTGGTTT'}

    expected_result = pd.DataFrame([['id_ex', 'MA0002.3', 0, 9, 'CTGTGGTTT', 'CTGTGGTTT', True, 'Runx1', 1.612, 1.612, 1.61, 1.61],
    ['id_ex', 'MA0002.3', 9, 18, 'TTGTGGTTT', 'CTGTGGTTT', False, 'Runx1', 1.533, 1.533, 1.61, 1.61]],
    columns = ['sequence_ID', 'jaspar_ID', 'start', 'end', 'subseq', 'consensus', 'is_consensus', 'TF', 'score', 'norm_score', 'TFBS_score', 'norm_TFBS_score'])

    result = tff.find_matches(test_stats, test_seqs)
    assert all(result == expected_result)

    test_seqs = {'id_ex':'CTGTGGTTWTTGTGGTTT'} # test case where non A C G T letters are provided
    expected_result = pd.DataFrame(columns = ['sequence_ID', 'jaspar_ID', 'start', 'end', 'subseq', 'consensus',
    'is_consensus', 'TF', 'score', 'norm_score', 'TFBS_score', 'norm_TFBS_score'])
    result = tff.find_matches(test_stats, test_seqs)
    assert all(result == expected_result)

def test_compare_matches():
    test_stats = {}
    test_seqs = {}
    with pytest.raises(SystemExit):
        tff.compare_matches(test_stats, test_seqs)

    test_stats = {'jaspar_id': ['MA0002.3'],'TF': ['Runx1'],'length': [9],'consensus': ['CTGTGGTTT'],
    'pssm': [{'A':[-2.02, -3.13,  -np.Inf, -2.52,  -np.Inf, -4.88, -5.64, -1.93,  0.0],
    'C':[ 1.1 ,  -np.Inf, -2.74, -1.98,  -np.Inf, -3.57, -0.32, -0.11, -1.66],
    'G':[-1.75, -6.16,  1.9 , -2.84,  1.99,  1.89, -0.99, -2.63, -0.79],
    'T':[ 0.39,  1.95, -3.24,  1.78, -5.27, -2.43,  1.42,  1.41,  1.07]}],
    'mean': [1.61],'dev': np.array([0.0])}

    test_seqs = {'id_ex':['CTGTGGTtT','TTgtgGTTTAGC']} # test use of lowercase letters and sequences of different lengths

    expected_result = pd.DataFrame([['id_ex', 'MA0002.3', 'CTGTGGTTT', 'CTGTAGTTT', 'CTGTGGTTT', True, False, 'Runx1', 1.612, -np.Inf, 1.612, -np.Inf, np.Inf, -0.79, 1.61]],
    columns = ['sequence_ID', 'jaspar_ID', 'subseq1', 'subseq2', 'consensus', 'is_consensus1', 'is_consensus2', 'TF', 'score1', 'score2',
       'norm_score1', 'norm_score2', 'diff_in_scores', 'diff_base_score', 'TFBS_score'])

    result = tff.compare_matches(test_stats, test_seqs)
    assert all(result == expected_result)

    test_seqs = {'id_ex':['CTGTGGTTT','NTGTGGTTT']} # test use of non A C T G letters
    expected_result = pd.DataFrame(columns = ['sequence_ID', 'jaspar_ID', 'subseq1', 'subseq2' , 'consensus',
                                         'is_consensus1', 'is_consensus2', 'TF', 'score1', 'score2',
                                         'norm_score1', 'norm_score2', 'diff_in_scores', 'diff_base_score', 'TFBS_score'])
    result = tff.compare_matches(test_stats, test_seqs)
    assert all(result == expected_result)

    test_seqs = {'id_ex':['CTGTGGTTT','CTGTGGTTT']} # test case when both sequences provided are the same
    result = tff.compare_matches(test_stats, test_seqs)
    assert all(result == expected_result)

def test_summary_tabel():
    test_scoring_df = pd.DataFrame([['id_ex', 'MA0002.3', 0, 9, 'CTGTGGTTT', 'CTGTGGTTT', True, 'Runx1', 1.612, 1.612, 1.61, 1.61],
    ['id_ex', 'MA0002.3', 9, 18, 'TTGTGGTTT', 'CTGTGGTTT', False, 'Runx1', 1.533, 1.533, 1.61, 1.61]],
    columns = ['sequence_ID', 'jaspar_ID', 'start', 'end', 'subseq', 'consensus', 'is_consensus', 'TF', 'score', 'norm_score', 'TFBS_score', 'norm_TFBS_score'])

    expected_result = pd.DataFrame([['Runx1', 2, 3.145, 1.5725]], columns = ['TF', 'counts', 'sum_scores', 'mean_scores'])
    result = tff.summary_table(test_scoring_df)
    assert all(result == expected_result)




















