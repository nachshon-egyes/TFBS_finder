# TFBS finder
### Locate transcription factor binding sites (TFBSs) in a set of sequences

## Quick example
`python tfbs_finder.py ex_compare_data.csv`

## Background

It has been established that much of the phenotypic differences we observe between closely-related species are due to changes in [*cis*-regulatory elements](https://en.wikipedia.org/wiki/Cis-regulatory_element) (CREs). To better understand the molecular outcomes of these changes, it is crucial to find the various [transcription factors](https://en.wikipedia.org/wiki/Transcription_factor) (TFs) that bind these CREs, and uncover the effect these changes have on binding affinity.  
  
My script aims to do so by 2 separate modes:
- Finding the TFs that bind a given set of sequences and their relative propnesity (regular mode).
- Find the effect a single nucleotide substitution mutation has on the affinity of a TF to its corresponding binding site (compare mode).

## How to use the code
1. Make sure all dependencies are installed. *requirements.txt* lists them all.
1. In terminal, clone the repository and change current directory (`cd`) to it's path.
1. After typing `python`, the first argument should be the main script file - `tfbs_finder.py`
1. For the next argument, choose a file containing sequence data you would like to have analyzed. If you wish to use 'compare mode', ensure the file provided is of type *csv* and has columns with 3 of the following names - **id, sequence1, sequence2**. The sequences in each pair should be identicle, excluding the nucleotide substitution.  
If you wish to use 'regular mode', the csv file should contain **id** and **sequence** columns or alternatively use a *fasta* file holding all sequences.
1. Lastly, provide a directory with JASPAR files where each one holds relevant TFBS data. **If left blank, the default directory is used.**  

#### Examples:
`python tfbs_finder.py ex_data.csv tfbs`  
`python tfbs_finder.py ex_data.fasta`  
`python tfbs_finder.py ex_data.csv`

## Output
### 4 different files
- **TFBS_data.csv**

| jaspar_id | TF | length | consensus | degen_consensus | matrix | pwm | pssm | max_scores | mean | dev |  
| --------- | -- | ------ | --------- | --------------- | ------ | --- | ---- | ---------- | ---- | --- |
| TFBS ID | transcription factor name | length of sequence the TF binds to | the sequence the TF is most likely to bind | the degenerate TFBS - special letters are used for positions with more than 1 possible nucleotide | number of times each nucleotide was shown to appear at each position of the TFBS | position weighted matrix - relative frequency of each nucleotide at each position | position-specific scoring matrix - score for " " " " " | list of all maximal scores for each position | mean of list | deviation of mean from the grand mean (means of all TFBSs) |

- **matches.csv** / **matches_comp.csv**

| sequence_ID |	jaspar_ID |	start |	end	| subseq | consensus | is_consensus |	TF | score | norm_score | TFBS_score | norm_TFBS_score |
| ----------- | --------- | ----- | --- | ------ | --------- | ------------ | -- | ----- | ---------- | ---------- | --------------- |
| ID of sequence tested | TFBS ID that was found within the sequence | start index of TFBS found within the test sequence | end index | the sub-sequence found | the TFBS consensus sequence | boolean indicating if the sub-sequence found is the consensus sequence | naem of TF that binds the TFBS | score of sub-sequence as calculatd using the pssm | score normalized to the grand mean (for comparison to other scores) | maximal score for that TFBS | the grand mean |

- **summary.csv**  
summary of all matching results

| TF | counts | sum_scores | mean_scores | mean_difference | abs_mean_diff |
| -- | ------ | ---------- | ----------- | --------------- | ------------- |
| | number of times the TF was found within the sequences | in regular mode - sum of all normalized scores for that TF | " " " - mean " " " " " " " | in compare mode - mean difference in scores for all sequence pairs with that TF | |  


Note - in compare mode, infinity values were set to 2 for enabling data visualization of score difference while still keeping the large difference in binding affinity that's observed.
 

- **TF_enrichment.png**  
For regular mode: a bar plot displaying the number of times the top 10 TFBSs were found within the sequences.  
For compare mode: a bar plot displaying the mean difference in scores for top 10 TFs found to bind the sequence pair.  


## Additional information
### Runtime:
Applying the code on ~20 reads of 600 nucleotides each should take a couple minutes of computation. Compare mode is much quicker.

### Testing:
Make sure pytest is installed and run `pytest` .

### Biological data:
JASPAR files from the JASPAR transcription factor binding site databse were used to search for TFBSs.  
Promoter DNA sequences from the Eukaryotic Promoter Database (EPD) were used for the example files.
