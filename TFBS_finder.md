# TFBS finder
### A quick way to locate transcription factor binding sites in a set of sequences

It has been established that much of the phenotypic differences we observe between closely-related species are due to changes in [*cis*-regulatory elements](https://en.wikipedia.org/wiki/Cis-regulatory_element) (CREs). To better understand the molecular outcomes of these changes, it is crucial to find the various [transcription factors](https://en.wikipedia.org/wiki/Transcription_factor) (TFs) that bind these CREs, and uncover the effect these changes have on binding affinity.

#### In my project, I aim to build a tool that provides insight into these effects and will do so as follows:
- It will receive a table with a column of DNA sequences.
- Search the sequences for likely transcription factor binding sites (TFBSs).
- Return putative TFBSs found within the sequences (many TFs can tolerate slight changes to their canonical TFBSs), including noteworthy information regarding the TFBS.
- Perform an enrichment analysis of TFBS, examining whether certain TFs underwent substantial regulatory changes.
- In the event that the sequence changes will also be specified by the user, a postulated change in binding affinity by the potential TF will also be returned.


![DNA-bound_TF](image.png)