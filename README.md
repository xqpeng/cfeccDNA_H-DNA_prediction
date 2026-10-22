# Introduction
cfeccDNA_H-DNA_prediction contains scripts for linearized cell-free eccDNA (cf-eccDNA) fragments analysis and H-DNA prediction models based certain cf-eccDNA fragments.


## SplitAnalysis
Scripts for mapping cfDNA sequences, extracting discordant and split alignments, and analyzing split alignments.

## Identical6merMotif_analysis
This fold contains scripts for evaluating the enrichment of identical 6-mer end motifs in Same-Strand 5' Ends fragments.

## Triplex_tool_analysis
This fold contians a pipeline for identification and analysis of H-DNA (triplex DNA) sites
from Same-Strand 5' Ends fragments with endpoint distance less than 5 bp by using the Triplex package. 
The pipeline covers the entire workflow from format conversion, H-DNA structural searching, and multi-dimensional statistical visualization.

## Negative_sample_generation
The script is used to generate negative samples for training dataset, valiation dataset, and test dataset.

## Triplex_prediction_models
Models for H-DNA Sequence Identification.
This fold contains the official PyTorch & Scikit-Learn implementations for benchmarking deep learning architectures and traditional machine learning baselines on high-throughput H-DNA sequence prediction datasets. 

