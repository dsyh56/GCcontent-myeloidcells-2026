This repo contains the scripts used to generate the data for analysing GC-content and other splicing-related properties in myeloid cell gene expression datasets:

1_prepRefFileFromGTF:
Main script for generating a core file annotated with exon-level properties. Used for annotating and retrieval of data in most downstream steps.

2_HPAanalysis:
Scripts for collating immune cell-type enriched/enhanced genes from the Human Protein Atlas database, filtering/grouping the subsets, and retrieving gene-level parameters for analysis.
Also used as a base for further exploration of the data (e.g. gene enrichment)

3_mappingUniprotToExons:
Scripts for mapping Uniprot protein features to exon coding sequences, as well as verifying their reading frames. Used for creating annotation files for downstream protein-relevant analyses.

4_splicingFunctionalEffects:
Scripts for simulating the effect of alternative splicing events on protein coding and/or NMD status of a transcript. Used for annotation of differential alternative splicing events.

5_rMATSfilters:
General script for filtering rMATS results based on FDR, deltaPSI, minimum junction counts, and even coverage of 3' and 5' ends of exons (when relevant). Also for calculating various exon properties for the dataset and assigning functional annotations for each splicing event.
