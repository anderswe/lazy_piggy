The other file in this folder:

for_reviewers/processed_data/wgs_lp_tumours/all_annotations_20211201.csv.gz

is a table you can decompress using the Archive utility on Mac, or directly by running:

gunzip <filename> in the terminal, once you have downloaded it.

From there, it is a regular .csv that can be opened in Microsoft Excel, or any R/Python/etc etc environment.

Each row in the data indicates the read support for specific insertion sites for each mouse + library combination.

There are 8 columns in the data:

1. Sample:
e.g. 01-07-16M_PT_IRL
where "01-07-16M" indicates the mouse ID, which corresponds to metadata stored in for_reviewers/processed_data/metadata/LP_sample_metadata.txt.
"PT" indicates this sample is a Primary Tumour, and "IRL" indicates the library, e.g. "Inverted Repeats, Left"...PBR = "piggyBac, Right"...JXR = "Junction library, Right"
2. Gene
3. Gene Location:
Indicates where in the gene the insertion is.
4. Functional effect: 
Predicted effect of the insertion on gene function. E.g. disrupt_CDS = disrupt coding sequence. E.g. 28kb = distance from the gene body.
5. chr = Chromosome.
6. loc = location (coordinate) on the chromosome.
7. Count = number of supporting reads.
8. Gene Orientation = orientation of the insertion relative to the annotated gene.
