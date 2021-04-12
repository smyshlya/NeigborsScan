# README #

This repository contains Perl script neighbors_scan.pl and an example input file: example_mapping_table.

neighbors_scan.pl implements the search of conserved proteins in the genome found within user defined proximity from proteins of interest. Proteins of interest are defined in the input mapping table file ('uniprot accession number' -> 'GI', such as 'example_mapping_table'). This file can be produced [here](http://www.uniprot.org/uploadlists/) using a list of Uniprot accessions. Conserved proteins to be scanned for are defined in the input hmm library.
 
The script 
* finds accession numbers of nucleotide sequence, where the input uniprot sequence is located 
* downloads all proteins within a specified range from the beginning of the input uniprot sequence
* runs hmmscan analysis against the input hmm library using the downloaded proteins as queries 
* parses the results into output file ('example_mapping_table.out'). The output file can be uploaded as a dataset at [iTol](http://itol.embl.de/upload.cgi)
 
![1](https://user-images.githubusercontent.com/58728948/113307335-5085d600-9305-11eb-8911-938a1a35cd86.png)
 
The script produces many additional output files:
*annotation files for nucleotide sequences (.ft)
*protein files (.prot)
*results of hmmscan (.hmm_res)
*files for convertion from uniprot accession numbers to genbank accession numbers and for retrieving corresponding coordinates (input_mapping_table.uniprot_to_an, ~.coordinates, ~.elink, ~.long)

# Citation #
Smyshlyaev G, Barabas O., Bateman A. [Sequence analysis allows functional annotation of tyrosine recombinases in prokaryotic genomes.](https://www.biorxiv.org/content/10.1101/542381v1) 2019. BioRxiv.
# Setting up #
* Dependencies:

Bioperl, hmmscan, hmmpress

* How to run tests:

For identifying conjugation proteins associated with the integrase proteins:
[Download](https://github.com/gem-pasteur/Macsyfinder_models/tree/master/models/Conjugation/profiles) CONJscan database. [This](https://www.ncbi.nlm.nih.gov/pubmed/31584169) is the paper where the database is described. Specify input files, folder for saving the results, range for the scan  and your [API key] (https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/#:~:text=To%20create%20the%20key%2C%20go,and%20copy%20the%20resulting%20key.) to access NCBI in the script. Then run:
```bat
mkdir results/
mkdir hmm/
cd hmm/
hmmpress input_hmm_library 
cd ..
perl neigbors_scan.pl
```
