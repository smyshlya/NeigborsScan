# README #

This repository contains Perl script neighbors_scan.pl and an example input file: example_mapping_table.

neighbors_scan.pl implements the search of conserved proteins in the genome found within user defined proximity from proteins of interest. Proteins of interest are defined in the input mapping table file ('uniprot accession number' -> 'GI', such as 'example_mapping_table'). This file can be produced [here](http://www.uniprot.org/uploadlists/) using a list of Uniprot accessions. Conserved proteins to be scanned for are defined in the input hmm library.
 
The script 
* finds accession numbers of nucleotide sequence, where the input uniprot sequence is located 
* downloads all proteins within a specified range from the beginning of the input uniprot sequence
* runs hmmscan analysis against the input hmm library using the downloaded proteins as queries 
* parses the results into output file ('example_mapping_table.out'). The output file can be uploaded as a dataset at [iTol](http://itol.embl.de/upload.cgi)
 
 ![1](https://user-images.githubusercontent.com/58728948/113287575-0c3b0b80-92ee-11eb-9db4-79a0556cfe75.png "Title")
 
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

[Download](https://github.com/gem-pasteur/Macsyfinder_models/tree/master/models/Conjugation/profiles) CONJscan database. [This](https://www.ncbi.nlm.nih.gov/pubmed/31584169) is the paper where the database is described. Specify input files, folder for saving the results and range for the scan in the script. Then run:
```bat
hmmpress input_hmm_library 
perl neigbors_scan.pl
```
