# README #

### What is this repository for? ###

This repository contains Perl script neighbors_scan.pl and two example input files: example_mapping_table and CONJscan.hmm3.

neighbors_scan.pl implements the search of conserved proteins in the genome found within user defined proximity from proteins of interest. Proteins of interest are defined in the input mapping table file ('uniprot accession number' -> 'GI', such as 'example_mapping_table'). This file can be retrieved from http://www.uniprot.org/uploadlists/. Conserved proteins to be discovered are defined in the input hmm library (example library is 'CONJscan.hmm3').
 
The script 1) finds accession numbers of nucleotide sequence, where the input uniprot sequence is located; 2) downloads all proteins within a specified range from the beginning of the input uniprot sequence; 3) runs hmmscan analysis against the input hmm library using the downloaded proteins as queries; 4) parses the results into output file ('example_mapping_table.out'). The output file can be uploaded as a dataset at http://itol.embl.de/upload.cgi
 
The script produces many additional output files:

-annotation files for nucleotide sequences (.ft)

-protein files (.prot)

-results of hmmscan (.hmm_res)

-files for convertion from uniprot accession numbers to genbank accession numbers and for retrieving corresponding coordinates (input_mapping_table.uniprot_to_an, ~.coordinates, ~.elink, ~.long)

### How do I get set up? ###

* Dependencies:

Bioperl, hmmscan, hmmpress

* How to run tests

Specify input files, folder for saving the results and range for the scan in the script. Then run:

hmmpress input_hmm_library 

cd script_folder

perl neigbors_scan.pl

* Deployment instructions

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact