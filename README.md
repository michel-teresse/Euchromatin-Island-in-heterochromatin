# Euchomatin-Island-in-heterochromatin


This package contains the following scripts :

search_island.pl

Identification_of_genes_into_island.py

## search_island.pl

search_island.pl is meant to identify euchromatic island within heterochromatin. It is based on the chromatin state of each base. The script detect the typical eucrhromatic states that are surrounded by a certain length (determined by the user) of state 9 and 8.

### What can search_island.pl do ?

Generate the list of coordinates (chr, start, stop) of the euchromatic islands within heterochromatin as a bed file.

### How to use search_island.pl ?
 
    search_island.pl [-v] chrom_state_header_less MIN_LEN_STATE_8 MIN_LEN_STATE_9
    
With : 

- chrom_state_header_less : The file contening the topography of all chromatin states for Arabidopsis whole genome as determined by Sequeira-Mendes et al with no header. The file `chrom_states_At_Sequeira-Mendes`, available at https://github.com/seeterlab/Chromatin-state-analysis can be used, after removing the header.

- MIN_LEN_STATE_8 : The minamal number of bases covered by state 8 to be considered

- MIN_LEN_STATE_9 : The minimal number of bases covered by state 9 to be considered

## Identification_of_genes_into_island.py

Identification_of_genes_into_island.py is meant to identify genes that are into the euchromatic islands.

### What can Identification_of_genes_into_island.py do ?
Generate a list of the genes included into the euchromatic island and also genes that go over right or left or both borders of the island. At least gene_min_len bases being included in the island.

- The genes that go over the right or left border are tagged with a `*`

- The genes that go over both borders are tagged with a `**`

### How to use  Identification_of_genes_into_island.py ?

    Identification_of_genes_into_island.py genes_ref_file interval_file
    
genes_ref_file should look like this, at least for the first 4 columns:

No_chromosome   gene_id   start   end

Exemple:

    1	AT1G01010	3631	5899	+ 	 26
    1	AT1G01020	6788	9130	- 	 42137
    1	AT1G01030	11649	13714	- 	 25
    1	AT1G01040	23121	31227	+ 	 137


interval_file should look like this (no header):

No_chromosome   start   end

Exemple:

    1 30001 32000
    1 52001 54000
    1 166001 168000

The file interval_file has to be sorted by chromosome number and by coordinates.
Chromosome number in the the file interval_file has to be only 1 to 5.



As a result,Identification_of_genes_into_island.py print interval_file with 2 supplementary columns as follow:
Number_genes_per_intervals  gene_name_list_comma_separated_no_space

Exemple:

    1	12685051	12711000 	 4 	 AT1G34630,AT1G34640,AT1G34650,AT1G34670
    1	12735601	12765300 	 5 	 AT1G34750,AT1G34760,AT1G34770,AT1G34780,AT1G34790*
    1	12839701	12845700 	 0

A star (*) is added when a part of the gene is not included in the interval but at least 500 bp (GENE_MIN_LEN) are.



