#!/usr/bin/env python3

# Identification of gene that are into the euchromatic islands.
#
# Copyright (c) 2017 Michel TERESE
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

CHROMOSOME_COUNT = 5
# Minimal length of the gene into the interval to be considered
GENE_MIN_LEN = 500

import sys
from bisect import bisect

# Check arguments
if len(sys.argv) -1 != 2 :
    print( "syntax: ", sys.argv[0],
'''genes_ref_file interval_file

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



Print interval_file with 2 supplementary columns as follow:
Number_genes_per_intervals  gene_name_list_comma_separated_no_space

Exemple:
1	12685051	12711000 	 4 	 AT1G34630,AT1G34640,AT1G34650,AT1G34670
1	12735601	12765300 	 5 	 AT1G34750,AT1G34760,AT1G34770,AT1G34780,AT1G34790*
1	12839701	12845700 	 0

A star (*) is added when a part of the gene is not included in the interval but at least 500 bp (GENE_MIN_LEN) are.
''' )
    sys.exit(1)

genes_ref_file = sys.argv[1]
interval_file = sys.argv[2]

# The tables below will be indexed with chromosome number (from 1 to n)
# Every elements in the tables will be a list of position for the concerned chromosome
# The index number will be initialized to None because not used
a_gene_ID = []
a_gene_start = []
a_gene_end = []

# Initialize the 2 tables
for i in range(0, CHROMOSOME_COUNT +1) :
    a_gene_ID.append(i)
    a_gene_ID[i] = [] if i > 0 else None
    a_gene_start.append(i)
    a_gene_start[i] = [] if i > 0 else None
    a_gene_end.append(i)
    a_gene_end[i] = [] if i > 0 else None

# Read the file genes_ref_file and initialize the table a_gene_end
with open( genes_ref_file ) as f:  
    for line in f :  
        # Skip the empty lines
        if not line.strip() :
            continue

        # Split the line
        (chrom, ID, gene_start, gene_end) = line.split(maxsplit=4)[0:4]
        chrom = int(chrom)
        gene_start = int(float(gene_start))
        gene_end = int(float(gene_end))

        a_gene_ID[chrom].append( ID )
        a_gene_start[chrom].append( gene_start )
        a_gene_end[chrom].append( gene_end )

# DEBUG: print the tables
# for chrom in range(1, CHROMOSOME_COUNT+1) :
#     for i in range(0, len(a_gene_end[chrom])) :
#          print( chrom, a_gene_end[chrom][i] )

# Read the file interval_file
with open( interval_file ) as f:  
    for line in f :  
        # Skip the empty lines
        if not line.strip() :
            continue

        # Split the line
        (chrom, island_start, island_end) = line.split()  
        chrom = int(chrom)
        island_start = int(float(island_start))
        island_end = int(float(island_end))
        island_len = island_end - island_start +1
        gene_min_len = min(island_len , GENE_MIN_LEN)

        genes_list = []

        #================================================================================
        # Search for genes included in the island some of which can go over the right border
        # of the island. At least gene_min_len bases are included in the island 
        # (or lenght of the island if it is shorter than gene_min_len)
        #================================================================================
        # Search for the closer "island_start" on the right in genes_ref_file
        i = bisect(a_gene_start[chrom], island_start)
        # if it exists
        if i < len(a_gene_start[chrom]) and a_gene_start[chrom][i] < island_end - gene_min_len :

            gene_name = a_gene_ID[chrom][i]

            if a_gene_end[chrom][i] > island_end :
                gene_name += '*'
                
                genes_list.append( gene_name )

            else :
                genes_list.append( gene_name )

                # Search for genes included in the interval
                while True :
                    i += 1
                    if a_gene_start[chrom][i] <= island_end - gene_min_len :
                        gene_name = a_gene_ID[chrom][i]

                        if a_gene_end[chrom][i] > island_end :
                            gene_name += '*'
                            genes_list.append( gene_name )
                            break

                        genes_list.append( gene_name )

                    else :
                        break

        #================================================================================
        # Search for genes included in the island which can go over both right and left 
        # island borders
        #================================================================================
        # Search for the closer "island_end" on the right in genes_ref_file
        i = bisect(a_gene_end[chrom], island_end)
#         if i >= 0 and a_gene_start[chrom][i] <= island_end - gene_min_len :
        if i >= 0 and a_gene_start[chrom][i] < island_start :
            # the found gene goes over the 2 borders, 2 stars are added
            gene_name = a_gene_ID[chrom][i] + '**'

            # Check that the found gene is not already in the list
            # Normaly not necessary but twice better than one
            if gene_name not in genes_list :
                genes_list.append( gene_name )


        #================================================================================
        # Search for genes included in the island some of which can go over the left border
        # of the island. At least gene_min_len bases are included in the island 
        # (or lenght of the island if it is shorter than gene_min_len)
        #================================================================================
        # search for the closer "island_end" on the right in genes_ref_file
#         i = bisect(a_gene_end[chrom], island_end) -1
        i -= 1
        # if the number of bases included in the island is >= gene_min_len
        if i >= 0 and a_gene_end[chrom][i] - island_start >= gene_min_len :

            gene_name = a_gene_ID[chrom][i]

            if a_gene_start[chrom][i] < island_start :
                gene_name += '*'
                
                # Check that the found gene is not already in the list
                if gene_name not in genes_list :
                    genes_list.append( gene_name )

            else :
            
                # Check that the found gene is not already in the list
                if gene_name not in genes_list :
                    genes_list.append( gene_name )

                # Search for the genes included in the interval
                while True :
                    i -= 1
                    if i < 0 :
                        break

                    gene_name = a_gene_ID[chrom][i]
                    if gene_name not in genes_list :
                        if a_gene_end[chrom][i] - island_start >= gene_min_len :

                            if a_gene_start[chrom][i] < island_start :
                                gene_name += '*'
                                genes_list.append( gene_name )
                                break

                            genes_list.append( gene_name )


        #================================================================================
        # Result
        #================================================================================
        # Print the line  adding 2 suplementary columns
        genes_list.sort()
        print( line.rstrip(), "\t", len(genes_list), "\t", ','.join(gene for gene in genes_list) )


