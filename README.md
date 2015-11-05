Introduction
------------

SWeBLAST (Sliding Window wEb-based BLAST) is a cross-platform command line program written in Python.
**SWeBLAST.py** splits one or more sequences into sub-sequences using a sliding window and blasts them against a local or online (NCBI) databases.
Sub-sequences are saved in *seq.fa* and BLAST results are saved into blast.txt. The *blast.txt* file is parsed by **SWeBLAST_parser.py**, creating some summary files in comma-separated files.

By default **SWeBLAST.py** use BLAST+ and uses the `-remote` option in order to query an NCBI database.

Requirements
------------
 * SWeBLAST is a set Python scripts and requires the [**requests**](https://github.com/kennethreitz/requests) module. These scripts were developed under Python 2.7.
 * [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) should be installed locally unless you want to you use the old BLAST API.

Quick start
-----------
Download the SWeBLAST.zip file and extract in a single folder, add the FASTA file containing sequence(s).
Open a command line window and navigate to the folder containing the SWeBLAST scripts.

To look at the available arguments or options, one can run the program using the option `--help`

````
python SWeBLAST.py --help
usage: SWeBLAST.py [-h] -i FILE -o FOLDER
                   [-p {blastn,blastp,blastx,megablast,rpsblast}]
                   [-d DATABASE] [-w WINDOW] [-s STEP] [-t THRESHOLD] [-l]
                   [-W] [-r] [-e EXPECT] [-z ENTREZ]

Cut sequences in a series of windows and BLAST them using BLAST+.

optional arguments:
  -h, --help   show this help message and exit
  -i FILE, --input FILE   an input sequence file in FASTA format.
  -o FOLDER, --output FOLDER   an output folder.
  -p {blastn,blastp,blastx,megablast,rpsblast}, --program {blastn,blastp,blastx,megablast,rpsblast}   a blast program.
  -d DATABASE, --database DATABASE   a blast database.
  -w WINDOW, --window WINDOW   window length.
  -s STEP, --step STEP   step length.
  -t THRESHOLD, --threshold THRESHOLD   hits with an e-value greater than the threshold will be discarded.
  -l, --local   use a local database.
  -W, --WWW   use old BLAST API instead of BLAST+.
  -r, --parse   call SWeBLAST_parser.
  -e EXPECT, --expect EXPECT   cutoff for the expected value.
  -z ENTREZ, --entrez ENTREZ   ENTREZ query.
````

Example:
````
python SWeBLAST.py -p blastn -d nr -i dna.fasta -o myRESULT -w 100 -s 50 -r -e 1.0
````

This will automatically start the parser and produce a *blast.txt* and *seq.fa* files in the folder *myRESULT*.
The parser may be re-run using a smaller cutoff with the following command:

````
python SWeBLAST_parser.py -i myRESULT/blast.txt -e 0.01
````

Note: The Entrez query must be written between double quotes as well as any arguments containing spaces such as the input file path.

Output files
------------

 * *seq.fas*: A Fasta format file of the sub-sequences that were sent to BLAST as query sequences.

 * *blast.txt*: A file containing a concatenation of all the BLAST output files

 * *evalue.csv*: A table of the e-values (exponential term only) for all the reported matches found by the BLAST program. Each row of this table corresponds to a sequence in the Genbank database found by BLAST to have matched with any of the submitted subsequences. Column A lists the matched sequences, and the other columns in the table correspond to the successive slices of the query sequence. In the body of the table are the reported e-values.

 * *rank.csv*: This table is like *evalue.csv*, but records the rank position of each match in the BLAST results for each sub-sequence submitted to BLAST.

 * *stat.csv*: The rows in this table also correspond to each database sequence found to match the query sequence, and these are in the same order as in the ‘evalue.csv’ and ‘rank.csv’ files. The first column of the file again records the names of the matched sequences, the second column records the number of matches obtained with all slices of the query sequence, and subsequent columns contain a series of statistics of the evalues and rankings obtained by the database sequence with every slice of the query sequence; their minima, maxima, ranges, medians, means and standard deviations.

 * *summary.csv*: This table differs from the others in that each row contains the results for one subsequence of the query sequence. Column A gives the name and position of the slice, and successive columns then give the names of the sequences with which it matched and the evalue of that match.

 * *frame.csv*: Generated only when blastx is used. A table of the frames for all the reported matches found by the BLAST program. Each row of this table corresponds to a sequence in the Genbank database found by BLAST to have matched with any of the submitted subsequences. Column A lists the matched sequences, and the other columns in the table correspond to the successive slices of the query sequence. In the body of the table are the reported frames.

Reference
---------
Fourment M, Gibbs AJ and Gibbs MJ. SWeBLAST: A Sliding Window Web-based BLAST tool for recombinant analysis. [J Virol Methods](http://www.sciencedirect.com/science/article/pii/S0166093408002164), 2008, 152, 98­101.
