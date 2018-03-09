# Scaff10X
Pipeline for scaffolding and breaking a genome assembly using 10x genomics linked-reads

Pipeline steps:

    1 Barcoded tags are extracted from 10Xg raw sequencing reads and appended to read names for further processing
    2 The reads are mapped to the draft assembly using either BWA or SMALT
    3 Barcodes are sorted together with contigs as well as mapping coordinates
    4 A relation matrix is built to record the shared barcodes among the contigs which may be linked
    5 Order and orientation of linked contigs are determined after nearest neighbours are found. 


### Download and Compile:
Requirements for compiling: gcc

    $ git clone  https://github.com/wtsi-hpag/Scaff10X.git 
    $ cd Scaff10X
    $ ./install.sh
		
If everything compiled successfully you must see the final comment: 
		"Congrats: installation successful!"		

(Tested with gcc-4.9.2)


#### External packages
The genome aligner BWA (http://bio-bwa.sourceforge.net) and SMALT (http://www.sanger.ac.uk/science/tools/smalt-0) are downloaded and compiled by Scaff10X.


#### Run:
           $ /full/path/to/Scaff10X/src/scaff10X -nodes <nodes> -align <aligner> -score <score> -matrix <matrix_size> -reads <min_reads> -longread <aggressive> -gap <gap_size> -edge <edge_len> -link <n_links> -block <block>  draft-asssembly.fasta read-BC_1.fastq read-BC_2.fastq output_scaffolds.fasta
           

	       parameters:
             nodes:    number of CPUs requested  [ default = 30 ]
             score: averaged mapping score on each barcode fragment [ default = 20 ]
             aligner:  sequence aligner: bwa or smalt [ default = bwa ]
             matrix_size:   relation matrix size [ default = 2000 ]
             min_reads:  minimum number of reads per barcode [ default = 10 ]
             edge_len e length of mapped reads to consider for scaffolding [ default = 50000 ]
             link:     - minimum number of shared barcodes [ default = 8 ]
             aggressive: 1 - aggressively mapping filtering on small PacBio/ONT contigs; 0 - no aggressive for short read assembly  [ default = 1 ]
             block  (2500)  - length to determine for nearest neighbours [ default = 2500 ]
             gap    (100)   - gap size in building scaffold [ default = 100 ]
