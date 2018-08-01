# Scaff10X
Pipeline for scaffolding and breaking a genome assembly using 10x genomics linked-reads.

Pipeline steps:
        
    Scaffolding with scaff10x:
      1 Barcoded tags are extracted from 10Xg raw sequencing reads and appended 
          to read names for further processing
      2 The reads are mapped to the draft assembly using either BWA or SMALT
      3 Barcodes are sorted together with contigs as well as mapping coordinates
      4 A relation matrix is built to record the shared barcodes among the contigs which may be linked
      5 Order and orientation of linked contigs are determined after nearest neighbours are found. 
      
    Breaking scaffolds (or contigs) with break10x:
      Steps 1, 2 and 3 are the same as for scaff10x
      4 The barcode depth is monitored for each scaffold, if this is reduced by more than 15% with 
      	respect to the scaffold average barcode depth, the point of minimum depth is considered to 
		be an assembly error and the scaffold is broken at that point. Exception are points close to 
		long N-gaps unless the N-Gaps was added by scaff10x.
      5 The location of breaking points are printed out in a file and the broken scaffolds are 
      	printed in a fasta file
      

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

### Run the pipelines

#### Prepare the fastq files:
#### (1) Make a file with read files
#### Say you have a number of R1, R2 files, you make a file named as file.dat

q1=/lustre/scratch116/vr/projects/Tes1_S1_L008_R1_001.fastq.gz

q2=/lustre/scratch116/vr/projects/Tes1_S1_L008_R2_001.fastq.gz

q1=/lustre/scratch116/vr/projects/Tes1_S2_L008_R1_001.fastq.gz

q2=/lustre/scratch116/vr/projects/Tes1_S2_L008_R2_001.fastq.gz

q1=/lustre/scratch116/vr/projects/Tes1_S3_L008_R1_001.fastq.gz

q2=/lustre/scratch116/vr/projects/Tes1_S3_L008_R2_001.fastq.gz

q1=/lustre/scratch116/vr/projects/Tes1_S4_L008_R1_001.fastq.gz

q2=/lustre/scratch116/vr/projects/Tes1_S4_L008_R2_001.fastq.gz

	$ /full/path/to/Scaff10X/scaff_reads file.dat reads_1.fastq reads_2.fastq > try.out

        where: file.dat is the file pointing read locations; reads_1.fastq reads_2.fastq are output read files with the barcodes cut out of the reads and attached to the names.

#### Or (2) Process files one by one

	$ /full/path/to/Scaff10X/src/scaff_BC-reads-1 read_1.fastq read-BC_1.fastq read-BC_1.name 
	$ /full/path/to/Scaff10X/src/scaff_BC-reads-2 read-BC_1.name read_2.fastq read-BC_2.fastq 
	
	where:
		read_1.fastq and read_2.fastq:   the raw 10Xg read fastq files
		read-BC_1.fastq and read-BC_2.fastq: the output fastq files to use with scaff10x/break10x.

Warning: for multiple runs prepare a read-BC_1.fastq and read-BC_2.fastq pair for each run and then cat them together.
This is to avoid using scaff_BC-reads-1 and scaff_BC-reads-2 on files with too many reads.
		

#### Run scaff10x:
           $ /full/path/to/Scaff10X/src/scaff10X -nodes <nodes> -align <aligner> -score <score> \
	   	 -matrix <matrix_size> -reads <min_reads> -longread <aggressive> -gap <gap_size> \
		 -edge <edge_len> -link <n_links> -block <block>  \
		 [ -sam input.sam ] \
		 draft-asssembly.fasta read-BC_1.fastq read-BC_2.fastq output_scaffolds.fasta
           

	       Parameters:
             nodes:    number of CPUs requested  [ default = 30 ]
             score: averaged mapping score on each barcode fragment [ default = 20 ]
             aligner:  sequence aligner: bwa or smalt [ default = bwa ]
             matrix_size:   relation matrix size [ default = 2000 ]
             min_reads:  minimum number of reads per barcode [ default = 10 ]
             edge_len:   length of mapped reads to consider for scaffolding [ default = 50000 ]
             n_links:      minimum number of shared barcodes [ default = 8 ]
             aggressive:   1 - aggressively mapping filtering on small PacBio/ONT contigs; 
	     		       0 - no aggressive for short read assembly  [ default = 1 ]
             block:    length to determine for nearest neighbours [ default = 2500 ]
             gap:     gap size in building scaffold [ default = 100 ]
	     
	       Files
	        input.sam:   input a sam file if it already exists, 
				and skip the mapping (Optional, please provde full path)
	        draft-asssembly.fasta:   initial draft assembly to scaffold (full path or local)
	        read-BC_1.fastq read-BC_2.fastq:  10Xg reads with barcode appended 
						 to read names, prepared as shown above (full path or local)
	        output_scaffolds.fasta:   name for the output scaffolded assembly (local)

Some notes and suggestions:
            
	a. SMALT is notably slower than BWA. So we suggest to try BWA first;
	b. The block value is very important. The default value of 2500 is very conservative
	   and you may increase this value to say 5000 or 10000 to improve the length of scaffolds; 
	c. The default numbers of -reads and -link are based on 30X read depth. 
	   These values should be increased if the read deapth is higher
	d. To get the best results we suggest an iterative approach, 
	   where scaff10x scaffolds (output_scaffolds.fasta) are scaffolded again using scaff10x.
	e. Alignments with mapping score < score are filtered out to reduce linking errors;
	f. By using the option of "-longread 1", the pipeline performs an aggressive 
	    mapping score filtering on small PacBio/ONT contigs.  


#### Run break10x:
           
	   $ /full/path/to/Scaff10X/src/break10x -nodes <nodes>  -score <score> -reads <min_reads> \
		-gap <gap_size> -cover <cover> -ratio <ratio> \
		scaffolds.fasta read-BC_1.fastq read-BC_2.fastq scaffolds-break.fasta scaffolds-break.name	     
	    

	       Parameters:
             nodes:    number of CPUs requested  [ default = 30 ]
             score:    minimum average mapping score on an area covered by reads with 
	     		the same barcode [ default = 20 ]
             min_reads:  minimum number of reads per barcode [ default = 5 ]
	         cover: minimum barcode coverage at the breakpoint [ default = 50 ]
	         gap:  scaffold gap size added by scaff10x (if used). 
	     	   If a breakpoint is close to a gap region, break10x checks if the gap was added by  
		   scaff10x when joining two contigs (using the same 10X data). 
		   If it was, the scaffold is broken. If the scaffolding was done by other 
		   means (different length from scaff10x gap added), it will not be broken as the 
		   gap could be very big and 10X barcodes might not cross over. [ default = 100 ]
        
##### break10x output
The command break10x will output the assembly after scaffold breaking (scaffolds-break.fasta) and
the list of breaking points.

The list of breaking point  in scaffolds-break.name has the format:

	Break2: Scaffold_name  Scaffold_index  Position Scaffold_Length Minimum_Coverage Average_Coverage  Gap_Length
	
	Scaffold_name   name of scaffold broken by break10x
	Scaffold_index  index of scaffold broken by break10x
	Position        Break-point position
	Scaffold_length  length of scaffold
	Minimum_Coverage  minimum barcode coverage in the breakpoint region
	Average_Coverage  scaffold average barcode coverage
	Gap_Length        gap length ( 0 indicates a contig break)

Example:
	
	Break2: Scaff10x_0 0 1527414 15892392 1 114 100
	Break2: Scaff10x_0 0 7246257 15892392 0 114 100
	Break2: Scaff10x_2 2 17232850 17798733 21 176 100
	Break2: Scaff10x_9 9 3709193 3904076 37 451 100
	Break1: Scaff10x_10 10 3713358 5150262 0 131 0
