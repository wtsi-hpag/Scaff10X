#tracking 

# Scaff10X v5.0
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
Requirements for compiling: gcc gcc-4.9.2 or late:

If you see this message,
cc1: error: unrecognised command line option ‘-std=c11’
make: *** [break10x.o] Error 1

you need a higher version of gcc
CC= /software/gcc-4.9.2/bin/gcc in the makefile


    $ git clone  https://github.com/wtsi-hpag/Scaff10X.git 
    $ cd Scaff10X
    $ ./install.sh
		
If everything compiled successfully you must see the final comment: 
		"Congrats: installation successful!"		


#### External packages
The genome aligner BWA (http://bio-bwa.sourceforge.net) and SMALT (http://www.sanger.ac.uk/science/tools/smalt-0) are downloaded and compiled by Scaff10X.

### Run the pipelines

#### Run scaff10x:
           $ /full/path/to/Scaff10X/src/scaff10X -nodes <nodes> -align <aligner> -score <score> \
	   	 -matrix <matrix_size> -read-s1 <min_reads_s1> -read-s2 <min_reads_s2> -longread <aggressive> -gap <gap_size> \
		 -edge <edge_len> -link-s1 <n_links_s1> -link-s2 <n_links_s2> -block <block>  \
		 [ -data input.dat ] [ -sam input.sam ] [ -bam input.bam ] -cram input.cram   \
                 [ -htag ema ] [ -plot barcode-length.png ] \
		 draft-assembly.fasta output_scaffolds.fasta
           

	       Parameters:
             nodes:        number of CPUs requested  [ default = 30 ]
             score:        averaged mapping score on each barcode fragment [ default = 20 ]
             aligner:      sequence aligner: bwa or smalt [ default = bwa ]
             matrix_size:  relation matrix size [ default = 2000 ]
             min_reads_s1: step 1: minimum number of reads per barcode [ default = 10 ]
             min_reads_s2: step 2: minimum number of reads per barcode [ default = 10 ]
             edge_len:     length of mapped reads to consider for scaffolding [ default = 50000 ]
             n_links_s1:   step 1: minimum number of shared barcodes [ default = 8 ]
             n_links_s2:   step 2: minimum number of shared barcodes [ default = 8 ]
             aggressive:   1 - aggressively mapping filtering on small PacBio/ONT contigs; 
	     		   0 - no aggressive for short read assembly  [ default = 1 ]
             block:        length to determine for nearest neighbours [ default = 50000 ]
             gap:          gap size in building scaffold [ default = 100 ]
             htag:         ema - cram file from ema linked-read alignment
                           longranger - cram file from longranger linked-read alignment 
             plot:         output image file with barcode length distributions and coverage stats 
	     
	       Files
                ==========
	        input.dat:   input a text file to point the locations of the reads in paired files

q1=/lustre/scratch116/vr/projects/Tes1_S1_L008_R1_001.fastq.gz \
q2=/lustre/scratch116/vr/projects/Tes1_S1_L008_R2_001.fastq.gz \
q1=/lustre/scratch116/vr/projects/Tes1_S2_L008_R1_001.fastq.gz \
q2=/lustre/scratch116/vr/projects/Tes1_S2_L008_R2_001.fastq.gz \
q1=/lustre/scratch116/vr/projects/Tes1_S3_L008_R1_001.fastq.gz \
q2=/lustre/scratch116/vr/projects/Tes1_S3_L008_R2_001.fastq.gz \
q1=/lustre/scratch116/vr/projects/Tes1_S4_L008_R1_001.fastq.gz \
q2=/lustre/scratch116/vr/projects/Tes1_S4_L008_R2_001.fastq.gz \
 
		             The scaff10x pipeline will read the gzipped files, trim the barcodes and pipe to bwa for alignment	
                             There will be no sam file anymore in this way. The input.dat can be local or with full path

                ==========
	        input.sam:   input a sam file if it already exists, 
				and skip the mapping (Optional, please provide full path)

                ==========
	        input.bam:   input a bam file which had been produced by using lariat in longranger, 
				and skip the mapping (Optional, please provide full path)
				(a). rename the assembly file (Optional):
				$ /full/path/to/Scaff10X/src/scaff-bin/scaff_rename your_assembly.fa longrang_refasm.fa
				(b). generate reference assembly file using longranger
				$ longranger mkref longrang_refasm.fa 
				(c). align 10x reads using lariat
				$ longranger align --fastq="reads_10x" --sample=fTakRub1 --reference="longrang_refasm" --localcores=50 --id=10x-align 
				Note: please provide full path
				(d). run scaff10x 
				$ scaff10X -bam /lustre/scratch117/possorted_bam.bam draft-assembly.fasta output_scaffolds.fasta 
	        draft-assembly.fasta:   initial draft assembly to scaffold (full path or local)
	        output_scaffolds.fasta:   name for the output scaffolded assembly (local)

Some notes and suggestions:
            
	a. SMALT is notably slower than BWA. So we suggest to try BWA first;
	b. The block value is very important. The default value of 2500 is very conservative
	   and you may increase this value to say 5000 or 10000 to improve the length of scaffolds; 
	c. The default numbers of -reads and -link are based on 30X read depth. 
	   These values should be increased if the read depth is higher
	d. Alignments with mapping score < score are filtered out to reduce linking errors;
	e. By using the option of "-longread 1", the pipeline performs an aggressive 
	   mapping score filtering on small PacBio/ONT contigs.
	f. An image png file is generated showing distributions of barcode length with "-plot " option 
	   It also compares with Human, Hummingbird, fish fAnaTes1 and fish fSimDai1 
	   Human, Hummingbird and fish fAnaTes1 are in good quality, while fSimDai1 is a failed sample. 
	g. File cover.dat can be produced when you use "-plot " option 
	   This file provides barcode information (length and coverage) as well as sequence coverage 

#### Remember: you only need to run scaff10x once (previously we suggested two iterations)


#### Run break10x:
           
	   $ /full/path/to/Scaff10X/src/break10x -nodes <nodes>  -score <score> -reads <min_reads> \
		-gap <gap_size> -cover <cover> -ratio <ratio> -data <input.dat>\
		scaffolds.fasta scaffolds-break.fasta scaffolds-break.name	     
	    

	       Parameters:
	         input.dat:  input a text file to point the locations of the reads in paired files, see the file format for scaff10x
             nodes:      number of CPUs requested  [ default = 30 ]
             score:      minimum average mapping score on an area covered by reads with 
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

#### Other applications:
           
##### Process barcodes and output paired gzip files
    
	   $ /full/path/to/Scaff10X/src/scaff_reads -nodes <nodes>  input.dat genome-BC_1.fastq.gz genome-BC_2.fastq.gz  \
              Here input.dat is a text file to point the locations of the reads in paired files  \
              The paired reads can also be used for other applications such as genomeScope or jellyfish etc
	    
##### Run scaff10x in another way   

           $ /full/path/to/Scaff10X/src/scaff10X -nodes <nodes> -align <aligner> -score <score> \
	   	 -matrix <matrix_size> -read-s1 <min_reads_s1> -read-s2 <min_reads_s2> -longread <aggressive> -gap <gap_size> \
		 -edge <edge_len> -link-s1 <n_links_s1> -link-s2 <n_links_s2> -block <block>  \
		 draft-assembly.fasta genome-BC_1.fastq.gz genome-BC_2.fastq.gz output_scaffolds.fasta \
	    
##### Run break10x in another way   

	   $ /full/path/to/Scaff10X/src/break10x -nodes <nodes>  -score <score> -reads <min_reads> \
		-gap <gap_size> -cover <cover> -ratio <ratio> \
		scaffolds.fasta genome-BC_1.fastq.gz genome-BC_2.fastq.gz scaffolds-break.fasta scaffolds-break.name	     
	    
##### Haplotagging data   

	   A new application has been added to process haplotagging data:
 
	   $ /full/path/to/Scaff10X/src/scaff10x -cram /lustre/scratch117/sciops/team117/hpag/ema/htag.cram -htag ema -plot ema-length.png 
             ref.fa ema-assembly.fasta > try.out \

           htag.cram          - a cram file with longranger or ema linked-read alignment
           ema-length.png     - barcode length distribution image
           ref.fa             - the reference assembly used for alignment
           ema-assembly.fasta - output scaffolded assembly

	   Here barcode sequences with 24 bases have been converted into 16 bases
           longranger is the default aligner
           Please contact us when the cram file can not be read when reference is not properly set up.  

##### Scaff10X v5.0 release note   
	   In this new release, the scaffolding pipeline has been modified to cope with genomes > 40Gb. 
           A new input parameter "-size" has been added. The default genome size is 10Gb.
           With huge genomes, say ~30Gb, you need to use "-size 30.0" in order for the code to process the data more effectively.  

           $ /full/path/to/Scaff10X/src/scaff10X -nodes <nodes> -align <aligner> -score <score> \
		 -size 30.0 \
		 [ -data input.dat ] [ -sam input.sam ] [ -bam input.bam ] -cram input.cram   \
                 [ -htag ema ] [ -plot barcode-length.png ] \
		 draft-assembly.fasta output_scaffolds.fasta
           

