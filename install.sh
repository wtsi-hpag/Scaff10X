#!/bin/bash


projdir=`pwd`

bindir=$projdir/src/scaff-bin/
mkdir -p $bindir
mkdir -p $projdir/src/log/

errs=0

##### Download and install BWA ######

echo "Downloading and installing BWA"
if [[ ! -s $bindir/bwa ]]; then

    if [[ ! -d $projdir/src/bwa ]]; then
	cd $projdir/src/
	git clone https://github.com/lh3/bwa.git &> $projdir/src/log/bwa_cloning.log
    fi

    if [[ ! -s $projdir/src/bwa/bwa ]]; then
	cd $projdir/src/bwa
	make &> $projdir/src/log/bwa_installation.log
    fi

    cp bwa $bindir
fi

if  [[ ! -s $bindir/bwa ]]; then
    echo " !! Error: bwa not installed properly!"; 
    echo "   Please check the log files:" 
    echo "   Check if bwa was downloaded properly:" $projdir/src/log/bwa_cloning.log 
    echo "   Check if the bwa was compiled properly:" $projdir/src/log/bwa_installation.log

    # Cleaning up
    cd $projdir/src
    rm -rf $projdir/src/bwa/bwa $bindir/bwa 
    
    errs=$(($errs+1))
else
    echo " BWA succesfully installed!"
    rm -rf $projdir/src/bwa/
fi


##### Download and install SMALT ######
echo; echo "Downloading and installing Smalt"
if [[ ! -s $bindir/smalt ]]; then
   
    if [[ ! -d $projdir/src/smalt-0.7.4 ]]; then
	cd $projdir/src/
	wget ftp://ftp.sanger.ac.uk/pub/resources/software/smalt/smalt-0.7.4.tgz &> $projdir/src/log/smalt_wget.log
	tar -xvzf smalt-0.7.4.tgz &> $projdir/src/log/smalt_untar.log
	rm -f smalt-0.7.4.tgz
    fi

    cp $projdir/src/smalt-0.7.4/smalt_x86_64 $bindir/smalt
fi
if  [[ ! -s $bindir/smalt ]]; then 
    echo " !! Error: smalt not installed properly!"; 
    echo "   Please check the log files:" 
    echo "   Check if smalt was downloaded properly:" $projdir/src/log/smalt_wget.log 
    echo "   Check if the folder was uncompressed properly:" $projdir/src/log/smalt_untar.log

    # Cleaning up	
    rm -rf $projdir/src/smalt-0.7.4/ $bindir/smalt

    errs=$(($errs+1))
else
    echo " Smalt succesfully installed!"
    rm -rf $projdir/src/smalt-0.7.4/
fi


##### Download and install pigz ######

echo "Downloading and installing pigz"
if [[ ! -s $bindir/pigz ]]; then

    if [[ ! -d $projdir/src/pigz ]]; then
	cd $projdir/src/
        wget -r -np -nd https://zlib.net/pigz/pigz-2.4.tar.gz &> $projdir/src/log/pigz_wget.log
        tar -xvzf pigz-2.4.tar.gz &> $projdir/src/log/pigz_untar.log
        rm -f pigz-2.4.tar.gz
    fi

    if [[ ! -s $projdir/src/pigz/pigz ]]; then
	cd $projdir/src/pigz-2.4
	make &> $projdir/src/log/pigz_installation.log
    fi

    cp pigz $bindir
fi

if  [[ ! -s $bindir/pigz ]]; then
    echo " !! Error: pigz not installed properly!"; 
    echo "   Please check the log files:" 
    echo "   Check if bwa was downloaded properly:" $projdir/src/log/pigz_cloning.log 
    echo "   Check if the bwa was compiled properly:" $projdir/src/log/pigz_installation.log

    # Cleaning up
    cd $projdir/src
    rm -rf $projdir/src/pigz/pigz $bindir/pigz 
    
    errs=$(($errs+1))
else
    echo " pigz succesfully installed!"
    rm -rf $projdir/src/pigz/
fi

###### Compile Scaff10x sources ######

echo; echo "Compiling scaff10X sources"

srcs=( break10x scaff_barcode-cover scaff_barcode-sort scaff_break-clean scaff_contigs-sort scaff_mapping-clean scaff_outbreak scaff_PCRdup scaff_samout scaff_BC-reads-1 scaff_break-names scaff_fastq scaff_mapping-sort scaff_outbreak-seq scaff_reads scaff_samprocess scaff10x scaff_barcode-screen scaff_BC-reads-2 scaff_bwa scaff_length scaff_matrix scaff_output scaff_rename scaff_ctgloci-sort scaff_lengthdis scaff_agp2agp scaff_RDplace scaff_superAGP scaff_FilePreProcess)

cd $projdir/src
make &> $projdir/src/log/sources_compilation.log

echo; echo "Checking installation:"
for src in "${srcs[@]}"; do
    if [[ ! -s $bindir/$src ]]; then 
        echo " !! Error: executable $src missing in $bindir"
	echo "    Please check for errors the log file:" $projdir/src/log/sources_*	
        errs=$(($errs+1))
    fi
done

if [  $errs -gt 0 ]; then echo; echo " ****  Errors occurred! **** "; echo; exit; 
else echo " Congrats: installation successful!"; fi




