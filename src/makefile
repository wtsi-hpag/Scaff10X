# Makefile for scaff10x 
CC= gcc
CFLAGS= -O2
LFLAGS= -lm

SOURCES=break10x.c  scaff_barcode-cover.c  scaff_barcode-length.c  scaff_barcode-sort.c  scaff_break-clean.c  scaff_contigs-sort.c  scaff_mapping-clean.c  scaff_outbreak.c scaff_PCRdup.c  scaff_samout.c  scaff_BC-reads-1.c scaff_break-names.c  scaff_fastq.c scaff_mapping-sort.c scaff_outbreak-seq.c  scaff_reads.c scaff_samprocess.c scaff10x.c  scaff_barcode-screen.c scaff_BC-reads-2.c scaff_bwa.c scaff_bwa-barcode.c scaff_length.c scaff_matrix.c scaff_output.c scaff_rename.c scaff_BClen-dis.c scaff_BCnum-dis.c scaff_superAGP.c


OBJS = $(patsubst %.c,%.o,$(SOURCES)) fast.o
EXECS = $(patsubst %.c,%,$(SOURCES))
EXECS_BIN = $(patsubst %.c,scaff-bin/%,$(SOURCES))
COMPILE = $(CC) $(CFLAGS)


all:  cleanall iprint $(OBJS) executables clean oprint

executables:
	for exe in $(EXECS);  do $(COMPILE) -o $$exe $$exe.o fast.o; cp $$exe scaff-bin/.; done

%.o: %.c fasta.h
	$(CC) $(CFLAGS)  -c $< -o $@

iprint:
	@echo '+++ Compiling All ... '

oprint:
	@echo 'All Done '


clean: 
	@echo '+++ Cleaning Up ... '
	@rm -f $(EXECS)
	@rm -f $(OBJS)
	@cp scaff-bin/scaff10x .
	@cp scaff-bin/break10x .

cleanall: 
	@echo '+++ Cleaning All ... '
	@rm -f $(EXECS)
	@rm -f $(OBJS) fast.o
	@rm -f $(EXECS_BIN)
