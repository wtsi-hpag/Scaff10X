/***********************************************************************\
 *                                                                     * 
 *                      PROJECT        PHUSION                         *
 *                                                                     * 
 *---------------------------------------------------------------------*
 *                                                                     * 
 *             A Whole Genome Shotgun Clustering Algorithm             *
 *                                                                     *
 *                                By                                   *
 *                                                                     *
 *                    Zemin Ning & James C. Mullikin                   *
 *                                                                     *
 *          Copyright (C) 2001-2002 by Genome Research Limited         *
 *                         All rights reserved                         *
 *                                                                     *
 *---------------------------------------------------------------------*
 #######################################################################
 # This software has been created by Genome Research Limited (GRL).    # 
 # GRL hereby grants permission to use, copy, modify and distribute    # 
 # this software and its documentation for non-commercial purposes     # 
 # without fee at the user's own risk on the basis set out below.      #
 # GRL neither undertakes nor accepts any duty whether contractual or  # 
 # otherwise in connection with the software, its use or the use of    # 
 # any derivative, and makes no representations or warranties, express #
 # or implied, concerning the software, its suitability, fitness for   #
 # a particular purpose or non-infringement.                           #
 # In no event shall the authors of the software or GRL be responsible # 
 # or liable for any loss or damage whatsoever arising in any way      # 
 # directly or indirectly out of the use of this software or its       # 
 # derivatives, even if advised of the possibility of such damage.     #
 # Our software can be freely distributed under the conditions set out # 
 # above, and must contain this copyright notice.                      #
 #######################################################################
 *---------------------------------------------------------------------*/


/****************************************************************************/

 
#include <math.h>
#include <values.h>
#include <stdio.h>
#include <netinet/in.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <ctype.h>
#include "fasta.h"
#include "ssaha.h"


static int *list,*readlength,*out_list,*ctg_list,*rd_group;
static long *hist;
static long n_Entry;
static int K2_LEN=17;
static int mhist=10;
static int mhist2=2;
static int breakmod=0;
static int NEDGE=250;
static int nmatch=30;
static int nmatch2=30;
static int N_SET=2000;
static int N2_SET=12000;
static char *SCGname=NULL;//SingleCopyGenome
static long baseSize;
static int n_grouped;
static int fastq_flag=1;
static int n_maxblock;
static int n_blockreads;
static int n_block;
static int Max_N_NameBase=60;

static void hpsortl(long n, unsigned long *ra)
{
	long l,j,ir,i;
	unsigned long rra;

	if(n < 2)
          return;
        l=(n >> 1)+1;
	ir=n;
	for (;;) {
		if (l > 1) {
			rra=ra[--l];
		} else {
			rra=ra[ir];
			ra[ir]=ra[1];
			if (--ir == 1) {
				ra[1]=rra;
				return;
			}
		}
		i=l;
		j=l << 1;
		while (j <= ir)	{
			if (j < ir && ra[j] < ra[j+1]) ++j;
			if (rra < ra[j]) {
				ra[i]=ra[j];
				j += (i=j);
			}
			else j=ir+1;
		}
		ra[i]=rra;
	}
}


/*   Subroutine to sort the DNA sequences into a matrix table   */
/* =============================================  */
void fastaSort(char **argv, int argc, int SEG_LEN)
/* =============================================  */
{
     long i,j,k,n_Sbase=SEG_LEN;
     long ps=0,ns=0,IntSeg=0,IntBase=0,lii;
     fasta *seq;
     fast  *seq2;
     int nSeq;
     long totalBases;
     long mhistc = 0;
     long mhistcc = 0;
     unsigned long *sarr, *sarrp;
     int nsorts = 1024;
     int nshift = (n_Sbase<<1)-10;//2^10=nsorts
     int qthresh=23;
//     int tb = 0;
     int args,ac;
     int m_st,m_ed,step_len,pmod;
     int nseq = 0,seqc,n_patch;
     int sshift,rshift=6;
     unsigned long nmask;
     int id_read[64],num_sect,n_reads,max_edge; //id_read will never exceed mhist
     int rd_st,rd_ed,n_input;
     unsigned int **imatrix(long nrl,long nrh,long ncl,long nch);
     char line[500] = {0},outName[60]={0},outFast[60]={0},syscmd[200]={0};
     char *ptr,base[20],zero[20]={0},line2[500]={0};
     FILE *fpMate,*fpOutname,*fpOutfast;
     int rd,read_pair[200];
     long kmask = (1L<<(n_Sbase<<1))-1;
     double Qerr[100];
     long gtBases=0,nclip=0;
     Qerr[0] = 1.0;
     void Phusion_Stage(int brr,int len,int crr,int m1,int m2,int step,int pmod, char **argv, int argc, int args);

     for(i=1;i<100;i++) Qerr[i] = pow((double)10.,(double)-.1*i);

/*   sort all the names of genes or name entries   */
     printf("Input data starts ...\n");
     fflush(stdout);
     system("ps aux | grep super2contig; date");

     args=1;
     for(i=1;i<argc;i++)
     {
       if(!strcmp(argv[i],"-kmer"))
       {
	 i++;
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-kmer2"))
       {
	 sscanf(argv[++i],"%d",&K2_LEN); 
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-depth"))
       {
         sscanf(argv[++i],"%d",&mhist);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-depth2"))
       {
         sscanf(argv[++i],"%d",&mhist2);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-match"))
       {
         sscanf(argv[++i],"%d",&nmatch);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-match2"))
       {
         sscanf(argv[++i],"%d",&nmatch2);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-fastq"))
       {
         sscanf(argv[++i],"%d",&fastq_flag);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-break"))
       {
         sscanf(argv[++i],"%d",&breakmod);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-matrix"))
       {
         sscanf(argv[++i],"%d",&NEDGE);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-set"))
       {
         sscanf(argv[++i],"%d",&N_SET);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-set2"))
       {
         sscanf(argv[++i],"%d",&N2_SET);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-SCG"))
       {
	 SCGname = argv[++i];
         args=args+2;
       } 
     }
     for(ac=0;ac<argc;ac++)
	printf("%s ",argv[ac]);
     printf("\n");
     if(SCGname != NULL) printf("SCG name:          %s\n",SCGname);
     fflush(stdout);
     n_Entry=1L<<(n_Sbase<<1);

     printf("Input file: %s ",argv[1]);
     m_st=nmatch;
     m_ed=nmatch;
     step_len=2;
     pmod=1;

     nseq = 0;
     Phusion_Stage(nseq,SEG_LEN,mhist,m_st,m_ed,step_len,pmod,argv,argc,args);


     fflush(stdout);

     printf("All jobs finished\n");
     fflush(stdout);
     system("ps aux | grep super2contig; date");
     return;
}


int ssaha_init( char **argv, int argc, int SEG_LEN)
{

    fastaSort(argv,argc,SEG_LEN);
    return(1);   
}

/*   Subroutine to sort the DNA sequences into a matrix table   */
/* =============================================  */
void Phusion_Stage(int nRead,int SEG_LEN,int n_depth,int mat_st,int mat_ed,int s_len,int pmod, char **argv, int argc, int args)
/* =============================================  */
{
     long i,j,k,iseq,n_Sbase=SEG_LEN;
     fasta *seqp;
     fasta *seq;
     fast  *seq2;
     long totalBases,total_len;
     int nSeq,num_Ns,i_contig;
     int qthresh=23;
     int ac,stopflag,n_contig,offset,n_Ns;
     int *remove,*readIndex,*head,*list;
     char outName[2000],name_tag[10],namep_tag[10];
     void ArraySort_String(int n,char Pair_Name[][Max_N_NameBase],int *brr);
     char **cmatrix(long nrl,long nrh,long ncl,long nch);
     FILE *fpOutfast,*fpOutfast2,*namef;


     printf("Relation matrix finished: %d\n",args);
     

     printf("names: %s %d %d\n",argv[1],ac,argc);
     seq = decodeFastq ( argv[1], &nSeq, &totalBases, qthresh);
     if(seq == NULL)
     {
     	printf("ERROR pileup: no data found\n");
     	exit(1);
     }
     fastaUC(seq,nSeq);
     if((readIndex= (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("Error Contig_Merge: calloc - readIndex\n");
       exit(1);
     } 
     if((head= (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("Error Contig_Merge: calloc - head\n");
       exit(1);
     }
     if((list= (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("Error Contig_Merge: calloc - list\n");
       exit(1);
     }
                 

     if((fpOutfast = fopen(argv[2],"w")) == NULL)
     {
       printf("Unable to open file for fastq out\n");
       exit(1);
     }
     if((fpOutfast2 = fopen(argv[3],"w")) == NULL)
     {
       printf("Unable to open file for fastq out\n");
       exit(1);
     }

     n_contig = 0;
     seqp=seq;
     offset=0;
     total_len = 0;
     memset(namep_tag,'\0',10);
     strncpy(namep_tag,"NNNNNN",5);
     memset(name_tag,'\0',10);
     strncpy(name_tag,seqp->name,5);
     for(iseq=0;iseq<nSeq;iseq++)
     {
        int kk,rc,slength,start=0,nline=0,outlen=0,n_base,olen=60,offset2=0;
        char *st;

        seqp=seq+iseq;
        slength=seqp->length;
	memset(name_tag,'\0',10);
	strncpy(name_tag,seqp->name,5);
        if(strcmp(name_tag,namep_tag) != 0)
	{
	  n_contig = 0;
          printf("%s %s %d\n",namep_tag,name_tag,n_contig);
          memset(namep_tag,'\0',10);
	  strncpy(namep_tag,name_tag,5);
	}

        i_contig = 1;
        n_base=0;
        n_Ns=0;
        for(j=0;j<slength;j++) 
        {
           kk=j+1;
           n_Ns=1;
           while((seqp->data[kk]=='N')&&(seqp->data[j]=='N')&&(kk<slength))
           {
             kk++;
             n_Ns++;
           }
           if((j==0)&&(n_Ns>=3))
           {
//             printf("offset: %d %d %d\n",j,n_Ns);
             offset2=n_Ns;
           }
           if((n_Ns>=3))
           {
             offset=kk-1;
             sprintf(outName,"%s %s",seqp->name,seqp->name2);
//             printf("%scontig_%09d\n",name_tag,n_contig);
             outlen=offset-n_Ns+1-start-offset2;
//             fprintf(fpOutfast2,"%20s %10d %10d %5d   W    %20s %10d %10d    +\n",seqp->name,start+1,kk-n_Ns,i_contig,outName,1,outlen);
             fprintf(fpOutfast2,"%s\t%d\t%d\t%d\t%c\t%s\t%d\t%d\t%c\n",seqp->name,start+1,kk-n_Ns,i_contig,'W',outName,1,outlen,'+');
             i_contig++;
             fprintf(fpOutfast2,"%s\t%d\t%d\t%d\t%c\t%d\t%s\t%s\n",seqp->name,kk+1-n_Ns,kk,i_contig,'N',n_Ns,"fragment","yes");
//               fprintf(fpOutfast2,"contig %s\n",outName);
//               fprintf(fpOutfast2,"gap %d * * *\n",n_Ns);
//               fprintf(fpOutfast,">%s\n",outName);
               nline=outlen/olen;
               total_len = total_len+outlen;
               st=seqp->data+start+offset2;
             i_contig++;
             if(fastq_flag)
             {
               fprintf(fpOutfast,"@%s\n",outName);
               for(i=0;i<outlen;i++)
                  fprintf(fpOutfast,"%c",seqp->data[start+offset2+i]);
               fprintf(fpOutfast,"\n");
               fprintf(fpOutfast,"+\n");
               putc(0+041,fpOutfast);
               for(i=1;i<outlen;i++)
               {
                  putc((seqp->qual[start+offset2+i]+041),fpOutfast);
               }
               fprintf(fpOutfast,"\n");
             }
             else
             {
//               fprintf(fpOutfast,">%s %s %d %d %d %d %d\n",outName,seqp->name,iseq,offset,outlen,start+offset2,n_Ns);
               for(i=0;i<nline;i++)
               {
                  for(k=0;k<olen;k++,st++)
                     fprintf(fpOutfast,"%c",*st);
                  fprintf(fpOutfast,"\n");
               }
               for(k=0;k<(outlen-(nline*olen));k++,st++)
                  fprintf(fpOutfast,"%c",*st);
               if(outlen%olen!=0)
                 fprintf(fpOutfast,"\n");
             }
             start=offset+1;
             n_Ns=0;
             n_base=0;
             offset2=0;
             n_contig++;
           }
           j=kk-1; 
        }

        num_Ns=n_Ns-1;
        sprintf(outName,"%s %s",seqp->name,seqp->name2);
//        outlen=slength-start;
        outlen=slength-start-n_Ns+1;
        if(outlen >=10)
        {
//          fprintf(fpOutfast2,"%20s %10d %10d %5d   W    %20s %10d %10d    +\n",seqp->name,start+1,kk-n_Ns+1,i_contig,outName,1,outlen);
          fprintf(fpOutfast2,"%s\t%d\t%d\t%d\t%c\t%s\t%d\t%d\t%c\n",seqp->name,start+1,kk-n_Ns+1,i_contig,'W',outName,1,outlen,'+');
          i_contig++;
//          fprintf(fpOutfast2,"contig %s\n",outName);
//          fprintf(fpOutfast,">%s\n",outName);
          nline=outlen/olen;
          total_len = total_len+outlen;
//          printf("offset: %d %d %d %ld\n",start,outlen,n_Ns,total_len);
          if(fastq_flag)
          {
            fprintf(fpOutfast,"@%s\n",outName);
            for(i=0;i<outlen;i++)
               fprintf(fpOutfast,"%c",seqp->data[start+i]);
            fprintf(fpOutfast,"\n");
            fprintf(fpOutfast,"+\n");
            putc(0+041,fpOutfast);
            for(i=1;i<outlen;i++)
               putc(seqp->qual[start+i]+041,fpOutfast);
            fprintf(fpOutfast,"\n");
          }
          else
          {
            fprintf(fpOutfast,">%s %s %d %d %d %d %d\n",outName,seqp->name,iseq,offset,outlen,start,n_Ns-1);
            st=seqp->data+start;
            for(i=0;i<nline;i++)
            {
               for(k=0;k<olen;k++,st++)
                  fprintf(fpOutfast,"%c",*st);
               fprintf(fpOutfast,"\n");
            }
            for(k=0;k<(outlen-(nline*olen));k++,st++)
               fprintf(fpOutfast,"%c",*st);
            if(outlen%olen!=0)
              fprintf(fpOutfast,"\n");
          }
          n_contig++;
        }
//        fprintf(fpOutfast2,"\n");
     }
     fclose(fpOutfast);
     fclose(fpOutfast2);
     if(seq){
         free(seq->name);
         free(seq);
         seq = NULL;
     }

}


/* creat an int matrix with subscript ange m[nrl...nrh][ncl...nch]  */
unsigned int     **imatrix(long nrl,long nrh,long ncl,long nch)
{
        long i, nrow=((nrh-nrl+1)/100 + 1)*100,ncol=nch-ncl+1;
        unsigned int  **m;
	long nri,nrn=nrow/100;

        /* allocate pointers to rows        */
        if((m=(unsigned int **)calloc(nrow,sizeof(int*)))==NULL)
        {
           printf("error imatrix: calloc error No. 1 \n");
           return(NULL);
        }

        /* allocate rows and set pointers to them        */
	/* allocate in 100 batches to use freed memory */
	nrl = 0;
	for(nri=0;nri<100;nri++,nrl+=nrn) {
           if((m[nrl]=(unsigned int *)calloc(nrn*ncol,sizeof(int)))==NULL)
           {
              printf("error imatrix: calloc error No. 2 \n");
              return(NULL);
           }

           for(i=1;i<nrn;i++)
              m[i+nrl]=m[i+nrl-1]+ncol;
	}
       /* return pointer to array of pointers to rows   */
        return m;
}

#define SWAP(a,b) temp=(a);(a)=b;(b)=temp;

/*   Subroutine to sort an array arr[0,...,n-1] into ascending order while
     making the corresponding reaarangement of the array brr[0,...,n-1]
     by the use of Quicksort (Sedgwick, R. 1978, Communications o fthe ACM,
     vol. 21, pp. 847-857) also see Numerical Recipes in C                  */  



/*   to swap the string arrays           */
/* ============================================= */
void s_swap(char Pair_Name[][Max_N_NameBase], int i, int j)
/* ============================================= */
{
     char temp[Max_N_NameBase];

     strcpy(temp,Pair_Name[i]);
     strcpy(Pair_Name[i],Pair_Name[j]);
     strcpy(Pair_Name[j],temp);
}


/*   to sort the string array in order          */
/* ============================================= */
void ArraySort_String(int n, char Pair_Name[][Max_N_NameBase], int *brr)
/* ============================================= */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     int temp,MIN=7;
     char p[Max_N_NameBase];

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             strcpy(p,Pair_Name[j]);
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(strcmp(Pair_Name[i],p)<=0) break;
                strcpy(Pair_Name[i+1],Pair_Name[i]);
                brr[i+1]=brr[i];
             }
             strcpy(Pair_Name[i+1],p);
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          s_swap(Pair_Name,k,m+1);
          SWAP(brr[k],brr[m+1]);

          if(strcmp(Pair_Name[m],Pair_Name[ir])>0)
          {
            s_swap(Pair_Name,m,ir);
            SWAP(brr[m],brr[ir]);
          }

          if(strcmp(Pair_Name[m+1],Pair_Name[ir])>0)
          {
            s_swap(Pair_Name,m+1,ir);
            SWAP(brr[m+1],brr[ir]);
          }

          if(strcmp(Pair_Name[m],Pair_Name[m+1])>0)
          {
            s_swap(Pair_Name,m,m+1);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          strcpy(p,Pair_Name[m+1]);
          b=brr[m+1];
          for(;;)
          {
             do i++; while (strcmp(Pair_Name[i],p)<0);
             do j--; while (strcmp(Pair_Name[j],p)>0);
             if(j<i) break;
             s_swap(Pair_Name,i,j);
             SWAP(brr[i],brr[j]);
          }
          strcpy(Pair_Name[m+1],Pair_Name[j]);
          strcpy(Pair_Name[j],p);
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}


/* creat a char matrix with subscript ange c[nrl...nrh][ncl...nch]  */
char     **cmatrix(long nrl,long nrh,long ncl,long nch)
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        char  **m;

        /* allocate pointers to rows        */
        if((m=(char **)calloc(nrow,sizeof(char*)))==NULL)
        {
           printf("error cmatrix: calloc error No. 1 \n");
           return(NULL);
        }
        m+=0;
        m-=nrl;

        /* allocate rows and set pointers to them        */
        if((m[nrl]=(char *)calloc(nrow*ncol,sizeof(char)))==NULL)
        {
           printf("error cmatrix: calloc error No. 2 \n");
           return(NULL);
        }
        m[nrl]+=0;
        m[nrl]-=nrl;

        for(i=nrl+1;i<=nrh;i++)
           m[i]=m[i-1]+ncol;
        /* return pointer to array of pointers to rows   */
        return m;
}





