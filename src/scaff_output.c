/****************************************************************************
 ****************************************************************************
 *                                                                          *
 *  Copyright (C) 2018  Genome Research Ltd.                                *
 *                                                                          *
 *  Author: Zemin Ning (zn1@sanger.ac.uk)                                   *
 *                                                                          *
 *  This file is part of Scaff10X pipeline.                                 *
 *                                                                          *
 *  Scaff10X is a free software: you can redistribute it and/or modify it   *
 *  under the terms of the GNU General Public License as published by the   *
 *  Free Software Foundation, either version 3 of the License, or (at your  *
 *  option) any later version.                                              *
 *                                                                          *
 *  This program is distributed in the hope that it will be useful, but     *
 *  WITHOUT ANY WARRANTY; without even the implied warranty of              *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU        *
 *  General Public License for more details.                                *
 *                                                                          *
 *  You should have received a copy of the GNU General Public License along *
 *  with this program.  If not, see <http://www.gnu.org/licenses/>.         *
 *                                                                          *
 ****************************************************************************
 ****************************************************************************/
/****************************************************************************/


#include <math.h>
#include <values.h>
#include <stdio.h>
#include <netinet/in.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h> 
#include <sys/socket.h>
#include <sys/wait.h>
#include <sys/signal.h>
#include <errno.h>
#include "fasta.h"

#define MAXLINE 4096
#define ENDS_EXTRA 0
#define PADCHAR '-'
#define Max_N_NameBase 50
#define Max_N_Pair 100
static int *ctg_index,*ctg_mask;
static int *ctg_rcdex,*ctg_list,*ctg_head;
static int *ctg_sfdex,*ctg_length,*ctg_mscore,*ctg_nlinks;
static int IMOD = 10;
static int Max_Gap = 200;
static int longread_flag =1;
static B64_long *cigar_head,sBase;

/* SSAS default parameters   */

typedef struct
{
       int foffset;
       int fsindex;
} SIO;

fasta *sub;


int main(int argc, char **argv)
{
    FILE *fp,*namef,*fpOutfast,*fpOutfast2;
    int i,j,nSeq,args,i_contig,idt,stopflag,num_hit,n_scaff,rcdex;
    char *st,*ed;
    char line[2000]={0},ctgnameout[Max_N_NameBase],tmptext[Max_N_NameBase],tmptext2[Max_N_NameBase];
    char **cmatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch);
    fasta *seq,*seq2; 
    void ArraySort_Mix(int n, B64_long *arr, int *brr);
    void Read_Pairs(char **argv,int args,int nLib,int nSeq);
    void Align_Process(char **argv,int args,int nRead);
    void Cigar_Filter(char **argv,int args,int nRead);
    fasta *segg,*seqp,*seqp2;
    B64_long Size_q_pdata;
    int num_seqque,size_range[5],*ctg_list,*ctg_head,*ctg_loci;
    char *pdata;


    if(argc < 2)
    {
      printf("Usage: %s <input_scaffold_fastq> <10X_scaffold-structure_file> <output_scaffold_fastq> <output_scaffold_agp>\n",argv[0]);
      exit(1);
    }

    nSeq=0;
    args=1;
    for(i=1;i<argc;i++)
    {
       if(!strcmp(argv[i],"-mod"))
       {
         sscanf(argv[++i],"%d",&IMOD); 
         args=args+2;
       }
       else if(!strcmp(argv[i],"-longread"))
       {
         sscanf(argv[++i],"%d",&longread_flag);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-gap"))
       {
         sscanf(argv[++i],"%d",&Max_Gap);
         args=args+2;
       }
    }

    if((fp=fopen(argv[args],"rb"))==NULL) printf("Cannot open file\n");
      fseek(fp, 0, SEEK_END);
    Size_q_pdata = ftell(fp) + 1;
    fclose(fp);
    if((pdata=(char*)calloc(Size_q_pdata,sizeof(char)))==NULL)
      printf("calloc pdata\n");
    num_seqque = extractFastq(argv[args],pdata,Size_q_pdata);
    if((segg=(fasta*)calloc((num_seqque),sizeof(fasta)))==NULL)
      printf("calloc segg\n");
    if((seq=decodeFastq(argv[args],&num_seqque,&sBase,pdata,Size_q_pdata,segg))==NULL)
      printf("no query data found.\n");
    nSeq = num_seqque;
    fastaUC(seq,nSeq);

    if((namef = fopen(argv[args+1],"r")) == NULL)
    {
      printf("ERROR main:: reads group file \n");
      exit(1);
    }
    i_contig = 0;
    while(!feof(namef))
    {
      if(fgets(line,2000,namef) == NULL)
      {
//        printf("fgets command error:\n);
      }
      if(feof(namef)) break;
      i_contig++;
    }
    fclose(namef);

      printf("num: %d %s\n",nSeq,argv[args+1]);
    i_contig = i_contig + nSeq;
    if((ctg_list= (int *)calloc(i_contig,sizeof(int))) == NULL)
    {
      printf("ERROR Memory_Allocate: Align_Process - ctg_list\n");
      exit(1);
    }
    if((ctg_head= (int *)calloc(i_contig,sizeof(int))) == NULL)
    {
      printf("ERROR Memory_Allocate: Align_Process - ctg_head\n");
      exit(1);
    }
    if((ctg_mask= (int *)calloc(i_contig,sizeof(int))) == NULL)
    {
      printf("ERROR Memory_Allocate: Align_Process - ctg_mask\n");
      exit(1);
    }
    if((ctg_mscore= (int *)calloc(i_contig,sizeof(int))) == NULL)
    {
      printf("ERROR Memory_Allocate: Align_Process - ctg_length\n");
      exit(1);
    }
    if((ctg_nlinks= (int *)calloc(i_contig,sizeof(int))) == NULL)
    {
      printf("ERROR Memory_Allocate: Align_Process - ctg_nlinks\n");
      exit(1);
    }
    if((ctg_length= (int *)calloc(i_contig,sizeof(int))) == NULL)
    {
      printf("ERROR Memory_Allocate: Align_Process - ctg_length\n");
      exit(1);
    }
    if((ctg_rcdex= (int *)calloc(i_contig,sizeof(int))) == NULL)
    {
      printf("ERROR Memory_Allocate: Align_Process - ctg_rcdex\n");
      exit(1);
    }
    if((ctg_index= (int *)calloc(i_contig,sizeof(int))) == NULL)
    {
      printf("ERROR Memory_Allocate: Align_Process - ctg_index\n");
      exit(1);
    }
    if((ctg_sfdex= (int *)calloc(i_contig,sizeof(int))) == NULL)
    {
      printf("ERROR Memory_Allocate: Align_Process - ctg_sfdex\n");
      exit(1);
    }
    if((namef = fopen(argv[args+1],"r")) == NULL)
    {
      printf("ERROR main:: reads group file \n");
      exit(1);
    }
    i = 0;

    while(fscanf(namef,"%s %d %d %d %d",tmptext,&ctg_sfdex[i],&ctg_index[i],&ctg_length[i],&rcdex)!=EOF)
    {
      st = strrchr(tmptext,':');
      ed = strrchr(tmptext,'_');
      ctg_mscore[ctg_index[i]] = atoi(st+1);
      ctg_nlinks[ctg_index[i]] = atoi(ed+1);
      ctg_mask[ctg_index[i]] = 1; 
      ctg_rcdex[ctg_index[i]] = rcdex;
      i++;
    }
    fclose(namef);

    i_contig = i;
    num_hit = 0;
    for(i=0;i<(i_contig);i++)
    {
       stopflag=0;
       j=i+1;
       while((j<nSeq)&&(stopflag==0))
       {
         if(ctg_sfdex[j]==ctg_sfdex[i])
         {
           j++;
         }
         else
           stopflag=1;
       }

       if((longread_flag)&&((j-i)>=2)&&(ctg_mscore[ctg_index[j-1]] <= 50)&&(ctg_length[j-1]<=60000))
       {
         ctg_list[num_hit] = j-i-1;
      printf("num: %d %d %d %d\n",j-1,ctg_index[j-1],ctg_mscore[ctg_index[j-1]],ctg_length[j-1]);
         ctg_rcdex[j-1] = 0;
         num_hit++;
         ctg_list[num_hit] = 1;
       }
       else
         ctg_list[num_hit] = j-i;
       num_hit++;
       i=j-1;
    }

      printf("num: %d %d\n",num_hit,i_contig);
    ctg_head[0] = 0;
    for(i=1;i<=num_hit;i++)
    {
       ctg_head[i] = ctg_head[i-1] + ctg_list[i-1];
//      printf("hist: %d %d %d\n",i-1,ctg_list[i-1],ctg_sfdex[ctg_head[i-1]]);
    }

/*  input read alignment info line   */
    if((fpOutfast = fopen(argv[args+2],"w")) == NULL)
    {
      printf("ERROR main:: reads group file: %s \n",argv[args]);
      exit(1);
    }
    if((fpOutfast2 = fopen(argv[args+3],"w")) == NULL)
    {
      printf("ERROR main:: reads group file: %s \n",argv[args]);
      exit(1);
    }

    n_scaff = 0;
    for(i=0;i<num_hit;i++)
    {
       int seq_st,seq_ed,rc,seq_len,st2,ed2,sup_len;
       int trash_flag = 0;
       char *dpp;
       seqp = seq + ctg_index[ctg_head[i]];
       fprintf(fpOutfast,"@scaff10x_%d\n",n_scaff);
//       printf("@%s\n",seqp->name);
       sup_len = 0;
       for(j=0;j<ctg_list[i];j++)
       {
          int idd = ctg_head[i]+j;
          seqp = seq + ctg_index[idd];
          seq_st = 0;
          seq_ed = seqp->length;
          if(ctg_rcdex[ctg_index[idd]] == 0)
          {
            for(rc=seq_st;rc<seq_ed;rc++)
               fprintf(fpOutfast,"%c",seqp->data[rc]);
          }
          else
          {
            dpp = seqp->data + seq_ed-1;
            for(rc=0;rc<seq_ed;rc++)
             {
                if(*dpp == 'A') fprintf(fpOutfast,"%c",'T');
                else if(*dpp == 'C') fprintf(fpOutfast,"%c",'G');
                else if(*dpp == 'G') fprintf(fpOutfast,"%c",'C');
                else if(*dpp == 'T') fprintf(fpOutfast,"%c",'A');
                else                 fprintf(fpOutfast,"%c",*dpp);
                dpp--;
             }
          }
          fprintf(fpOutfast2,"scaff10x_%d %d %d %d %d %d %d %d\n",n_scaff,ctg_index[idd],seqp->length,sup_len,Max_Gap,ctg_rcdex[ctg_index[idd]],ctg_mscore[ctg_index[idd]],ctg_nlinks[ctg_index[idd]]);
          if(j<(ctg_list[i]-1))
          {
            for(rc=0;rc<Max_Gap;rc++)
               fprintf(fpOutfast,"%c",'N');
            sup_len = sup_len + Max_Gap;
          }
          sup_len = sup_len + seqp->length;
       }
       fprintf(fpOutfast,"\n");
       fprintf(fpOutfast,"+\n");
       for(j=0;j<ctg_list[i];j++)
       {
          int idd = ctg_head[i]+j;
          seqp = seq + ctg_index[idd];
          seq_st = 0;
          seq_ed = seqp->length;
          for(rc=0;rc<seq_ed;rc++)
             putc(40+041,fpOutfast);
          if(j<(ctg_list[i]-1))
          {
            for(rc=0;rc<Max_Gap;rc++)
               putc(0+041,fpOutfast); 
          }
       }
       fprintf(fpOutfast,"\n");
       n_scaff++;
    }

    for(i=0;i<nSeq;i++)
    {
       int seq_st,seq_ed,rc;
       if(ctg_mask[i] == 0)
       {
         seqp = seq + i;
         seq_st = 0;
         seq_ed = seqp->length;
         fprintf(fpOutfast,"@scaff10x_%d\n",n_scaff);
//         fprintf(fpOutfast,"@%s\n",seqp->name);
         for(rc=seq_st;rc<seq_ed;rc++)
            fprintf(fpOutfast,"%c",seqp->data[rc]);
         fprintf(fpOutfast,"\n");
         fprintf(fpOutfast,"+\n");
         for(rc=0;rc<seq_ed;rc++)
             putc(40+041,fpOutfast);
         fprintf(fpOutfast,"\n");
         fprintf(fpOutfast2,"scaff10x_%d %d %d %d %d %d %d %d\n",n_scaff,i,seqp->length,0,0,Max_Gap,6000,0);
         n_scaff++;
       }
    }
    fclose(fpOutfast);
    fclose(fpOutfast2);
    return EXIT_SUCCESS;

}
/* end of the main */

#define SWAP(a,b) temp=(a);(a)=b;(b)=temp;

/*   Subroutine to sort an array arr[0,...,n-1] into ascending order while
     making the corresponding reaarangement of the array brr[0,...,n-1]
     by the use of Quicksort (Sedgwick, R. 1978, Communications o fthe ACM,
     vol. 21, pp. 847-857) also see Numerical Recipes in C                  */  

/* =============================== */
void ArraySort_Long(int n, B64_long *arr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,NSTACK=50,istack[NSTACK];
     B64_long a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
             }
             arr[i+1]=a;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
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


/* =============================== */
void ArraySort_Int(int n, int *arr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,NSTACK=50,istack[NSTACK];
     int a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
             }
             arr[i+1]=a;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
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


/* =============================== */
void ArraySort_Mix(int n, B64_long *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     B64_long a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
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

/* =============================== */
void ArraySort_Int2(int n, int *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     int a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
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

/*   function to sort an array into a decreasing order:  a>b>c>....    */  
/* =============================== */
void ArraySort2_Int2(int n, int *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     int a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]>=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);

          if(arr[m]<arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
          }

          if(arr[m+1]<arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
          }

          if(arr[m]<arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          for(;;)
          {
             do i++; while (arr[i]>a);
             do j--; while (arr[j]<a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
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

/* =============================== */
void ArraySort_Mix3(int n, B64_long *arr, int *brr, int *crr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,c,NSTACK=50,istack[NSTACK];
     B64_long a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             c=crr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
                crr[i+1]=crr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
             crr[i+1]=c;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);
          SWAP(crr[k],crr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
            SWAP(crr[m],crr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
            SWAP(crr[m+1],crr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
            SWAP(crr[m],crr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          c=crr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
             SWAP(crr[i],crr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          brr[m+1]=brr[j];
          brr[j]=b;
          crr[m+1]=crr[j];
          crr[j]=c;
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


/*   to swap the string arrays           */
/* ============================================= */
void s_swap(char **Pair_Name, int i, int j)
/* ============================================= */
{
     char temp[Max_N_NameBase];

     strcpy(temp,Pair_Name[i]);
     strcpy(Pair_Name[i],Pair_Name[j]);
     strcpy(Pair_Name[j],temp);
}


/*   to sort the string array in order          */
/* ============================================= */
void ArraySort_String(int n, char **Pair_Name, int *brr)
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


/* creat an int matrix with subscript ange m[nrl...nrh][ncl...nch]  */
int     **mmatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch)
{
        B64_long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        int  **m;

        
        /* allocate pointers to rows        */
        if((m=(int **)calloc(nrow,sizeof(int*)))==NULL)
        {
           printf("error imatrix: calloc error No. 1 \n");
           return(NULL);
        }
        m+=0;
        m-=nrl;

        /* allocate rows and set pointers to them        */
        if((m[nrl]=(int *)calloc(nrow*ncol,sizeof(int)))==NULL)
        {
           printf("error imatrix: calloc error No. 2 \n");
           return(NULL);
        }
        m[nrl]+=0;
        m[nrl]-=nrl;

        for(i=nrl+1;i<=nrh;i++)
           m[i]=m[i-1]+ncol;
        /* return pointer to array of pointers to rows   */
        return m;
}

/* creat an int matrix with subscript ange m[nrl...nrh][ncl...nch]  */
int     **imatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch)
{
        B64_long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        int  **m;

        /* allocate pointers to rows        */
        if((m=(int **)calloc(nrow,sizeof(int*)))==NULL)
        {
           printf("error imatrix: calloc error No. 1 \n");
           return(NULL);
        }
        m+=0;
        m-=nrl;

        /* allocate rows and set pointers to them        */
        if((m[nrl]=(int *)calloc(nrow*ncol,sizeof(int)))==NULL)
        {
           printf("error imatrix: calloc error No. 2 \n");
           return(NULL);
        }
        m[nrl]+=0;
        m[nrl]-=nrl;

        for(i=nrl+1;i<=nrh;i++)
           m[i]=m[i-1]+ncol;
        /* return pointer to array of pointers to rows   */
        return m;
}

/* creat char matrix with subscript ange cm[nrl...nrh][ncl...nch]  */
char    **cmatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch)
{
        B64_long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        char **cm;

        /* allocate pointers to rows        */
        if((cm=(char **)calloc(nrow,sizeof(char*)))==NULL)
        {
           printf("error cmatrix: calloc error No. 1 \n");
           return(NULL);
        }
        cm+=0;
        cm-=nrl;

        /* allocate rows and set pointers to them        */
        if((cm[nrl]=(char *)calloc(nrow*ncol,sizeof(char)))==NULL)
        {
           printf("error cmatrix: calloc error No. 2 \n");
           return(NULL);
        }
        cm[nrl]+=0;
        cm[nrl]-=nrl;

        for(i=nrl+1;i<=nrh;i++)
           cm[i]=cm[i-1]+ncol;
        /* return pointer to array of pointers to rows   */
        return cm;
}

