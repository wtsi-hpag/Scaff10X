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
static char **S_Name;
static int *ctg_index,*ctg_rcdex,*ctg_list,*ctg_head,*ctg_id2id;
static int *ctg_sfdex,*ctg_length,*ctg_mscore,*ctg_nlinks,*ctg_gaplen;
static int *ctg_index2,*ctg_rcdex2,*ctg_list2,*ctg_head2;
static int *ctg_sfdex2,*ctg_length2,*ctg_mscore2,*ctg_nlinks2,*ctg_gaplen2;
static int IMOD = 10;
static int Max_Gap = 100;
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
    int i,j,k,nSeq,args,i_contig,idt,stopflag,num_hit,n_scaff,rcdex;
    char *st,*ed,RC[2]={0};
    char line[2000]={0},ctgnameout[Max_N_NameBase],tmptext[Max_N_NameBase],temp2[Max_N_NameBase];
    char **cmatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch);
    fasta *seq,*seq2; 
    void ArraySort_Mix(int n, B64_long *arr, int *brr);
    void Read_Pairs(char **argv,int args,int nLib,int nSeq);
    void Align_Process(char **argv,int args,int nRead);
    void Cigar_Filter(char **argv,int args,int nRead);
    int n_contig,n_contig1,n_contig2,offset,start,nRead;

    if(argc < 2)
    {
      printf("Usage: %s <AGP_stage1> <AGP_stage2> <Tag_file> <Mgered_AGP>\n",argv[0]);
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

    if((namef = fopen(argv[args],"r")) == NULL)
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

    n_contig1 = i_contig;
    printf("num: %d %s\n",i_contig,argv[args]);
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
    if((ctg_gaplen= (int *)calloc(i_contig,sizeof(int))) == NULL)
    {
      printf("ERROR Memory_Allocate: Align_Process - ctg_gaplen\n");
      exit(1);
    }
    if((ctg_length= (int *)calloc(i_contig,sizeof(int))) == NULL)
    {
      printf("ERROR Memory_Allocate: Align_Process - ctg_length\n");
      exit(1);
    }
    if((ctg_id2id= (int *)calloc(i_contig,sizeof(int))) == NULL)
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

    if((namef = fopen(argv[args],"r")) == NULL)
    {
      printf("ERROR main:: reads group file \n");
      exit(1);
    }
    i = 0;

    while(fscanf(namef,"%s %d %d %s %d %d %d %d",tmptext,&ctg_index[i],&ctg_length[i],temp2,&ctg_gaplen[i],&ctg_rcdex[i],&ctg_mscore[i],&ctg_nlinks[i])!=EOF)
    {
      st = strrchr(tmptext,'_');
      ctg_sfdex[i] = atoi(st+1);
      ctg_id2id[i] = ctg_index[i];
      i++;
    }
    fclose(namef);

    n_contig1 = i;

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

    n_contig2 = i_contig;

    if((ctg_list2= (int *)calloc(i_contig,sizeof(int))) == NULL)
    {
      printf("ERROR Memory_Allocate: Align_Process - ctg_list\n");
      exit(1);
    }
    if((ctg_head2= (int *)calloc(i_contig,sizeof(int))) == NULL)
    {
      printf("ERROR Memory_Allocate: Align_Process - ctg_head\n");
      exit(1);
    }
    if((ctg_mscore2= (int *)calloc(i_contig,sizeof(int))) == NULL)
    {
      printf("ERROR Memory_Allocate: Align_Process - ctg_length\n");
      exit(1);
    }
    if((ctg_nlinks2= (int *)calloc(i_contig,sizeof(int))) == NULL)
    {
      printf("ERROR Memory_Allocate: Align_Process - ctg_nlinks\n");
      exit(1);
    }
    if((ctg_gaplen2= (int *)calloc(i_contig,sizeof(int))) == NULL)
    {
      printf("ERROR Memory_Allocate: Align_Process - ctg_gaplen2\n");
      exit(1);
    }
    if((ctg_length2= (int *)calloc(i_contig,sizeof(int))) == NULL)
    {
      printf("ERROR Memory_Allocate: Align_Process - ctg_length\n");
      exit(1);
    }
    if((ctg_rcdex2= (int *)calloc(i_contig,sizeof(int))) == NULL)
    {
      printf("ERROR Memory_Allocate: Align_Process - ctg_rcdex\n");
      exit(1);
    }
    if((ctg_index2= (int *)calloc(i_contig,sizeof(int))) == NULL)
    {
      printf("ERROR Memory_Allocate: Align_Process - ctg_index\n");
      exit(1);
    }
    if((ctg_sfdex2= (int *)calloc(i_contig,sizeof(int))) == NULL)
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

    while(fscanf(namef,"%s %d %d %s %d %d %d %d",tmptext,&ctg_index2[i],&ctg_length2[i],temp2,&ctg_gaplen2[i],&ctg_rcdex2[i],&ctg_mscore2[i],&ctg_nlinks2[i])!=EOF)
    {
      st = strrchr(tmptext,'_');
      ctg_sfdex2[i] = atoi(st+1);
      i++;
    }
    fclose(namef);

    num_hit = 0;
    for(i=0;i<n_contig1;i++)
    {
       stopflag=0;
       j=i+1;
       while((j<n_contig1)&&(stopflag==0))
       {
         if(ctg_sfdex[j]==ctg_sfdex[i])
         {
           j++;
         }
         else
           stopflag=1;
       }

       ctg_list[num_hit] = j-i;
       num_hit++;
       i=j-1;
    }

    printf("num: %d %d\n",num_hit,n_contig1);
    ctg_head[0] = 0;
    for(i=1;i<=num_hit;i++)
    {
       ctg_head[i] = ctg_head[i-1] + ctg_list[i-1];
    }

/*
    for(j=0;j<num_hit;j++)
    {
       for(i=0;i<ctg_list[j];i++)
       {
          int idd = ctg_head[j]+i;
          int idk = ctg_index[idd];
          printf("hist: %d %d %d %d %d\n",i,j,ctg_list[j],idd,idk); 
       }
    }   */

    num_hit = 0;
    for(i=0;i<n_contig2;i++)
    {
       stopflag=0;
       j=i+1;
       while((j<n_contig2)&&(stopflag==0))
       {
         if(ctg_sfdex2[j]==ctg_sfdex2[i])
         {
           j++;
         }
         else
           stopflag=1;
       }

       ctg_list2[num_hit] = j-i;
       num_hit++;
       i=j-1;
    }

    printf("num: %d %d\n",num_hit,n_contig2);
    ctg_head2[0] = 0;
    for(i=1;i<=num_hit;i++)
    {
       ctg_head2[i] = ctg_head2[i-1] + ctg_list2[i-1];
//      printf("hist: %d %d %d\n",i-1,ctg_list[i-1],ctg_sfdex[ctg_head[i-1]]);
    }

    nRead=0;
    if((namef = fopen(argv[args+2],"r")) == NULL)
    {
      printf("ERROR main:: args+1 \n");
      exit(1);
    }
    while(!feof(namef))
    {
      fgets(line,2000,namef);
      if(feof(namef)) break;
      nRead++;
    }
    fclose(namef);

    S_Name=cmatrix(0,nRead+10,0,Max_N_NameBase);

    if((namef = fopen(argv[args+2],"r")) == NULL)
    {
      printf("ERROR main:: reads group file \n");
      exit(1);
    }

/*  read the alignment files         */
    i=0;
    while(fscanf(namef,"%s %s %s %s",temp2,temp2,temp2,S_Name[i])!=EOF)
    {
        i++;
    }
    fclose(namef);


/*  input read alignment info line   */
    if((fpOutfast = fopen(argv[args+3],"w")) == NULL)
    {
      printf("ERROR main:: reads group file: %s \n",argv[args]);
      exit(1);
    }

    n_scaff = 0;
    n_contig = 0;
    for(i=0;i<num_hit;i++)
    {
       start = 0;
       offset = 0;
       i_contig = 1;
       for(j=0;j<ctg_list2[i];j++)
       {
          int idd = ctg_head2[i]+j;
          int idk = ctg_index2[idd];
//             printf("ddd_%d %d %d %d %d\n",n_scaff,idd,idk,ctg_list2[i],ctg_list[idk]);
          for(k=0;k<ctg_list[idk];k++)
          {
             int idd1,idk1,rcdex,idt;
             if(ctg_rcdex2[idd] == 0)
             {
               idd1 = ctg_head[idk]+k;
               idk1 = ctg_index[idd1];
               idt = idd1+1;
               rcdex = 0;
             }
             else
             {
               idd1 = ctg_head[idk]-k-1+ctg_list[idk];
               idk1 = ctg_index[idd1];
               idt = idd1;
               rcdex = 0;
             }
             offset = offset+ctg_length[idd1];
             rcdex = ctg_rcdex[idd1]^ctg_rcdex2[idd];
             if(rcdex == 0)
               RC[0] = '+';
             else
               RC[0] = '-'; 
//             printf("scaff10x_%d\t%d\t%d\t%d\t%c\t%s\t%d\t%d\t%s\n",n_scaff,start+1,offset,i_contig,'W',S_Name[ctg_id2id[idd1]],1,ctg_length[idd1],RC);
             fprintf(fpOutfast,"scaff10x_%d\t%d\t%d\t%d\t%c\t%s\t%d\t%d\t%s\n",n_scaff,start+1,offset,i_contig,'W',S_Name[ctg_id2id[idd1]],1,ctg_length[idd1],RC);
             i_contig++;
             if(k < (ctg_list[idk]-1))
             {
//               printf("scaff10x_%d\t\%d\t%d\t%d\t%c\t%d\t%d\t%d\t%s\n",n_scaff,offset+1,offset+100,i_contig,'N',ctg_gaplen[idd1],ctg_mscore[idd1],ctg_nlinks[idd1],"10x");
               fprintf(fpOutfast,"scaff10x_%d\t\%d\t%d\t%d\t%c\t%d\t%d\t%d\t%s\n",n_scaff,offset+1,offset+100,i_contig,'N',ctg_gaplen[idd1],ctg_mscore[idt],ctg_nlinks[idt],"10x");
               offset = offset + 100;
             }
             start = offset;
             i_contig++;
             n_contig++;
          }
          if(j < (ctg_list2[i]-1))
          {
            fprintf(fpOutfast,"scaff10x_%d\t\%d\t%d\t%d\t%c\t%d\t%d\t%d\t%s\n",n_scaff,offset+1,offset+100,i_contig,'N',ctg_gaplen2[idd],ctg_mscore2[idd+1],ctg_nlinks2[idd+1],"10x");
            offset = offset + 100;
            start = offset;
          }
       }
       n_scaff++;
    }

    fclose(fpOutfast);
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

