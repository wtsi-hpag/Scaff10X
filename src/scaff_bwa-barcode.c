/***********************************************************************\
 *                                                                     * 
 *                     PROJECT   ssaha_Assembly                        *
 *                                                                     * 
 *---------------------------------------------------------------------*
 *                  Assembly Program using SSAHA                       *
 *                                                                     * 
 *      Sequence Search and Alignment by the Hashing Algorithm         *
 *                                                                     *
 *                                By                                   *
 *                                                                     *
 *                    Zemin Ning & James C. Mullikin                   *
 *                                                                     *
 *          Copyright (C) 1999-2000 by Genome Research Limited         *
 *                       All rights reserved                           *
 *                                                                     *
 *---------------------------------------------------------------------*
 #######################################################################
 #                                                                     #
 #             <------   LICENCE NOTICE   ------>                      #
 #                                                                     #
 # This is a licensed software developed by Genome Research Limited    #
 # (GRL) for genomic sequence assembling. For both commercial and non- # 
 # commercial purposes, a licence must be obtained from GRL before     #
 # use. Under no circumstances, the users should copy, modify and      #
 # redistribut the software as a whole or even part of it without      #
 # permission from GRL. For more information about the software and    #
 # its ducumentation particularly, please contact either one of the    # 
 # authors or GRL.                                                     #
 #######################################################################
 *---------------------------------------------------------------------*/

#include <math.h>
#include <values.h>
#include <stdio.h>
#include <netinet/in.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <ctype.h>
#include "fasta.h"

#define GT '>'
#define GT4 (((((GT<<8)+GT)<<8)+GT)<<8)+GT

#define ENDS_EXTRA 0
#define PADCHAR '-'
#define MAX_N_BRG 50000 
#define MAX_N_ROW 50000 
#define Max_N_NameBase 20
#define Max_N_Pair 100
static char **S_Name,**R_Name1,**R_Name2,**T_Name,**cellname;
static int *hit_index,*ctg_list,*ctg_head,*hit_length,*ctg_index,*hit_mask;

/* SSAS default parameters   */
static int IMOD=0;
static int n_type=0;
static int barreads=10;
static int file_flag=2;
static int tiles_flag=0;
static int edge_set=2;
static int edge_flag=0;
static int nContig=0;
static int max_len = 100000;
typedef struct
{
       int foffset;
       int fsindex;
} SIO;

fasta *expt;

int main(int argc, char **argv)
{
    FILE *namef,*namef2;
    int i,j,nSeq,args;
    int n_contig,nseq,n_reads;
    fasta *seq;
    void decodeReadpair(int nSeq);
    void Mapping_Process(char **argv,int args,int nSeq,int nRead);
    char temp1[100],temp2[100],temp3[100],temp4[100],temp5[100],temp6[100],rdname[100];
    char line[2000]={0},cc[60],RC[2],*st,*ed;
    char **cmatrix(long nrl,long nrh,long ncl,long nch);

    seq=NULL;
    if(argc < 2)
    {
      printf("Usage: %s <Assembly_tag> <align.dat> <Processed file>\n",argv[0]);

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
       else if(!strcmp(argv[i],"-type"))
       {
         sscanf(argv[++i],"%d",&n_type); 
         args=args+2;
       }
       else if(!strcmp(argv[i],"-edge"))
       {
         sscanf(argv[++i],"%d",&edge_set);
         edge_flag=1;
         args=args+2;
       }
       else if(!strcmp(argv[i],"-tile"))
       {
         sscanf(argv[++i],"%d",&tiles_flag);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-reads"))
       {
         sscanf(argv[++i],"%d",&barreads);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-max"))
       {
         sscanf(argv[++i],"%d",&max_len);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-file"))
       {
         sscanf(argv[++i],"%d",&file_flag);
         args=args+2;
       }
    }

    nseq=0;
    n_reads=0;

    if((namef = fopen(argv[args],"r")) == NULL)
    {
      printf("ERROR main:: reads group file \n");
      exit(1);
    }

    while(!feof(namef))
    {
      fgets(line,2000,namef);
      if(feof(namef)) break;
      nseq++;
    }
    fclose(namef);

    R_Name1=cmatrix(0,nseq,0,Max_N_NameBase);
    R_Name2=cmatrix(0,nseq,0,Max_N_NameBase);
    T_Name=cmatrix(0,nseq+nseq,0,Max_N_NameBase);

    if((namef = fopen(argv[args],"r")) == NULL)
    {
      printf("ERROR main:: reads group file \n");
      exit(1);
    }

/*  read the alignment files         */
    i=0;
    while(fscanf(namef,"%s %s %s %s",R_Name1[i],temp2,temp3,R_Name2[i])!=EOF)
    {
        i++;
    }
    fclose(namef);
    
    if((namef = fopen(argv[args+1],"r")) == NULL)
    {
      printf("ERROR main:: reads group file \n");
      exit(1);
    }

    while(!feof(namef))
    {
      fgets(line,2000,namef);
      if(feof(namef)) break;
      n_reads++;
    }
    fclose(namef);

    if((hit_index = (int *)calloc(n_reads,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - ctg_index\n");
      exit(1);
    }
    if((ctg_index = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - ctg_index\n");
      exit(1);
    }
    if((ctg_list = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - ctg_list\n");
      exit(1);
    }
    if((ctg_head = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - ctg_head\n");
      exit(1);
    }

    S_Name=cmatrix(0,n_reads,0,Max_N_NameBase);

    if((namef = fopen(argv[args+1],"r")) == NULL)
    {
      printf("ERROR main:: reads group file \n");
      exit(1);
    }

/*  read the alignment files         */
    i=0;
    while(fscanf(namef,"%s %s %s %s %s %s",temp1,temp2,temp3,S_Name[i],temp5,temp6)!=EOF)
    {
        i++;
    }
    fclose(namef);

    n_reads = i;
    printf("%d %d\n",nseq,n_reads);
    Mapping_Process(argv,args,nseq,n_reads);

    if((namef = fopen(argv[args+1],"r")) == NULL)
    {
      printf("ERROR main:: reads group file \n");
      exit(1);
    }
    i=0;
    if((namef2 = fopen(argv[args+2],"w")) == NULL)
    {
      printf("ERROR main:: args+1 \n");
      exit(1);
    }

    while(fscanf(namef,"%s %s %s %s %s %s",temp1,temp2,temp3,temp4,temp5,temp6)!=EOF)
    {
        int idt,idd;
        memset(rdname,'\0',100);
        st = temp2;
        ed= strrchr(temp2,'Z');
        strcpy(rdname,ed+2);
        idt = hit_index[i];
        idd = ctg_index[idt]; 
        fprintf(namef2,"%s_%s %s %s %s %s\n",temp1,rdname,temp3,R_Name1[idd],temp5,temp6);
        i++;
    }
    fclose(namef);
    fclose(namef2);

//    printf("Job finished for %d reads!\n",nSeq);
    return EXIT_SUCCESS;

}
/* end of the main */

/*   subroutine to sort out read pairs    */
/* =============================== */
void Mapping_Process(char **argv,int args,int nSeq,int nRead)
/* =============================== */
{
     FILE *namef,*namef2;
     int i,j,k,m,n,idd,num_rd_find;
     int num_hits,n_contigs,n_scaffs;
     int stopflag,*readIndex;
     int offset;
     void ArraySort_Mix(int n, long *arr, int *brr);
     char **DBname,RC[2];
     void ArraySort_Int2(int n, int *arr, int *brr);
     char **cmatrix(long nrl,long nrh,long ncl,long nch);
     void ArraySort_String(int n,char **Pair_Name,int *brr);

     n_contigs = 3*nSeq;          
     if((readIndex= (int *)calloc(n_contigs,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - hit_maps\n");
       exit(1);
     }
     DBname=cmatrix(0,n_contigs,0,Max_N_NameBase);

     num_hits =0;
     k = 0;
     offset = 0;
     for(i=0;i<nRead;i++)
     {
        stopflag=0;
        j=i+1;
        while((j<nRead)&&(stopflag==0))
        {
          if(strcmp(S_Name[i],S_Name[j])==0)
          {
            j++;
          }
          else
            stopflag=1;
        }
        if((j-i)>=1) 
        {
	  num_hits = j-i;
          for(n=i;n<j;n++)
             hit_index[n] = k;
          strcpy(T_Name[k],S_Name[i]);
          ctg_list[k] = num_hits;
//            printf("www: %s %s %d\n",S_Name[i],R_Name1[k],j-i);
          k++;
        }
        else
        {
//            printf("yyy: %s %s %d %d %d\n",S_Name[i],R_Name[i],hit_locus[i],ctg_head[i],hit_length[i]);
        }
        i=j-1;
     }

     n_contigs = k;
     printf("Numbers of contigs: %d %d\n",n_contigs,nSeq);

     if(n_contigs != nSeq)
     {
        printf("Numbers of contigs are difference, please check reference assembly! %d %d\n",n_contigs,nSeq);
        exit(1);
     }
     ctg_head[0] = 0;
     for(i=1;i<n_contigs;i++)
     {
        ctg_head[i] = ctg_list[i-1]+ctg_head[i-1];
     }

/*   find out the read name match   */
     for(j=0;j<nSeq;j++)
     {
        strcpy(DBname[j],R_Name2[j]);
        readIndex[j]=j;
     }
     for(j=0;j<n_contigs;j++)
     {
        strcpy(DBname[j+nSeq],T_Name[j]);
        readIndex[j+nSeq]=j+nSeq;
     }

     n_scaffs=nSeq+n_contigs;
     ArraySort_String(n_scaffs,DBname,readIndex);

     num_rd_find=0;
     for(i=0;i<n_scaffs;i++)
     {
/*      search reads with an index < i     */
/*      search reads with an index > i     */
        stopflag=0;
        j=i+1;
        if(readIndex[i]>=nSeq)
          idd = readIndex[i];
        while((j<n_scaffs)&&(stopflag==0))
        {
          if(strcmp(DBname[j],DBname[i])==0)
          {
            if(readIndex[j]>=nSeq)
              idd = readIndex[j];
//            num_rd_find++;
            j++;
          }
          else
            stopflag=1;
        }
        if((j-i)>=2)
        {
          for(k=i;k<j;k++)
          {
             if(readIndex[k]<nSeq)
             {
               ctg_index[readIndex[k]] = idd-nSeq;
//         printf("name: %d %s\n",k,DBname[k]);
               num_rd_find++;
             }
          }
        }
        i=j-1;
     }

//     for(i=0;i<nSeq;i++)
//        printf("Name: %s %d %d %d %d\n",R_Name1[i],i,nSeq,num_rd_find,ctg_index[i]);

/*
     if((namef2 = fopen(argv[args+2],"w")) == NULL)
     {
       printf("ERROR main:: args+1 \n");
       exit(1);
     }

     fclose(namef2);   */
     
}


#define SWAP(a,b) temp=(a);(a)=b;(b)=temp;

/*   Subroutine to sort an array arr[0,...,n-1] into ascending order while
     making the corresponding reaarangement of the array brr[0,...,n-1]
     by the use of Quicksort (Sedgwick, R. 1978, Communications o fthe ACM,
     vol. 21, pp. 847-857) also see Numerical Recipes in C                  */  

/* =============================== */
void ArraySort_Long(int n, long *arr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,NSTACK=50,istack[NSTACK];
     long a,temp,MIN=7;

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
void ArraySort_Mix(int n, long *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     long a,temp,MIN=7;

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
void ArraySort_Mix3(int n, long *arr, int *brr, int *crr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,c,NSTACK=50,istack[NSTACK];
     long a,temp,MIN=7;

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
int     **imatrix(long nrl,long nrh,long ncl,long nch)
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
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
char    **cmatrix(long nrl,long nrh,long ncl,long nch)
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
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

