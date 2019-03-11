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

#define PADCHAR '-'
#define Max_N_NameBase 60
static char **R_Name;
static int *hit_rddex,*hit_score,*hit_rcdex,*hit_locus1,*superlength,*hit_matlocu1,*hit_matlocu2,*hit_matindex;
static int *ctg_length,*hit_index,*hit_masks,*hit_linkdex1,*hit_linkdex2;;
static int *scf_offset,*scf_length,*ctg_index,*ctg_scaffs,*ctg_offset,*ctg_mscore;
static int *scf_index,*ctg_gaplen,*ctg_scaffs,*ctg_rcdex,*ctg_offset,*ctg_mscore;

/* SSAS default parameters   */
static int IMOD=0;
static int nContig=0;
static int file_flag = 1;
static int min_len = 100000;
static int max_ctg = 1000000;
static int Gap_len = 200;
static int read_len = 150;
static int n_depth = 60;
static float m_score = 400.0;
static int i_getindex = 2;

int main(int argc, char **argv)
{
    FILE *namef,*namef2;
    int i,j,k,nSeq,args,idt,rcdex;
    int n_contig,n_reads,nseq,iscaf1,idd1,offset1;
    int stopflag;
    char *st,*ed;
    char line[2000]={0},tempc1[60],rdname[60];
    char **cmatrix(long nrl,long nrh,long ncl,long nch);

    if(argc < 2)
    {
      printf("Usage: %s <Input_10X-alignment> <input_scaff_agp> <Output_10X-alignment>\n",argv[0]);

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
       else if(!strcmp(argv[i],"-depth"))
       {
         sscanf(argv[++i],"%d",&n_depth);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-len"))
       {
         sscanf(argv[++i],"%d",&Gap_len);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-score"))
       {
         sscanf(argv[++i],"%f",&m_score);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-index"))
       {
         sscanf(argv[++i],"%d",&i_getindex);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-file"))
       {
         sscanf(argv[++i],"%d",&file_flag);
         args=args+2;
       }
    }


    fflush(stdout);
    if(system("ps aux | grep scaff10X_pmatrix; date") == -1)
    {
//        printf("System command error:\n);
    }

    nseq=0;
    if((namef = fopen(argv[args],"r")) == NULL)
    {
      printf("ERROR main:: args \n");
      exit(1);
    }
    while(!feof(namef))
    {
      if(fgets(line,2000,namef) == NULL)
      {
//        printf("fgets command error:\n);
      }
      if(feof(namef)) break;
      nseq++;
    }
    fclose(namef); 

    n_contig=0;
    if((namef = fopen(argv[args+1],"r")) == NULL)
    {
      printf("ERROR main:: args \n");
      exit(1);
    }
    while(!feof(namef))
    {
      if(fgets(line,2000,namef) == NULL)
      {
//        printf("fgets command error:\n);
      }
      if(feof(namef)) break;
      n_contig++;
    }
    fclose(namef); 



    if((hit_score = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - hit_score\n");
      exit(1);
    }
    if((hit_locus1 = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - hit_locus1\n");
      exit(1);
    }
    if((superlength = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - superlength\n");
      exit(1);
    }
    if((hit_index = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - hit_index\n");
      exit(1);
    }

    nSeq=nseq;
    R_Name=cmatrix(0,nseq+10,0,Max_N_NameBase);

    printf("reads: %d %s\n",nseq,argv[args]);
    if((namef = fopen(argv[args],"r")) == NULL)
    {
      printf("ERROR main:: reads group file \n");
      exit(1);
    }

/*  read the alignment files         */
    i=0;
    max_ctg = 0;
    while(fscanf(namef,"%s %s %d %d %d",R_Name[i],rdname,&hit_locus1[i],&hit_score[i],&superlength[i])!=EOF)
    {
        st = rdname;
        ed = strrchr(rdname,'_');
        idt = atoi(ed+1);
        hit_index[i] = idt;
        i++;
    }
    fclose(namef);
    n_reads = i;

    printf("reads: %d %d %s\n",nseq,n_reads,argv[args]);
    if((ctg_index = (int *)calloc(n_contig,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - ctg_index\n");
      exit(1);
    }
    if((scf_index = (int *)calloc(n_contig,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - scf_index\n");
      exit(1);
    }
    if((scf_length = (int *)calloc(n_contig,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - scf_index\n");
      exit(1);
    }
    if((ctg_offset = (int *)calloc(n_contig,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - ctg_offset\n");
      exit(1);
    }
    if((ctg_length = (int *)calloc(n_contig,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - ctg_offset\n");
      exit(1);
    }
    if((ctg_gaplen = (int *)calloc(n_contig,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - ctg_gaplen\n");
      exit(1);
    }
    if((ctg_rcdex = (int *)calloc(n_contig,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - ctg_rcdex\n");
      exit(1);
    }
    if((ctg_mscore = (int *)calloc(n_contig,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - ctg_mscore\n");
      exit(1);
    }
    if((ctg_scaffs = (int *)calloc(n_contig,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - ctg_scaffs\n");
      exit(1);
    }
    if((scf_offset = (int *)calloc(n_contig,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - ctg_scaffs\n");
      exit(1);
    }

    if((namef = fopen(argv[args+1],"r")) == NULL)
    {
      printf("ERROR main:: reads group file \n");
      exit(1);
    }

/*  read the alignment files         */
    i=0;
    while(fscanf(namef,"%s %d %d %d %d %d %d %s",rdname,&ctg_index[i],&ctg_length[i],&ctg_offset[i],&ctg_gaplen[i],&rcdex,&ctg_mscore[i],tempc1)!=EOF)
    {
        st = rdname;
        ed = strrchr(rdname,'_');
        idt = atoi(ed+1);
        scf_index[ctg_index[i]] = idt;
        scf_length[idt] = scf_length[idt]+ctg_length[i]+ctg_gaplen[i];
        scf_offset[ctg_index[i]] = ctg_offset[i];
        ctg_scaffs[ctg_index[i]] = idt;
        ctg_rcdex[ctg_index[i]] = rcdex;
        i++;
    }
    fclose(namef);
    
    nseq=0;
    if((namef = fopen(argv[args+2],"w")) == NULL)
    {
      printf("ERROR main:: args \n");
      exit(1);
    }


    for(i=0;i<(n_reads-1);i++)
    {
       stopflag=0;
       j=i+1;
       while((j<n_reads)&&(stopflag==0))
       {
         if(strcmp(R_Name[i],R_Name[j])==0)
         {
           j++;
         }
         else
           stopflag=1;
       }
       if((j-i)>=5) 
       {
         for(k=i;k<j;k++)
         {
            int iscaf = scf_index[hit_index[k]];
            int offset  = 0;
            int idd = hit_index[k];

            if(ctg_rcdex[idd] == 0)
              offset = hit_locus1[k]+scf_offset[idd]; 
            else
              offset = superlength[k] - hit_locus1[k]+scf_offset[idd]-read_len;
            fprintf(namef,"%s tarseq_%d %d %d %d\n",R_Name[k],iscaf,offset,hit_score[k],scf_length[iscaf]-Gap_len);
//            printf("%s tarseq_%d %d %d %d\n",R_Name[k],iscaf,offset,hit_score[k],scf_length[iscaf]-Gap_len);
         }
//      if((scf_length[iscaf1]-Gap_len) < 0)
//        printf("%s tarseq_%d %d %d %d %d %d %d %d\n",R_Name[i+1],iscaf1,offset1,ctg_rcdex[idd],ctg_rcdex[idd1],scf_length[iscaf1]-Gap_len,scf_length[iscaf1],hit_locus1[i],scf_offset[idd]); 
       }
       else
       {
//         printf("www: %s %d %d\n",R_Name[i],hit_read2[i],superlength[i]);
       }
       i=j-1;
    }
    fclose(namef); 

    printf("Job finished for %d reads!\n",n_reads);
    return EXIT_SUCCESS;

}
/* end of the main */

/*   subroutine to sort out read pairs    */
/* =============================== */
void PCRdup_Process(char **argv,int args,int nSeq,int offset1,int offset2)
/* =============================== */
{
     int i,j,k,m,n,n_scaff,n_blocks;
     int BAR,stopflag,nstep,base,num_steps;
     float rate;
     int *heap_list1,*heap_list2,*heap_steps,num_heap1,num_heap2;

     if((heap_list1 = (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_score2\n");
       exit(1);
     }
     if((heap_list2 = (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_score2\n");
       exit(1);
     }
     if((heap_steps = (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_score2\n");
       exit(1);
     }

     num_heap1 = 0;
     num_heap2 = 0;
     BAR = 5000;
     nstep = 20000;
     num_steps = 0;

//     for(k=0;k<nSeq;k++)
//        printf("kkk: %d %d %d %s\n",k,hit_matlocu1[offset1+k],hit_matlocu1[offset2+k],R_Name[hit_linkdex2[offset2+k]]);                
     for(k=0;k<nSeq;k++)
     {
/*      search reads with an index < i     */
/*      search reads with an index > i     */
        stopflag=0;
        j=k+1;
        base = hit_matlocu1[offset1+k];; 
        while((j<nSeq)&&(stopflag==0))
        {
          if((hit_matlocu1[offset1+j]<(BAR+nstep))&&(hit_matlocu1[offset1+k]>=BAR))
          {
            base = base + hit_matlocu1[offset1+j];
            j++;
          }
          else
            stopflag=1;
        }
        if((j-k)>1)
        {
          rate = (j-k)*100;
          rate = rate/nSeq;
//          printf("frequency1:%d %d %f\n",j-k,BAR,rate);
          BAR = BAR+nstep;
          heap_list1[num_heap1] = j-k;
          heap_steps[num_heap1] = k; 
          num_heap1++;
          num_steps++;
        }
        else if((j-k)==1)
        {
          rate = 100;
          rate = rate/nSeq;
//          printf("frequency11: %d %f\n",BAR,rate);
          BAR = hit_matlocu1[offset1+k];
          num_steps++;
        }
        k = j-1;
     }

     for(i=0;i<num_heap1;i++)
     {
        int rat1,rat2;

        if(i==0)
        {
          if(heap_list1[i+1] == 0)
            rat1 = 10;
          else
            rat1 = 3*heap_list1[i]/heap_list1[i+1];
          rat2 = rat1;
        }
        else if(i==(num_heap1-1))
        {
          if(heap_list1[i-1] == 0)
            rat1 = 10;
          else
            rat1 = 3*heap_list1[i]/heap_list1[i-1];
          rat2 = rat1;
        }
        else
        {
          rat1 = heap_list1[i]/heap_list1[i-1];
          rat2 = heap_list1[i]/heap_list1[i+1];
        }
        if((rat1 >= 3)&&(rat2 >=3))
        {
          int num_dups = heap_list1[i];
          int idk = heap_steps[i];
          for(j=0;j<num_dups;j++)
          {
             int idd = hit_linkdex1[offset1+idk+j];
             if((idd%2) == 0)
             {
               hit_masks[idd] = 1;
               hit_masks[idd+1] = 1;
             }
             else
             {
               hit_masks[idd] = 1;
               hit_masks[idd-1] = 1;
             }
          }
        }
     }
 //    printf("hit1: %d %d\n",num_heap1,nSeq);


     num_heap1 = 0;
     num_heap2 = 0;
     BAR = 5000;
     nstep = 20000;
     num_steps = 0;
     base = 0;
     memset(heap_list1,0,4*nSeq);
     memset(heap_steps,0,4*nSeq);

//     for(k=0;k<nSeq;k++)
//        printf("kkk: %d %d %d %s\n",k,hit_matlocu1[offset1+k],hit_matlocu2[offset2+k],R_Name[hit_linkdex2[offset2+k]]);

     for(k=0;k<nSeq;k++)
     {
/*      search reads with an index < i     */
/*      search reads with an index > i     */
        stopflag=0;
        j=k+1;
        base = hit_matlocu2[offset2+k];; 
        while((j<nSeq)&&(stopflag==0))
        {
          if((hit_matlocu2[offset2+j]<(BAR+nstep))&&(hit_matlocu2[offset2+k]>=BAR))
          {
            base = base + hit_matlocu2[offset2+j];
            j++;
          }
          else
            stopflag=1;
        }
        if((j-k)>1)
        {
          rate = (j-k)*100;
          rate = rate/nSeq;
//          printf("frequency2:%d %d %f\n",j-k,BAR,rate);
          BAR = BAR+nstep;
          heap_list1[num_heap1] = j-k;
          heap_steps[num_heap1] = k; 
          num_heap1++;
          num_steps++;
        }
        else if((j-k)==1)
        {
          rate = 100;
          rate = rate/nSeq;
//          printf("frequency22: %d %f\n",BAR,rate);
          BAR = hit_matlocu2[offset2+k];
          num_steps++;
        }
        k = j-1;
     }

     for(i=0;i<num_heap1;i++)
     {
        int rat1,rat2;

        if(i==0)
        {
          if(heap_list1[i+1] == 0)
            rat1 = 10;
          else
            rat1 = 3*heap_list1[i]/heap_list1[i+1];
          rat2 = rat1;
        }
        else if(i==(num_heap1-1))
        {
          if(heap_list1[i-1] == 0)
            rat1 = 10; 
          else
            rat1 = 3*heap_list1[i]/heap_list1[i-1];
          rat2 = rat1;
        }
        else
        {
          rat1 = heap_list1[i]/heap_list1[i-1];
          rat2 = heap_list1[i]/heap_list1[i+1];
        }
        if((rat1 >= 3)&&(rat2 >=3))
        {
          int num_dups = heap_list1[i];
          int idk = heap_steps[i];
          for(j=0;j<num_dups;j++)
          {
             int idd = hit_linkdex2[offset2+idk+j];
             if((idd%2) == 0)
             {
               hit_masks[idd] = 1;
               hit_masks[idd+1] = 1;
             }
             else
             {
               hit_masks[idd] = 1;
               hit_masks[idd-1] = 1;
             }
          }
        }
     }
//     printf("hit2: %d %d\n",num_heap1,nSeq);
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

/* =============================== */
void ArraySort_float(int n, float *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,b,jstack=0,NSTACK=50,istack[NSTACK],MIN=7;
     float a,temp;

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

/*   function to sort an array into a decreasing order:  a>b>c>....    */  
/* =============================== */
void ArraySort_float2(int n, float *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,b,jstack=0,NSTACK=50,istack[NSTACK],MIN=7;
     float a,temp;

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

/* creat char matrix with subscript ange fm[nrl...nrh][ncl...nch]  */
float   **fmatrix(long nrl,long nrh,long ncl,long nch)
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        float **fm;

        /* allocate pointers to rows        */
        if((fm=(float **)calloc(nrow,sizeof(float*)))==NULL)
        {
           printf("error fmatrix: calloc error No. 1 \n");
           return(NULL);
        }
        fm+=0;
        fm-=nrl;

        /* allocate rows and set pointers to them        */
        if((fm[nrl]=(float *)calloc(nrow*ncol,sizeof(float)))==NULL)
        {
           printf("error fmatrix: calloc error No. 2 \n");
           return(NULL);
        }
        fm[nrl]+=0;
        fm[nrl]-=nrl;

        for(i=nrl+1;i<=nrh;i++)
           fm[i]=fm[i-1]+ncol;
        /* return pointer to array of pointers to rows   */
        return fm;
}


