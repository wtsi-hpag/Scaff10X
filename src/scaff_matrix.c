/****************************************************************************
 ****************************************************************************
 *                                                                          *
 *  Copyright (C) 2017  Genome Research Ltd.                                *
 *                                                                          *
 *  Author: Zemin Ning (zn1@sanger.ac.uk)                                   *
 *                                                                          *
 *  This file is part of scaff10x pipeline.                                 *
 *                                                                          *
 *  Scaff10x is a free software: you can redistribute it and/or modify it   *
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
#include "fasta.h"

#define GT '>'
#define GT4 (((((GT<<8)+GT)<<8)+GT)<<8)+GT

#define ENDS_EXTRA 0
#define PADCHAR '-'
#define MAX_N_BRG 50000 
#define MAX_N_ROW 50000 
#define Max_N_NameBase 60
#define Max_N_Pair 100
static char **S_Name,**R_Name,**R_Name2,**T_Name,**cellname;
static int *ctg_length,*ctg_mscore,*ctg_pairmp,*ctg_rcdex1,*ctg_rcdex2,*ctg_index,*ctg_dcode;
static int *ctg_used,*ctg_mask,*ctg_outp,*ctg_mapp1,*ctg_mapp2,*ctg_idex1,*ctg_idex2;

/* SSAS default parameters   */
static int IMOD=0;
static int n_length=5000;
static int i_getindex=0;
static int file_flag=0;
static int n_matrix=2000;
static int n_links=10;
static int n_PacBio=0;
static int uplinks=50;
static int i_max=0;
static int max_barcodes=1000;

typedef struct
{
       int foffset;
       int fsindex;
} SIO;

fasta *expt;

static char rc_char[500000];
static char rc_sub[5000];

int ReverseComplement(int seqdex)
{
        int i,len;
        char *tp,*dp;
        fasta *seqp;

        seqp=expt+seqdex;
        len=seqp->length;
        memset(rc_sub,'\0',5000);
        dp=rc_sub;      
        tp = seqp->data+len;
        for(i=len;--i>=0;)
        {
                int tmp = *--tp;
                if     (tmp == 't') *dp++ = 'a';
                else if(tmp == 'g') *dp++ = 'c';
                else if(tmp == 'c') *dp++ = 'g';
                else if(tmp == 'a') *dp++ = 't';
                else                *dp++ = tmp;
        }
        return(0);
}


int Reverse_Complement_Contig(char c_array[],int num_len)
{
        int i,len;
        char *tp,*dp;

        len=num_len;
        dp=rc_char;
        tp = c_array+len;
        for(i=len;--i>=0;)
        {
                int tmp = *--tp;
                if     (tmp == 't') *dp++ = 'a';
                else if(tmp == 'g') *dp++ = 'c';
                else if(tmp == 'c') *dp++ = 'g';
                else if(tmp == 'a') *dp++ = 't';
                else                *dp++ = tmp;
        }
        return(0);
}


int main(int argc, char **argv)
{
    FILE *namef;
    int i,j,nSeq,args,hitlen,nRead,n_ctgs;
    int n_contig,n_reads,n_readsMaxctg,nseq;
    fasta *seq;
    void decodeReadpair(int nSeq);
    void HashFasta_Head(int i, int nSeq);
    void HashFasta_Table(int i, int nSeq);
    void Search_SM(fasta *seq,int nSeq);
    void Assemble_SM(int arr,int brr);
    void Readname_match(fasta *seq,char **argv,int args,int nSeq,int nRead);
    void Barcoding_Process(char **argv,int args,int nSeq);
    void Memory_Allocate(int arr);
    char line[2000]={0},tempc1[60],cc[60],RC[2],readname[60],*st,*ed;
    char **cmatrix(long nrl,long nrh,long ncl,long nch);
    void Read_Pairs(char **argv,int args,fasta *seq,int nSeq);

    seq=NULL;
    if(argc < 2)
    {
      printf("Usage: %s [-file 0] [-matrix 2000] [-link 10] <barcode-sort.clean file> <tag_file> <contig-cluster_file>\n",argv[0]);

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
       else if(!strcmp(argv[i],"-len"))
       {
         sscanf(argv[++i],"%d",&n_length); 
         args=args+2;
       }
       else if(!strcmp(argv[i],"-matrix"))
       {
         sscanf(argv[++i],"%d",&n_matrix);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-link"))
       {
         sscanf(argv[++i],"%d",&n_links);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-longread"))
       {
         sscanf(argv[++i],"%d",&n_PacBio);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-uplink"))
       {
         sscanf(argv[++i],"%d",&uplinks);
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
    system("ps aux | grep scaff_matrix; date");

    nseq=0;
    if((namef = fopen(argv[args],"r")) == NULL)
    {
      printf("ERROR main:: args \n");
      exit(1);
    }
    while(!feof(namef))
    {
      fgets(line,2000,namef);
      if(feof(namef)) break;
      nseq++;
    }
    fclose(namef); 
   
    nRead=0;
    if((namef = fopen(argv[args+1],"r")) == NULL)
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

    n_ctgs = nseq+nRead;
    if((ctg_index = (int *)calloc(n_ctgs,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - ctg_index\n");
      exit(1);
    }
    if((ctg_mscore = (int *)calloc(n_ctgs,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - ctg_index\n");
      exit(1);
    }
    if((ctg_pairmp = (int *)calloc(n_ctgs,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - ctg_pairmp\n");
      exit(1);
    }
    if((ctg_dcode = (int *)calloc(n_ctgs,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - ctg_dcode\n");
      exit(1);
    }
    if((ctg_used = (int *)calloc(n_ctgs,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - ctg_used\n");
      exit(1);
    }
    if((ctg_length = (int *)calloc(n_ctgs,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - ctg_length\n");
      exit(1);
    }

    nSeq=nseq;
    S_Name=cmatrix(0,nseq+10,0,Max_N_NameBase);
    n_readsMaxctg=0;
    n_contig=0;
    n_reads=0;
    hitlen = 0;
    if((namef = fopen(argv[args],"r")) == NULL)
    {
      printf("ERROR main:: reads group file \n");
      exit(1);
    }

/*  read the alignment files         */
    i=0;
    i_max = 0;
    while(fscanf(namef,"%s %d %s %s %d %s %s %s %d %d",cc,&ctg_mscore[i],readname,S_Name[i],&ctg_index[i],cc,cc,cc,&ctg_dcode[i],&hitlen)!=EOF)
    {
        int idt;
//        if(ctg_index[i]>i_max)
//          i_max = ctg_index[i];
//        ctg_length[ctg_index[i]] = hitlen;
//    printf("%d %s %d\n",i,S_Name[i],i_max);
        i++;
    }
    fclose(namef);

    i_max = nRead;
    printf("Contigs: %d %d\n",i,i_max);
    n_reads=i;

    if((namef = fopen(argv[args+1],"r")) == NULL)
    {
      printf("ERROR main:: reads group file \n");
      exit(1);
    }

/*  read the alignment files         */
    i=0;
    while(fscanf(namef,"%s %s %d %s",cc,tempc1,&ctg_length[i],cc)!=EOF)
    {
        i++;
    }
    fclose(namef);

//    Readname_match(seq,argv,args,n_reads,nRead);
    Barcoding_Process(argv,args,n_reads);
//    Read_Pairs(argv,args,seq,n_reads);

//    printf("Job finished for %d reads!\n",nSeq);
    return EXIT_SUCCESS;

}
/* end of the main */

/*   subroutine to sort out read pairs    */
/* =============================== */
void Barcoding_Process(char **argv,int args,int nSeq)
/* =============================== */
{
     int i,j,k,m,n,n_scaff,n_uplinks;
     FILE *namef;
     int num_hits,hit_ray[5000];
     int stopflag,*n_list,*ctg_links,*ctg_hitnum;
     int offset,**r_matrix,**b_matrix,**m1_matrix,**m2_matrix,**bc_links;
     int **l0_matrix,**l1_matrix,**l2_matrix,**r0_matrix,**r1_matrix,**r2_matrix,**h12_matrix,**h21_matrix,**h22_matrix;
     void ArraySort_Mix(int n, long *arr, int *brr);
     char **DBname,RC[2];
     void ArraySort_Int2(int n, int *arr, int *brr);
     char **cmatrix(long nrl,long nrh,long ncl,long nch);
     void ArraySort_String(int n,char **Pair_Name,int *brr);
     int **imatrix(long nrl,long nrh,long ncl,long nch);
          
     if((n_list = (int *)calloc(i_max,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - n_list\n");
       exit(1);
     }
    if((ctg_idex1 = (int *)calloc(i_max,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - ctg_idex1\n");
      exit(1);
    }
    if((ctg_idex2 = (int *)calloc(i_max,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - ctg_idex2\n");
      exit(1);
    }
    if((ctg_mapp1 = (int *)calloc(i_max,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - ctg_mapp1\n");
      exit(1);
    }
    if((ctg_mapp2 = (int *)calloc(i_max,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - ctg_mapp2\n");
      exit(1);
    }
    if((ctg_rcdex1 = (int *)calloc(i_max,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - ctg_rcdex1\n");
      exit(1);
    }
    if((ctg_rcdex2 = (int *)calloc(i_max,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - ctg_rcdex2\n");
      exit(1);
    }
    if((ctg_links = (int *)calloc(i_max,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - ctg_links\n");
      exit(1);
    }
    if((ctg_hitnum = (int *)calloc(i_max,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - ctg_hitnum\n");
      exit(1);
    }
    if((ctg_mask = (int *)calloc(i_max,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - ctg_mask\n");
      exit(1);
    }
    if((ctg_outp = (int *)calloc(i_max,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - ctg_outp\n");
      exit(1);
    }
     num_hits =0;
     k = 0;
     offset = 0;
     r_matrix=imatrix(0,i_max,0,n_matrix);
     b_matrix=imatrix(0,i_max,0,n_matrix);
     m1_matrix=imatrix(0,i_max,0,n_matrix);
     m2_matrix=imatrix(0,i_max,0,n_matrix);
     l0_matrix=imatrix(0,i_max,0,n_matrix);
     l1_matrix=imatrix(0,i_max,0,n_matrix);
     l2_matrix=imatrix(0,i_max,0,n_matrix);
     r0_matrix=imatrix(0,i_max,0,n_matrix);
     r1_matrix=imatrix(0,i_max,0,n_matrix);
     r2_matrix=imatrix(0,i_max,0,n_matrix);
     bc_links=imatrix(0,i_max,0,i_max);

/*   Ourput the cigar line file   */
     for(i=0;i<(nSeq-1);i++)
     {
        stopflag=0;
        j=i+1;
        while((j<nSeq)&&(stopflag==0))
        {
          if(strcmp(S_Name[i],S_Name[j])==0)
          {
            j++;
          }
          else
            stopflag=1;
        }
        num_hits = j-i;
        if((j-i)>=2) 
        {
          int ii,id;
          for(k=0;k<num_hits;k++)
          {
             int idk = ctg_index[i+k];
             int cidk = n_list[idk];
             for(m=0;m<num_hits;m++)
             {
                int idi = ctg_index[i+m];
                if((k!=m)&&(n_list[idk]<= (n_matrix-1))&&(idk!=idi))
                {
                  int match=0;

//          if(idk == 137)
//            printf("xxx: %d %d %d %d %s %d %d\n",k,idk,idi,n_list[idk],S_Name[i],ctg_mscore[i+k],ctg_mscore[i+m]);
                  for(n=0;n<cidk;n++)
                  {
                     if(b_matrix[idk][n]==idi)
                     {
                       r_matrix[idk][n]++;
                       m1_matrix[idk][n] = m1_matrix[idk][n] + ctg_mscore[i+k];
                       m2_matrix[idk][n] = m2_matrix[idk][n] + ctg_mscore[i+m];
//             printf("yyy: %d %d %d %d %d %s\n",cidk,idk,idi,r_matrix[idk][n],b_matrix[idk][n],S_Name[i]);
//                       n_list[idk]++;
                       if(ctg_dcode[i+k] == 0)
                         l0_matrix[idk][n]++; 
                       if(ctg_dcode[i+k] == 1)
                         l1_matrix[idk][n]++; 
                       if(ctg_dcode[i+k] == 2)
                         l2_matrix[idk][n]++; 
                       if(ctg_dcode[i+m] == 0)
                         r0_matrix[idk][n]++; 
                       if(ctg_dcode[i+m] == 1)
                         r1_matrix[idk][n]++; 
                       if(ctg_dcode[i+m] == 2)
                         r2_matrix[idk][n]++; 
                       match=1;
                       break;
                     } 
                  }
                  if(match==0)
                  {
                    r_matrix[idk][cidk] = 1;
                    b_matrix[idk][cidk] = idi;
                    m1_matrix[idk][cidk] = ctg_mscore[i+k];
                    m2_matrix[idk][cidk] = ctg_mscore[i+m];
                    if(ctg_dcode[i+k] == 0)
                      l0_matrix[idk][cidk]++; 
                    if(ctg_dcode[i+k] == 1)
                      l1_matrix[idk][cidk]++; 
                    if(ctg_dcode[i+k] == 2)
                      l2_matrix[idk][cidk]++; 
                    if(ctg_dcode[i+m] == 0)
                      r0_matrix[idk][cidk]++; 
                    if(ctg_dcode[i+m] == 1)
                      r1_matrix[idk][cidk]++; 
                    if(ctg_dcode[i+m] == 2)
                      r2_matrix[idk][cidk]++; 
                    cidk++;
                    n_list[idk]++;
                  }
                }        
             } 
          }
          ii = n_list[ctg_index[i]];
          id = ctg_index[i];
//          for(k=0;k<ii;k++)
//             printf("www: %d %d %d %d %d %d %d %s\n",k,i,ctg_index[i],ii,id,r_matrix[id][k],b_matrix[id][k],S_Name[i]);
        }
	num_hits = j-i;
	offset = offset+num_hits;
        i=j-1;
     }

//        printf("www-uuu: %d\n",i_max);
     for(k=0;k<i_max;k++)
     {
        ctg_mapp1[k] = -1;
        ctg_mapp2[k] = -1;
        if(n_list[k]>=(n_matrix-1))
        {
          printf("www: %d %d %d\n",k,n_list[k],ctg_length[k]);
          ctg_mask[k] = 1;
        } 
     }
//     ctg_mask[133] = 1;
          for(k=0;k<i_max;k++)
          {
             if((n_list[k]>=1)&&(ctg_length[k]>=n_length))
             {
               int num_edge = 0;
               int hitmax1 = 0;
               int hitmax2 = 0;
               int ctgmax1 = -1;
               int ctgmax2 = -1;
               int rcdex = 0;
               int rcdex2 = 0;

              if(file_flag)
               printf("scaffold: %d %d %d\n",k,n_list[k],ctg_length[k]);
               if((k==i_getindex))
               {
//                 for(i=0;i<n_list[k];i++)
//                    printf("yyy: %s %d %d %d %d\n",S_Name[i],i,k,r_matrix[k][i],b_matrix[k][i]);
               }
               n_uplinks = 0;
               for(i=0;i<n_list[k];i++)
               {
                  int r_gap = 0;
                  int l_gap = 0;

                  if(l0_matrix[k][i] > l1_matrix[k][i])
                    l_gap = l0_matrix[k][i] - l1_matrix[k][i];
                  else
                    l_gap = l1_matrix[k][i] - l0_matrix[k][i];         
                  if(r0_matrix[k][i] > r1_matrix[k][i])
                    r_gap = r0_matrix[k][i] - r1_matrix[k][i];
                  else
                    r_gap = r1_matrix[k][i] - r0_matrix[k][i];         
//                  if((r_matrix[k][i]>=n_links)&&(r_matrix[k][i]<=max_barcodes)&&(l0_matrix[k][i]!=l1_matrix[k][i])&&(r0_matrix[k][i]!=r1_matrix[k][i]))
                  if((r_matrix[k][i]>=n_links)&&(r_matrix[k][i]<=max_barcodes)&&(l_gap >= 2)&&(r_gap >= 2))
                  {
                    n_uplinks++;
                    if(r_matrix[k][i] > hitmax1)
                    {
                      if(l0_matrix[k][i] > l1_matrix[k][i])
                        rcdex = 1;
                      else
                        rcdex = 2;;
                      hitmax1 = r_matrix[k][i];
                      ctgmax1 = i; 
                    }
                  }
               }
               printf("uplinks: %d %d %d\n",k,n_list[k],n_uplinks);
               for(i=0;i<n_list[k];i++)
               {
                  int r_gap = 0;
                  int l_gap = 0;

                  if(l0_matrix[k][i] > l1_matrix[k][i])
                    l_gap = l0_matrix[k][i] - l1_matrix[k][i];
                  else
                    l_gap = l1_matrix[k][i] - l0_matrix[k][i];         
                  if(r0_matrix[k][i] > r1_matrix[k][i])
                    r_gap = r0_matrix[k][i] - r1_matrix[k][i];
                  else
                    r_gap = r1_matrix[k][i] - r0_matrix[k][i];         
                  if((r_matrix[k][i]>=n_links)&&(r_matrix[k][i]<=max_barcodes)&&(ctg_mask[b_matrix[k][i]]==0))
                  {
                    if((r_matrix[k][i] >= hitmax2)&&(r_matrix[k][i]!=hitmax1)&&(l_gap >= 2)&&(r_gap >= 2))
                    {
                      if(l0_matrix[k][i] > l1_matrix[k][i])
                        rcdex2 = 1;
                      else
                        rcdex2 = 2;
                      if(rcdex != rcdex2)
                      {
                        hitmax2 = r_matrix[k][i];
                        ctgmax2 = i;
                      }
                    }
                    num_edge++;
                    if((k==i_getindex))
                      printf("xxx: %s %d %d %d %d %d\n",S_Name[i],hitmax1,hitmax2,b_matrix[k][i],ctgmax1,ctgmax2);
                   if(file_flag)
                    printf(" %d %d %d | %d %d %d %d %d %d || ",r_matrix[k][i],b_matrix[k][i],ctg_length[b_matrix[k][i]],l0_matrix[k][i],l1_matrix[k][i],l2_matrix[k][i],r0_matrix[k][i],r1_matrix[k][i],r2_matrix[k][i]);
                  }
               }
                   if(file_flag)
               printf("\n");
//               if((num_edge>0)&&(num_edge<=2))
               if((n_list[k]< (n_matrix-10))&&(n_uplinks < uplinks)) 
               {             
                 for(i=0;i<n_list[k];i++)
                 {
                    int idd = b_matrix[k][i];
                    int nn;
                    int len1 = ctg_length[k];
                    int len2 = ctg_length[idd];
                    int r_gap = 0;
                    int l_gap = 0;
                    int mps1 = idd;
                    int mps2 = idd;
                    int PacBio_filter = 1;
                    if(idd > 0)
                    {
                      if(len1 > 60000)
                        mps1 = 60;
                      else
                        mps1 = m1_matrix[k][i]/idd;
                      if(len2 > 60000)
                         mps2 = 60;
                      else
                        mps2 = m2_matrix[k][i]/idd;
                    }
                    else
                    {
                      mps1 = 0;
                      mps2 = 0;
                      PacBio_filter = 1;
                    }
                    if(n_PacBio)
                    {
                      if((mps1 >= 50)&&(mps2 >= 50))
                        PacBio_filter = 1;
                      else
                        PacBio_filter = 0;
                    }
                    else
                      PacBio_filter = 1; 
                    if(l0_matrix[k][i] > l1_matrix[k][i])
                      l_gap = l0_matrix[k][i] - l1_matrix[k][i];
                    else
                      l_gap = l1_matrix[k][i] - l0_matrix[k][i];         
                    if(r0_matrix[k][i] > r1_matrix[k][i])
                      r_gap = r0_matrix[k][i] - r1_matrix[k][i];
                    else
                      r_gap = r1_matrix[k][i] - r0_matrix[k][i];         
//                    if((r_matrix[k][i]>=n_links)&&(ctg_mask[b_matrix[k][i]]==0)&&(r_matrix[k][i]<=max_barcodes)&&(l_gap >= 2)&&(r_gap >= 2))
                    if((PacBio_filter)&&(r_matrix[k][i]>=n_links)&&(ctg_mask[b_matrix[k][i]]==0)&&(r_matrix[k][i]<=max_barcodes)&&(l_gap >= 2)&&(r_gap >= 2))
                    {

                      if(i==ctgmax1)
                      {
                        ctg_mapp1[k] = b_matrix[k][i];
                        ctg_idex1[k] = i;
                      }
                      if(i==ctgmax2)
                      {
                        ctg_mapp2[k] = b_matrix[k][i]; 
                        ctg_idex2[k] = i;
                      }
//                      if(((i==ctgmax1)||(i==ctgmax2))&&(r_matrix[k][i]>=n_links))
                      if(((i==ctgmax1)||(i==ctgmax2))&&(r_matrix[k][i]>=n_links)&&(ctg_mask[b_matrix[k][i]]==0))
                      {
                        if(i==ctgmax1)
                        {
                          if(l0_matrix[k][i] > l1_matrix[k][i])
                          {
                            if(r0_matrix[k][i] > r1_matrix[k][i])
                              ctg_rcdex1[k] = 0;
                            else if(r0_matrix[k][i] != r1_matrix[k][i])
                              ctg_rcdex1[k] = 1;
                            else
                              ctg_rcdex1[k] = 4;
                          }
                          else if(l0_matrix[k][i] != l1_matrix[k][i])
                          {
                            if(r0_matrix[k][i] > r1_matrix[k][i])
                              ctg_rcdex1[k] = 2;
                            else if(r0_matrix[k][i] != r1_matrix[k][i])
                              ctg_rcdex1[k] = 3;
                            else 
                              ctg_rcdex1[k] = 4;
                          }
                          else
                            ctg_rcdex1[k] = 4;
                          if(ctg_rcdex1[k] != 4)
                          {
                            int m1_score = r_matrix[k][i];
                            int m2_score = r_matrix[k][i];
                            if(m1_score > 0)
                              m1_score = m1_matrix[k][i]/m1_score;
                            else
                              m1_score = 0;
                            if(m2_score > 0)
                              m2_score = m2_matrix[k][i]/m2_score;
                            else
                              m2_score = 0;
                            ctg_hitnum[k]++;
                            if(ctg_pairmp[b_matrix[k][i]] == 0)
                              ctg_pairmp[b_matrix[k][i]] = m2_score;
                              bc_links[k][b_matrix[k][i]] = r_matrix[k][i];;
                    if((b_matrix[k][i]==i_getindex))
                      printf("xxx: %s %d %d %d %d %d %d %d\n",S_Name[i],hitmax1,hitmax2,b_matrix[k][i],ctgmax1,ctgmax2,m1_score,m2_score);
                            if(file_flag)
                            printf("mapped: 1 %d %d %d %d %d || %d %d || %d %d %d %d %d %d || %d\n",i,r_matrix[k][i],b_matrix[k][i],hitmax1,hitmax2,m1_score,m2_score,l0_matrix[k][i],l1_matrix[k][i],l2_matrix[k][i],r0_matrix[k][i],r1_matrix[k][i],r2_matrix[k][i],ctg_rcdex1[k]);

                          }
                        }
                        if((i==ctgmax2)&&(i!=ctgmax1))
                        {
                          if(l0_matrix[k][i] > l1_matrix[k][i])
                          {
                            if(r0_matrix[k][i] > r1_matrix[k][i])
                              ctg_rcdex2[k] = 0;
                            else if(r0_matrix[k][i] != r1_matrix[k][i])
                              ctg_rcdex2[k] = 1;
                            else
                              ctg_rcdex2[k] = 4;
                          }
                          else if(l0_matrix[k][i] != l1_matrix[k][i])
                          {
                            if(r0_matrix[k][i] > r1_matrix[k][i])
                              ctg_rcdex2[k] = 2;
                            else if(r0_matrix[k][i] != r1_matrix[k][i])
                              ctg_rcdex2[k] = 3;
                            else
                              ctg_rcdex2[k] = 4;
                          }
                          else
                            ctg_rcdex2[k] = 4;
                          if(ctg_rcdex2[k] != 4)
                          {
                            int m1_score = r_matrix[k][i];
                            int m2_score = r_matrix[k][i];
                            if(m1_score > 0)
                              m1_score = m1_matrix[k][i]/m1_score;
                            else
                              m1_score = 0;
                            if(m2_score > 0)
                              m2_score = m2_matrix[k][i]/m2_score;
                            else
                              m2_score = 0;
                            ctg_hitnum[k]++;
                            if(ctg_pairmp[b_matrix[k][i]] == 0)
                              ctg_pairmp[b_matrix[k][i]] = m2_score;
                              bc_links[k][b_matrix[k][i]] = r_matrix[k][i];;
                    if((b_matrix[k][i]==i_getindex))
                      printf("xxx: %s %d %d %d %d %d %d %d\n",S_Name[i],hitmax1,hitmax2,b_matrix[k][i],ctgmax1,ctgmax2,m1_score,m2_score);
                            if(file_flag)
                            printf("mapped: 2 %d %d %d %d %d || %d %d || %d %d %d %d %d %d || %d\n",i,r_matrix[k][i],b_matrix[k][i],hitmax1,hitmax2,m1_score,m2_score,l0_matrix[k][i],l1_matrix[k][i],l2_matrix[k][i],r0_matrix[k][i],r1_matrix[k][i],r2_matrix[k][i],ctg_rcdex2[k]);
                          }
                        }
                        ctg_links[k]++;
                      }
                    }
                 }
               }
             }
          }

          
     if((namef = fopen(argv[args+2],"w")) == NULL)
     {
       printf("ERROR main:: alignment file 2 \n");
       exit(1);
     }
          n_scaff = 0;
          for(k=0;k<i_max;k++)
          {
            int tbase = 0;
//             if(k==i_getindex)
//               printf("hhh: %d %d \n",k,ctg_hitnum[k]);
            if((ctg_hitnum[k]==0)&&(ctg_outp[k]==0))
            {
              printf("supercontig: tarseq_%d %d\n",k,k);
              fprintf(namef,"contig-1:%d_%d %d %d %d 0\n",ctg_pairmp[k],0,n_scaff,k,ctg_length[k]);
              printf("contig-1: %d %d %d 0\n",n_scaff,k,ctg_length[k]);
              tbase = ctg_length[k];
              printf("bases: %d %d %d\n",k,n_scaff,tbase);
              ctg_outp[k] = 1;
              n_scaff++;
            }
            else if(ctg_hitnum[k]==1)
            {
             int idd = k;
             int idk = k;
             int stopflag=0;
             int num_loops=0;

//             if(k==i_getindex)
//               printf("hhh: %d %d %d %d %d %d %d %d\n",k,idd,idk,ctg_hitnum[k],ctg_used[k],ctg_used[idd],ctg_mapp1[idd],ctg_mapp2[idd]);
             while((idd >= 0)&&(ctg_mapp1[idd] >= 0)||(ctg_mapp2[idd] >= 0)&&(stopflag == 0)&&(ctg_links[idk]>0)&&(ctg_mask[idk]==0))
             {
               int rc_idk = -1;
               int rc_idi = 0;
               if(ctg_used[idk] == 0)
               {
//             printf("hits: %d %d %d %d %d %d %d %d\n",k,idd,idk,ctg_hitnum[k],ctg_used[k],ctg_used[idd],ctg_mapp1[idd],ctg_mapp2[idd]);
                 if(ctg_used[ctg_mapp1[idk]] == 0)
                 {
                   idk = ctg_mapp1[idd];
                   rc_idk = ctg_rcdex1[idd];
                 }
                 else if(ctg_used[ctg_mapp2[idk]] == 0)
                 {
                   idk = ctg_mapp2[idd];
                   rc_idk = ctg_rcdex2[idd];      
                 }
                 if((ctg_hitnum[k]==1)&&(ctg_used[k]==0)&&(idk>=0))
                 {
                   printf("supercontig: tarseq_%d %d %d %d %d %d\n",k,idd,idk,rc_idk,ctg_rcdex1[k],ctg_hitnum[k]);
                   if(rc_idk==0)
                   {
                     fprintf(namef,"contigg1:%d_%d %d %d %d 1\n",ctg_pairmp[k],bc_links[idd][idk],n_scaff,k,ctg_length[k]);
                     printf("contigg1: %d %d %d 1\n",n_scaff,k,ctg_length[k]);
                     ctg_outp[k] = 1;
                     if(idk!=k)
                     {
                       if(ctg_outp[idk] == 0)
                       {
                         fprintf(namef,"contigg2:%d_%d %d %d %d 0\n",ctg_pairmp[idk],bc_links[idd][idk],n_scaff,idk,ctg_length[idk]);
                         printf("contigg2: %d %d %d 0\n",n_scaff,idk,ctg_length[idk]);
                       }
                       ctg_outp[idk] = 1;
                     }
                   }
                   else if(rc_idk==1)
                   {
                     fprintf(namef,"contigg1:%d_%d %d %d %d 1\n",ctg_pairmp[k],bc_links[idd][idk],n_scaff,k,ctg_length[k]);
                     printf("contigg1: %d %d %d 1\n",n_scaff,k,ctg_length[k]);
                     ctg_outp[k] = 1;
                     if(idk!=k)
                     {
                       if(ctg_outp[idk] == 0)
                       {
                         fprintf(namef,"contigg2:%d_%d %d %d %d 1\n",ctg_pairmp[idk],bc_links[idd][idk],n_scaff,idk,ctg_length[idk]);
                         printf("contigg2: %d %d %d 1\n",n_scaff,idk,ctg_length[idk]);
                       }
                       ctg_outp[idk] = 1;
                     }
                   }
                   else if(rc_idk==2)
                   {
                     fprintf(namef,"contigg1:%d_%d %d %d %d 0\n",ctg_pairmp[k],bc_links[idd][idk],n_scaff,k,ctg_length[k]);
                     printf("contigg1: %d %d %d 0\n",n_scaff,k,ctg_length[k]);
                     ctg_outp[k] = 1;
                     if(idk!=k)
                     {
                       if(ctg_outp[idk] == 0)
                       {
                         fprintf(namef,"contigg2:%d_%d %d %d %d 0\n",ctg_pairmp[idk],bc_links[idd][idk],n_scaff,idk,ctg_length[idk]);
                         printf("contigg2: %d %d %d 0\n",n_scaff,idk,ctg_length[idk]);
                       }
                       ctg_outp[idk] = 1;
                     }
                   }
                   else if(rc_idk==3)
                   {
                     fprintf(namef,"contigg1:%d_%d %d %d %d 0\n",ctg_pairmp[k],bc_links[idd][idk],n_scaff,k,ctg_length[k]);
                     printf("contigg1: %d %d %d 0\n",n_scaff,k,ctg_length[k]);
                     ctg_outp[k] = 1;
                     if(idk!=k)
                     {
                       if(ctg_outp[idk] == 0)
                       {
                         fprintf(namef,"contigg2:%d_%d %d %d %d 1\n",ctg_pairmp[idk],bc_links[idd][idk],n_scaff,idk,ctg_length[idk]);
                         printf("contigg2: %d %d %d 1\n",n_scaff,idk,ctg_length[idk]);
                       }
                       ctg_outp[idk] = 1;
                     }
                   }
                   else
                   {
                     fprintf(namef,"contigg1:%d_%d %d %d %d 0\n",ctg_pairmp[k],bc_links[idd][idk],n_scaff,k,ctg_length[k]);
                     printf("contigg1: %d %d %d 0\n",n_scaff,k,ctg_length[k]);
                     ctg_outp[k] = 1;
                   }
                   tbase = tbase + ctg_length[k];
                   if(idk!=k)
                     tbase = tbase + ctg_length[idk];
                   ctg_used[k] = 0;
                   num_loops++;
                 }
                 else if((idd>=0)&&(ctg_mapp1[idd]>=0)&&(ctg_mapp2[idd]>=0)&&(idd!=idk))
                 {
                   int rc_idd = 0;
                   int rc_ide = 0;

//                   if(idk>=0)
                   tbase = tbase + ctg_length[idk];
                   if(rc_idk==0)
                     rc_idd = 0;
                   else if(rc_idk==1)
                     rc_idd = 1;
                   else if(rc_idk==2)
                     rc_idd = 0;
                   else if(rc_idk==3)
                     rc_idd = 1;
                   if(rc_idd!=rc_idi)
                     rc_ide = 1;
                   else
                     rc_ide = 0;
//             printf("hhh: %d %d %d %d %d %d\n",k,idd,idk,ctg_hitnum[k],ctg_mapp1[idd],ctg_mapp2[idd]);
                   if(ctg_outp[idk] == 0)
                   {
                     fprintf(namef,"contig-0:%d_%d %d %d %d %d\n",ctg_pairmp[idk],bc_links[idd][idk],n_scaff,idk,ctg_length[idk],rc_ide);
//                     fprintf(namef,"contig-0:%d_%d %d %d %d %d\n",ctg_pairmp[idk],bc_links[idd][idk],n_scaff,idk,ctg_length[idk],rc_ide);
                     printf("contig-0: %d %d %d %d %d %d\n",n_scaff,idk,ctg_length[idk],rc_ide,idd,idk);
                   }
                   ctg_outp[idk] = 1;
                   rc_idi = rc_ide;       
                   num_loops++;
                 }
                 else if(ctg_mapp1[idd]>0)
                 {
//             printf("hhh: %d %d %d %d %d %d\n",k,idd,idk,ctg_hitnum[k],ctg_mapp1[idd],ctg_mapp2[idd]);
                 }
//             if(idd==i_getindex)
//               printf("hhh3: %d %d %d %d\n",k,idd,ctg_hitnum[idd],ctg_used[idd]);
                 ctg_used[idd] = 1;
                 idd = idk;
               }
               else
                 stopflag=1;
               if(stopflag == 1)
                 break;
             }
             if(tbase == 0)
               tbase = ctg_length[idk];
             if(num_loops != 0)
             {
               printf("bases: %d %d %d\n",idk,n_scaff,tbase);
               n_scaff++;
             }
            }
          }
          for(k=0;k<i_max;k++)
          {
             int tbase = 0;
             if(ctg_outp[k] == 0)
             {
               printf("supercontig: tarseq_%d %d\n",k,k);
               fprintf(namef,"contig-n:%d_%d %d %d %d 0\n",ctg_pairmp[k],bc_links[k][k],n_scaff,k,ctg_length[k]);
               printf("contig-n: %d %d %d 0\n",n_scaff,k,ctg_length[k]);
               tbase = ctg_length[k];
               printf("bases: %d %d %d\n",k,n_scaff,tbase);
               n_scaff++;
             }
          }
//     fclose(namef); 
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

