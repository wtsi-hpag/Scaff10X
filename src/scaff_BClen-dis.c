#include <values.h>
#include <stdio.h>
#include <netinet/in.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <ctype.h>

static int *nn;

int main(int argc, char **argv)
{
  int i=0,j=0,k,len=0,num_steps,nSeq,num_base;
  int *s_len,BAR = 0,nstep = 0,stopflag,base;
  char line[100],tempc1[100];
  FILE *namef;
  float rate;

    if((s_len= (int *)calloc(100000000,sizeof(int))) == NULL)
    {
      printf("ERROR Memory_Allocate: Align_Process - hit_rcdex\n");
      exit(1);
    }

    if((namef = fopen(argv[1],"r")) == NULL)
    {
      printf("ERROR main:: reads group file \n");
      exit(1);
    }

/*  read the alignment files         */
    i=0;
    num_base = 0;
    while(fscanf(namef,"%s %d",tempc1,&s_len[i])!=EOF)
    {
//         printf("s: %d %d %s %s\n",i,s_len[i],tempc1,argv[1]);
      num_base = num_base + s_len[i];
      i++;
    }
    fclose(namef);

    nSeq = i;
    num_steps = 0;
    BAR = 100;
    nstep = 100;
    for(i=0;i<nSeq;i++)
    {
/*     search reads with an index < i     */
/*     search reads with an index > i     */
       stopflag=0;
       j=i+1;
       base = s_len[i]; 
//         printf("www: %d %d %d %d\n",i,s_len[i],s_len[j],BAR);
       while((j<nSeq)&&(stopflag==0))
       {
         if((s_len[j]<(BAR+nstep))&&(s_len[i]>=BAR))
         {
           base = base + s_len[j];
           j++;
         }
         else
           stopflag=1;
       }
       if((j-i)>1)
       {
         rate = (j-i)*100;
//         rate = rate/num_base;
         rate = rate/nSeq;
         printf("frequency: %d %f\n",BAR,rate);
         BAR = BAR+nstep;
         num_steps++;
       }
       else if((j-i)==1)
       {
         rate = 100;
         rate = rate/nSeq;
         printf("frequency2: %d %f\n",BAR,rate);
         BAR = s_len[i];
         num_steps++;
       }
       i=j-1;
     }
     return EXIT_SUCCESS;
}


