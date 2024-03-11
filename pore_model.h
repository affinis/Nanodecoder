#ifndef PORE_MODEL_H
#define PORE_MODEL_H

#ifndef level_t
#define level_t float
#endif


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define M_SIZE 4096
#define MER_LENGTH 6
typedef struct {
    level_t l_mean,l_sd,s_mean,s_sd;
} entry_t;

static inline entry_t *read_pore_model(const char *f, const level_t mean_to_subtract){
    FILE *file;
    char line[1024];
    const char *delim="\t";
    int i=0;
    entry_t *model= (entry_t *)calloc(M_SIZE, sizeof(entry_t));

    file = fopen(f, "r");
    if (file == NULL) {
        fprintf(stderr, "read_pore_model error: Unable to open the file %s\n", f);
        return NULL;
    }    
    while(fgets(line, sizeof(line), file)){
        line[strcspn(line, "\n")]=0; // remove newline
        char *token=strtok(line,delim);
        if(strcmp(token, "kmer") == 0 && i==0)
            continue;

        if(token != NULL){
            token=strtok(NULL, delim); //strtok looks like a hack, with NULL
            //model[i].l_mean = (strtod(token, NULL));
            model[i].l_mean = (strtod(token, NULL))-mean_to_subtract;
            token=strtok(NULL, delim);
            model[i].l_sd = (strtod(token, NULL));
            token=strtok(NULL, delim);
            model[i].s_mean = (strtod(token, NULL));
            token=strtok(NULL, delim);
            model[i].s_sd = (strtod(token, NULL));
        }
        i+=1;
    }
    fprintf(stderr,"read_pore_model: %i entries read.\n",i);
    return model;
}

static inline level_t mer2level(const entry_t *m, const char *mer){
  unsigned int x=0;
  if(strlen(mer) < MER_LENGTH){
    printf("mer size %lu\n", strlen(mer));
    return -1;
  }
  for(int i=0; i<MER_LENGTH; i++){
    switch(mer[i]){
    case 'a':
      x=(x<<2);break;
    case 'c':
      x=(x<<2)+1;break;
    case 'g':
      x=(x<<2)+2;break;
    case 't':
      x=(x<<2)+3;break;
    case 'A':
      x=(x<<2);break;
    case 'C':
      x=(x<<2)+1;break;
    case 'G':
      x=(x<<2)+2;break;
    case 'T':
      x=(x<<2)+3;break;
    default:
      fprintf(stderr,"base %c\n", mer[i]); return -2;
    }
  }
  return m[x].l_mean;
}

static inline void seq2level(const entry_t *m, char *seq, level_t *ret, long max_size, long append_inf){
  unsigned long l=strlen(seq);
  if(l>max_size+MER_LENGTH-1-append_inf)
    l=max_size+MER_LENGTH-1-append_inf;
  if(l < MER_LENGTH){
    printf("seq2level error: seq size %lu\n", l);
    return;
  }
  char *ptr=seq;
  for(long i=0; i<(l-MER_LENGTH+1); i++){
    ret[i]=mer2level(m, ptr);
    ptr++;
  }
  for(long i=l-MER_LENGTH+1; i< l-MER_LENGTH+1+append_inf; i++){
    ret[i]=10000;   ///hack, impossibly large level
  }
}

static inline void seq2merqual(char *seq, char *ret, long max_size, long append){
  unsigned long l=strlen(seq);
  if(l>max_size+MER_LENGTH-1-append-1)   // one fewer than level due to ending \0
    l=max_size+MER_LENGTH-1-append-1;
  if(l < MER_LENGTH){
    printf("seq2merqual error: seq size %lu\n", l);
    return;
  }
  ///this version simply returns the last qual of the kmer, assuming left to right base called
  for(long i=0; i<(l-MER_LENGTH+1); i++){
    ret[i]=seq[i+MER_LENGTH-1];
  }
  for(long i=l-MER_LENGTH+1; i< l-MER_LENGTH+1+append; i++){
    ret[i]=35;   ///hack, equals '#'
  }
  ret[l-MER_LENGTH+1+append]=0;
}

static inline void seq2merqual2(char *seq, char *ret, long max_size, long append){
  unsigned long l=strlen(seq);
  if(l>max_size+MER_LENGTH-1-append-1)   // one fewer than level due to ending \0
    l=max_size+MER_LENGTH-1-append-1;
  if(l < MER_LENGTH){
    printf("seq2merqual error: seq size %lu\n", l);
    return;
  }
  ///this version returns the mer average
  int tot=0;
  for(long i=0; i<MER_LENGTH; i++){
  	tot+=seq[i];
  }
  ret[0]=tot/MER_LENGTH;
  for(long i=1; i<(l-MER_LENGTH+1); i++){
    tot-=seq[i-1];
    tot+=seq[i+MER_LENGTH-1];
    ret[i]=tot/MER_LENGTH;
  }
  for(long i=l-MER_LENGTH+1; i< l-MER_LENGTH+1+append; i++){
    ret[i]=35;   ///hack, equals '#'
  }
  ret[l-MER_LENGTH+1+append]=0;
}


static inline void printLevels(const level_t *l, const size_t length){
  for(int i=0; i<length; i++){
    printf(" %4.2f", (float) l[i]);
  }
  printf("\n");
}

static void revcom(char *s){
  int i;
  char tmp;
  int len=strlen(s);
  for(i=0; i<floor(len/2);i++){
    tmp=s[i];
    s[i]=s[len-i-1];
    s[len-i-1]=tmp;
  }
  for(i=0; i<len; i++){
  switch(s[i]){
    case 'a':
      s[i]='T';break;
    case 'c':
      s[i]='G';break;
    case 'g':
      s[i]='C';break;
    case 't':
      s[i]='A';break;
    case 'A':
      s[i]='T';break;
    case 'C':
      s[i]='G';break;
    case 'G':
      s[i]='C';break;
    case 'T':
      s[i]='A';break;
    default:
      s[i]='N';
  }
  }
}

#endif
