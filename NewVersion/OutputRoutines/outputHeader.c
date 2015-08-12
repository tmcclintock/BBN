#include "outputHeader.h"

#include <stdlib.h>
#include <stdio.h>

void outputHeader(FILE *hp,int NUCS,int REACS, int STEPS){
  fwrite(&NUCS,sizeof(int),1,hp);
  fwrite(&REACS,sizeof(int),1,hp);
  fwrite(&STEPS,sizeof(int),1,hp);
  printf("Header file filled:\tnnuc=%i\tnreac=%i\tnprints=%i\n",NUCS,REACS,STEPS);
}
