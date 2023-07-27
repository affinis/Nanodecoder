#include <unistd.h>
#include <stdlib.h>
#include <cstring>
#include <string>
#include <fstream>
#include <cstring>
#include <vector>
#include <unordered_map>
#include <iostream>

#include <htslib/sam.h>

using std::ifstream;
using std::cout;
using std::endl;
using std::string;

#define bam_is_umapped(b) (((b)->core.flag&BAM_FUNMAP) != 0)

int main(int argc, char *argv[]){
  int option;
  const char* Bamfilename;
  while((option=getopt(argc,argv,"f:"))!=-1){
    switch(option){
    case 'f':
      Bamfilename=optarg;
      break;
    } // end switch
  }// end while
  int r=0;
  htsFile *testfile=hts_open(Bamfilename,"r");
  bam1_t *b = NULL;
  bam_hdr_t *h = NULL;
  h = sam_hdr_read(testfile);
  if (h == NULL) {
      cout << "Cound not read header" << endl;
      return 0;
  }
  b=bam_init1();
  while ((r = sam_read1(testfile, h, b)) >= 0) {
    string flag;
    if(bam_is_rev(b)){
      flag="reverse";
    }else if(bam_is_umapped(b)){
      flag="*";
    }else{
      flag="forward";
    }
    cout << bam_get_qname(b) << "\t" << flag << endl;
  }
  bam_destroy1(b);
  bam_hdr_destroy(h);
}

