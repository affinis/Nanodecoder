#include <string>
#include <fstream>
#include <processSeq.hpp>
#include <unistd.h>


using std::string;
using std::fstream;

int main(int argc, char *argv[]){
  int option;
  const char* barcodefilename;
  const char* currenttablefilename;
  while((option=getopt(argc,argv,"f:t:"))!=-1){
    switch(option){
    case 'f':
      barcodefilename=optarg;
      break;
    case 't':
      currenttablefilename=optarg;
      break;
    } // end switch
  }// end while
  fstream barcodefile;
  barcodefile.open(barcodefilename,fstream::in);
  while(!barcodefile.eof()){
    string barcode;
  }
}