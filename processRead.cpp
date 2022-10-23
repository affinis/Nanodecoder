#include <unistd.h>
#include "fastq.hpp"

using std::string;


int main(int argc, char *argv[]){
  int option;
  const char* Readfilename;
  const char* Barcoedefilename;
  int MaxMismatch=1;
  int MaxTSOMismatch=2;
  int MinR1Overlap=10;
  vector<int> WordSize;
  int K=16;
  int BarcodeLength=16;
  int barcodeRange=100;
  int tech=5;
  while((option=getopt(argc,argv,"f:b:m:k:r:l:t:"))!=-1){
    switch(option){
    case 'f':
      Readfilename=optarg;
      break;
    case 'b':
      Barcoedefilename=optarg;
      break;
    case 'm':
      MaxMismatch=atoi(optarg);
      break;
    case 'k':
      K=atoi(optarg);
      break;
    case 'r':
      barcodeRange=atoi(optarg);
      break;
    case 'l':
      BarcodeLength=atoi(optarg);
      break;
    case 't':
      tech=atoi(optarg);
      break;
    }
  }
  if(K>barcodeRange){
    cout << "Barcode detect range should be larger than k" << endl;
    return 0;
  }
  if(MaxMismatch>BarcodeLength/2){
    cout << "Mismatch should not larger than half of barcode length" << endl;
    return 0;
    }
  int minSegment=MaxMismatch+1;
  vector<int> SegmentsLengths=getMaxComplexitySegments(BarcodeLength,minSegment);
  cout << "Loading barcodes and building dicts..." << endl;
  BarcodeFile barcodeFile(Barcoedefilename,SegmentsLengths,BarcodeLength);
  cout << "Dicts done, " << barcodeFile.barcodes.size() << " barcodes loaded." << endl;
  ReadFile InputFile(Readfilename,barcodeFile,SegmentsLengths,barcodeRange,MaxMismatch,tech);
  return 0;
}
