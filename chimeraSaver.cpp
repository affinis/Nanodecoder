#include "fastq.hpp"

int main(int argc, char *argv[]){
  int option;
  const char* chimeraFastqFilename;
  const char* chimeraMappingStatFilename;
  const char* chimeraReadSummaryFilename;
  const char* barcodeFileName;
  const char* outFileName;
  int mismatch=1;
  int barcoderange=100;
  int tech;
  while((option=getopt(argc,argv,"f:m:s:d:l:t:b:"))!=-1){
    switch(option){
    case 'f':
      chimeraFastqFilename=optarg;
      break;
    case 'm':
      chimeraMappingStatFilename=optarg;
      break;
    case 's':
      chimeraReadSummaryFilename=optarg;
      break;
    case 'd':
      mismatch=atoi(optarg);
      break;
    case 'l':
      barcoderange=atoi(optarg);
      break;
    case 't':
      tech=atoi(optarg);
      break;
    case 'b':
      barcodeFileName=optarg;
      break;
    } // end switch
  }// end while
  int minSegment=mismatch+1;
  string firstline;
  fstream testfile;
  testfile.open(barcodeFileName,fstream::in);
  getline(testfile,firstline,'\n');
  int BarcodeLength=firstline.length();
  ReadFile chimeraReads(chimeraFastqFilename,chimeraMappingStatFilename,chimeraReadSummaryFilename);
  vector<int> SegmentsLengths=getMaxComplexitySegments(BarcodeLength,minSegment);
  BarcodeFile barcodeFile(barcodeFileName,SegmentsLengths,BarcodeLength);
  chimeraReads.identifyBarcodes(barcodeFile,SegmentsLengths,barcoderange,mismatch,tech,"chimeric_filtered.fastq.gz");
}
