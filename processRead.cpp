#include <unistd.h>
#include "fastq.mod.hpp"

using std::string;

int main(int argc, char *argv[]){
  int option;
  const char* Readfilename;
  const char* Barcoedefilename;
  const char* Modelfilename;
  const char* Alignmentsummaryfilename;
  const char* previousFile=nullptr;
  string latest_read;
  int MaxMismatch=1;
  int MaxTSOMismatch=2;
  int MinR1Overlap=10;
  vector<int> WordSize;
  int BarcodeLength=16;
  int barcodeRange=100;
  int tech=5;
  int threads=1;
  bool use_dtw=false;
  bool use_anno=false;
  bool restart=false;
  while((option=getopt(argc,argv,"hf:b:m:a:r:l:t:d:c:@:"))!=-1){
    switch(option){
    case 'h':
      printf("nanodecoder, identify barcode in long reads, Lin Lyu, 2023\n\n\
              \t-h\tprint help\n\
              \t-f\tinput file\n\
              \t-b\tBC whitelist file\n\
              \t-m\tmismatch allowed\n\
              \t-a\talignmrnt summary file\n\
              \t-r\tBC search range\n\
              \t-t\ttech used 5' or 3'\n\
              \t-d\tusing DTW to rescue reads\n\
              \t-c\tcontinue from a broken run, give the filename.tsv produced before\n\
              \t-@\tthreads used\n\n");
      return 0;
    case 'f':
      Readfilename=optarg;
      break;
    case 'b':
      Barcoedefilename=optarg;
      break;
    case 'm':
      MaxMismatch=atoi(optarg);
      break;
    case 'a':
      fprintf(stderr,"Using read annotation.\n");
      use_anno=true;
      Alignmentsummaryfilename=optarg;
      break;
    case 'r':
      barcodeRange=atoi(optarg);
      break;
    case 't':
      tech=atoi(optarg);
      break;
    case 'd':
      fprintf(stderr,"Using DTW as supplementary method, this will be much slower.\n");
      use_dtw=true;
      Modelfilename=optarg;
      break;
    case 'c':
      restart=true;
      previousFile=optarg;
      break;
    case '@':
      threads=atoi(optarg);
      break;
    } // end switch
  }// end while
  
  if(MaxMismatch>BarcodeLength/2){
    cout << "Mismatch should not larger than half of barcode length" << endl;
    return 0;
    }
  int minSegment=MaxMismatch+1;
  vector<int> SegmentsLengths=getMaxComplexitySegments(BarcodeLength,minSegment);
  int r=3;
  WL_DB* wl_db=new WL_DB[60000];
  entry_t* model=new entry_t[4096];
  wl_db=init_WL_DB(60000, 1024);
  level_t mean_to_subtract=90;
  char qbase[100];
  strcpy(qbase,TENX_R1.c_str());
  //load whitelist for classical alignment based method
  BarcodeFile barcodeFile(Barcoedefilename,SegmentsLengths,BarcodeLength);
  unordered_map<string, ReadAnno> AnnoMap;
  fstream oldfile;
  fstream appendfile;
  if(restart){
    fprintf(stderr,"Restarting from previous run.\n");
    string lastline;
    oldfile.open(previousFile);
    while(!oldfile.eof()){
      string line;
      string this_read;
      getline(oldfile,line,'\n');
      if(line==""){
        continue;
      }else{
        this_read=tokenize(line,'\t')[0];
        if(latest_read==""){
          latest_read=this_read;
          continue;
        }else{
          if(this_read[this_read.length()-2]=='-'){
            continue;
          }else{
            latest_read=this_read;
          }
        }
      }
    }
    fprintf(stderr,"Last read of previous run is: %s\n",latest_read.c_str());
    oldfile.close();
    appendfile.open(previousFile,ios::app);
  }
  if(use_dtw){
    //load whitelist and current model for DTW based method
    model=read_pore_model(Modelfilename, mean_to_subtract);
    long wl_size=read_WL(Barcoedefilename,model, qbase, 0, r, wl_db);
    if(use_anno){
      readAlignmentSummaryFile(AnnoMap,Alignmentsummaryfilename);
    }
    ReadFile InputFile(Readfilename,barcodeFile,SegmentsLengths,barcodeRange,MaxMismatch,
                       tech,threads,AnnoMap,appendfile,model,wl_db,use_dtw,false,latest_read);
  }else{
    if(use_anno){
      readAlignmentSummaryFile(AnnoMap,Alignmentsummaryfilename);
    }
    ReadFile InputFile(Readfilename,barcodeFile,SegmentsLengths,barcodeRange,MaxMismatch,
                       tech,threads,AnnoMap,appendfile,nullptr,nullptr,use_dtw,false,latest_read);
    
  }
  if(restart){
    appendfile.close();
  }
  delete [] wl_db;
  wl_db=nullptr;
  delete [] model;
  model=nullptr;
  return 0;
}
