#include "processMapping.hpp" 

int main(int argc, char *argv[]){
  int option;
  int tech;
  const char* Bamfilename;
  const char* GTFfilename;
  bool loose_transcript=false;
  bool debug=false;
  cout.precision(2);
  while((option=getopt(argc,argv,"b:g:t:lx"))!=-1){
    switch(option){
    case 'b':
      Bamfilename=optarg;
      break;
    case 'g':
      GTFfilename=optarg;
      break;
    case 't':
      tech=atoi(optarg);
      break;
    case 'l': //in loose transcript mode, all transcripts from GTF file will be recorded even if they were not well supported 
      loose_transcript=true;
      break;
    case 'x':
      debug=true;
      break;
    } // end switch
  }// end while
  
  ofstream File;
  File.open("read_annotation.tsv", ios::out);
  
  unordered_map<string, gene_feature> gene_features;
  unordered_map<string, mapping_info> ReadMapping;
  unordered_map<string, map<int,vector<seq_feature>>> sorted_gene_features;
  
  readFeaturesFromGTF(GTFfilename,gene_features,loose_transcript);
  for(auto feature:gene_features){
    feature.second.debugPrintInfo(true);
  }
  sortFeaturesByCoord(gene_features,sorted_gene_features);
  readMappingFile(Bamfilename,ReadMapping,sorted_gene_features,gene_features, tech, debug,File);
  
  File.close();
}

