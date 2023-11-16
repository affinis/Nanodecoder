#include "processMapping.hpp" 

int main(int argc, char *argv[]){
  int option;
  const char* Bamfilename;
  const char* GTFfilename;
  cout.precision(2);
  while((option=getopt(argc,argv,"b:g:"))!=-1){
    switch(option){
    case 'b':
      Bamfilename=optarg;
      break;
    case 'g':
      GTFfilename=optarg;
      break;
    } // end switch
  }// end while
  
  ofstream File;
  File.open("read_annotation.tsv", ios::out);
  
  unordered_map<string, gene_feature> gene_features;
  unordered_map<string, mapping_info> ReadMapping;
  unordered_map<string, map<int,pair<int,string>>> sorted_gene_features;
  
  readFeaturesFromGTF(GTFfilename,gene_features);
  sortFeaturesByCoord(gene_features,sorted_gene_features);
  readMappingFile(Bamfilename,ReadMapping,sorted_gene_features,gene_features,false,File);
  
  File.close();
}

