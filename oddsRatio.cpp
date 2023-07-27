#include "fastq.hpp"
#include <regex>

int main(int argc, char *argv[]){
  int option;
  const char* filename;
  while((option=getopt(argc,argv,"f:"))!=-1){
    switch(option){
    case 'f':
      filename=optarg;
      break;
    } // end switch
  }// end while
  fstream outFile;
  outFile.open("read_status.tsv",fstream::out|fstream::app);
  outFile << "read_id\tlength\tnum_polyA\tnum_polyT\tnum_R1\tR1_indices\tR1_strand\tnum_TSO\tTSO_indices\tTSO_strand" << endl;
  ReadFile Reads=ReadFile(filename);
  int readCount=0;
  std::regex patternPolyA("[GCT]A{13,13}");
  std::regex patternPolyT("[GCA]T{13,13}");
  for(Read read:Reads.Reads){
//    cout << read.ID << endl;
    // count polyA/T number
    int countPolyA = 0;
    int countPolyT = 0;
    for (std::regex_iterator<std::string::iterator> i(read.Seq.begin(), read.Seq.end(), patternPolyA); i != std::regex_iterator<std::string::iterator>(); ++i) {
      countPolyA++;
    }
    for (std::regex_iterator<std::string::iterator> i(read.Seq.begin(), read.Seq.end(), patternPolyT); i != std::regex_iterator<std::string::iterator>(); ++i) {
      countPolyT++;
    }
    
    //count TSO number, strand, dists
    vector<pair<int,int>> TSO_info;
    string TSO_rc=reverse_complement(TENX_TSO);
    int TSO_len=TENX_TSO.length();
    vector<string> kmers4tso=genrateKmerFromSeq(read.Seq,TSO_len);
    int kmer4tsoIndex=0;
    while(kmer4tsoIndex<kmers4tso.size()){
      if(editDistance(kmers4tso[kmer4tsoIndex],TENX_TSO,1)<=1){
        pair<int, int> kmer_info_particular={kmer4tsoIndex+1,1};
        TSO_info.push_back(kmer_info_particular);
        kmer4tsoIndex=kmer4tsoIndex+13;
      }else if(editDistance(kmers4tso[kmer4tsoIndex],TSO_rc,1)<=1){
        pair<int, int> kmer_info_particular={kmer4tsoIndex+1,-1};
        TSO_info.push_back(kmer_info_particular);
        kmer4tsoIndex=kmer4tsoIndex+13;
      }else{
        kmer4tsoIndex++;
      }
    }
    
    //count adapter number, strand, dists
    vector<pair<int,int>> R1_info;
    string TENX_R1_rc=reverse_complement(TENX_R1);
    int R1_len=TENX_R1.length();
    vector<string> kmers4r1=genrateKmerFromSeq(read.Seq,R1_len);
    int kmer4r1Index=0;
    while(kmer4r1Index<kmers4r1.size()){
      if(editDistance(kmers4r1[kmer4r1Index],TENX_R1,2)<=2){
        pair<int, int> kmer_info_particular={kmer4r1Index+1,1};
        R1_info.push_back(kmer_info_particular);
        kmer4r1Index=kmer4r1Index+22;
      }else if(editDistance(kmers4r1[kmer4r1Index],TENX_R1_rc,2)<=2){
        pair<int, int> kmer_info_particular={kmer4r1Index+1,-1};
        R1_info.push_back(kmer_info_particular);
        kmer4r1Index=kmer4r1Index+22;
      }else{
        kmer4r1Index++;
      }
    }
    
    string r1_indices_str;
    string r1_strand_str;
    for(pair<int, int> r1:R1_info){
      r1_indices_str=r1_indices_str+to_string(r1.first)+",";
      r1_strand_str=r1_strand_str+to_string(r1.second)+",";
    }
    if(R1_info.size()==0){
      r1_indices_str="*";
      r1_strand_str="*";
    }
    
    
    string tso_indices_str;
    string tso_strand_str;
    for(pair<int, int> tso:TSO_info){
      tso_indices_str=tso_indices_str+to_string(tso.first)+",";
      tso_strand_str=tso_strand_str+to_string(tso.second)+",";
    }
    if(TSO_info.size()==0){
      tso_indices_str="*";
      tso_strand_str="*";
    }
    
    outFile << read.ID << "\t" << read.Seq.length() << "\t" << countPolyA << "\t" << countPolyT << "\t";
    outFile << R1_info.size() << "\t" << r1_indices_str << "\t" << r1_strand_str << "\t";
    outFile << TSO_info.size() << "\t" << tso_indices_str << "\t" << tso_strand_str << endl;
  }
  outFile.close();
}