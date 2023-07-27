#include <string>
#include <string_view>
#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <errno.h>
#include <unordered_map>
#include <regex>
#include <thread>
#include <omp.h>

#include "kseq.h"
#include "prosessSeq.hpp"

using std::unordered_map;
using std::ifstream;
using std::pair;
using std::to_string;
using std::fstream;
using std::regex;
using std::smatch;
using std::regex_search;
using std::thread;

//firstly sample from reads, and find TSO and r1 to determine barcode range

KSEQ_INIT(gzFile, gzread)

unordered_map<int, string> num_to_strand{
  {0,"original"},{1,"reverse"},{2,"complement"},{3,"reverse_complement"},{-1,"*"}
  };

string TENX_TSO="TTTCTTATATGGG";
string TENX_R1="CTACACGACGCTCTTCCGATCT";
string POLYT="TTTTTTTTTTTTT";

class BarcodeFile{
public:
  ifstream barcodeFile;
  vector<string> barcodes;
  unordered_map<string, vector<pair<int, int>>> barcode_dict;
  unordered_map<string, string> barcode_ori;
  unordered_map<string, string> barcode_rc;
  
  BarcodeFile(){
    
  }
  
  BarcodeFile(const BarcodeFile& source){
    barcodes=source.barcodes;
    barcode_dict=source.barcode_dict;
    barcode_ori=source.barcode_ori;
    barcode_rc=source.barcode_rc;
  }
  
  BarcodeFile(const char* filename, vector<int> word_sizes, int barcode_length){
    barcodeFile.open(filename);
    while(!barcodeFile.eof()){
      string barcode;
      getline(barcodeFile,barcode,'\n');
      if(barcode==""){
        continue;
      }
      if(barcode.length()!=barcode_length){
        cout << "Barcode " << barcodes.size()+1 << ", " << barcode << ", length incorrect" << endl;
        cout << "Barcode length not correct, the default length is 16, if length of barcodes in the white list provided is not 16, provide actual length using -l" << endl;
        barcodeFile.close();
        exit(0);
      }
      string barcode_rc_seq=reverse_complement(barcode);
      barcodes.push_back(barcode);
      
      barcode_ori[barcode]=barcode;
      barcode_rc[barcode_rc_seq]=barcode;
      
      vector<string> barcode_segments=getMaxComplexitySegments(barcode,word_sizes);
      vector<string> barcode_segments_rc=getMaxComplexitySegments(barcode_rc_seq,word_sizes);
      for(int j=0;j<barcode_segments.size();j++){
        string key;
        pair<int, int> barcode_value;
        key=barcode_segments[j]+to_string(j)+to_string(0);
        barcode_value.first=barcodes.size();
        barcode_value.second=0;
        barcode_dict[key].push_back(barcode_value);
      }
      for(int j=0;j<barcode_segments_rc.size();j++){
        string key;
        pair<int, int> barcode_value;
        key=barcode_segments[j]+to_string(j)+to_string(3);
        barcode_value.first=barcodes.size();
        barcode_value.second=3;
        barcode_dict[key].push_back(barcode_value);
      }
    }
    barcodeFile.close();
  }
  
  /*  
   void printWordToBarcode(string word){
   cout << "Printing barcodes containing " << word << " on oringial strand.." << endl;
   for(int i=0;i<barcode_dict[word].size();i++){
   cout << barcode_dict[word][i] << " ";
   }
   cout << endl;
   }
   */
};


class Read{
public:
  string ID;
  string Seq;
  string Quality;
  string UMI="*";
  string Seq_trimmed;
  string Quality_trimmed;
  int OUTER;
  int INNER;
  vector<string> Kmers;
  int barcodeStart;
  string barcode;
  string barcodeStrand;
  int fusionpoint;
  int barcodeMismatch;
  string stat;
  int polyAT_containing=0;
  int barcode_in_tolerance=0;
  vector<int> barcode_intolerance_indices;
  vector<string> barcode_intolerance_detail;
  vector<string> matched_kmer_quality;
  vector<int> barcode_mismatches;
  
  Read(){
    ID="defalut";
    Seq="NNNNN";
    Quality="NNNNN";
  }
  
  Read(const Read* source){
    ID=source->ID;
    Seq=source->Seq;
    Quality=source->Quality;
  }
  
  void printBarcodeInfo(){
    string barcode_intolerance_detail_str="";
    string barcode_intolerance_indices_str="";
    string barcode_mismatches_str="";
    if(barcode_intolerance_detail.size()==0){
      barcode_intolerance_detail_str="*";
      barcode_intolerance_indices_str="*";
      string barcode_mismatches_str="";
    }else{
      for(int i=0;i<barcode_intolerance_detail.size();++i){
        barcode_intolerance_detail_str=barcode_intolerance_detail_str+barcode_intolerance_detail[i]+",";
        barcode_intolerance_indices_str=barcode_intolerance_indices_str+std::to_string(barcode_intolerance_indices[i])+",";
        barcode_mismatches_str=barcode_mismatches_str+std::to_string(barcode_mismatches[i])+",";
      }
    }
    cout << ID << "\t" << stat << "\t" << OUTER << "\t";
    cout << barcode << "\t" << barcodeStart << "\t" << INNER;
    cout << "\t" << barcodeStrand << "\t" << barcodeMismatch << "\t" << polyAT_containing << "\t" << barcode_in_tolerance;
    cout << "\t" << barcode_intolerance_detail_str << "\t" << barcode_intolerance_indices_str << "\t" << barcode_mismatches_str << endl;
  }
  
  void printBarcodeInfo(std::stringstream &stream){
    string barcode_intolerance_detail_str="";
    string barcode_intolerance_indices_str="";
    string barcode_mismatches_str="";
    if(barcode_intolerance_detail.size()==0){
      barcode_intolerance_detail_str="*";
      barcode_intolerance_indices_str="*";
      barcode_mismatches_str="*";
    }else{
      for(int i=0;i<barcode_intolerance_detail.size();++i){
        barcode_intolerance_detail_str=barcode_intolerance_detail_str+barcode_intolerance_detail[i]+",";
        barcode_intolerance_indices_str=barcode_intolerance_indices_str+std::to_string(barcode_intolerance_indices[i])+",";
        barcode_mismatches_str=barcode_mismatches_str+std::to_string(barcode_mismatches[i])+",";
      }
    }
    //1. read_id 2. status 3. OCS_start 4. barcode 5. barcode_start 6. ICS_start 7. barcode_strand 8. barcode_mismatch 9. polyAT_containing
    //10 num_barcode 11. barcodes 12. barcode_starts 13. barcode_mismatches
    stream << ID << "\t" << stat << "\t" << OUTER << "\t";
    stream << barcode << "\t" << barcodeStart << "\t" << INNER;
    stream << "\t" << barcodeStrand << "\t" << barcodeMismatch << "\t" << polyAT_containing << "\t" << barcode_in_tolerance;
    stream << "\t" << barcode_intolerance_detail_str << "\t" << barcode_intolerance_indices_str << "\t" << barcode_mismatches_str << endl;
  }
  
  void fq_gz_write(gzFile out_file, bool trim=false) {
    std::stringstream stream;
    if(trim){
      string name=this->ID+"_"+this->barcode+"_"+this->UMI;
      string seq=this->Seq_trimmed;
      string qual=this->Quality_trimmed;
      stream << "@" << name << "\n" << 
        seq << "\n" << 
          "+" << "\n" << 
            qual << "\n";
    }else{
      string name=this->ID+"_"+this->barcode;
      stream << "@" << name << "\n" << 
      this->Seq << "\n" << 
        "+" << "\n" << 
          this->Quality << "\n";
    }
    gzputs(out_file, stream.str().c_str());
  }
  
  void fq_gz_write(std::stringstream &stream, bool trim=false) {
    if(trim){
      string name=this->ID+"_"+this->barcode+"_"+this->UMI;
      string seq=this->Seq_trimmed;
      string qual=this->Quality_trimmed;
      stream << "@" << name << "\n" << 
        seq << "\n" << 
          "+" << "\n" << 
            qual << "\n";
    }else{
      string name=this->ID+"_"+this->barcode;
      stream << "@" << name << "\n" << 
        this->Seq << "\n" << 
          "+" << "\n" << 
            this->Quality << "\n";
    }
  }
  
  int checkPolyAT(){
    string testseq=Seq;
    int ployA=testseq.find("AAAAAAAAAA");
    int polyT=testseq.find("TTTTTTTTTT");
    if(ployA>=0&&polyT<0){
      return 0;
    }else if(ployA<0&&polyT>=0){
      return 3;
    }else{
      return -1;
    }
  } //
  
  int localAlign(string& query, string& subject, int max_mismatch, char preference){
    int alignStart=0;
    int min_dist_align_start=0;
    int qlen=query.length();
    pair<int, int> best_align={-1, qlen};
    int max_indel_num=max_mismatch;
    vector<string> kmers=genrateKmerFromSeq(subject, qlen);
    for(string kmer:kmers){
      int dist=editDistance(query,kmer,max_indel_num);
      if(dist<=max_mismatch&&dist<=best_align.second){
        best_align.first=alignStart;
        best_align.second=dist;
        if(best_align.second==0&&preference=='l'){
          return best_align.first;
          break;
        }else{
          alignStart++;
          continue;
        }
      }else{
        alignStart++;
        continue;
      }
    }
    return best_align.first;
  }
  
  pair<int, int> isConstantSequenceContaining(string& constantseq, int testRange, int max_testseq_mismatch, int strand){
    string testForward;
    string testseq;
    string testReverse;
    string testseq_rc;
    pair<int, int> res;
//    vector<string> testsegF;
//    vector<string> testsegR;
//    int segments=max_testseq_mismatch+1;
    int index1=-1;
    int index2=-1;
    if(strand==0){
      testForward=this->Seq.substr(0,testRange);
      testseq=constantseq;
      index1=localAlign(testseq,testForward,max_testseq_mismatch,'l');
      /*
      testsegF=getMaxComplexitySegments(testseq,getMaxComplexitySegments(testseq.length(),segments));
      int i=0;
      while(i<testsegF.size()&&index1<0&&index2<0){
        index1=testForward.find(testsegF[i]);
        i++;
      } //while
       */
    }else if(strand==3){
      testReverse=this->Seq.substr(this->Seq.length()-testRange,testRange);
      testseq_rc=reverse_complement(constantseq);
      int align_result=localAlign(testseq_rc,testReverse,max_testseq_mismatch,'r');
      if(align_result<0){
        index2=align_result;
      }else{
        index2=align_result+this->Seq.length()-testRange;
      }
      /*
      testsegR=getMaxComplexitySegments(testseq_rc,getMaxComplexitySegments(testseq_rc.length(),segments));
      int i=0;
      while(i<testsegR.size()&&index1<0&&index2<0){
        index2=testReverse.find(testsegR[i])+this->Seq.length()-testRange;
        i++;
      } //while
       */
    }else{
      testForward=this->Seq.substr(0,testRange);
      testseq=constantseq;
//      testsegF=getMaxComplexitySegments(testseq,getMaxComplexitySegments(testseq.length(),segments));
      testReverse=this->Seq.substr(this->Seq.length()-testRange,testRange);
      testseq_rc=reverse_complement(constantseq);
      index1=localAlign(testseq,testForward,max_testseq_mismatch,'l');
      int rv_align_result=localAlign(testseq_rc,testReverse,max_testseq_mismatch,'r');
      if(rv_align_result<0){
        index2=rv_align_result;
      }else{
        index2=rv_align_result+this->Seq.length()-testRange;
      }
//      testsegR=getMaxComplexitySegments(testseq_rc,getMaxComplexitySegments(testseq_rc.length(),segments));
//      int j=0;
//      while(j<testsegF.size()&&index1<0){
//        index1=testForward.find(testsegF[j]);
//        j++;
//      } //while
//      int i=0;
//      while(i<testsegR.size()&&index2<0){
//        index2=testReverse.find(testsegR[i]);
//        i++;
//      } //while
//      if(index2>=0){
//        index2=index2+Seq.length()-testRange;
//      }
    }
    res.first=index1;
    res.second=index2;
    return res;
  } //isConstantSequenceContaining
  
  void genrateKmer(int k){
    int kmerSize=k;
      const int len=this->Seq.length();
      for(int j=0;j+k<=len;j++){
        string kmer=this->Seq.substr(j,k);
        //        cout << kmer << endl;
        Kmers.push_back(kmer);
      }
    }
  
  int identifyBarcodes(BarcodeFile& barcodes, vector<int> segments, int maxDistToEnds, int maxMismatch, int tech){
    int barcode_length=barcodes.barcodes[0].length();
    int kmer_index_start=-1;
    int kmer_index_end=-1;
    int constant_seq_strand;
    int umiLength;
    pair<int, int> inner_stat;
    pair<int, int> outer_stat;
    string strand_by_polyAT;
    
    if(tech==3){
      umiLength=12;
    }else{
      umiLength=10;
    }
    
    //if read shorter than barcode search range * 2, discard
    if(this->Seq.length()<=maxDistToEnds*2){
      this->INNER=-1;
      this->OUTER=-1;
      this->barcode="*";
      this->barcodeStrand="*";
      this->barcodeStart=-1;
      this->barcodeMismatch=-1;
      this->stat="Read_too_short";
//      this->printBarcodeInfo();
      return 0;
    }
    if(this->Kmers.size()==0){
      this->genrateKmer(barcode_length);
    }
    //if (find polyA and 5')||(find polyT and 3'), means the barcode should be at left end
    if((checkPolyAT()==0&&tech==5)||(checkPolyAT()==3&&tech==3)){
      constant_seq_strand=0;
      this->polyAT_containing=1;
      if(tech==5){
        inner_stat=isConstantSequenceContaining(TENX_TSO,maxDistToEnds,2,0);
        outer_stat=isConstantSequenceContaining(TENX_R1,maxDistToEnds,3,0);
      }else if(tech==3){
        inner_stat=isConstantSequenceContaining(POLYT,maxDistToEnds,0,0);
        outer_stat=isConstantSequenceContaining(TENX_R1,maxDistToEnds,3,0);
      }else{
        cout << "Failed, 5' or 3' not specified or specified a wrong argument";
      }
      if(inner_stat.first<0&&outer_stat.first<0){
        this->INNER=-1;
        this->OUTER=-1;
        kmer_index_start=0;
        kmer_index_end=maxDistToEnds;
      }else if(inner_stat.first>=0&&outer_stat.first<0){
        this->INNER=inner_stat.first+1;
        this->OUTER=-1;
        kmer_index_start=0;
        kmer_index_end=inner_stat.first;
      }else if(inner_stat.first<0&&outer_stat.first>=0){
        this->INNER=-1;
        this->OUTER=outer_stat.first+1;
        kmer_index_start=outer_stat.first;
        kmer_index_end=maxDistToEnds;
      }else{
        this->INNER=inner_stat.first+1;
        this->OUTER=outer_stat.first+1;
        kmer_index_start=outer_stat.first;
        kmer_index_end=inner_stat.first;
      }
      //if (find polyT and 5')||(find polyA and 3'), means the barcode should be at right end
    }else if((checkPolyAT()==3&&tech==5)||(checkPolyAT()==0&&tech==3)){
      this->polyAT_containing=1;
      constant_seq_strand=3;
      if(tech==5){
        inner_stat=isConstantSequenceContaining(TENX_TSO,maxDistToEnds,2,3);
        outer_stat=isConstantSequenceContaining(TENX_R1,maxDistToEnds,3,3);
      }else if(tech==3){
        inner_stat=isConstantSequenceContaining(POLYT,maxDistToEnds,0,3);
        outer_stat=isConstantSequenceContaining(TENX_R1,maxDistToEnds,3,3);
      }else{
        cout << "Failed, 5' or 3' not specified or specified a wrong argument";
      }
      if(inner_stat.second<0&&outer_stat.second<0){
        this->INNER=-1;
        this->OUTER=-1;
        kmer_index_start=this->Seq.length()-maxDistToEnds;
        kmer_index_end=this->Seq.length()-barcode_length;
      }else if(inner_stat.second>=0&&outer_stat.second<0){
        this->INNER=inner_stat.second+1;
        this->OUTER=-1;
        kmer_index_start=inner_stat.second;
        kmer_index_end=this->Seq.length()-barcode_length;
      }else if(inner_stat.second<0&&outer_stat.second>=0){
        this->INNER=-1;
        this->OUTER=outer_stat.second+1;
        kmer_index_start=this->Seq.length()-maxDistToEnds;
        kmer_index_end=outer_stat.second;
      }else{
        this->INNER=inner_stat.second+1;
        this->OUTER=outer_stat.second+1;
        kmer_index_start=inner_stat.second;
        kmer_index_end=outer_stat.second;
      }
    }else{
      this->polyAT_containing=0;
      if(tech==5){
        inner_stat=isConstantSequenceContaining(TENX_TSO,maxDistToEnds,2,-1);
        outer_stat=isConstantSequenceContaining(TENX_R1,maxDistToEnds,3,-1);
      }else if(tech==3){
        inner_stat=isConstantSequenceContaining(POLYT,maxDistToEnds,0,-1);
        outer_stat=isConstantSequenceContaining(TENX_R1,maxDistToEnds,3,-1);
      }else{
        cout << "Failed, 5' or 3' not specified or specified a wrong argument";
      }
      if(inner_stat.first<0&&outer_stat.first<0&&inner_stat.second<0&&outer_stat.second<0){
        this->INNER=-1;
        this->OUTER=-1;
        this->barcode="*";
        this->barcodeStrand="*";
        this->barcodeStart=-1;
        this->barcodeMismatch=-1;
        this->stat="ICS_R1_missing";
//        this->printBarcodeInfo();
        return 0;
      }else if(inner_stat.first<0&&outer_stat.first>=0&&inner_stat.second<0&&outer_stat.second<0){
        this->INNER=-1;
        this->OUTER=outer_stat.first+1;
        kmer_index_start=outer_stat.first;
        kmer_index_end=maxDistToEnds;
        constant_seq_strand=0;
      }else if(inner_stat.first<0&&outer_stat.first<0&&inner_stat.second<0&&outer_stat.second>=0){
        this->INNER=-1;
        this->OUTER=outer_stat.second+1;
        kmer_index_start=this->Seq.length()-maxDistToEnds;
        kmer_index_end=outer_stat.second;
        constant_seq_strand=3;
      }else if(inner_stat.first<0&&outer_stat.first>=0&&inner_stat.second<0&&outer_stat.second>=0){
        this->INNER=-1;
        this->OUTER=-2;
        this->barcode="*";
        this->barcodeStrand="*";
        this->barcodeStart=-1;
        this->barcodeMismatch=-1;
        this->stat="ICS_missing_R1_both_ends";
//        this->printBarcodeInfo();
        return 0;
      }else if(inner_stat.first>=0&&outer_stat.first<0&&inner_stat.second<0&&outer_stat.second<0){
        this->INNER=inner_stat.first+1;
        this->OUTER=-1;
        kmer_index_start=0;
        kmer_index_end=inner_stat.first;
        constant_seq_strand=0;
      }else if(inner_stat.first>=0&&outer_stat.first>=0&&inner_stat.second<0&&outer_stat.second<0){
        this->INNER=inner_stat.first+1;
        this->OUTER=outer_stat.first+1;
        kmer_index_start=outer_stat.first;
        kmer_index_end=inner_stat.first;
        constant_seq_strand=0;
      }else if(inner_stat.first>=0&&outer_stat.first<0&&inner_stat.second<0&&outer_stat.second>=0){
        this->INNER=inner_stat.first+1;
        this->OUTER=outer_stat.second+1;
        this->barcode="*";
        this->barcodeStrand="*";
        this->barcodeStart=-1;
        this->barcodeMismatch=-1;
        this->stat="ICS_R1_different_end";
//        this->printBarcodeInfo();
        return 0;
      }else if(inner_stat.first>=0&&outer_stat.first>=0&&inner_stat.second<0&&outer_stat.second>=0){
        this->INNER=inner_stat.first+1;
        this->OUTER=outer_stat.first+1;
        kmer_index_start=outer_stat.first;
        kmer_index_end=inner_stat.first;
        constant_seq_strand=0;
      }else if(inner_stat.first<0&&outer_stat.first<0&&inner_stat.second>=0&&outer_stat.second<0){
        this->INNER=inner_stat.second+1;
        this->OUTER=-1;
        kmer_index_start=inner_stat.second;
        kmer_index_end=this->Seq.length()-barcode_length;
        constant_seq_strand=3;
      }else if(inner_stat.first<0&&outer_stat.first>=0&&inner_stat.second>=0&&outer_stat.second<0){
        this->INNER=inner_stat.second+1;
        this->OUTER=outer_stat.first+1;
        this->barcode="*";
        this->barcodeStrand="*";
        this->barcodeStart=-1;
        this->barcodeMismatch=-1;
        this->stat="ICS_R1_different_end";
//        this->printBarcodeInfo();
        return 0;
      }else if(inner_stat.first<0&&outer_stat.first<0&&inner_stat.second>=0&&outer_stat.second>=0){
        this->INNER=inner_stat.second+1;
        this->OUTER=outer_stat.second+1;
        kmer_index_start=inner_stat.second;
        kmer_index_end=outer_stat.second;
        constant_seq_strand=3;
      }else if(inner_stat.first<0&&outer_stat.first>=0&&inner_stat.second>=0&&outer_stat.second>=0){
        this->INNER=inner_stat.second+1;
        this->OUTER=outer_stat.second+1;
        kmer_index_start=inner_stat.second;
        kmer_index_end=outer_stat.second;
        constant_seq_strand=3;
      }else if(inner_stat.first>=0&&outer_stat.first<0&&inner_stat.second>=0&&outer_stat.second<0){
        this->INNER=-2;
        this->OUTER=-1;
        this->barcode="*";
        this->barcodeStrand="*";
        this->barcodeStart=-1;
        this->barcodeMismatch=-1;
        this->stat="ICS_both_ends_R1_missing";
//        this->printBarcodeInfo();
        return 0;
      }else if(inner_stat.first>=0&&outer_stat.first>=0&&inner_stat.second>=0&&outer_stat.second<0){
        this->INNER=inner_stat.first+1;
        this->OUTER=outer_stat.first+1;
        kmer_index_start=outer_stat.first;
        kmer_index_end=inner_stat.first;
        constant_seq_strand=0;
      }else if(inner_stat.first>=0&&outer_stat.first<0&&inner_stat.second>=0&&outer_stat.second>=0){
        this->INNER=inner_stat.second+1;
        this->OUTER=outer_stat.second+1;
        kmer_index_start=inner_stat.second;
        kmer_index_end=outer_stat.second;
        constant_seq_strand=3;
      }else{
        if(inner_stat.first>outer_stat.first&&inner_stat.second<outer_stat.second){
          this->INNER=-2;
          this->OUTER=-2;
          this->barcode="*";
          this->barcodeStrand="*";
          this->barcodeStart=-1;
          this->barcodeMismatch=-1;
          this->stat="ICS_R1_both_ends";
//          this->printBarcodeInfo();
          return 0;
        }else if(inner_stat.first<outer_stat.first&&inner_stat.second<outer_stat.second){
          this->INNER=inner_stat.second+1;
          this->OUTER=outer_stat.second+1;
          kmer_index_start=inner_stat.second;
          kmer_index_end=outer_stat.second;
          constant_seq_strand=3;
        }else if(inner_stat.first>outer_stat.first&&inner_stat.second>outer_stat.second){
          this->INNER=inner_stat.first+1;
          this->OUTER=outer_stat.first+1;
          kmer_index_start=outer_stat.first;
          kmer_index_end=inner_stat.first;
          constant_seq_strand=0;
        }else{
          this->INNER=-1;
          this->OUTER=-1;
          this->barcode="*";
          this->barcodeStrand="*";
          this->barcodeStart=-1;
          this->barcodeMismatch=-1;
          this->stat="ICS_R1_wrongly_oriented";
//          this->printBarcodeInfo();
          return 0;
        }
      }
    }
    if(kmer_index_start>kmer_index_end){
      this->barcode="*";
      this->barcodeStrand="*";
      this->barcodeStart=-1;
      this->barcodeMismatch=-1;
      this->stat="ICS_R1_wrongly_oriented";
//      this->printBarcodeInfo();
      return 0;
    }
//    cout << this->ID << "," << this->OUTER << "," << this->INNER << endl;
//    cout << "determine kmers by inner and outer constant seq start index" << endl;
    std::vector<string> kmers;
    //determine kmers by inner and outer constant seq start index
    
    //1. both CS found at left end
    if(this->INNER>0&&this->OUTER>0&&this->INNER<Seq.length()/2){
      kmer_index_start=this->OUTER-1+TENX_R1.length()-4;
      kmers.insert(kmers.end(),this->Kmers.begin()+this->OUTER-1+TENX_R1.length()-4,this->Kmers.begin()+this->OUTER-1+TENX_R1.length()+4);
      
    //2. both CS found at right end
    }else if(this->INNER>0&&this->OUTER>0&&this->INNER>Seq.length()/2){
      kmer_index_start=this->OUTER-1-barcode_length-4;
      kmers.insert(kmers.end(),this->Kmers.begin()+this->OUTER-1-barcode_length-4,this->Kmers.begin()+this->OUTER-1-barcode_length+4);
      
    //3. only ICS found at left end
    }else if(this->INNER>0&&this->OUTER<0&&this->INNER<Seq.length()/2){
      if(this->INNER-1-barcode_length-umiLength-4<0){
        this->barcode="*";
        this->barcodeStrand="*";
        this->barcodeStart=-1;
        this->barcodeMismatch=-1;
        this->stat="CS_too_close_to_end_or_to_each_other";
        return 0;
      }else{
        kmer_index_start=this->INNER-1-barcode_length-umiLength-4;
        kmers.insert(kmers.end(),this->Kmers.begin()+this->INNER-1-barcode_length-umiLength-4,this->Kmers.begin()+this->INNER-1-barcode_length-umiLength+4);
      }
      
    //4. only OCS found at left end
    }else if(this->INNER<0&&this->OUTER>0&&this->OUTER<Seq.length()/2){
      kmer_index_start=this->OUTER-1+TENX_R1.length()-4;
      kmers.insert(kmers.end(),this->Kmers.begin()+this->OUTER-1+TENX_R1.length()-4,this->Kmers.begin()+this->OUTER-1+TENX_R1.length()+4);
      
    //5. only ICS found at right end
    }else if(this->INNER>0&&this->OUTER<0&&this->INNER>Seq.length()/2){
      if(this->INNER+TENX_TSO.length()+umiLength+barcode_length+4-1>this->Seq.length()){
        this->barcode="*";
        this->barcodeStrand="*";
        this->barcodeStart=-1;
        this->barcodeMismatch=-1;
        this->stat="CS_too_close_to_end_or_to_each_other";
        return 0;
      }else{
        kmer_index_start=this->INNER-1+TENX_TSO.length()+umiLength-4;
        kmers.insert(kmers.end(),this->Kmers.begin()+this->INNER-1+TENX_TSO.length()+umiLength-4,this->Kmers.begin()+this->INNER-1+TENX_TSO.length()+umiLength+4);
      }
      
    //6. only OCS found at right end
    }else if(this->INNER<0&&this->OUTER>0&&this->OUTER>Seq.length()/2){
      kmer_index_start=this->OUTER-1-barcode_length-4;
      kmers.insert(kmers.end(),this->Kmers.begin()+this->OUTER-1-barcode_length-4,this->Kmers.begin()+this->OUTER-1-barcode_length+4);
      
    //7. both CS not found
    }else{
      this->barcode="*";
      this->barcodeStrand="*";
      this->barcodeStart=-1;
      this->barcodeMismatch=-1;
      this->stat="ICS_R1_missing";
      return 0;
    }
    
    //if no valid kmer can be used
    if(kmers.size()==0){
      this->barcode="*";
      this->barcodeStrand="*";
      this->barcodeStart=-1;
      this->barcodeMismatch=-1;
      this->stat="CS_too_close_to_end_or_to_each_other";
      return 0;
    }
    
/*
    if((kmer_index_end+barcode_length-1>this->Seq.length())&&(kmer_index_start+barcode_length-1<this->Seq.length())){
      kmers.insert(kmers.end(),this->Kmers.begin()+kmer_index_start,this->Kmers.end());
    }else if(kmer_index_start+barcode_length-1>this->Seq.length()){
      this->INNER=-1;
      this->OUTER=1;
      this->barcode="*";
      this->barcodeStrand="*";
      this->barcodeStart=-1;
      this->barcodeMismatch=-1;
      this->stat="R1_wrongly_oriented";
//      this->printBarcodeInfo();
      return 0;
    }else{
      //if constant seqs at left end
      kmers.insert(kmers.end(),this->Kmers.begin()+kmer_index_start+TENX_R1.length()-2,this->Kmers.begin()+kmer_index_end-barcode_length-TENX_TSO.length()+2);
    }
 */
    //cout << "ok" << endl;
    //for kmers, try to find exact match in barcode list, else to find ambiguous match
    for(int j=0;j<kmers.size();j++){
      int index=kmer_index_start+j+1;
      /*
       if(j<=kmers.size()/2){
       index=j+1;
       }else{
       index=Reads[i].Seq.length()-maxDistToEnds+j-(kmers.size()/2)+1;
       }
       */
      if(barcodes.barcode_ori[kmers[j]].length()!=0&&constant_seq_strand==0){
        this->barcode=kmers[j];
        this->barcodeStrand="original";
        this->barcodeStart=index;
        this->barcodeMismatch=0;
        this->stat="fine";
        this->barcode_in_tolerance++;
        this->barcode_intolerance_detail.push_back(kmers[j]);
        this->barcode_intolerance_indices.push_back(index);
        this->barcode_mismatches.push_back(0);
//        this->matched_kmer_quality=this->Quality.substr(index,16);
        break;
      }else if(barcodes.barcode_rc[kmers[j]].length()!=0&&constant_seq_strand==3){
        this->barcode=barcodes.barcode_rc[kmers[j]];
        this->barcodeStrand="reverse_complement";
        this->barcodeStart=index;
        this->barcodeMismatch=0;
        this->stat="fine";
        this->barcode_in_tolerance++;
        this->barcode_intolerance_indices.push_back(index);
        this->barcode_intolerance_detail.push_back(barcodes.barcode_rc[kmers[j]]);
        this->barcode_mismatches.push_back(0);
//        this->matched_kmer_quality=this->Quality.substr(index,16);
        break;
      }
    } //for(int j=0;j<kmers.size();j++)
    //if exact match found, stop function.
    if(this->barcode.length()!=0){
//      substring to generate a seq downstream barcode for 10nt as UMI, 
//      and modify reverse_com seq to original one, cut off adapter to UMI
      this->guessUMI(tech);
      return 0;
      //start to find ambiguous match
    }else{
      unordered_map<string, string> barcode_record;
      //cout << "start to find ambiguous hits.." << endl;
      for(int j=0;j<kmers.size();j++){
        int index=kmer_index_start+j+1;
        vector<pair<int, int>> candidates;
        vector<string> kmer_segment=getMaxComplexitySegments(kmers[j],segments);
        for(int n=0;n<kmer_segment.size();n++){
          string query=kmer_segment[n]+to_string(n)+to_string(constant_seq_strand);
          candidates.insert(candidates.end(),barcodes.barcode_dict[query].begin(),barcodes.barcode_dict[query].end());
        }
        if(candidates.size()==0){
          continue;
        }else{
          for(int m=0;m<candidates.size();m++){
            int barcode_dict_index=candidates[m].first-1;
            string candidate_ori=barcodes.barcodes[barcode_dict_index];
            if(barcode_record.count(candidate_ori)){
              continue;
            }
            barcode_record[candidate_ori]=candidate_ori;
            string candidate=barcodes.barcodes[barcode_dict_index];
            int strand=candidates[m].second;
            if(strand==1){
              candidate=reverse(candidate);
            }else if(strand==2){
              candidate=complement(candidate);
            }else if(strand==3){
              candidate=reverse_complement(candidate);
            }else{
              if(strand!=0){
                cout << "Unsupported strand coding, Read: " << this->ID << ", kmer: " << kmers[j] << endl;
              }
            }
            int mismatch=editDistance(candidate,kmers[j],maxMismatch);
            //take care here!!!! once kmer meets the threshould, it will be considered to be the barcode, 
            //even there may be better ones after it.
            if(mismatch>maxMismatch){
              continue;
            }else{
              this->barcode=candidate_ori;
              this->barcodeStrand=num_to_strand[strand];
              this->barcodeStart=index;
              this->barcodeMismatch=mismatch;
              this->stat="fine";
              this->barcode_in_tolerance++;
              this->barcode_intolerance_detail.push_back(candidate_ori);
              this->barcode_intolerance_indices.push_back(index);
              this->barcode_mismatches.push_back(mismatch);
//              this->matched_kmer_quality=this->Quality.substr(index,16);
              continue;
            }
          } //for(int m=0;m<candidates.size();m++)
        }// if(candidates.size()==0) else
      }//for(int j=0;j<kmers.size();j++)
    }// if(tmpRead.barcode.length()!=0) else
    //check again if barcodes found, if not, fill the field with '*'.
    //if multiple barcode found, get the barcode with min mismatch
    if(barcode_mismatches.size()>1){
      int minPosition=min_element(barcode_mismatches.begin(),barcode_mismatches.end()) - barcode_mismatches.begin();
      this->barcode=barcode_intolerance_detail[minPosition];
      this->barcodeStart=barcode_intolerance_indices[minPosition];
      this->barcodeMismatch=barcode_mismatches[minPosition];
    }
    
    if(this->barcode.length()!=0){
    //substring to generate a seq downstream barcode for 10nt as UMI, 
    //and modify reverse_com seq to original one, cut off adapter to UMI
        this->guessUMI(tech);
    }else{
      this->barcode="*";
      this->barcodeStrand="*";
      this->barcodeStart=-1;
      this->barcodeMismatch=-1;
      this->stat="barcode_missing";
      return 0;
    }
  } //identifyBarcodes

private:
  void guessUMI(int tech){
    if(tech==3){
      if(this->barcode.find("*")!=0){
        if(this->barcodeStrand.find("original")==0){
          int UMI_start=this->barcodeStart+16-1;
          int UMI_end=UMI_start+12-1;
          this->UMI=this->Seq.substr(UMI_start,12);
          this->Seq_trimmed=reverse_complement(this->Seq.substr(UMI_end+1));
          this->Quality_trimmed=reverse(this->Quality.substr(UMI_end+1));
        }else{
          int UMI_start=this->barcodeStart-12-1;
          int UMI_end=UMI_start+12-1;
          this->UMI=reverse_complement(this->Seq.substr(UMI_start,12));
          this->Seq_trimmed=this->Seq.substr(0,UMI_start);
          this->Quality_trimmed=this->Quality.substr(0,UMI_start);
        }
      }else{
        return;
      }
    }else if(tech==5){
      if(this->barcode.find("*")!=0){
        if(this->barcodeStrand.find("original")==0){
          int UMI_start=this->barcodeStart+16-1;
          int UMI_end=UMI_start+10-1;
          this->UMI=this->Seq.substr(UMI_start,10);
          this->Seq_trimmed=this->Seq.substr(UMI_end+1);
          this->Quality_trimmed=this->Quality.substr(UMI_end+1);
        }else{
          int UMI_start=this->barcodeStart-10-1;
          int UMI_end=UMI_start+10-1;
          this->UMI=reverse_complement(this->Seq.substr(UMI_start,10));
          this->Seq_trimmed=reverse_complement(this->Seq.substr(0,UMI_start));
          this->Quality_trimmed=reverse(this->Quality.substr(0,UMI_start));
        }
      }else{
        return;
      }
    }else{
      return;
    }
  } //guessUMI
};

struct AmbiguousBarcode{
  string seq;
  int strand;
  int mismatch;
  };

vector<string> tokenize(string const &str, const char delim){
  size_t start;
  size_t end = 0;
  vector<string> splits;
  while ((start = str.find_first_not_of(delim, end)) != string::npos)
  {
    end = str.find(delim, start);
    splits.push_back(str.substr(start, end - start));
  }
  return splits;
}

class ReadFile{
private:
  //return 0-based index of constant sequence find in specified read part. 
  pair<int, int> isConstantSequenceContaining(Read& read, string& constantseq, int testRange, int max_testseq_mismatch, int strand){
    string testForward;
    string testseq;
    string testReverse;
    string testseq_rc;
    pair<int, int> res;
    vector<string> testsegF;
    vector<string> testsegR;
    int segments=max_testseq_mismatch+1;
    int index1=-1;
    int index2=-1;
    if(strand==0){
      testForward=read.Seq.substr(0,testRange);
      testseq=constantseq;
      testsegF=getMaxComplexitySegments(testseq,getMaxComplexitySegments(testseq.length(),segments));
      int i=0;
      while(i<testsegF.size()&&index1<0&&index2<0){
        index1=testForward.find(testsegF[i]);
        i++;
      } //while
    }else if(strand==3){
      testReverse=read.Seq.substr(read.Seq.length()-testRange,testRange);
      testseq_rc=reverse_complement(constantseq);
      testsegR=getMaxComplexitySegments(testseq_rc,getMaxComplexitySegments(testseq_rc.length(),segments));
      int i=0;
      while(i<testsegR.size()&&index1<0&&index2<0){
        index2=testReverse.find(testsegR[i])+read.Seq.length()-testRange;
        i++;
      } //while
    }else{
      testForward=read.Seq.substr(0,testRange);
      testseq=constantseq;
      testsegF=getMaxComplexitySegments(testseq,getMaxComplexitySegments(testseq.length(),segments));
      testReverse=read.Seq.substr(read.Seq.length()-testRange,testRange);
      testseq_rc=reverse_complement(constantseq);
      testsegR=getMaxComplexitySegments(testseq_rc,getMaxComplexitySegments(testseq_rc.length(),segments));
      int j=0;
      while(j<testsegF.size()&&index1<0){
        index1=testForward.find(testsegF[j]);
        j++;
      } //while
      int i=0;
      while(i<testsegR.size()&&index2<0){
        index2=testReverse.find(testsegR[i]);
        i++;
      } //while
      if(index2>=0){
        index2=index2+read.Seq.length()-testRange;
      }
    }
    res.first=index1;
    res.second=index2;
    return res;
  } //isConstantSequenceContaining
  
  int checkPolyAT(Read& read){
    string testseq=read.Seq;
    int ployA=testseq.find("AAAAAAAAAA");
    int polyT=testseq.find("TTTTTTTTTT");
    if(ployA>=0&&polyT<0){
       return 0;
    }else if(ployA<0&&polyT>=0){
       return 3;
    }else{
       return -1;
      }
    } //
public:
  int ReadCount;
  gzFile File;
  kseq_t *seq;
  vector<struct Read> Reads;
  int MaxTSOMismatch;
  int MinR1Overlap;
  unordered_map<string, AmbiguousBarcode> ambiguous_barcode;
  int kmerSize;
  
  //basic constructor
  //constructor reading whole fastq into vector
  ReadFile(const char* filename){
    ReadCount=0;
    File=gzopen(filename, "r");
    seq = kseq_init(File);
    while (kseq_read(seq) >= 0){
      ReadCount++;
      Read tmpRead;
      tmpRead.ID=seq->name.s;
      tmpRead.Seq=seq->seq.s;
      tmpRead.Quality=seq->qual.s;
      Reads.push_back(tmpRead);
      }
    kseq_destroy(seq);
    gzclose(File);
  };
  
  //constructor removing chimera by identifying multiple polyT/A
  ReadFile(const char* filename, const char* read_summary_filename, int tech){
    fstream mapping_stat;
    fstream read_summary;
    fstream read_split_report;
    gzFile RawReadFile=gzopen("chimeric_raw_split_fastq.gz", "wb2");
    ReadCount=0;
    unordered_map<string,pair<int, string>> read2barcodeStar_strand;
    cout << "loading chimera read summary..." << endl;
    read_summary.open(read_summary_filename,fstream::in);
    //    int line_counter=0;
    while(!read_summary.eof()){
      string record;
      vector<string> fields;
      string read_id;
      pair<int, string> barcodeStar_strand;
      getline(read_summary,record,'\n');
      //      cout << line << endl;
      //      line_counter++;
      //      if(line_counter%100000==0){
      //        cout << "processed " << line_counter << " line" << endl;
      //      }
      if(record==""){
        continue;
      }
      fields=tokenize(record,'\t');
      read_id=fields[0];
      barcodeStar_strand.first=stoi(fields[4]);
      barcodeStar_strand.second=fields[6];
      read2barcodeStar_strand[read_id]=barcodeStar_strand;
    }
    read_summary.close();
    cout << "start read processing" << endl;
    read_split_report.open("read_split_report.tsv",fstream::in|fstream::out|fstream::app);
    File=gzopen(filename, "r");
    seq = kseq_init(File);
    while (kseq_read(seq) >= 0){
      ReadCount++;
      Read tmpRead;
      int split_position;
      int read_length;
      tmpRead.ID=seq->name.s;
      tmpRead.ID=tokenize(tmpRead.ID,'_')[0];
      tmpRead.Seq=seq->seq.s;
      tmpRead.Quality=seq->qual.s;
      tmpRead.barcodeStrand=read2barcodeStar_strand[tmpRead.ID].second;
      tmpRead.barcodeStart=read2barcodeStar_strand[tmpRead.ID].first;
      read_length=tmpRead.Seq.length();
      if(read_length<=200){
        read_split_report << tmpRead.ID << "\t" << tmpRead.ID << "\t" << tmpRead.barcodeStart << "\t" << -1 << "\t" << "read_too_short" << endl;
        fq_gz_write(RawReadFile,tmpRead);
        continue;
      }
      if(tmpRead.barcodeStrand=="original"){
          Read splitRead1;
          Read splitRead2;
          int first_read_length;
          if(tech==5){
            vector<string> polyAs;
            regex regexp("A{10,10}[^A]");
            smatch polyAs_match;
            string testseq=tmpRead.Seq;
            while(regex_search(testseq, polyAs_match, regexp)){
              //              cout << polyTs[0] << endl;
              polyAs.push_back(polyAs_match[0]);
              testseq=polyAs_match.suffix().str();
            }
            if(polyAs.empty()){
              read_split_report << tmpRead.ID << "\t" << tmpRead.ID << "\t" << tmpRead.barcodeStart << "\t" << -1 << "\t" << "fusion_gene_like" << endl;
              fq_gz_write(RawReadFile,tmpRead);
              continue;
            }else{
              int polyAIndex;
              polyAIndex=tmpRead.Seq.find(polyAs[0]);
//              cout << polyAs[0] << endl;
//              cout << polyAIndex << endl;
              if(polyAIndex+10>read_length-200){
                read_split_report << tmpRead.ID << "\t" << tmpRead.ID << "\t" << tmpRead.barcodeStart << "\t" << -1 << "\t" << "fusion_gene_like" << endl;
                fq_gz_write(RawReadFile,tmpRead);
                continue;
              }else{
                first_read_length=polyAIndex+10;
              }
            }//if(polyAs.empty())
          }else if(tech==3){
            vector<string> polyTs;
            regex regexp("[^T]T{10,10}");
            smatch polyTs_match;
            string testseq=tmpRead.Seq;
            while(regex_search(testseq, polyTs_match, regexp)){
              //              cout << polyTs[0] << endl;
              polyTs.push_back(polyTs_match[0]);
              testseq=polyTs_match.suffix().str();
            }
            if(polyTs.empty()){
              read_split_report << tmpRead.ID << "\t" << tmpRead.ID << "\t" << tmpRead.barcodeStart << "\t" << -1 << "\t" << "fusion_gene_like" << endl;
              fq_gz_write(RawReadFile,tmpRead);
              continue;
            }else if(polyTs.size()==1){
              int polyTIndex;
              polyTIndex=tmpRead.Seq.find(polyTs[0]);
              if(polyTIndex+1<=70){
                read_split_report << tmpRead.ID << "\t" << tmpRead.ID << "\t" << tmpRead.barcodeStart << "\t" << -1 << "\t" << "fusion_gene_like" << endl;
                fq_gz_write(RawReadFile,tmpRead);
                continue;
              }else{
                first_read_length=polyTIndex+1-70;
              }
            }else{
              int polyTIndex;
              if(polyTs[0]==polyTs[1]){
                polyTIndex=tmpRead.Seq.find(polyTs[0],tmpRead.Seq.find(polyTs[0])+1);
              }else{
                polyTIndex=tmpRead.Seq.find(polyTs[1]);
              }
              first_read_length=polyTIndex+1-70;
            }
          }else{
            cout << "Invalid tech while splitting reads" << endl;
            exit(0);
          }
          splitRead1.ID=tmpRead.ID + "-1";
          splitRead1.Seq=tmpRead.Seq.substr(0,first_read_length);
          splitRead1.Quality=tmpRead.Quality.substr(0,first_read_length);
          splitRead2.ID=tmpRead.ID + "-2";
          splitRead2.Seq=tmpRead.Seq.substr(first_read_length,read_length-first_read_length);
          splitRead2.Quality=tmpRead.Quality.substr(first_read_length,read_length-first_read_length);
          read_split_report << tmpRead.ID << "\t" << splitRead1.ID+",-2" << "\t" << tmpRead.barcodeStart << "\t" << first_read_length << "\t" << "splited" << endl;
          fq_gz_write(RawReadFile,splitRead1);
          fq_gz_write(RawReadFile,splitRead2);
          Reads.push_back(splitRead1);
          Reads.push_back(splitRead2);
      }else{
        //reverse complement read
          Read splitRead1;
          Read splitRead2;
          int first_read_length;
          if(tech==5){
            vector<string> polyTs;
            regex regexp("[^T]T{10,10}");
            smatch polyTs_match;
            string testseq=tmpRead.Seq;
            while(regex_search(testseq, polyTs_match, regexp)){
//              cout << polyTs[0] << endl;
              polyTs.push_back(polyTs_match[0]);
              testseq=polyTs_match.suffix().str();
            }
            if(polyTs.empty()){
//              cout << "no polyT found" << endl;
              read_split_report << tmpRead.ID << "\t" << tmpRead.ID << "\t" << tmpRead.barcodeStart << "\t" << -1 << "\t" << "fusion_gene_like" << endl;
              fq_gz_write(RawReadFile,tmpRead);
              continue;
            }else{
              int polyTIndex;
              int sameTail=0;
              for(int i=0;i<polyTs.size()-1;i++){
                if(polyTs[i]==polyTs[polyTs.size()-1]){
                  sameTail++;
                }
              }
//              cout << sameTail << endl;
              if(sameTail==0){
//                cout << "unique tail" << endl;
                polyTIndex=tmpRead.Seq.find(polyTs[polyTs.size()-1]);
              }else{
 //               cout << "iterating tails" << endl;
                polyTIndex=-1;
                for(int j=0;j<=sameTail;j++){
                  polyTIndex=tmpRead.Seq.find(polyTs[polyTs.size()-1],polyTIndex+1);
                }
              }
              if(polyTIndex+1<200){
//                cout << "polyT too close to left end" << endl;
                read_split_report << tmpRead.ID << "\t" << tmpRead.ID << "\t" << tmpRead.barcodeStart << "\t" << -1 << "\t" << "fusion_gene_like" << endl;
                fq_gz_write(RawReadFile,tmpRead);
                continue;
              }else{
                first_read_length=read_length-(polyTIndex+1);
              }
            }
          }else if(tech==3){
            vector<string> polyAs;
            regex regexp("A{10,10}[^A]");
            smatch polyAs_match;
            string testseq=tmpRead.Seq;
            while(regex_search(testseq, polyAs_match, regexp)){
              //              cout << polyTs[0] << endl;
              polyAs.push_back(polyAs_match[0]);
              testseq=polyAs_match.suffix().str();
            }
            if(polyAs.empty()){
              read_split_report << tmpRead.ID << "\t" << tmpRead.ID << "\t" << tmpRead.barcodeStart << "\t" << -1 << "\t" << "fusion_gene_like" << endl;
              fq_gz_write(RawReadFile,tmpRead);
              continue;
            }else if(polyAs.size()==1){
              int polyAIndex;
              polyAIndex=tmpRead.Seq.find(polyAs[0]);
              if(polyAIndex+10>read_length-70){
                read_split_report << tmpRead.ID << "\t" << tmpRead.ID << "\t" << tmpRead.barcodeStart << "\t" << -1 << "\t" << "fusion_gene_like" << endl;
                fq_gz_write(RawReadFile,tmpRead);
                continue;
              }else{
                first_read_length=read_length-(polyAIndex+10+70);
              }
            }else{
              int polyAIndex;
              polyAIndex=tmpRead.Seq.find(polyAs[polyAs.size()-2]);
              first_read_length=read_length-(polyAIndex+10+70);
            }
          }else{
            cout << "Invalid tech while splitting reads" << endl;
            exit(0);
          }
          splitRead1.ID=tmpRead.ID+"-1";
          splitRead1.Seq=tmpRead.Seq.substr(read_length-first_read_length,first_read_length);
          splitRead1.Quality=tmpRead.Quality.substr(read_length-first_read_length,first_read_length);
          splitRead2.ID=tmpRead.ID+"-2";
          splitRead2.Seq=tmpRead.Seq.substr(0,read_length-first_read_length);
          splitRead2.Quality=tmpRead.Quality.substr(0,read_length-first_read_length);
          read_split_report << tmpRead.ID << "\t" << splitRead1.ID+",-2" << "\t" << tmpRead.barcodeStart << "\t" << read_length-first_read_length << "\t" << "splited" << endl;
          fq_gz_write(RawReadFile,splitRead1);
          fq_gz_write(RawReadFile,splitRead2);
          Reads.push_back(splitRead1);
          Reads.push_back(splitRead2);
        }
    }
    kseq_destroy(seq);
    gzclose(File);
    gzclose(RawReadFile);
    read_split_report.close();
  };
  
  
  //constructor aiming at processing chimera reads with mapping stat
  ReadFile(const char* filename, const char* mapping_stat_filename, const char* read_summary_filename){
    fstream mapping_stat;
    fstream read_summary;
    fstream read_split_report;
    ReadCount=0;
    unordered_map<string,vector<int>> read2clippos;
    unordered_map<string,pair<int, string>> read2barcodeStar_strand;
    cout << "loading chimera mapping stat..." << endl;
    mapping_stat.open(mapping_stat_filename,fstream::in);
    while(!mapping_stat.eof()){
      string line;
      string read_id;
      vector<int> clippos;
      vector<string> fields;
      getline(mapping_stat,line,'\n');
      if(line==""){
        continue;
      }
      fields=tokenize(line,'\t');
      read_id=tokenize(fields[0],'_')[0];
      if(stoi(fields[1])==-1){
        clippos.push_back(0);
      }else{
        clippos.push_back(stoi(fields[1]));
      }
      if(stoi(fields[3])==-1){
        clippos.push_back(0);
      }else{
        clippos.push_back(stoi(fields[3]));
      }
      if(stoi(fields[4])==-1){
        clippos.push_back(0);
      }else{
        clippos.push_back(stoi(fields[4]));
      }
      if(stoi(fields[6])==-1){
        clippos.push_back(0);
      }else{
        clippos.push_back(stoi(fields[6]));
      }
      read2clippos[read_id]=clippos;
    }
    mapping_stat.close();
    cout << "loading chimera read summary..." << endl;
    read_summary.open(read_summary_filename,fstream::in);
//    int line_counter=0;
    while(!read_summary.eof()){
      string record;
      vector<string> fields;
      string read_id;
      pair<int, string> barcodeStar_strand;
      getline(read_summary,record,'\n');
//      cout << line << endl;
//      line_counter++;
//      if(line_counter%100000==0){
//        cout << "processed " << line_counter << " line" << endl;
//      }
      if(record==""){
        continue;
      }
      fields=tokenize(record,'\t');
      read_id=fields[0];
      barcodeStar_strand.first=stoi(fields[4]);
      barcodeStar_strand.second=fields[6];
      read2barcodeStar_strand[read_id]=barcodeStar_strand;
    }
    read_summary.close();
    cout << "start read processing" << endl;
    read_split_report.open("read_split_report.tsv",fstream::in|fstream::out|fstream::app);
    File=gzopen(filename, "r");
    seq = kseq_init(File);
    while (kseq_read(seq) >= 0){
      ReadCount++;
      Read tmpRead;
      int mapping_start;
      int split_position;
      int read_length;
      tmpRead.ID=seq->name.s;
      tmpRead.ID=tokenize(tmpRead.ID,'_')[0];
      tmpRead.Seq=seq->seq.s;
      tmpRead.Quality=seq->qual.s;
      tmpRead.barcodeStrand=read2barcodeStar_strand[tmpRead.ID].second;
      tmpRead.barcodeStart=read2barcodeStar_strand[tmpRead.ID].first;
      read_length=tmpRead.Seq.length();
      if(tmpRead.barcodeStrand=="original"){
        mapping_start=min(read2clippos[tmpRead.ID][0],read2clippos[tmpRead.ID][2]);
        if(mapping_start-tmpRead.barcodeStart>=200){
          read_split_report << tmpRead.ID << "\t" << tmpRead.ID << "\t" << tmpRead.barcodeStart << "\t" << mapping_start << "\t" << -1 << "\t" << "complex_fusion" << endl;
          continue;
        }else{
          Read splitRead1;
          Read splitRead2;
          int second_read_length;
          if(mapping_start==read2clippos[tmpRead.ID][0]){
            second_read_length=read2clippos[tmpRead.ID][1];
          }else{
            second_read_length=read2clippos[tmpRead.ID][3];
          }
          split_position=tmpRead.Seq.length()-second_read_length;
          splitRead1.ID=tmpRead.ID+"-1";
          splitRead1.Seq=tmpRead.Seq.substr(0,split_position);
          splitRead1.Quality=tmpRead.Quality.substr(0,split_position);
          splitRead2.ID=tmpRead.ID+"-2";
          splitRead2.Seq=tmpRead.Seq.substr(split_position,second_read_length);
          splitRead2.Quality=tmpRead.Quality.substr(split_position,second_read_length);
          read_split_report << tmpRead.ID << "\t" << splitRead1.ID+",-2" << "\t" << tmpRead.barcodeStart << "\t" << mapping_start << "\t" << split_position << "\t" << "splited" << endl;
          Reads.push_back(splitRead1);
          Reads.push_back(splitRead2);
        }
      }else{
        mapping_start=tmpRead.Seq.length()-min(read2clippos[tmpRead.ID][1],read2clippos[tmpRead.ID][3])-1;
        if(tmpRead.barcodeStart-mapping_start>=200){
          read_split_report << tmpRead.ID << "\t" << tmpRead.ID << "\t" << tmpRead.barcodeStart << "\t" << mapping_start << "\t" << -1 << "\t" << "complex_fusion" << endl;
          continue;
        }else{
          Read splitRead1;
          Read splitRead2;
          int second_read_length;
          if(mapping_start==(tmpRead.Seq.length()-read2clippos[tmpRead.ID][1])-1){
            second_read_length=read2clippos[tmpRead.ID][0];
          }else{
            second_read_length=read2clippos[tmpRead.ID][2];
          }
          split_position=second_read_length;
          splitRead1.ID=tmpRead.ID+"-1";
          splitRead1.Seq=tmpRead.Seq.substr(0,second_read_length);
          splitRead1.Quality=tmpRead.Quality.substr(0,second_read_length);
          splitRead2.ID=tmpRead.ID+"-2";
          splitRead2.Seq=tmpRead.Seq.substr(second_read_length,tmpRead.Seq.length()-second_read_length);
          splitRead2.Quality=tmpRead.Quality.substr(second_read_length,tmpRead.Seq.length()-second_read_length);
          read_split_report << tmpRead.ID << "\t" << splitRead1.ID+",-2" << "\t" << tmpRead.barcodeStart << "\t" << mapping_start << "\t" << split_position << "\t" << "splited" << endl;
          Reads.push_back(splitRead1);
          Reads.push_back(splitRead2);
        }
      }
    }
    kseq_destroy(seq);
    gzclose(File);
    read_split_report.close();
  };
 
  //function used in ReadFile(const char* filename, const char* read_summary_filename, int tech)
  void fq_gz_write(gzFile out_file, Read& read) {
    std::stringstream stream;
    string name=read.ID+"_"+read.barcode;
    stream << "@" << name << "\n" << 
      read.Seq << "\n" << 
        "+" << "\n" << 
      read.Quality << "\n";
    gzputs(out_file, stream.str().c_str());
  }

  
  //active constructor
  //constructor generate output while reading reads
  ReadFile(const char* filename, BarcodeFile& barcodes, vector<int> segments, int maxDistToEnds, int maxMismatch, int tech, int threadnum){
    string output_filename="filtered.fastq.gz";
    string trimmed_output_filename="trimmed_filtered.fastq.gz";
    gzFile outFile=gzopen(output_filename.c_str(), "wb2");
    gzFile trimmed_outFile=gzopen(trimmed_output_filename.c_str(), "wb2");
    ReadCount=0;
    vector<Read> Reads;
    vector<vector<Read>> Readmatrix;
    std::stringstream streamFQ;
    std::stringstream streamFQ_trim;
    std::stringstream streambarcodeinfo;
    int barcode_length=barcodes.barcodes[0].length();
    File=gzopen(filename, "r");
    seq=kseq_init(File);
    int threadcount=0;
    while(true){
      if(kseq_read(seq)>=0){
        Read tmpRead;
        tmpRead.ID=seq->name.s;
        tmpRead.Seq=seq->seq.s;
        tmpRead.Quality=seq->qual.s;
        Reads.push_back(tmpRead);
        if(Reads.size()==100){
          Readmatrix.push_back(Reads);
          threadcount++;
//          cout << "thread count: " << threadcount << endl;
          Reads.clear();
          if(Readmatrix.size()==threadnum){
//            cout << "excuting in if(threadcount==threadnum)" << endl;
            #pragma omp parallel for
            for(int i=0;i<Readmatrix.size();++i){
//              cout << "thread: " << i << endl;
              BarcodeFile threadbarcode=barcodes;
              vector<int> threadsegments=segments;
              int threadmaxDistToEnds=maxDistToEnds;
              int threadmaxMismatch=maxMismatch;
              int threadtech=tech;
              identifyBarcodesinBlock(Readmatrix[i],threadbarcode,threadsegments,threadmaxDistToEnds,threadmaxMismatch,threadtech);
            } //for(int i=0;i<threadnum;++i)
            
            writeReadData(Readmatrix,streambarcodeinfo,streamFQ,streamFQ_trim);
            threadcount=0;
            Readmatrix.clear();
            cout << streambarcodeinfo.str().c_str() << endl;
            gzputs(outFile,streamFQ.str().c_str());
            gzputs(trimmed_outFile,streamFQ_trim.str().c_str());
            streambarcodeinfo.str("");
            streamFQ.str("");
            streamFQ_trim.str("");
          }
        }else{
          continue;
        } //if(Reads.size()==10000)
      }else{
        Readmatrix.push_back(Reads);
        threadcount++;
        Reads.clear();
//        cout << "excuting in if(threadcount==threadnum)else" << endl;
        #pragma omp parallel for
        for(int i=0;i<Readmatrix.size();++i){
          BarcodeFile threadbarcode=barcodes;
          vector<int> threadsegments=segments;
          int threadmaxDistToEnds=maxDistToEnds;
          int threadmaxMismatch=maxMismatch;
          int threadtech=tech;
          identifyBarcodesinBlock(Readmatrix[i],threadbarcode,threadsegments,threadmaxDistToEnds,threadmaxMismatch,threadtech);
        } //for(int i=0;i<Readmatrix.size();++i)

        writeReadData(Readmatrix,streambarcodeinfo,streamFQ,streamFQ_trim);
        threadcount=0;
        Readmatrix.clear();
        cout << streambarcodeinfo.str().c_str() << endl;
        gzputs(outFile,streamFQ.str().c_str());
        gzputs(trimmed_outFile,streamFQ_trim.str().c_str());
        streambarcodeinfo.str("");
        streamFQ.str("");
        streamFQ_trim.str("");
        break;
      } //if(kseq_read(seq)>=0)else{}
    } //while(true)
  };

/*  
  //active constructor for multithreads
  //constructor generate output while reading reads
  ReadFile(const char* filename, BarcodeFile& barcodes, vector<int> segments, int maxDistToEnds, int maxMismatch, int tech, int threads){
    ReadCount=0;
    string output_filename="filtered.fastq.gz";
    gzFile outFile=gzopen(output_filename.c_str(), "wb2");
    int barcode_length=barcodes.barcodes[0].length();
    File=gzopen(filename, "r");
    seq=kseq_init(File);
    while (kseq_read(seq) >= 0){
      ReadCount++;
      Read* tmpRead=NULL;
      tmpRead=new Read;
      tmpRead->ID=seq->name.s;
      tmpRead->Seq=seq->seq.s;
      tmpRead->Quality=seq->qual.s;
      thread t1(&Read::identifyBarcodes,tmpRead,barcodes,segments,maxDistToEnds,maxMismatch,tech,outFile);
      t1.join();
      delete tmpRead;
    } //while (kseq_read(seq) >= 0)
    kseq_destroy(seq);
    gzclose(outFile);
    gzclose(File);
  };
   */
  
  void printReadIDs(){
    for(int i=0;i<Reads.size();++i){
      cout << Reads[i].ID << endl;
      }
    }

  void genrateKmerFromRead(int k){
    kmerSize=k;
    for(int i=0;i<Reads.size();i++){
      const int len=Reads[i].Seq.length();
      for(int j=0;j+k<=len;j++){
        string kmer=Reads[i].Seq.substr(j,k);
//        cout << kmer << endl;
        Reads[i].Kmers.push_back(kmer);
      }
    }
  }
  
  void identifyBarcodesinBlock(vector<Read>& blockreads,BarcodeFile& barcode,vector<int> segments,int maxDistToEnds,int maxMismatch,int tech){
    for(int i=0;i<blockreads.size();++i){
      blockreads[i].identifyBarcodes(barcode,segments,maxDistToEnds,maxMismatch,tech);
    }
  }
  
  void writeReadData(vector<vector<Read>>& readmatrix,std::stringstream &readinfostream, std::stringstream &filteredfastqstream, std::stringstream &trimfilteredfastqstream){
    for(int i=0;i<readmatrix.size();++i){
      for(int j=0;j<readmatrix[i].size();++j){
        readmatrix[i][j].printBarcodeInfo(readinfostream);
        if(readmatrix[i][j].barcode!="*"){
          readmatrix[i][j].fq_gz_write(filteredfastqstream,false);
          readmatrix[i][j].fq_gz_write(trimfilteredfastqstream,true);
        }
      }
    }
  }
}; //class ReadFile

