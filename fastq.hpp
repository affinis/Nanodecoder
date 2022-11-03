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

#include "kseq.h"
#include "prosessSeq.hpp"

using std::unordered_map;
using std::ifstream;
using std::pair;
using std::to_string;
using std::fstream;

//firstly sample from reads, and find TSO and r1 to determine barcode range

KSEQ_INIT(gzFile, gzread)

unordered_map<int, string> num_to_strand{
  {0,"original"},{1,"reverse"},{2,"complement"},{3,"reverse_complement"},{-1,"*"}
  };

string TENX_TSO="TTTCTTATATGGG";
string TENX_R1="CTACACGACGCTCTTCCGATCT";
string POLYT="TTTTTTTTTTT";

class BarcodeFile{
public:
  ifstream barcodeFile;
  vector<string> barcodes;
  unordered_map<string, vector<pair<int, int>>> barcode_dict;
  unordered_map<string, string> barcode_ori;
  unordered_map<string, string> barcode_rc;
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


struct Read{
  string ID;
  string Seq;
  string Quality;
  int OUTER;
  int INNER;
  vector<string> Kmers;
  int barcodeStart;
  string barcode;
  string barcodeStrand;
  int fusionpoint;
  int barcodeMismatch;
  string stat;
  int polyAT_containing;
  
  void printBarcodeInfo(){
    cout << ID << "\t" << stat << "\t" << OUTER << "\t";
    cout << barcode << "\t" << barcodeStart << "\t" << INNER;
    cout << "\t" << barcodeStrand << "\t" << barcodeMismatch << "\t" << polyAT_containing << endl;
  }
  
  void genrateKmer(int k){
    int kmerSize=k;
      const int len=Seq.length();
      for(int j=0;j+k<=len;j++){
        string kmer=Seq.substr(j,k);
        //        cout << kmer << endl;
        Kmers.push_back(kmer);
      }
    }
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
  
  //constructor aiming at processing chimera reads
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
    int line_counter=0;
    while(!read_summary.eof()){
      string record;
      vector<string> fields;
      string read_id;
      pair<int, string> barcodeStar_strand;
      getline(read_summary,record,'\n');
//      cout << line << endl;
      line_counter++;
      if(line_counter%100000==0){
//        cout << "processed " << line_counter << " line" << endl;
      }
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
      tmpRead.ID=seq->name.s;
      tmpRead.ID=tokenize(tmpRead.ID,'_')[0];
      tmpRead.Seq=seq->seq.s;
      tmpRead.Quality=seq->qual.s;
      tmpRead.barcodeStrand=read2barcodeStar_strand[tmpRead.ID].second;
      tmpRead.barcodeStart=read2barcodeStar_strand[tmpRead.ID].first;
      if(tmpRead.barcodeStrand=="original"){
        mapping_start=min(read2clippos[tmpRead.ID][0],read2clippos[tmpRead.ID][2]);
        if(mapping_start-tmpRead.barcodeStart>=100){
          read_split_report << tmpRead.ID << "\t" << tmpRead.ID << "\t" << tmpRead.barcodeStart << "\t" << mapping_start << "\t" << -1 << "\t" << "complex_fusion" << endl;
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
        if(tmpRead.barcodeStart-mapping_start>84){
          read_split_report << tmpRead.ID << "\t" << tmpRead.ID << "\t" << tmpRead.barcodeStart << "\t" << mapping_start << "\t" << -1 << "\t" << "complex_fusion" << endl;
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
  ReadFile(const char* filename, BarcodeFile& barcodes, vector<int> segments, int maxDistToEnds, int maxMismatch, int tech){
    ReadCount=0;
    string output_filename="filtered.fastq.gz";
    gzFile outFile=gzopen(output_filename.c_str(), "wb2");
    int barcode_length=barcodes.barcodes[0].length();
    File=gzopen(filename, "r");
    seq = kseq_init(File);
    while (kseq_read(seq) >= 0){
      ReadCount++;
      Read tmpRead;
      tmpRead.ID=seq->name.s;
      tmpRead.Seq=seq->seq.s;
      tmpRead.Quality=seq->qual.s;
      int barcode_length=barcodes.barcodes[0].length();
      int kmer_index_start=-1;
      int kmer_index_end=-1;
      int constant_seq_strand;
      pair<int, int> inner_stat;
      pair<int, int> outer_stat;
      string strand_by_polyAT;
      //if read shorter than barcode, discard
      if(tmpRead.Seq.length()<=maxDistToEnds){
        tmpRead.INNER=-1;
        tmpRead.OUTER=-1;
        tmpRead.barcode="*";
        tmpRead.barcodeStrand="*";
        tmpRead.barcodeStart=-1;
        tmpRead.barcodeMismatch=-1;
        tmpRead.stat="Read_too_short";
        tmpRead.printBarcodeInfo();
        continue;
      }
      if(tmpRead.Kmers.size()==0){
        tmpRead.genrateKmer(barcode_length);
      }

      //if (find polyA and 5')||(find polyT and 3'), means the barcode should be at left end
      if((checkPolyAT(tmpRead)==0&&tech==5)||(checkPolyAT(tmpRead)==3&&tech==3)){
        constant_seq_strand=0;
        tmpRead.polyAT_containing=1;
        if(tech==5){
          inner_stat=isConstantSequenceContaining(tmpRead,TENX_TSO,maxDistToEnds,1,0);
          outer_stat=isConstantSequenceContaining(tmpRead,TENX_R1,maxDistToEnds,2,0);
        }else if(tech==3){
          inner_stat=isConstantSequenceContaining(tmpRead,POLYT,maxDistToEnds,0,0);
          outer_stat=isConstantSequenceContaining(tmpRead,TENX_R1,maxDistToEnds,2,0);
        }else{
          cout << "Failed, 5' or 3' not specified or specified a wrong argument";
        }
        if(inner_stat.first<0&&outer_stat.first<0){
          tmpRead.INNER=-1;
          tmpRead.OUTER=-1;
          kmer_index_start=0;
          kmer_index_end=maxDistToEnds;
        }else if(inner_stat.first>=0&&outer_stat.first<0){
          tmpRead.INNER=inner_stat.first+1;
          tmpRead.OUTER=-1;
          kmer_index_start=0;
          kmer_index_end=inner_stat.first;
        }else if(inner_stat.first<0&&outer_stat.first>=0){
          tmpRead.INNER=-1;
          tmpRead.OUTER=outer_stat.first+1;
          kmer_index_start=outer_stat.first;
          kmer_index_end=maxDistToEnds;
        }else{
          tmpRead.INNER=inner_stat.first+1;
          tmpRead.OUTER=outer_stat.first+1;
          kmer_index_start=outer_stat.first;
          kmer_index_end=inner_stat.first;
        }
      //if (find polyT and 5')||(find polyA and 3'), means the barcode should be at right end
      }else if((checkPolyAT(tmpRead)==3&&tech==5)||(checkPolyAT(tmpRead)==0&&tech==3)){
        tmpRead.polyAT_containing=1;
//        kmer_index_start=tmpRead.Seq.length()-maxDistToEnds;
//        kmer_index_end=tmpRead.Seq.length()-barcode_length;
        constant_seq_strand=3;
        if(tech==5){
          inner_stat=isConstantSequenceContaining(tmpRead,TENX_TSO,maxDistToEnds,1,3);
          outer_stat=isConstantSequenceContaining(tmpRead,TENX_R1,maxDistToEnds,2,3);
        }else if(tech==3){
          inner_stat=isConstantSequenceContaining(tmpRead,POLYT,maxDistToEnds,0,3);
          outer_stat=isConstantSequenceContaining(tmpRead,TENX_R1,maxDistToEnds,2,3);
        }else{
          cout << "Failed, 5' or 3' not specified or specified a wrong argument";
        }
        if(inner_stat.second<0&&outer_stat.second<0){
          tmpRead.INNER=-1;
          tmpRead.OUTER=-1;
          kmer_index_start=tmpRead.Seq.length()-maxDistToEnds;
          kmer_index_end=tmpRead.Seq.length()-barcode_length;
        }else if(inner_stat.second>=0&&outer_stat.second<0){
          tmpRead.INNER=inner_stat.second+1;
          tmpRead.OUTER=-1;
          kmer_index_start=inner_stat.second;
          kmer_index_end=tmpRead.Seq.length()-barcode_length;
        }else if(inner_stat.second<0&&outer_stat.second>=0){
          tmpRead.INNER=-1;
          tmpRead.OUTER=outer_stat.second+1;
          kmer_index_start=tmpRead.Seq.length()-maxDistToEnds;
          kmer_index_end=outer_stat.second;
        }else{
          tmpRead.INNER=inner_stat.second+1;
          tmpRead.OUTER=outer_stat.second+1;
          kmer_index_start=inner_stat.second;
          kmer_index_end=outer_stat.second;
        }
      }else{
        tmpRead.polyAT_containing=0;
        if(tech==5){
          inner_stat=isConstantSequenceContaining(tmpRead,TENX_TSO,maxDistToEnds,1,-1);
          outer_stat=isConstantSequenceContaining(tmpRead,TENX_R1,maxDistToEnds,2,-1);
        }else if(tech==3){
          inner_stat=isConstantSequenceContaining(tmpRead,POLYT,maxDistToEnds,0,-1);
          outer_stat=isConstantSequenceContaining(tmpRead,TENX_R1,maxDistToEnds,2,-1);
        }else{
          cout << "Failed, 5' or 3' not specified or specified a wrong argument";
        }
        if(inner_stat.first<0&&outer_stat.first<0&&inner_stat.second<0&&outer_stat.second<0){
          tmpRead.INNER=-1;
          tmpRead.OUTER=-1;
          tmpRead.barcode="*";
          tmpRead.barcodeStrand="*";
          tmpRead.barcodeStart=-1;
          tmpRead.barcodeMismatch=-1;
          tmpRead.stat="ICS_R1_missing";
          tmpRead.printBarcodeInfo();
          continue;
        }else if(inner_stat.first<0&&outer_stat.first>=0&&inner_stat.second<0&&outer_stat.second<0){
          tmpRead.INNER=-1;
          tmpRead.OUTER=outer_stat.first+1;
          kmer_index_start=outer_stat.first;
          kmer_index_end=maxDistToEnds;
          constant_seq_strand=0;
        }else if(inner_stat.first<0&&outer_stat.first<0&&inner_stat.second<0&&outer_stat.second>=0){
          tmpRead.INNER=-1;
          tmpRead.OUTER=outer_stat.second+1;
          kmer_index_start=tmpRead.Seq.length()-maxDistToEnds;
          kmer_index_end=outer_stat.second;
          constant_seq_strand=3;
        }else if(inner_stat.first<0&&outer_stat.first>=0&&inner_stat.second<0&&outer_stat.second>=0){
          tmpRead.INNER=-1;
          tmpRead.OUTER=-2;
          tmpRead.barcode="*";
          tmpRead.barcodeStrand="*";
          tmpRead.barcodeStart=-1;
          tmpRead.barcodeMismatch=-1;
          tmpRead.stat="ICS_missing_R1_both_ends";
          tmpRead.printBarcodeInfo();
          continue;
        }else if(inner_stat.first>=0&&outer_stat.first<0&&inner_stat.second<0&&outer_stat.second<0){
          tmpRead.INNER=inner_stat.first+1;
          tmpRead.OUTER=-1;
          kmer_index_start=0;
          kmer_index_end=inner_stat.first;
          constant_seq_strand=0;
        }else if(inner_stat.first>=0&&outer_stat.first>=0&&inner_stat.second<0&&outer_stat.second<0){
          tmpRead.INNER=inner_stat.first+1;
          tmpRead.OUTER=outer_stat.first+1;
          kmer_index_start=outer_stat.first;
          kmer_index_end=inner_stat.first;
          constant_seq_strand=0;
        }else if(inner_stat.first>=0&&outer_stat.first<0&&inner_stat.second<0&&outer_stat.second>=0){
          tmpRead.INNER=inner_stat.first+1;
          tmpRead.OUTER=outer_stat.second+1;
          tmpRead.barcode="*";
          tmpRead.barcodeStrand="*";
          tmpRead.barcodeStart=-1;
          tmpRead.barcodeMismatch=-1;
          tmpRead.stat="ICS_R1_different_end";
          tmpRead.printBarcodeInfo();
          continue;
        }else if(inner_stat.first>=0&&outer_stat.first>=0&&inner_stat.second<0&&outer_stat.second>=0){
          tmpRead.INNER=inner_stat.first+1;
          tmpRead.OUTER=outer_stat.first+1;
          kmer_index_start=outer_stat.first;
          kmer_index_end=inner_stat.first;
          constant_seq_strand=0;
        }else if(inner_stat.first<0&&outer_stat.first<0&&inner_stat.second>=0&&outer_stat.second<0){
          tmpRead.INNER=inner_stat.second+1;
          tmpRead.OUTER=-1;
          kmer_index_start=inner_stat.second;
          kmer_index_end=tmpRead.Seq.length()-barcode_length;
          constant_seq_strand=3;
        }else if(inner_stat.first<0&&outer_stat.first>=0&&inner_stat.second>=0&&outer_stat.second<0){
          tmpRead.INNER=inner_stat.second+1;
          tmpRead.OUTER=outer_stat.first+1;
          tmpRead.barcode="*";
          tmpRead.barcodeStrand="*";
          tmpRead.barcodeStart=-1;
          tmpRead.barcodeMismatch=-1;
          tmpRead.stat="ICS_R1_different_end";
          tmpRead.printBarcodeInfo();
          continue;
        }else if(inner_stat.first<0&&outer_stat.first<0&&inner_stat.second>=0&&outer_stat.second>=0){
          tmpRead.INNER=inner_stat.second+1;
          tmpRead.OUTER=outer_stat.second+1;
          kmer_index_start=inner_stat.second;
          kmer_index_end=outer_stat.second;
          constant_seq_strand=3;
        }else if(inner_stat.first<0&&outer_stat.first>=0&&inner_stat.second>=0&&outer_stat.second>=0){
          tmpRead.INNER=inner_stat.second+1;
          tmpRead.OUTER=outer_stat.second+1;
          kmer_index_start=inner_stat.second;
          kmer_index_end=outer_stat.second;
          constant_seq_strand=3;
        }else if(inner_stat.first>=0&&outer_stat.first<0&&inner_stat.second>=0&&outer_stat.second<0){
          tmpRead.INNER=-2;
          tmpRead.OUTER=-1;
          tmpRead.barcode="*";
          tmpRead.barcodeStrand="*";
          tmpRead.barcodeStart=-1;
          tmpRead.barcodeMismatch=-1;
          tmpRead.stat="ICS_both_ends_R1_missing";
          tmpRead.printBarcodeInfo();
          continue;
        }else if(inner_stat.first>=0&&outer_stat.first>=0&&inner_stat.second>=0&&outer_stat.second<0){
          tmpRead.INNER=inner_stat.first+1;
          tmpRead.OUTER=outer_stat.first+1;
          kmer_index_start=outer_stat.first;
          kmer_index_end=inner_stat.first;
          constant_seq_strand=0;
        }else if(inner_stat.first>=0&&outer_stat.first<0&&inner_stat.second>=0&&outer_stat.second>=0){
          tmpRead.INNER=inner_stat.second+1;
          tmpRead.OUTER=outer_stat.second+1;
          kmer_index_start=inner_stat.second;
          kmer_index_end=outer_stat.second;
          constant_seq_strand=3;
        }else{
          if(inner_stat.first>outer_stat.first&&inner_stat.second<outer_stat.second){
            tmpRead.INNER=-2;
            tmpRead.OUTER=-2;
            tmpRead.barcode="*";
            tmpRead.barcodeStrand="*";
            tmpRead.barcodeStart=-1;
            tmpRead.barcodeMismatch=-1;
            tmpRead.stat="ICS_R1_both_ends";
            tmpRead.printBarcodeInfo();
            continue;
          }else if(inner_stat.first<outer_stat.first&&inner_stat.second<outer_stat.second){
            tmpRead.INNER=inner_stat.second+1;
            tmpRead.OUTER=outer_stat.second+1;
            kmer_index_start=inner_stat.second;
            kmer_index_end=outer_stat.second;
            constant_seq_strand=3;
          }else if(inner_stat.first>outer_stat.first&&inner_stat.second>outer_stat.second){
            tmpRead.INNER=inner_stat.first+1;
            tmpRead.OUTER=outer_stat.first+1;
            kmer_index_start=outer_stat.first;
            kmer_index_end=inner_stat.first;
            constant_seq_strand=0;
          }else{
            tmpRead.INNER=-1;
            tmpRead.OUTER=-1;
            tmpRead.barcode="*";
            tmpRead.barcodeStrand="*";
            tmpRead.barcodeStart=-1;
            tmpRead.barcodeMismatch=-1;
            tmpRead.stat="ICS_R1_wrongly_oriented";
            tmpRead.printBarcodeInfo();
            continue;
          }
        }
      }
      if(kmer_index_start>kmer_index_end){
        tmpRead.barcode="*";
        tmpRead.barcodeStrand="*";
        tmpRead.barcodeStart=-1;
        tmpRead.barcodeMismatch=-1;
        tmpRead.stat="ICS_R1_wrongly_oriented";
        tmpRead.printBarcodeInfo();
        continue;
      }
      std::vector<string> kmers;
      //determine kmers by inner and outer constant seq start index
      //cout << "determine kmers by inner and outer constant seq start index" << endl;
      //cout << kmer_index_start << " " << kmer_index_end << " " << barcode_length << " " << tmpRead.Seq.length() << endl;
      if((kmer_index_end+barcode_length-1>tmpRead.Seq.length())&&(kmer_index_start+barcode_length-1<tmpRead.Seq.length())){
        kmers.insert(kmers.end(),tmpRead.Kmers.begin()+kmer_index_start,tmpRead.Kmers.end());
      }else if(kmer_index_start+barcode_length-1>tmpRead.Seq.length()){
        tmpRead.INNER=-1;
        tmpRead.OUTER=1;
        tmpRead.barcode="*";
        tmpRead.barcodeStrand="*";
        tmpRead.barcodeStart=-1;
        tmpRead.barcodeMismatch=-1;
        tmpRead.stat="R1_wrongly_oriented";
        tmpRead.printBarcodeInfo();
        continue;
      }else{
        kmers.insert(kmers.end(),tmpRead.Kmers.begin()+kmer_index_start,tmpRead.Kmers.begin()+kmer_index_end);
      }
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
          tmpRead.barcode=kmers[j];
          tmpRead.barcodeStrand="original";
          tmpRead.barcodeStart=index;
          tmpRead.barcodeMismatch=0;
          tmpRead.stat="fine";
          break;
        }else if(barcodes.barcode_rc[kmers[j]].length()!=0&&constant_seq_strand==3){
          tmpRead.barcode=barcodes.barcode_rc[kmers[j]];
          tmpRead.barcodeStrand="reverse_complement";
          tmpRead.barcodeStart=index;
          tmpRead.barcodeMismatch=0;
          tmpRead.stat="fine";
          break;
        }
      }
      //if exact match found, next read.
      if(tmpRead.barcode.length()!=0){
        //        cout << "exact hits found" << endl;
        tmpRead.printBarcodeInfo();
        fq_gz_write(outFile,tmpRead);
        continue;
        //start to find ambiguous match
      }else{
        //cout << "start to find ambiguous hits.." << endl;
        for(int j=0;j<kmers.size();j++){
          int index=kmer_index_start+j+1;
          //calculate position of this kmer
          /*
           if(j<=kmers.size()/2){
           index=j+1;
           }else{
           index=Reads[i].Seq.length()-maxDistToEnds+j-(kmers.size()/2)+1;
           }
           */
          if(ambiguous_barcode[kmers[j]].seq.length()!=0&&ambiguous_barcode[kmers[j]].strand==constant_seq_strand){
            tmpRead.barcode=ambiguous_barcode[kmers[j]].seq;
            tmpRead.barcodeStrand=num_to_strand[constant_seq_strand];
            tmpRead.barcodeStart=index;
            tmpRead.barcodeMismatch=ambiguous_barcode[kmers[j]].mismatch;
            tmpRead.stat="fine";
            break;
          }
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
                  cout << "Unsupported strand coding, Read: " << tmpRead.ID << ", kmer: " << kmers[j] << endl;
                }
              }
              int mismatch=editDistance(candidate,kmers[j]);
              //take care here!!!! once kmer meets the threshould, it will be considered to be the barcode, 
              //even there may be better ones after it.
              if(mismatch>maxMismatch){
                continue;
              }else{
                tmpRead.barcode=candidate_ori;
                tmpRead.barcodeStrand=num_to_strand[strand];
                tmpRead.barcodeStart=index;
                tmpRead.barcodeMismatch=mismatch;
                tmpRead.stat="fine";
                AmbiguousBarcode ambiguous_match;
                ambiguous_match.seq=candidate_ori;
                ambiguous_match.strand=strand;
                ambiguous_match.mismatch=mismatch;
                ambiguous_barcode[kmers[j]]=ambiguous_match;
                break;
              }
            } //for(int m=0;m<candidates.size();m++)
          }// if(candidates.size()==0) else
        }//for(int j=0;j<kmers.size();j++)
      }// if(tmpRead.barcode.length()!=0) else
      //check again if barcodes found, if not, fill the field with '*'.
      if(tmpRead.barcode.length()!=0){
        tmpRead.printBarcodeInfo();
        fq_gz_write(outFile,tmpRead);
        continue;
      }else{
        tmpRead.barcode="*";
        tmpRead.barcodeStrand="*";
        tmpRead.barcodeStart=-1;
        tmpRead.barcodeMismatch=-1;
        tmpRead.stat="barcode_missing";
        tmpRead.printBarcodeInfo();
        continue;
      }
    } //while (kseq_read(seq) >= 0)
    kseq_destroy(seq);
    gzclose(outFile);
    gzclose(File);
  };
  
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
  
  void identifyBarcodes(BarcodeFile& barcodes, vector<int> segments, int maxDistToEnds, int maxMismatch, int tech, const char* outputfilename){
    int numReads=Reads.size();
    int barcode_length=barcodes.barcodes[0].length();
    gzFile outFile=gzopen(outputfilename, "wb2");
    cout << "Identifying barcodes from " << numReads << " reads..." << endl;
    for(int i=0;i<numReads;i++){
      Read tmpRead=Reads[i];
      int barcode_length=barcodes.barcodes[0].length();
      int kmer_index_start=-1;
      int kmer_index_end=-1;
      int constant_seq_strand;
      pair<int, int> inner_stat;
      pair<int, int> outer_stat;
      string strand_by_polyAT;
      //if read shorter than barcode, discard
      if(tmpRead.Seq.length()<=maxDistToEnds){
        tmpRead.INNER=-1;
        tmpRead.OUTER=-1;
        tmpRead.barcode="*";
        tmpRead.barcodeStrand="*";
        tmpRead.barcodeStart=-1;
        tmpRead.barcodeMismatch=-1;
        tmpRead.stat="Read_too_short";
        tmpRead.printBarcodeInfo();
        continue;
      }
      if(tmpRead.Kmers.size()==0){
        tmpRead.genrateKmer(barcode_length);
      }
      
      //if (find polyA and 5')||(find polyT and 3'), means the barcode should be at left end
      if((checkPolyAT(tmpRead)==0&&tech==5)||(checkPolyAT(tmpRead)==3&&tech==3)){
        constant_seq_strand=0;
        tmpRead.polyAT_containing=1;
        if(tech==5){
          inner_stat=isConstantSequenceContaining(tmpRead,TENX_TSO,maxDistToEnds,1,0);
          outer_stat=isConstantSequenceContaining(tmpRead,TENX_R1,maxDistToEnds,2,0);
        }else if(tech==3){
          inner_stat=isConstantSequenceContaining(tmpRead,POLYT,maxDistToEnds,0,0);
          outer_stat=isConstantSequenceContaining(tmpRead,TENX_R1,maxDistToEnds,2,0);
        }else{
          cout << "Failed, 5' or 3' not specified or specified a wrong argument";
        }
        if(inner_stat.first<0&&outer_stat.first<0){
          tmpRead.INNER=-1;
          tmpRead.OUTER=-1;
          kmer_index_start=0;
          kmer_index_end=maxDistToEnds;
        }else if(inner_stat.first>=0&&outer_stat.first<0){
          tmpRead.INNER=inner_stat.first+1;
          tmpRead.OUTER=-1;
          kmer_index_start=0;
          kmer_index_end=inner_stat.first;
        }else if(inner_stat.first<0&&outer_stat.first>=0){
          tmpRead.INNER=-1;
          tmpRead.OUTER=outer_stat.first+1;
          kmer_index_start=outer_stat.first;
          kmer_index_end=maxDistToEnds;
        }else{
          tmpRead.INNER=inner_stat.first+1;
          tmpRead.OUTER=outer_stat.first+1;
          kmer_index_start=outer_stat.first;
          kmer_index_end=inner_stat.first;
        }
        //if (find polyT and 5')||(find polyA and 3'), means the barcode should be at right end
      }else if((checkPolyAT(tmpRead)==3&&tech==5)||(checkPolyAT(tmpRead)==0&&tech==3)){
        tmpRead.polyAT_containing=1;
        //        kmer_index_start=tmpRead.Seq.length()-maxDistToEnds;
        //        kmer_index_end=tmpRead.Seq.length()-barcode_length;
        constant_seq_strand=3;
        if(tech==5){
          inner_stat=isConstantSequenceContaining(tmpRead,TENX_TSO,maxDistToEnds,1,3);
          outer_stat=isConstantSequenceContaining(tmpRead,TENX_R1,maxDistToEnds,2,3);
        }else if(tech==3){
          inner_stat=isConstantSequenceContaining(tmpRead,POLYT,maxDistToEnds,0,3);
          outer_stat=isConstantSequenceContaining(tmpRead,TENX_R1,maxDistToEnds,2,3);
        }else{
          cout << "Failed, 5' or 3' not specified or specified a wrong argument";
        }
        if(inner_stat.second<0&&outer_stat.second<0){
          tmpRead.INNER=-1;
          tmpRead.OUTER=-1;
          kmer_index_start=tmpRead.Seq.length()-maxDistToEnds;
          kmer_index_end=tmpRead.Seq.length()-barcode_length;
        }else if(inner_stat.second>=0&&outer_stat.second<0){
          tmpRead.INNER=inner_stat.second+1;
          tmpRead.OUTER=-1;
          kmer_index_start=inner_stat.second;
          kmer_index_end=tmpRead.Seq.length()-barcode_length;
        }else if(inner_stat.second<0&&outer_stat.second>=0){
          tmpRead.INNER=-1;
          tmpRead.OUTER=outer_stat.second+1;
          kmer_index_start=tmpRead.Seq.length()-maxDistToEnds;
          kmer_index_end=outer_stat.second;
        }else{
          tmpRead.INNER=inner_stat.second+1;
          tmpRead.OUTER=outer_stat.second+1;
          kmer_index_start=inner_stat.second;
          kmer_index_end=outer_stat.second;
        }
      }else{
        tmpRead.polyAT_containing=0;
        if(tech==5){
          inner_stat=isConstantSequenceContaining(tmpRead,TENX_TSO,maxDistToEnds,1,-1);
          outer_stat=isConstantSequenceContaining(tmpRead,TENX_R1,maxDistToEnds,2,-1);
        }else if(tech==3){
          inner_stat=isConstantSequenceContaining(tmpRead,POLYT,maxDistToEnds,0,-1);
          outer_stat=isConstantSequenceContaining(tmpRead,TENX_R1,maxDistToEnds,2,-1);
        }else{
          cout << "Failed, 5' or 3' not specified or specified a wrong argument";
        }
        if(inner_stat.first<0&&outer_stat.first<0&&inner_stat.second<0&&outer_stat.second<0){
          tmpRead.INNER=-1;
          tmpRead.OUTER=-1;
          tmpRead.barcode="*";
          tmpRead.barcodeStrand="*";
          tmpRead.barcodeStart=-1;
          tmpRead.barcodeMismatch=-1;
          tmpRead.stat="ICS_R1_missing";
          tmpRead.printBarcodeInfo();
          continue;
        }else if(inner_stat.first<0&&outer_stat.first>=0&&inner_stat.second<0&&outer_stat.second<0){
          tmpRead.INNER=-1;
          tmpRead.OUTER=outer_stat.first+1;
          kmer_index_start=outer_stat.first;
          kmer_index_end=maxDistToEnds;
          constant_seq_strand=0;
        }else if(inner_stat.first<0&&outer_stat.first<0&&inner_stat.second<0&&outer_stat.second>=0){
          tmpRead.INNER=-1;
          tmpRead.OUTER=outer_stat.second+1;
          kmer_index_start=tmpRead.Seq.length()-maxDistToEnds;
          kmer_index_end=outer_stat.second;
          constant_seq_strand=3;
        }else if(inner_stat.first<0&&outer_stat.first>=0&&inner_stat.second<0&&outer_stat.second>=0){
          tmpRead.INNER=-1;
          tmpRead.OUTER=-2;
          tmpRead.barcode="*";
          tmpRead.barcodeStrand="*";
          tmpRead.barcodeStart=-1;
          tmpRead.barcodeMismatch=-1;
          tmpRead.stat="ICS_missing_R1_both_ends";
          tmpRead.printBarcodeInfo();
          continue;
        }else if(inner_stat.first>=0&&outer_stat.first<0&&inner_stat.second<0&&outer_stat.second<0){
          tmpRead.INNER=inner_stat.first+1;
          tmpRead.OUTER=-1;
          kmer_index_start=0;
          kmer_index_end=inner_stat.first;
          constant_seq_strand=0;
        }else if(inner_stat.first>=0&&outer_stat.first>=0&&inner_stat.second<0&&outer_stat.second<0){
          tmpRead.INNER=inner_stat.first+1;
          tmpRead.OUTER=outer_stat.first+1;
          kmer_index_start=outer_stat.first;
          kmer_index_end=inner_stat.first;
          constant_seq_strand=0;
        }else if(inner_stat.first>=0&&outer_stat.first<0&&inner_stat.second<0&&outer_stat.second>=0){
          tmpRead.INNER=inner_stat.first+1;
          tmpRead.OUTER=outer_stat.second+1;
          tmpRead.barcode="*";
          tmpRead.barcodeStrand="*";
          tmpRead.barcodeStart=-1;
          tmpRead.barcodeMismatch=-1;
          tmpRead.stat="ICS_R1_different_end";
          tmpRead.printBarcodeInfo();
          continue;
        }else if(inner_stat.first>=0&&outer_stat.first>=0&&inner_stat.second<0&&outer_stat.second>=0){
          tmpRead.INNER=inner_stat.first+1;
          tmpRead.OUTER=outer_stat.first+1;
          kmer_index_start=outer_stat.first;
          kmer_index_end=inner_stat.first;
          constant_seq_strand=0;
        }else if(inner_stat.first<0&&outer_stat.first<0&&inner_stat.second>=0&&outer_stat.second<0){
          tmpRead.INNER=inner_stat.second+1;
          tmpRead.OUTER=-1;
          kmer_index_start=inner_stat.second;
          kmer_index_end=tmpRead.Seq.length()-barcode_length;
          constant_seq_strand=3;
        }else if(inner_stat.first<0&&outer_stat.first>=0&&inner_stat.second>=0&&outer_stat.second<0){
          tmpRead.INNER=inner_stat.second+1;
          tmpRead.OUTER=outer_stat.first+1;
          tmpRead.barcode="*";
          tmpRead.barcodeStrand="*";
          tmpRead.barcodeStart=-1;
          tmpRead.barcodeMismatch=-1;
          tmpRead.stat="ICS_R1_different_end";
          tmpRead.printBarcodeInfo();
          continue;
        }else if(inner_stat.first<0&&outer_stat.first<0&&inner_stat.second>=0&&outer_stat.second>=0){
          tmpRead.INNER=inner_stat.second+1;
          tmpRead.OUTER=outer_stat.second+1;
          kmer_index_start=inner_stat.second;
          kmer_index_end=outer_stat.second;
          constant_seq_strand=3;
        }else if(inner_stat.first<0&&outer_stat.first>=0&&inner_stat.second>=0&&outer_stat.second>=0){
          tmpRead.INNER=inner_stat.second+1;
          tmpRead.OUTER=outer_stat.second+1;
          kmer_index_start=inner_stat.second;
          kmer_index_end=outer_stat.second;
          constant_seq_strand=3;
        }else if(inner_stat.first>=0&&outer_stat.first<0&&inner_stat.second>=0&&outer_stat.second<0){
          tmpRead.INNER=-2;
          tmpRead.OUTER=-1;
          tmpRead.barcode="*";
          tmpRead.barcodeStrand="*";
          tmpRead.barcodeStart=-1;
          tmpRead.barcodeMismatch=-1;
          tmpRead.stat="ICS_both_ends_R1_missing";
          tmpRead.printBarcodeInfo();
          continue;
        }else if(inner_stat.first>=0&&outer_stat.first>=0&&inner_stat.second>=0&&outer_stat.second<0){
          tmpRead.INNER=inner_stat.first+1;
          tmpRead.OUTER=outer_stat.first+1;
          kmer_index_start=outer_stat.first;
          kmer_index_end=inner_stat.first;
          constant_seq_strand=0;
        }else if(inner_stat.first>=0&&outer_stat.first<0&&inner_stat.second>=0&&outer_stat.second>=0){
          tmpRead.INNER=inner_stat.second+1;
          tmpRead.OUTER=outer_stat.second+1;
          kmer_index_start=inner_stat.second;
          kmer_index_end=outer_stat.second;
          constant_seq_strand=3;
        }else{
          if(inner_stat.first>outer_stat.first&&inner_stat.second<outer_stat.second){
            tmpRead.INNER=-2;
            tmpRead.OUTER=-2;
            tmpRead.barcode="*";
            tmpRead.barcodeStrand="*";
            tmpRead.barcodeStart=-1;
            tmpRead.barcodeMismatch=-1;
            tmpRead.stat="ICS_R1_both_ends";
            tmpRead.printBarcodeInfo();
            continue;
          }else if(inner_stat.first<outer_stat.first&&inner_stat.second<outer_stat.second){
            tmpRead.INNER=inner_stat.second+1;
            tmpRead.OUTER=outer_stat.second+1;
            kmer_index_start=inner_stat.second;
            kmer_index_end=outer_stat.second;
            constant_seq_strand=3;
          }else if(inner_stat.first>outer_stat.first&&inner_stat.second>outer_stat.second){
            tmpRead.INNER=inner_stat.first+1;
            tmpRead.OUTER=outer_stat.first+1;
            kmer_index_start=outer_stat.first;
            kmer_index_end=inner_stat.first;
            constant_seq_strand=0;
          }else{
            tmpRead.INNER=-1;
            tmpRead.OUTER=-1;
            tmpRead.barcode="*";
            tmpRead.barcodeStrand="*";
            tmpRead.barcodeStart=-1;
            tmpRead.barcodeMismatch=-1;
            tmpRead.stat="ICS_R1_wrongly_oriented";
            tmpRead.printBarcodeInfo();
            continue;
          }
        }
      }
      if(kmer_index_start>kmer_index_end){
        tmpRead.barcode="*";
        tmpRead.barcodeStrand="*";
        tmpRead.barcodeStart=-1;
        tmpRead.barcodeMismatch=-1;
        tmpRead.stat="ICS_R1_wrongly_oriented";
        tmpRead.printBarcodeInfo();
        continue;
      }
      std::vector<string> kmers;
      //determine kmers by inner and outer constant seq start index
      //cout << "determine kmers by inner and outer constant seq start index" << endl;
      //cout << kmer_index_start << " " << kmer_index_end << " " << barcode_length << " " << tmpRead.Seq.length() << endl;
      if((kmer_index_end+barcode_length-1>tmpRead.Seq.length())&&(kmer_index_start+barcode_length-1<tmpRead.Seq.length())){
        kmers.insert(kmers.end(),tmpRead.Kmers.begin()+kmer_index_start,tmpRead.Kmers.end());
      }else if(kmer_index_start+barcode_length-1>tmpRead.Seq.length()){
        tmpRead.INNER=-1;
        tmpRead.OUTER=1;
        tmpRead.barcode="*";
        tmpRead.barcodeStrand="*";
        tmpRead.barcodeStart=-1;
        tmpRead.barcodeMismatch=-1;
        tmpRead.stat="R1_wrongly_oriented";
        tmpRead.printBarcodeInfo();
        continue;
      }else{
        kmers.insert(kmers.end(),tmpRead.Kmers.begin()+kmer_index_start,tmpRead.Kmers.begin()+kmer_index_end);
      }
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
          tmpRead.barcode=kmers[j];
          tmpRead.barcodeStrand="original";
          tmpRead.barcodeStart=index;
          tmpRead.barcodeMismatch=0;
          tmpRead.stat="fine";
          break;
        }else if(barcodes.barcode_rc[kmers[j]].length()!=0&&constant_seq_strand==3){
          tmpRead.barcode=barcodes.barcode_rc[kmers[j]];
          tmpRead.barcodeStrand="reverse_complement";
          tmpRead.barcodeStart=index;
          tmpRead.barcodeMismatch=0;
          tmpRead.stat="fine";
          break;
        }
      }
      //if exact match found, next read.
      if(tmpRead.barcode.length()!=0){
        //        cout << "exact hits found" << endl;
        tmpRead.printBarcodeInfo();
        fq_gz_write(outFile,tmpRead);
        continue;
        //start to find ambiguous match
      }else{
        //cout << "start to find ambiguous hits.." << endl;
        for(int j=0;j<kmers.size();j++){
          int index=kmer_index_start+j+1;
          //calculate position of this kmer
          /*
           if(j<=kmers.size()/2){
           index=j+1;
           }else{
           index=Reads[i].Seq.length()-maxDistToEnds+j-(kmers.size()/2)+1;
           }
           */
          if(ambiguous_barcode[kmers[j]].seq.length()!=0&&ambiguous_barcode[kmers[j]].strand==constant_seq_strand){
            tmpRead.barcode=ambiguous_barcode[kmers[j]].seq;
            tmpRead.barcodeStrand=num_to_strand[constant_seq_strand];
            tmpRead.barcodeStart=index;
            tmpRead.barcodeMismatch=ambiguous_barcode[kmers[j]].mismatch;
            tmpRead.stat="fine";
            break;
          }
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
                  cout << "Unsupported strand coding, Read: " << tmpRead.ID << ", kmer: " << kmers[j] << endl;
                }
              }
              int mismatch=editDistance(candidate,kmers[j]);
              //take care here!!!! once kmer meets the threshould, it will be considered to be the barcode, 
              //even there may be better ones after it.
              if(mismatch>maxMismatch){
                continue;
              }else{
                tmpRead.barcode=candidate_ori;
                tmpRead.barcodeStrand=num_to_strand[strand];
                tmpRead.barcodeStart=index;
                tmpRead.barcodeMismatch=mismatch;
                tmpRead.stat="fine";
                AmbiguousBarcode ambiguous_match;
                ambiguous_match.seq=candidate_ori;
                ambiguous_match.strand=strand;
                ambiguous_match.mismatch=mismatch;
                ambiguous_barcode[kmers[j]]=ambiguous_match;
                break;
              }
            } //for(int m=0;m<candidates.size();m++)
          }// if(candidates.size()==0) else
        }//for(int j=0;j<kmers.size();j++)
      }// if(tmpRead.barcode.length()!=0) else
      //check again if barcodes found, if not, fill the field with '*'.
      if(tmpRead.barcode.length()!=0){
        tmpRead.printBarcodeInfo();
        fq_gz_write(outFile,tmpRead);
        continue;
      }else{
        tmpRead.barcode="*";
        tmpRead.barcodeStrand="*";
        tmpRead.barcodeStart=-1;
        tmpRead.barcodeMismatch=-1;
        tmpRead.stat="barcode_missing";
        tmpRead.printBarcodeInfo();
        continue;
      }
    } //for(int i=0;i<numReads;i++){
    gzclose(outFile);
  } //void identifyBarcodes
}; //class ReadFile

