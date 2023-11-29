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

#include "prosessSeq.hpp"
#include "kseq.h"
#include "UCR_DTW2.h"
#include "pore_model.h"
#include "decoder.h"

using std::unordered_map;
using std::ifstream;
using std::pair;
using std::to_string;
using std::fstream;
using std::regex;
using std::smatch;
using std::regex_search;
using std::thread;

KSEQ_INIT(gzFile, gzread)

string TENX_TSO="TTTCTTATATGGG";
string TENX_R1="CTACACGACGCTCTTCCGATCT";
string POLYT="TTTTTTTTTTTTT";

struct ReadAnno{
  string read_id;
  vector<string> gene_id;
  vector<string> gene_name;
  vector<int> read_orientation;
  vector<pair<int,int>> free_regions;
  vector<pair<int,int>> qhits;
  vector<vector<int>> hit_overlaps;
  vector<pair<int,int>> occupied_regions;
  vector<float> qcovs;
  vector<float> rcovs;
  vector<int> is_spliced;
  int flag;
};

class BarcodeFile{
public:
  ifstream barcodeFile;
  vector<string> barcodes;
  string barcoedefilename;
  unordered_map<string, vector<pair<int, int>>> barcode_dict;
  unordered_map<string, string> barcode_ori;
  unordered_map<string, string> barcode_rc;
  
  BarcodeFile(){
    
  }
  
  BarcodeFile(const BarcodeFile& source){
    barcodes=source.barcodes;
    barcoedefilename=source.barcoedefilename;
    barcode_dict=source.barcode_dict;
    barcode_ori=source.barcode_ori;
    barcode_rc=source.barcode_rc;
  }
  
  ~BarcodeFile(){

  }
  
  BarcodeFile(const char* filename, vector<int> word_sizes, int barcode_length){
    barcodeFile.open(filename);
    barcoedefilename=(string)filename;
    fprintf(stderr,"BarcodeFile: filename %s\n",barcoedefilename.c_str());
    fprintf(stderr,"BarcodeFile: Loading barcodes and building dicts...\n");
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
        key=barcode_segments_rc[j]+to_string(j)+to_string(3);
        /*
        if(barcode=="AAACAACGAATAGTTC"){
          cout << key << endl;
        }
        */
        barcode_value.first=barcodes.size();
        barcode_value.second=3;
        barcode_dict[key].push_back(barcode_value);
      }
    }
    barcodeFile.close();
  }
  
};


class Read{
public:
  string ID="default";
  string Seq;
  string Quality;
  string UMI="*";
  string Seq_trimmed;
  string Quality_trimmed;
  vector<pair<int,int>> freeRegions;
  vector<pair<int,int>> qhits;
  vector<int> qhitOrient;
  vector<float> qhit_covs;
  vector<vector<int>> hit_overlaps;
  vector<pair<int,int>> occupiedRegions;
  vector<string> hit_names;
  vector<string> hit_ids;
  vector<float> hit_covs;
  vector<int> spliced;
  vector<int> R1Aps;
  vector<int> R1Ass;
  vector<int> R1Ams;
  vector<int> TSOps;
  vector<int> TSOss;
  vector<int> TSOms;
  int OUTER=-1;
  int INNER=-1;
  vector<string> Kmers;
  int barcodeStart=-1;
  string barcode="*";
  string barcodeStrand;
  int fusionPoint;
  int barcodeMismatch;
  string stat;
  int polyAT_containing=0;
  int barcode_in_tolerance=0;
  vector<int> barcode_intolerance_indices;
  vector<string> barcode_intolerance_detail;
  vector<string> matched_kmer_quality;
  vector<int> barcode_mismatches;
  string final_annotation="*";
  
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
    string R1a_position_str="";
    string R1a_strand_str="";
    vector<string> R1a_strands;
    string TSO_position_str="";
    string TSO_strand_str="";
    vector<string> TSO_strands;
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
    if(R1Ass.empty()){
      if(OUTER==-1){
        R1a_position_str="-1";
        R1a_strand_str="*";
      }else{
        R1a_position_str=to_string(OUTER);
        R1a_strand_str="*";
      }
    }else{
      for(int j:this->R1Ass){
        R1a_strands.push_back(num2rna[j]);
      }
      R1a_position_str=stringCat(this->R1Aps);
      R1a_strand_str=stringCat(R1a_strands);
    }
    
    if(TSOss.empty()){
      if(INNER==-1){
        TSO_position_str="-1";
        TSO_strand_str="*";
      }else{
        TSO_position_str=to_string(INNER);
        TSO_strand_str="*";
      }
    }else{
      for(int j:this->TSOss){
        TSO_strands.push_back(num2rna[j]);
      }
      TSO_position_str=stringCat(this->TSOps);
      TSO_strand_str=stringCat(TSO_strands);
    }
    string hit_regions_str=regions2string(this->qhits);
    string hit_strands_str=strands2string(this->qhitOrient);
    cout << ID << "\t" << stat << "\t" << R1a_position_str << "\t" << R1a_strand_str << "\t";
    cout << barcode << "\t" << barcodeStart << "\t" << barcodeStrand << "\t" << barcodeMismatch << "\t";
    cout << barcode_in_tolerance << "\t" << barcode_intolerance_detail_str << "\t" << barcode_intolerance_indices_str << "\t" << barcode_mismatches_str << "\t";
    cout << hit_regions_str << "\t" << hit_strands_str << "\t" << TSO_position_str << "\t" << TSO_strand_str << "\t" << this->final_annotation << endl;
  }
  
  void printBarcodeInfo(std::stringstream &stream){
    if(this->stat=="Read_unmapped"){
      stream << ID << "\t" << stat << endl;
      return;
    }
    string barcode_intolerance_detail_str="";
    string barcode_intolerance_indices_str="";
    string barcode_mismatches_str="";
    string R1a_position_str="";
    string R1a_strand_str="";
    vector<string> R1a_strands;
    string TSO_position_str="";
    string TSO_strand_str="";
    vector<string> TSO_strands;
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
    if(R1Ass.empty()){
      if(OUTER==-1){
        R1a_position_str="-1";
        R1a_strand_str="*";
      }else{
        R1a_position_str=to_string(OUTER);
        R1a_strand_str="*";
      }
    }else{
      for(int j:this->R1Ass){
        R1a_strands.push_back(num2rna[j]);
      }
      R1a_position_str=stringCat(this->R1Aps);
      R1a_strand_str=stringCat(R1a_strands);
    }
    
    if(TSOss.empty()){
      if(INNER==-1){
        TSO_position_str="-1";
        TSO_strand_str="*";
      }else{
        TSO_position_str=to_string(INNER);
        TSO_strand_str="*";
      }
    }else{
      for(int j:this->TSOss){
        TSO_strands.push_back(num2rna[j]);
      }
      TSO_position_str=stringCat(this->TSOps);
      TSO_strand_str=stringCat(TSO_strands);
    }
    
    string hit_regions_str=regions2string(this->qhits);
    string hit_strands_str=strands2string(this->qhitOrient);
    //     1. read_id    2. status          3. R1a_start              4. R1a_strand
    stream << ID << "\t" << stat << "\t" << R1a_position_str << "\t" << R1a_strand_str << "\t";
    //     5. barcode          6. barcode_start         7. barcode_strand        8. barcode_mismatch
    stream << barcode << "\t" << barcodeStart << "\t" << barcodeStrand << "\t" << barcodeMismatch << "\t";
    //         9. num_of barcodes               10. barcodes                             11. barcode_starts                     12. barcode_mismatches
    stream << barcode_in_tolerance << "\t" << barcode_intolerance_detail_str << "\t" << barcode_intolerance_indices_str << "\t" << barcode_mismatches_str << "\t";
    //          13. hit regions           14. hit orientations         15. TSO position           16. TSO strand          17. final annotation
    stream << hit_regions_str << "\t" << hit_strands_str << "\t" << TSO_position_str << "\t" << TSO_strand_str << "\t" << this->final_annotation << endl;
  }
  
  string choose_annotation(string id, string name){
    if(name=="*"){
      if(id=="*"){
        return("non-transcriptome");
      }else{
        return(id);
      }
    }else{
      return(name);
    }
  }
  
  int find_nearest_region(vector<pair<int,int>> regions, int position, int strand){
    if(regions.empty()){
      fprintf(stderr,"find_nearest_region: regions empty.");
      exit(0);
    }
    int min_dist=this->Seq.length();
    int min_dist_index=regions.size();
    for(int i=0;i<regions.size();i++){
      int this_dist;
      if(strand==0){
        this_dist=regions[i].first-position;
      }else if(strand==3){
        this_dist=position-regions[i].second;
      }else{
        fprintf(stderr,"find_nearest_region: invalid strand.");
        exit(0);
      }
      if(this_dist<min_dist){
        min_dist=this_dist;
        min_dist_index=i;
        continue;
      }
    }
    return(min_dist_index);
  }
  
  vector<int> filter_hit_by_strand(vector<int>& strands, int target_strand){
    vector<int> valid_hit_idx;
    bool all_non_transcript=true;
    for(int i=0;i<strands.size();i++){
      if(strands[i]==target_strand){
        valid_hit_idx.push_back(i);
        all_non_transcript=false;
        continue;
      }
    }
    if(all_non_transcript){
      valid_hit_idx.push_back(-1);
    }
    return(valid_hit_idx);
  }
  
  vector<int> filter_hit_by_strand_and_region(vector<pair<int,int>>& occupied_regions,
              vector<pair<int,int>>& hits, vector<int> hit_strands, int target_position, int target_strand){
    vector<int> res;
    pair<int,int> region_for_R1a=occupied_regions[find_nearest_region(occupied_regions,target_position,target_strand)];
    bool all_non_transcript=true;
    for(int idx=0;idx<hits.size();idx++){
      if(isContaining(hits[idx],region_for_R1a)&&hit_strands[idx]==target_strand){
          res.push_back(idx);
          all_non_transcript=false;
      }
    }
    if(all_non_transcript){
      res.push_back(-1);
      return(res);
    }
    return(res);
  }
  
  bool is_fusion(vector<int>& valid_hit,float threshould=0.4){
    for(int i:valid_hit){
      for(int j:valid_hit){
        if(i==j){
          continue;
        }
        pair<int,int> query=this->qhits[i];
        pair<int,int> subject=this->qhits[j];
        float cov=calculate_cov(query,subject);
        if(cov>0.4){
          return(false);
        }
      }
    }
    return(true);
  }
  
  void annotate_read_by_valid_hit(vector<int> valid_hit_idx){
    // very little case, all of the hit in occupied region have diff strand with R1a
    if(valid_hit_idx.empty()){
      this->final_annotation="Non-transcriptome";
      // ideal case, only one hit have same strand with R1a
    }else if(valid_hit_idx.size()==1){
      if(valid_hit_idx[0]!=-1){
        this->final_annotation=choose_annotation(this->hit_ids[valid_hit_idx[0]],this->hit_names[valid_hit_idx[0]]);
      }else{
        this->final_annotation="Non-transcriptome";
      }
      
      //multiple hits on same strand, fusion genes or ambiguous mapping
    }else{
      if(is_fusion(valid_hit_idx)){
        this->final_annotation="Fusion: ";
        for(int i:valid_hit_idx){
          this->final_annotation=this->final_annotation+choose_annotation(this->hit_ids[i],this->hit_names[i])+"-";
        }
        this->final_annotation.erase(this->final_annotation.length()-1);
      }else{
        int best_align_count=0;
        float best_align_score=0;
        int best_align_index=-1;
        for(int i:valid_hit_idx){
          float align_score=this->qhit_covs[i]*this->hit_covs[i];
          if(align_score>0&&align_score>best_align_score){
            best_align_count=1;
            best_align_score=align_score;
            best_align_index=i;
          }else if(align_score>0&&align_score==best_align_score){
            best_align_count++;
          }else{
            continue;
          }
        }
        if(best_align_count==1){
          //cout << best_align_index << endl;
          //cout << best_align_score << endl;
          this->final_annotation=choose_annotation(this->hit_ids[best_align_index],
                                                   this->hit_names[best_align_index]);
        }else if(best_align_count>1){
          this->final_annotation="Ambiguous_mapping";
        }else{
          this->final_annotation="Non-transcriptome";
        }
      }
    }
  }
  
  /*!
   * @abstract this function will only be called when this read has no more than ONE R1 adapter and mapped to genome, 
   *           in other case, read should be split.
   *           
   *           Annotation will be added according to strand and position of R1 adapter, strand of hits.
   *           
   *           Multiple overlapped hit will be considered ambiguous or transcripts of fusion 
   *           genes according to the proportion of overlap to hits.
   * 
   */
  void add_annotation(int tech){
    if(this->R1Aps.empty()){
      this->final_annotation="Undetermined";
      return;
    }
    int desired_hit_strand;
    if(tech==3){
      desired_hit_strand=reverse_strand[this->R1Ass[0]];
    }else{
      desired_hit_strand=this->R1Ass[0];
    }
    //most common condition, only one hit and one occupied region
    if(this->hit_ids.size()==1&&this->occupiedRegions.size()==1){
      if(this->qhitOrient[0]==desired_hit_strand){
        this->final_annotation=choose_annotation(this->hit_ids[0],this->hit_names[0]);
      }else{
        if(this->qhitOrient[0]==-1){
          this->final_annotation="Non-transcriptome";
        }else{
          this->final_annotation="Non-transcriptome";
        }
      }
    //this condition should not exist in real data as num of hits < occupied region, but I processed it here.
    }else if(this->hit_ids.size()==1&&this->occupiedRegions.size()>1){
      pair<int, int> region_for_R1a=occupiedRegions[find_nearest_region(this->occupiedRegions,this->R1Aps[0]-1,desired_hit_strand)];
      if(isContaining(this->qhits[0],region_for_R1a)){
        if(this->qhitOrient[0]==desired_hit_strand){
          this->final_annotation=choose_annotation(this->hit_ids[0],this->hit_names[0]);
        }else{
          if(this->qhitOrient[0]==-1){
            this->final_annotation="Non-transcriptome";
          }else{
            this->final_annotation="Non-transcriptome";
          }
        }
      }else{
        this->final_annotation="Odds0";
      }
      
    // often happen, multiple hits to one region, maybe pseudo- or homologous-genes (ambiguous mapping) and fusion genes as well.
    }else if(this->hit_ids.size()>1&&this->occupiedRegions.size()==1){
      vector<int> valid_hit=filter_hit_by_strand(this->qhitOrient,desired_hit_strand);
      //cout << "debug1.5.1" << endl;
      this->annotate_read_by_valid_hit(valid_hit);

    //happen some times
    }else if(this->hit_ids.size()>1&&this->occupiedRegions.size()>1){
      //cout << "debug1.5.2" << endl;
      vector<int> valid_hit=filter_hit_by_strand_and_region(this->occupiedRegions,this->qhits,this->qhitOrient,this->R1Aps[0],desired_hit_strand);
      this->annotate_read_by_valid_hit(valid_hit);
      
    //else condition may not exist
    }else{
      this->final_annotation="Odds_x";
    }
  }
  
  bool R1a_homo_orient(){
    for(int orien:this->R1Ass){
      if(orien!=this->R1Ass[0]){
        return false;
      }
    }
    return true;
  }
  
  void add_complex_structure_annotation(){
    if(this->R1Aps.size()==0||this->R1Aps.size()==1){
      return;
    }else if(this->R1Aps.size()==2){
        if(this->R1Ass[0]==0&&this->R1Ass[1]==0){
          if(this->R1Aps[1]-this->R1Aps[0]>50){
            this->final_annotation="Chimera";
            return;
          }else{
            this->final_annotation="R1a_duplication";
            return;
          }
        }else if(this->R1Ass[0]==0&&this->R1Ass[1]==3){
          if(this->R1Aps[1]-this->R1Aps[0]>100){
            this->final_annotation="Dual_header";
            return;
          }else{
            this->final_annotation="Others";
            return;
          }
        }else if(this->R1Ass[0]==3&&this->R1Ass[1]==0){
          if(this->R1Aps[1]-this->R1Aps[0]<100){
            this->final_annotation="Anti_dual_header";
            return;
          }else{
            this->final_annotation="Others";
            return;
          }
        }else{
          if(this->R1Aps[1]-this->R1Aps[0]>50){
            this->final_annotation="Chimera";
            return;
          }else{
            this->final_annotation="R1a_duplication";
            return;
          }
        }
    }else{
      if(this->R1a_homo_orient()){
        this->final_annotation="Chimera";
        return;
      }else{
        this->final_annotation="Others";
        return;
      }
    }
  }
  
  void add_flag(char flag){
    switch(flag){
      case 'U':
        this->stat="Read_umapped";
        this->barcode="*";
        this->barcodeStart=-1;
        this->barcodeStrand="*";
        this->barcodeMismatch=-1;
        break;
      case 'L':
        this->stat="Read_too_long";
        this->barcode="*";
        this->barcodeStart=-1;
        this->barcodeStrand="*";
        this->barcodeMismatch=-1;
        break;
      case 'S':
        this->stat="Read_too_short";
        this->barcode="*";
        this->barcodeStart=-1;
        this->barcodeStrand="*";
        this->barcodeMismatch=-1;
        break;
      case 'M':
        this->stat="Barcode_missing";
        this->barcode="*";
        this->barcodeStart=-1;
        this->barcodeStrand="*";
        this->barcodeMismatch=-1;
        break;
      case 'R':
        this->stat="R1a_missing";
        this->barcode="*";
        this->barcodeStart=-1;
        this->barcodeStrand="*";
        this->barcodeMismatch=-1;
        break;
      case 'C':
        this->stat="Complex_structure";
        this->barcode="*";
        this->barcodeStart=-1;
        this->barcodeStrand="*";
        this->barcodeMismatch=-1;
        break;
    }
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
  
  /*!
   @abstract Find certain TSO sequence on read.
   
   @param mismatch_max: max mismatch allowed.
   @param search_start: start position of search range, 1-based.
   @param search_end: end position of search range, 1-based.
   @param seq: target sequence for searching read1 adapter.
   
   @return null, add TSO data to TSOps (position, 1-based), TSOms (mismatch), TSOss (strand).
   */
  void findTSO(int mismatch_max, int search_start, int search_end, int tech){
    if(tech==3){
      TENX_TSO="CCCATGTACTCTGCGTTGATACCACTGCTT";
    }
    int minDist=TENX_TSO.length();
    int minDistIndex;
    int minDistStrand;
    if(search_end>this->Seq.length()){
      fprintf(stderr,"findTSO: search_end larger than read length\n");
      exit(0);
    }
    string testseq=this->Seq.substr(search_start-1,search_end-search_start+1);
    vector<pair<int, int>> res;
    vector<string> kmers=genrateKmerFromSeq(testseq,TENX_TSO.length());
    string TENX_TSO_REV=reverse_complement(TENX_TSO);
    int i=0;
    for(string kmer:kmers){
      int strand;
      int dist=editDistance(TENX_TSO, kmer, mismatch_max);
      strand=0;
      if(editDistance(TENX_TSO_REV, kmer, mismatch_max)<dist){
        dist=editDistance(TENX_TSO_REV, kmer, mismatch_max);
        strand=3;
      }
      if(dist>minDist || dist>mismatch_max){
        i++;
        continue;
      }else{
        minDist=dist;
        minDistIndex=i+search_start;
        minDistStrand=strand;
        //        cout << minDist << "\t" << minDistIndex << "\t" << minDistStrand << endl;
        if(this->TSOps.size()==0||minDistIndex-TSOps[TSOps.size()-1]>=TENX_TSO.length()){
          this->TSOps.push_back(minDistIndex);
          this->TSOss.push_back(minDistStrand);
          this->TSOms.push_back(minDist);
          i++;
          continue;
        }else if(minDistIndex-TSOps[TSOps.size()-1]<TENX_TSO.length()){
          this->TSOps.pop_back();
          this->TSOss.pop_back();
          this->TSOms.pop_back();
          this->TSOps.push_back(minDistIndex);
          this->TSOss.push_back(minDistStrand);
          this->TSOms.push_back(minDist);
          i++;
          continue;
        }
      }
    }
  }
  
  /*!
   @abstract Find Read1 adapter on sequence.
   
   @param mismatch_max: max mismatch allowed for Read1 adapter.
   @param search_start: start position of search range, 1-based.
   @param search_end: end position of search range, 1-based.
   @param seq: target sequence for searching read1 adapter.
   
   @return null, add read1 adapter data to R1Aps (position, 1-based), R1Ams (mismatch), R1Ass (strand).
   */
  void findR1a(int mismatch_max, int search_start, int search_end){
    int minDist=22;
    int minDistIndex;
    int minDistStrand;
    if(search_end>this->Seq.length()){
      fprintf(stderr,"findR1a: search_end larger than read length\n");
      exit(0);
    }
    string testseq=this->Seq.substr(search_start-1,search_end-search_start+1);
    vector<pair<int, int>> res;
    vector<string> kmers=genrateKmerFromSeq(testseq,TENX_R1.length());
    string TENX_R1_REV=reverse_complement(TENX_R1);
    int i=0;
    for(string kmer:kmers){
      int strand;
      int dist=editDistance(TENX_R1, kmer, mismatch_max);
      strand=0;
      if(editDistance(TENX_R1_REV, kmer, mismatch_max)<dist){
        dist=editDistance(TENX_R1_REV, kmer, mismatch_max);
        strand=3;
      }
      if(dist>minDist || dist>mismatch_max){
        i++;
        continue;
      }else{
        minDist=dist;
        minDistIndex=i+search_start;
        minDistStrand=strand;
//        cout << minDist << "\t" << minDistIndex << "\t" << minDistStrand << endl;
        if(this->R1Aps.size()==0||minDistIndex-R1Aps[R1Aps.size()-1]>=22){
          this->R1Aps.push_back(minDistIndex);
          this->R1Ass.push_back(minDistStrand);
          this->R1Ams.push_back(minDist);
          i++;
          continue;
        }else if(minDistIndex-R1Aps[R1Aps.size()-1]<22){
          this->R1Aps.pop_back();
          this->R1Ass.pop_back();
          this->R1Ams.pop_back();
          this->R1Aps.push_back(minDistIndex);
          this->R1Ass.push_back(minDistStrand);
          this->R1Ams.push_back(minDist);
          i++;
          continue;
        }
      }
    }
  }
  
  /*!
   * @abstract find R1a related hit region with R1a_position (1-based) and R1a_orientation.
   *           return start position and length of entire R1a-BC-UMI-transcript region. 
   *           
   * @param R1a_position: R1a start position, 1-based.
   * @param R1a_orientation
   * 
   * @return .first: 0-based start position; .second: length of entire region.
   */
  void find_complete_region_by_R1a(Read& split_container, int R1a_position, int R1a_orientation, int R1a_index){
    pair<int,int> R1a_related_region;
    pair<int,int> R1a_related_occ_region;
    pair<int,int> entire_region={-1,-1};
    bool R1a_related_occ_region_exist=false;
    int free_region_index=0;
    int occ_region_index=0;
    if(R1a_orientation==0&&R1a_index!=this->R1Aps.size()-1){
      if(this->R1Ass[R1a_index+1]==0&&this->R1Aps[R1a_index+1]-R1a_position<=50){
        return;
      }
    }
    if(R1a_orientation==3&&R1a_index!=0){
      if(this->R1Ass[R1a_index-1]==3&&R1a_position-this->R1Aps[R1a_index-1]<=50){
        return;
      }
    }
    for(pair<int,int> free_region:this->freeRegions){
      if(isContaining(R1a_position,free_region)){
        R1a_related_region=free_region;
        break;
      }
      free_region_index++;
    }
    for(pair<int,int> occ_region:this->occupiedRegions){
      if(R1a_orientation==0&&occ_region.first==R1a_related_region.second+1&&occ_region.first-R1a_position<100){
        R1a_related_occ_region=occ_region;
        R1a_related_occ_region_exist=true;
        break;
      }
      if(R1a_orientation==3&&occ_region.second==R1a_related_region.first-1&&R1a_position-occ_region.second<100){
        R1a_related_occ_region=occ_region;
        R1a_related_occ_region_exist=true;
        break;
      }
      occ_region_index++;
    }
    if(!R1a_related_occ_region_exist){
      return;
    }else{
      if(R1a_orientation==0){
        entire_region.first=R1a_related_region.first-1;
        entire_region.second=pair_length(R1a_related_region)+pair_length(R1a_related_occ_region);
      }else if(R1a_orientation==3){
        entire_region.first=R1a_related_occ_region.first-1;
        entire_region.second=pair_length(R1a_related_region)+pair_length(R1a_related_occ_region);
      }else{
        fprintf(stderr,"ERROR: find_complete_region_by_R1a: undefind orientation provided.");
        exit(0);
      }
    }
    
    pair<int,int> entire_region_1_based_start_end={entire_region.first+1,entire_region.first+entire_region.second};
    //fprintf(stderr,"%s: %i-%i\n",split_container.ID.c_str(),entire_region_1_based_start_end.first,entire_region_1_based_start_end.second);
    
    split_container.Seq=this->Seq.substr(entire_region.first,entire_region.second);
    split_container.Quality=this->Quality.substr(entire_region.first,entire_region.second);
    split_container.R1Aps.push_back(R1a_position-entire_region.first);
    split_container.R1Ams.push_back(this->R1Ams[R1a_index]);
    split_container.R1Ass.push_back(this->R1Ass[R1a_index]);
    split_container.freeRegions.push_back({this->freeRegions[free_region_index].first-entire_region.first,
                                          this->freeRegions[free_region_index].second-entire_region.first});
    split_container.occupiedRegions.push_back({this->occupiedRegions[occ_region_index].first-entire_region.first,
                                              this->occupiedRegions[occ_region_index].second-entire_region.first});
    for(int i=0; i<this->qhits.size();i++){
      if(isContaining(this->qhits[i],entire_region_1_based_start_end)){
        split_container.qhits.push_back({this->qhits[i].first-entire_region.first,this->qhits[i].second-entire_region.first});
        split_container.qhitOrient.push_back(this->qhitOrient[i]);
        split_container.hit_ids.push_back(this->hit_ids[i]);
        split_container.hit_names.push_back(this->hit_names[i]);
        split_container.qhit_covs.push_back(this->qhit_covs[i]);
        split_container.hit_covs.push_back(this->hit_covs[i]);
        split_container.spliced.push_back(this->spliced[i]);
      }
    }
    return;
  }
  
  string generate_R1a_summary(){
    string res="";
    if(this->R1Aps.size()==0){
      return("*");
    }
    for(int index:this->R1Aps){
      res=res+to_string(index)+',';
    }
    res.erase(res.begin()+res.length()-1);
    res=res+':';
    for(int strand:this->R1Ass){
      res=res+num2rna[strand]+',';
    }
    res.erase(res.begin()+res.length()-1);
    res=res+':';
    for(int mismatch:this->R1Ams){
      res=res+to_string(mismatch)+',';
    }
    res.erase(res.begin()+res.length()-1);
    return(res);
  }
  
  pair<int, int> isConstantSequenceContaining(string& constantseq, int testRange, int max_testseq_mismatch, int strand, bool loose=false){
    string testForward;
    string testseq;
    string testReverse;
    string testseq_rc;
    pair<int, int> res;
    //try to loose the standard
    if(loose){
      max_testseq_mismatch=max_testseq_mismatch+1;
    }
    int index1=-1;
    int index2=-1;
    if(strand==0){
      testForward=this->Seq.substr(0,testRange);
      testseq=constantseq;
      index1=localAlign(testseq,testForward,max_testseq_mismatch,'l');
    }else if(strand==3){
      testReverse=this->Seq.substr(this->Seq.length()-testRange,testRange);
      testseq_rc=reverse_complement(constantseq);
      int align_result=localAlign(testseq_rc,testReverse,max_testseq_mismatch,'r');
      if(align_result<0){
        index2=align_result;
      }else{
        index2=align_result+this->Seq.length()-testRange;
      }
    }else{
      testForward=this->Seq.substr(0,testRange);
      testseq=constantseq;
      testReverse=this->Seq.substr(this->Seq.length()-testRange,testRange);
      testseq_rc=reverse_complement(constantseq);
      index1=localAlign(testseq,testForward,max_testseq_mismatch,'l');
      int rv_align_result=localAlign(testseq_rc,testReverse,max_testseq_mismatch,'r');
      if(rv_align_result<0){
        index2=rv_align_result;
      }else{
        index2=rv_align_result+this->Seq.length()-testRange;
      }
    }
    res.first=index1;
    res.second=index2;
    return res;
  } //isConstantSequenceContaining
  
  /*!
   * @abstract genrate kmers from whole read sequence
   * 
   * @param k: kmer length
   * 
   * @return fill in this->Kmers
   */
  void genrateKmer(int k){
      const int len=this->Seq.length();
      for(int j=0;j+k<=len;j++){
        string kmer=this->Seq.substr(j,k);
        this->Kmers.push_back(kmer);
      }
    }
  
  /*!
   * @abstract genrate kmers by start and end position desired
   * 
   * @param k: kmer length
   * @param start: 1-based, do not care if end greater than length
   * @param end: 1-based, do not care if end greater than length
   * 
   * @return key: kmer start position, 1-based, value: kmer sequence
   * 
   */
  unordered_map<int, string> genrateKmer(int k, int start, int end){
    unordered_map<int, string> kmers;
    const int len=this->Seq.length();
    start--;
    if(start<0){
      start=0;
    }
    if(start>len-k){
      return(kmers);
    }
    if(end>len-k){
      end=len-k;
    }else{
      end=end-k;
    }
    for(int j=start;j<=end;j++){
      kmers[j]=this->Seq.substr(j,k);
    }
    return(kmers);
  }
  
  /*!
   * @abstract identify barcodes with read annotation
   * 
   */
  
  void identifyBarcodes(vector<Read>& parent_block, BarcodeFile& barcodes, vector<int> segments, int maxMismatch, int tech, unordered_map<string, ReadAnno>& annoMap,
   const entry_t *model=nullptr, WL_DB *wl_db=nullptr, bool dtw=false){
    int barcode_length=barcodes.barcodes[0].length();
    //if read too long, discard it
    if(this->Seq.length()>=100000){
      this->add_flag('L');
      return;
    }
    if(this->Kmers.size()==0){
      this->genrateKmer(barcode_length);
    }
    
    if(!annoMap[this->ID].read_id.empty()){
      //cout << "debug1.1" << endl;
      vector<pair<int,int>> freeRegions=annoMap[this->ID].free_regions;
      this->freeRegions=freeRegions;
      this->qhits=annoMap[this->ID].qhits;
      this->hit_overlaps=annoMap[this->ID].hit_overlaps;
      this->qhitOrient=annoMap[this->ID].read_orientation;
      this->occupiedRegions=annoMap[this->ID].occupied_regions;
      this->hit_ids=annoMap[this->ID].gene_id;
      this->hit_names=annoMap[this->ID].gene_name;
      this->hit_covs=annoMap[this->ID].rcovs;
      this->qhit_covs=annoMap[this->ID].qcovs;
      this->spliced=annoMap[this->ID].is_spliced;
      //cout << "debug1.2" << endl;
    }

    if(this->qhits.empty()){
      this->add_flag('U');
      return;
    }
    
    if(this->R1Aps.empty()){
      for(pair<int,int> region:freeRegions){
        if(region.second-region.first+1<TENX_R1.length()){
          continue;
        }else{
          this->findR1a(maxMismatch*3,region.first,region.second);
        }
      }
    }
    
    //cout << "debug1.3" << endl;
    if(this->TSOps.empty()){
      int max_tso_mismatch=maxMismatch;
      if(tech==3){
        max_tso_mismatch=max_tso_mismatch*3;
        TENX_TSO="CCCATGTACTCTGCGTTGATACCACTGCTT";
      }
      for(pair<int,int> region:freeRegions){
        if(region.second-region.first+1<TENX_TSO.length()){
          continue;
        }else{
          this->findTSO(max_tso_mismatch,region.first,region.second,tech);
        }
      }
    }
    
    //cout << "debug1.4" << endl;
    if(this->R1Aps.empty()){
      if(dtw){
        this->identifybarcodesDTW(wl_db,model);
        if(this->barcode=="*"){
          this->add_flag('M');
          return;
        }
      }else{
        this->add_flag('R');
        return;
      }
    }else if(this->R1Ams.size()==1){
      this->OUTER=R1Aps[0];
      unordered_map<int, string> kmers;
      if(this->R1Ass[0]==0){
        kmers=genrateKmer(barcode_length,this->R1Aps[0]+1+TENX_R1.length()-maxMismatch*2,this->R1Aps[0]+1+TENX_R1.length()+barcode_length+maxMismatch);
        for(auto kmer:kmers){
          findmismatchingbarcodebykmer(barcodes,segments,kmer.second,kmer.first+1,"forward",maxMismatch);
        }
      }else{
        kmers=genrateKmer(barcode_length,this->R1Aps[0]+1-barcode_length-maxMismatch,this->R1Aps[0]+1+maxMismatch*2);
        for(auto kmer:kmers){
          findmismatchingbarcodebykmer(barcodes,segments,kmer.second,kmer.first+1,"reverse",maxMismatch);
        }
      }
      
      //cout << "debug1.5" << endl;
      //if multiple barcode found, get the barcode with min mismatch
      if(this->barcode_mismatches.size()>1){
        int minPosition=min_element(this->barcode_mismatches.begin(),this->barcode_mismatches.end()) - this->barcode_mismatches.begin();
        this->barcode=barcode_intolerance_detail[minPosition];
        this->barcodeStart=barcode_intolerance_indices[minPosition];
        this->barcodeMismatch=barcode_mismatches[minPosition];
      }
      
      //if barcode can not be found through alignment based method, use DTW instead 
      if((this->barcode.length()==0||this->barcode.length()==1)&&dtw){
        this->identifybarcodesDTW(wl_db,model);
      }
      
      if(this->barcode!="*"){
        this->add_annotation(tech);
      }else{
        this->add_flag('M');
      }
    }else{
      //cout << "debug1.6" << endl;
      this->add_flag('C');
      this->add_complex_structure_annotation();
      for(int i=0;i<this->R1Aps.size();++i){
        Read split_read;
        split_read.ID=this->ID+"-"+to_string(i);
        //cout << "debug1.7" << endl;
        find_complete_region_by_R1a(split_read,this->R1Aps[i],this->R1Ass[i],i);
        //cout << "debug1.8" << endl;
        if(!split_read.Seq.empty()){
          parent_block.push_back(split_read);
          continue;
        }else{
          continue;
        }
      }
    }
  }
  
  /*!
   * @abstract identify barcodes without mapping info
   * 
   */
  void identifyBarcodes(BarcodeFile& barcodes, vector<int> segments, int maxDistToEnds, int maxMismatch, int tech, const entry_t *model=nullptr, WL_DB *wl_db=nullptr, bool dtw=false){
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
    
    //if read too long, discard it
    if(this->Seq.length()>=100000){
      this->add_flag('L');
      return;
    }
    
    //if read shorter than barcode search range * 2, discard
    if(this->Seq.length()<=maxDistToEnds*2){
      this->add_flag('S');
      return;
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
        outer_stat=isConstantSequenceContaining(TENX_R1,maxDistToEnds,6,0,true);
      }else if(tech==3){
        inner_stat=isConstantSequenceContaining(POLYT,maxDistToEnds,0,0);
        outer_stat=isConstantSequenceContaining(TENX_R1,maxDistToEnds,6,0,true);
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
        outer_stat=isConstantSequenceContaining(TENX_R1,maxDistToEnds,6,3,true);
      }else if(tech==3){
        inner_stat=isConstantSequenceContaining(POLYT,maxDistToEnds,0,3);
        outer_stat=isConstantSequenceContaining(TENX_R1,maxDistToEnds,6,3,true);
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
        outer_stat=isConstantSequenceContaining(TENX_R1,maxDistToEnds,6,-1,true);
      }else if(tech==3){
        inner_stat=isConstantSequenceContaining(POLYT,maxDistToEnds,0,-1);
        outer_stat=isConstantSequenceContaining(TENX_R1,maxDistToEnds,6,-1,true);
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
        return;
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
        this->stat="Dual_header";
        return;
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
        this->stat="Dual_header";
        return;
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
        this->stat="Dual_header";
        return;
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
        this->stat="Dual_header";
        return;
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
          this->stat="Dual_header";
          return;
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
          this->stat="Dual_header";
          return;
        }
      }
    }
    if(kmer_index_start>kmer_index_end){
      this->barcode="*";
      this->barcodeStrand="*";
      this->barcodeStart=-1;
      this->barcodeMismatch=-1;
      this->stat="ICS_R1_wrongly_oriented";
      return;
    }
    std::vector<string> kmers;
    //determine kmers by inner and outer constant seq start index
    
    //1. both CS found at left end
    if(this->INNER>0&&this->OUTER>0&&this->INNER<Seq.length()/2){
      kmer_index_start=this->OUTER-1+TENX_R1.length()-6;
      kmers.insert(kmers.end(),this->Kmers.begin()+this->OUTER-1+TENX_R1.length()-6,this->Kmers.begin()+this->OUTER-1+TENX_R1.length()+6);
      
    //2. both CS found at right end
    }else if(this->INNER>0&&this->OUTER>0&&this->INNER>Seq.length()/2){
      kmer_index_start=this->OUTER-1-barcode_length-6;
      kmers.insert(kmers.end(),this->Kmers.begin()+this->OUTER-1-barcode_length-6,this->Kmers.begin()+this->OUTER-1-barcode_length+6);
      
    //3. only ICS found at left end
    }else if(this->INNER>0&&this->OUTER<0&&this->INNER<Seq.length()/2){
      if(this->INNER-1-barcode_length-umiLength-6<0){
        this->barcode="*";
        this->barcodeStrand="*";
        this->barcodeStart=-1;
        this->barcodeMismatch=-1;
        this->stat="Read_too_short";
        return;
      }else{
        kmer_index_start=this->INNER-1-barcode_length-umiLength-6;
        kmers.insert(kmers.end(),this->Kmers.begin()+this->INNER-1-barcode_length-umiLength-6,this->Kmers.begin()+this->INNER-1-barcode_length-umiLength+6);
      }
      
    //4. only OCS found at left end
    }else if(this->INNER<0&&this->OUTER>0&&this->OUTER<Seq.length()/2){
      kmer_index_start=this->OUTER-1+TENX_R1.length()-6;
      kmers.insert(kmers.end(),this->Kmers.begin()+this->OUTER-1+TENX_R1.length()-6,this->Kmers.begin()+this->OUTER-1+TENX_R1.length()+6);
      
    //5. only ICS found at right end
    }else if(this->INNER>0&&this->OUTER<0&&this->INNER>Seq.length()/2){
      if(this->INNER+TENX_TSO.length()+umiLength+barcode_length+6-1>this->Seq.length()){
        this->barcode="*";
        this->barcodeStrand="*";
        this->barcodeStart=-1;
        this->barcodeMismatch=-1;
        this->stat="Read_too_short";
        return;
      }else{
        kmer_index_start=this->INNER-1+TENX_TSO.length()+umiLength-6;
        kmers.insert(kmers.end(),this->Kmers.begin()+this->INNER-1+TENX_TSO.length()+umiLength-6,this->Kmers.begin()+this->INNER-1+TENX_TSO.length()+umiLength+6);
      }
      
    //6. only OCS found at right end
    }else if(this->INNER<0&&this->OUTER>0&&this->OUTER>Seq.length()/2){
      kmer_index_start=this->OUTER-1-barcode_length-6;
      kmers.insert(kmers.end(),this->Kmers.begin()+this->OUTER-1-barcode_length-6,this->Kmers.begin()+this->OUTER-1-barcode_length+6);
      
    //7. both CS not found
    }else{
      this->barcode="*";
      this->barcodeStrand="*";
      this->barcodeStart=-1;
      this->barcodeMismatch=-1;
      this->stat="ICS_R1_missing";
      return;
    }
    
    //if no valid kmer can be used
    if(kmers.size()==0){
      this->barcode="*";
      this->barcodeStrand="*";
      this->barcodeStart=-1;
      this->barcodeMismatch=-1;
      this->stat="Read_too_short";
      return;
    }
    
    //for kmers, try to find exact match in barcode list, else to find ambiguous match
    for(int j=0;j<kmers.size();j++){
      int index=kmer_index_start+j+1;
      findexactbarcodebykmer(barcodes,kmers[j],index,num2rna[constant_seq_strand]);
    } //for(int j=0;j<kmers.size();j++)
    
    if(this->barcode.length()!=0){
//      substring to generate a seq downstream barcode for 10nt as UMI, 
//      and modify reverse_com seq to original one, cut off adapter to UMI
      //this->guessUMI(tech);
      //start to find ambiguous match
    }else{
      for(int j=0;j<kmers.size();j++){
        int index=kmer_index_start+j+1;
        findmismatchingbarcodebykmer(barcodes, segments, kmers[j], index, num2rna[constant_seq_strand], maxMismatch);
      }//for(int j=0;j<kmers.size();j++)
    }// if(tmpRead.barcode.length()!=0) else
    
    //if multiple barcode found, get the barcode with min mismatch
    if(barcode_mismatches.size()>1){
      int minPosition=min_element(barcode_mismatches.begin(),barcode_mismatches.end()) - barcode_mismatches.begin();
      this->barcode=barcode_intolerance_detail[minPosition];
      this->barcodeStart=barcode_intolerance_indices[minPosition];
      this->barcodeMismatch=barcode_mismatches[minPosition];
    }
    
    //if barcode can not be found through alignment based method, use DTW instead 
    if((this->barcode.length()==0||this->barcode.length()==1)&&dtw){
      this->identifybarcodesDTW(wl_db,model);
    }
    
    if(this->barcode!="*"){
      return;
    }else{
      this->add_flag('M');
      return;
    }
  } //identifyBarcodes

private:
  void guessUMI(int tech){
    if(tech==3){
      if(this->barcode.find("*")!=0){
        if(this->barcodeStrand.find("forward")==0){
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
        if(this->barcodeStrand.find("forward")==0){
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
  
  void identifybarcodesDTW(WL_DB *wl_db, const entry_t *model){
    long max_size=100000;
    level_t lvl[max_size];
    char qual[max_size],qual2[max_size];
    char qbase[100]="CTACACGACGCTCTTCCGATCT";
    level_t qlvl[100];
    char qbase2[100];
    strcpy(qbase2, qbase);
    level_t qlvl2[100];
    revcom(qbase2);
    int qlvl_len=strlen(qbase)-MER_LENGTH+1;
    seq2level(model, qbase, qlvl, 100, 0);
    seq2level(model, qbase2, qlvl2, 100, 0);
//    fp=gzopen(argv[2], "r");
//    seq=kseq_init(fp);
    result_t *resq=(result_t *) calloc(2, sizeof(result_t));
    result_t *wl_resq=(result_t *) calloc(2, sizeof(result_t));
    result_t the_res;   // not initialized but ....
    int ret;
    int r=3; // specify bandwidth here
    level_t bc[128];
    level_t u_d[128], l_d[128];
    char tmpstr[128];  //for printing purpose
    seq2level(model, Seq.c_str(), lvl, max_size, r);
    seq2merqual(Quality.c_str(), qual, max_size, r);
    seq2merqual2(Quality.c_str(), qual2, max_size, r);
    reset_resultq(resq);
    ret = ucrdtw(lvl, Seq.length() - MER_LENGTH+1+r, qlvl, qlvl_len, r, 0, resq, 200, 1);
    ret = ucrdtw(lvl, Seq.length() - MER_LENGTH+1+r, qlvl2, qlvl_len, r, 0, resq, add_result_cap2(resq, -1,INF,1,0,5), -1);
    the_res=pick_result_cap2(resq, Seq.length()-MER_LENGTH+1);
    reset_resultq(wl_resq);
    
    if(the_res.location!=-1){
//      cout << "debug2.1" << endl;
      
      if(the_res.flag==1){//forward
//        cout << "debug2.1.1" << endl;
        //dtw3_bt(lvl+the_res.location, qlvl, qual+the_res.location, qual2+the_res.location, qlvl_len, r);
        level_t best=500;
        long best_loc=-1;
        long j=0;
        long bc_len=33+3;
        for(long i=the_res.location; i< the_res.location + bc_len ; i++){
          bc[j++]=lvl[i];
        }
        //printLevels(bc, bc_len);
        strncpy(tmpstr, Seq.c_str()+the_res.location, bc_len+MER_LENGTH-1);
        tmpstr[bc_len+MER_LENGTH-1]=0;
//        printf("debugSeq %s\n", tmpstr);
        lower_upper_lemire(bc, bc_len, r, l_d, u_d);
        //cout << "debug2.1.1.1" << endl;
        //cout << "wl_db size: " << wl_db->size << endl;
//        cout << "debug2.1.1.1" << endl;
        for(long i=0; i<wl_db->size; i++){
          //printLevels(wll[i],33);
//          cout << i << endl;
          reset_resultq(resq);
          //cout << "debug2.1.1.1.1" << endl;
          ret=ucrdtw_checkWL(bc, bc_len, u_d, l_d, wl_db, i, 0, 0, resq, best, 1);
          //cout << "debug2.1.1.1.2" << endl;
          //debug
          //printf("tmpdebug %s\n", wl_db->seq[i]+22);
          //dtw3_bt(lvl+the_res.location, wl_db->lvl[i], qual+the_res.location, qual2+the_res.location, 33,r );
          //debug
          //printResult(resq);
          if(resq[0].distance<best){
            best=add_result_cap2(wl_resq, i, resq[0].distance, 1, 0,0);
          }
        }
//        cout << "debug2.1.2" << endl;
        if(wl_resq[0].location>=0){
//          cout << getSeq(wl_db,wl_resq[0].location,22) << endl;
//          cout << the_res.location << endl;
//          cout << wl_resq[0].location << endl;
          
          //dtw3_bt(lvl+the_res.location, wl_db->lvl[wl_resq[0].location], qual+the_res.location, qual2+the_res.location, 33,r ); //buggy
          //dtw3_bt(lvl+the_res.location, wl_db->lvl[wl_resq[1].location], qual+the_res.location, qual2+the_res.location, 33,r );
          this->stat="fine_DTW";
          this->barcode=getSeq(wl_db,wl_resq[0].location,22);
          this->barcodeStrand="forward";
          this->barcodeMismatch=(int)(wl_resq[0].distance-the_res.distance);
//          printf("final %s %.4f forward %s %.4f %s %.4f\n", ID.c_str(), the_res.distance,getSeq(wl_db,wl_resq[0].location,22),wl_resq[0].distance,getSeq(wl_db,wl_resq[1].location,22),wl_resq[1].distance);
        }else{
//          printf("debugX,notFound\nfinal,notFound\n");
        }
//      cout << "debug2.1.3" << endl;
      }else{         //reverse
        //dtw3_bt(lvl+the_res.location, qlvl2, qual+the_res.location,qual2+the_res.location,qlvl_len, r);
        level_t best=500;
        long best_loc=-1;
        long j=0;
        long bc_len=33+8;
        for(long i=the_res.location-19; i< the_res.location-19 + bc_len; i++){
          bc[j++]=lvl[i];
        }
        
        //printLevels(bc, bc_len);
//        cout << "debug2.1.3.1" << endl;
        strncpy(tmpstr, Seq.c_str()+the_res.location-19, bc_len+MER_LENGTH-1);
//        cout << "debug2.1.3.2" << endl;
        tmpstr[bc_len+MER_LENGTH-1]=0;
//        printf("debugSeq %s\n", tmpstr);
        lower_upper_lemire(bc, bc_len, r, l_d, u_d);
//        cout << "debug2.1.4" << endl;
        for(long i=0; i<wl_db->size; i++){
          //printLevels(wllr[i],33);
          reset_resultq(resq);
          ret=ucrdtw_checkWL(bc, bc_len, u_d, l_d, wl_db, i, 1, 0, resq, best, -1);			
          //printResult(resq);
          if(resq[0].distance<best){
            //best=add_result_cap2(wl_resq, i, resq[0].distance, -1, 0,0);
            best=add_result_cap2(wl_resq, i, resq[0].distance, resq[0].location, 0,0);
            //offset=resq[0].location;
          }
        }
        //cout << "debug2.1.5" << endl;
        if(wl_resq[0].location>=0){
//          cout << getSeq(wl_db,wl_resq[0].location,22) << endl;
//          printf("debugY,");
//          printResult(wl_resq);
          //printf(",%s,",wl[wl_resq[0].location]);
          //printf(",%s,",wl_db->seq[wl_resq[0].location]);
//          printf(",%s,",getSeq(wl_db, wl_resq[0].location,0));
//          printResult(&wl_resq[1]);
//          printf("\n");
          //dtw3_bt(lvl+the_res.location-19+wl_resq[0].flag, wl_db->lvl_r[wl_resq[0].location], qual+the_res.location-19+wl_resq[0].flag, qual2+the_res.location-19+wl_resq[0].flag, 33,r );
          this->stat="fine_DTW";
          this->barcode=getSeq(wl_db,wl_resq[0].location,22);
          this->barcodeStrand="reverse";
          this->barcodeMismatch=(int)(wl_resq[0].distance-the_res.distance);
          //dtw3_bt(lvl+the_res.location-19+wl_resq[1].flag, wl_db->lvl_r[wl_resq[1].location], qual+the_res.location-19+wl_resq[1].flag, qual2+the_res.location-19+wl_resq[1].flag, 33,r );
          //printf("final\t%s\t%.4f\treverse\t%s\t%.4f\t%s\t%.4f\n", ID.c_str(), the_res.distance,getSeq(wl_db,wl_resq[0].location,22),wl_resq[0].distance,getSeq(wl_db,wl_resq[1].location,22),wl_resq[1].distance);
          //printf("final\t%s\t%.4f\treverse\t%s\t%.4f\t%s\t%.4f\n", seq->name.s, the_res.distance,wl_db->seq[wl_resq[0].location]+22,wl_resq[0].distance,wl_db->seq[wl_resq[1].location]+22,wl_resq[1].distance);
        }else{
//          printf("debugY,notFound\nfinal,notFound\n");
        } 
        
      }
    
    }
  } //void identifybarcodesDTW
  
  /*!
   * @abstract find exact matched barcode with a kmer sequence
   * 
   * @param barcodes: barcode dicts, refer to class BarcodeFile
   * @param kmer: kmer sequence
   * @param kmer_pos: kmer start position on read, 1-based
   * @param kmer_strand: "forward"/"reverse"
   */
  void findexactbarcodebykmer(BarcodeFile& barcodes, string kmer, int kmer_pos, string kmer_strand){
    if(barcodes.barcode_ori[kmer].length()!=0&&kmer_strand=="forward"){
      this->barcode=kmer;
      this->barcodeStrand="forward";
      this->barcodeStart=kmer_pos;
      this->barcodeMismatch=0;
      this->stat="fine";
      this->barcode_in_tolerance++;
      this->barcode_intolerance_detail.push_back(kmer);
      this->barcode_intolerance_indices.push_back(kmer_pos);
      this->barcode_mismatches.push_back(0);
    }else if(barcodes.barcode_rc[kmer].length()!=0&&kmer_strand=="reverse"){
      this->barcode=barcodes.barcode_rc[kmer];
      this->barcodeStrand="reverse";
      this->barcodeStart=kmer_pos;
      this->barcodeMismatch=0;
      this->stat="fine";
      this->barcode_in_tolerance++;
      this->barcode_intolerance_indices.push_back(kmer_pos);
      this->barcode_intolerance_detail.push_back(barcodes.barcode_rc[kmer]);
      this->barcode_mismatches.push_back(0);
    }
  } // void findexactbarcodebykmer
  
  /*!
   * @abstract find barcode with mismatch with a kmer sequence
   * 
   * @param barcodes: barcode dicts, refer to class BarcodeFile
   * @param kmer: kmer sequence
   * @param kmer_pos: kmer start position on read, 1-based
   * @param kmer_strand: "forward"/"reverse"
   * @param max_mismatch
   */
  void findmismatchingbarcodebykmer(BarcodeFile& barcodes, vector<int> segments, string kmer, int kmer_pos, string kmer_strand, int max_mismatch){
    //first is index in barcodes.barcode_dict, second is strand
    vector<pair<int, int>> candidates;
    vector<string> kmer_segment=getMaxComplexitySegments(kmer,segments);
    for(int n=0;n<kmer_segment.size();n++){
      string query=kmer_segment[n]+to_string(n)+to_string(strand2num[kmer_strand]);
      candidates.insert(candidates.end(),barcodes.barcode_dict[query].begin(),barcodes.barcode_dict[query].end());
    }
    if(candidates.size()==0){
      return;
    }else{
      for(int m=0;m<candidates.size();m++){
        int barcode_dict_index=candidates[m].first-1;
        string candidate_ori=barcodes.barcodes[barcode_dict_index];
        if(std::count(this->barcode_intolerance_detail.begin(),this->barcode_intolerance_detail.end(),candidate_ori)){
          continue;
        }
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
            cout << "Unsupported strand coding, Read: " << this->ID << ", kmer: " << kmer << endl;
          }
        }
        int mismatch=editDistance(candidate,kmer,max_mismatch);
        //take care here!!!! once kmer meets the threshould, it will be considered to be the barcode, 
        //even there may be better ones after it.
        if(mismatch>max_mismatch){
          continue;
        }else{
          this->barcode=candidate_ori;
          this->barcodeStrand=num2rna[strand];
          this->barcodeStart=kmer_pos;
          this->barcodeMismatch=mismatch;
          this->stat="fine";
          this->barcode_in_tolerance++;
          this->barcode_intolerance_detail.push_back(candidate_ori);
          this->barcode_intolerance_indices.push_back(kmer_pos);
          this->barcode_mismatches.push_back(mismatch);
          //              this->matched_kmer_quality=this->Quality.substr(index,16);
          continue;
        }
      } //for(int m=0;m<candidates.size();m++)
    }// if(candidates.size()==0) else
  }
  
};


struct AmbiguousBarcode{
  string seq;
  int strand;
  int mismatch;
};

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
  
  //constructor identify all R1a on reads
  ReadFile(const char* filename, bool best=true){
    ReadCount=0;
    File=gzopen(filename, "r");
    seq = kseq_init(File);
    while (kseq_read(seq) >= 0){
      ReadCount++;
      Read tmpRead;
      tmpRead.ID=seq->name.s;
      tmpRead.Seq=seq->seq.s;
      tmpRead.Quality=seq->qual.s;
//      cout << "identifying R1a" << endl;
      tmpRead.findR1a(4,1,tmpRead.Seq.length());
      string R1asum=tmpRead.generate_R1a_summary();
      cout << tmpRead.ID << "\t" << R1asum << endl;
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
      if(tmpRead.barcodeStrand=="forward"){
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
      if(tmpRead.barcodeStrand=="forward"){
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
  ReadFile(const char* filename, BarcodeFile& barcodes, vector<int> segments, int maxDistToEnds,
           int maxMismatch, int tech, int threadnum, unordered_map<string, ReadAnno>& annoMap, 
           const entry_t *model=nullptr, WL_DB *wl_db=nullptr, bool dtw=false,bool debug=false){
    fprintf(stderr,"Start BC identification and read annotation\n");
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
        if(debug){
          //cout << tmpRead.ID << endl;
          vector<Read> debug_chimeras;
          //cout << "debug1" << endl;
          tmpRead.identifyBarcodes(debug_chimeras,barcodes,segments,maxMismatch,tech, annoMap, model, wl_db, dtw);
          //cout << "debug2" << endl;
          tmpRead.printBarcodeInfo();
          //cout << "debug3" << endl;
          if(!debug_chimeras.empty()){
            //cout << "debug4" << endl;
            for(Read read:debug_chimeras){
              read.identifyBarcodes(debug_chimeras,barcodes,segments,maxMismatch,tech, annoMap, model, wl_db, dtw);
              read.printBarcodeInfo();
            }
          }
          continue;
        }
        Reads.push_back(tmpRead);
        //cout << tmpRead.ID << endl;
        if(Reads.size()==50){
          //cout << "block full" << endl;
          Readmatrix.push_back(Reads);
          //cout << "number of block: " << Readmatrix.size() << endl;
          threadcount++;
          Reads.clear();
          if(Readmatrix.size()==threadnum){
            //cout << "Batch start" << endl;
            #pragma omp parallel for
            for(int i=0;i<Readmatrix.size();++i){
              BarcodeFile threadbarcode=barcodes;
              vector<int> threadsegments=segments;
              int threadmaxDistToEnds=maxDistToEnds;
              int threadmaxMismatch=maxMismatch;
              int threadtech=tech;
              if(annoMap.empty()){
                identifyBarcodesinBlock(Readmatrix[i],threadbarcode,threadsegments,threadmaxDistToEnds,threadmaxMismatch,threadtech, model, wl_db, dtw);
              }else{
                identifyBarcodesinBlockWithAnno(Readmatrix[i],threadbarcode,threadsegments,threadmaxMismatch,threadtech,annoMap, model, wl_db, dtw);
              }
            } //for(int i=0;i<threadnum;++i)
            //cout << "Batch end" << endl;
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
        #pragma omp parallel for
        for(int i=0;i<Readmatrix.size();++i){
          BarcodeFile threadbarcode=barcodes;
          vector<int> threadsegments=segments;
          int threadmaxDistToEnds=maxDistToEnds;
          int threadmaxMismatch=maxMismatch;
          int threadtech=tech;
          if(annoMap.empty()){
            identifyBarcodesinBlock(Readmatrix[i],threadbarcode,threadsegments,threadmaxDistToEnds,threadmaxMismatch,threadtech, model, wl_db, dtw);
          }else{
            identifyBarcodesinBlockWithAnno(Readmatrix[i],threadbarcode,threadsegments,threadmaxMismatch,threadtech,annoMap, model, wl_db, dtw);
          }
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
  
  void identifyBarcodesinBlock(vector<Read>& blockreads,BarcodeFile& barcode,vector<int> segments,
                               int maxDistToEnds,int maxMismatch,int tech,
                               const entry_t *model=nullptr, WL_DB *wl_db=nullptr, bool dtw=false){
    for(int i=0;i<blockreads.size();++i){
      blockreads[i].identifyBarcodes(barcode,segments,maxDistToEnds,maxMismatch, tech, model, wl_db, dtw);
    }
  }
  
  void identifyBarcodesinBlockWithAnno(vector<Read>& blockreads,BarcodeFile& barcode,vector<int> segments,
                                       int maxMismatch,int tech, unordered_map<string, ReadAnno>& annoMap,
                                       const entry_t *model=nullptr,WL_DB *wl_db=nullptr, bool dtw=false){
    for(int i=0;i<blockreads.size();++i){
      blockreads[i].identifyBarcodes(blockreads,barcode,segments,maxMismatch,tech, annoMap, model, wl_db, dtw);
    }
  }
  
  void writeReadData(vector<vector<Read>>& readmatrix,std::stringstream &readinfostream,
                     std::stringstream &filteredfastqstream, std::stringstream &trimfilteredfastqstream){
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

void readAlignmentSummaryFile(unordered_map<string, ReadAnno>& annoMap, string filename){
  fprintf(stderr,"Start to read annotation file: %s\n",filename.c_str());
  fstream align_summary_file;
  align_summary_file.open(filename,fstream::in);
  int i=0;
  while(!align_summary_file.eof()){
    string line;
    ReadAnno thisAnno;
    getline(align_summary_file,line,'\n');
    i++;
    if(line==""){
      continue;
    }
    vector<string> fields=tokenize(line,'\t');
    thisAnno.read_id=fields[0];
    thisAnno.gene_id=tokenize(fields[1],',');
    thisAnno.gene_name=tokenize(fields[2],',');
    thisAnno.read_orientation=string2strands(fields[3]);
    thisAnno.free_regions=string2regions(fields[4]);
    thisAnno.qhits=string2regions(fields[5]);
    thisAnno.qcovs=string2covs(fields[6]);
    thisAnno.rcovs=string2covs(fields[8]);
    thisAnno.is_spliced=string2splice(fields[11]);
    thisAnno.hit_overlaps=string2hit_overlaps(fields[9]);
    thisAnno.occupied_regions=string2regions(fields[10]);
    annoMap[thisAnno.read_id]=thisAnno;
  }
  fprintf(stderr,"Annotation file read.\n");
}

