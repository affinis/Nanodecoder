#include <unistd.h>
#include <stdlib.h>
#include <cstring>
#include <string>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <htslib/sam.h>
#include <regex>
#include <map>

#include "prosessSeq.hpp"

using std::ifstream;
using std::ofstream;
using std::ios;
using std::cout;
using std::endl;
using std::string;
using std::pair;
using std::vector;
using std::unordered_map;
using std::to_string;
using std::regex;
using std::regex_replace;
using std::map;
using std::setprecision;
using std::reverse;

#define bam_is_umapped(b) (((b)->core.flag&BAM_FUNMAP) != 0)

regex gene_regex("(gene_id \"| gene_type \"| gene_biotype \"| gene_name \"|\")");
regex transcript_regex("( transcript_id \"| transcript_type \"| transcript_name \"|\")");

bool isSpliced(uint32_t* test_cigar, int ncigar){
  for(int i=0;i<ncigar;i++){
    if(bam_cigar_opchr(test_cigar[i])=='N'){
      return(true);
    }
  }
  return(false);
}



/*!
 * @abstract get specific field from aux field of GTF file with key_word, eg: "transcript_id", "transcript_name" etc.
 * 
 * @param key_word: key word a field contain, must to be specific.
 * @param aux_data: vector of string from tokenization of auxiliary filed of GTF file by character ';', notice that 
 *                  each filed may contain space or other undesired character, further processing may be needed.
 *                  
 * @param afileds_info: a key-value map storing field number, the default value each key point to is -2, if the function
 *                      find field number is -2, it will try to find the exact position of desired field and update 
 *                      this key-value map. If the exact position can not be found, it will be updated to -1. 
 * @feature_type: feature type by 3rd field of GTF file, eg "gene", "exon", "transcript", this value can affect the position
 *                of the desired field in aux_data.
 *                
 * @return: raw string from aux_data, eg: 'gene_id "LOC116943895"', ' db_xref "GeneID:116943895"', etc.
 */  
string get_aux_field(const string& key_word,vector<string>& aux_data, unordered_map<string, int>& afileds_info, const string& feature_type){
  string res;
  string query_key=feature_type+"_"+key_word+"_"+"field";
  if(afileds_info[query_key]==-2){
    int j=0;
    for(string field:aux_data){
      if(field.find(key_word)!=-1){
        afileds_info[query_key]=j;
        break;
      }else{
        j++;
        continue;
      }
    }
    if(afileds_info[query_key]==-2){
      afileds_info[query_key]=-1;
    }
  }
  if(afileds_info[query_key]==-1){
    res="*";
  }else{
    res=aux_data[afileds_info[query_key]];
  }
  
  return(res);
}

/*!
 * @abstract sequence feature
 * 
 * @start 1-based
 * @end 1-based
 * @id
 * 
 */
struct seq_feature{
  int start;
  int end;
  int strand;
  string id;
};

/*!
 * @abstract exon feature
 * 
 * @pos 1-based start and end position of exon.
 * @id id of this exon to its respective gene.
 */
struct exon_feature{
  pair<int, int> pos;
  int id;
};

/*!
 * @abstract transcript feature
 * 
 * @start 1-based
 * @end 1-based
 * @transcript_id
 * @transcript_name
 * @exons pairs of 1-based start&end sites of this gene indicating exon positions
 */
struct transcript_feature{
  int start=-1;
  int end=-1;
  string transcript_id;
  string transcript_name;
  vector<int> exons;
  vector<int> loose_exons;
};

/*!
 * @abstract gene features read from GTF annotation file
 * 
 * @index 1-based number indicate the order this feature parsed from GTF
 * @start 1-based
 * @end 1-based
 * @strand 0=forward/+, 1=reverse/-
 * @gene_type
 * @chromosome not strictly defined, can be certain scaffold that haven't been integrated
 * @gene_id
 * @gene_name
 * @exons
 * @transcripts
 */
struct gene_feature{
  int index;
  int start;
  int end;
  int strand;
  string gene_type;
  string chromosome;
  string gene_id;
  string gene_name;
  map<int, exon_feature> exons;
  map<int, exon_feature> loose_exons;
  unordered_map<string,transcript_feature> transcripts;
  
  void add_new_exon(int exon_key, pair<int,int> pos){
    exon_feature new_exon;
    new_exon.id=this->exons.size();
    //fprintf(stderr,"%i: %i-%i\n",new_exon.id,pos.first,pos.second);
    new_exon.pos=pos;
    this->exons[exon_key]=new_exon;
  }
  
  void add_new_loose_exon(int exon_key, pair<int,int> pos){
    exon_feature new_exon;
    new_exon.id=this->loose_exons.size();
    //fprintf(stderr,"%i: %i-%i\n",new_exon.id,pos.first,pos.second);
    new_exon.pos=pos;
    this->loose_exons[exon_key]=new_exon;
  }
  
  void update_loose_exon_info(int exon_key, pair<int,int> new_pos_to_integrate){
    int new_exon_key;
    exon_feature new_exon;
    int old_exon_id=this->loose_exons[exon_key].id;
    new_exon.id=old_exon_id;
    int new_start=std::min(this->loose_exons[exon_key].pos.first,new_pos_to_integrate.first);
    int new_end=std::max(this->loose_exons[exon_key].pos.second,new_pos_to_integrate.second);
    new_exon.pos={new_start,new_end};
    new_exon_key=ceil((new_start+new_end)/2);
    this->loose_exons.erase(exon_key);
    this->loose_exons[new_exon_key]=new_exon;
  }
  
  /*!
   * @abstract fill in data from parsed GTF fields, notice when input data are from transcript/exon rather than gene, 
   *           related gene_id/transcript_id should already exists. 
   *            
   * @param core_data, tokenized (with '\t') fields, corresponding to first 8 column of GTF
   * @param aux_data, last (9th) column of GTF, tokenized with ';', the info in this fields depends on core_data[2] (feature type)
   */
  void fill_in_feature_info(vector<string>& core_data, vector<string>& aux_data, 
                            unordered_map<string, int>& afileds_info){
    string gene_type_raw;
    string gene_name_raw;
    string transcript_id_raw;
    string transcript_name_raw;
    if(core_data[2]=="gene"){
      this->chromosome=core_data[0];
      this->start=stoi(core_data[3]);
      this->end=stoi(core_data[4]);
      this->strand=strand2num[core_data[6]];
      string gene_id_raw=aux_data[0];
      gene_type_raw=get_aux_field("gene_type",aux_data,afileds_info,"gene");
      if(gene_type_raw=="*"){
        gene_type_raw=get_aux_field("gene_biotype",aux_data,afileds_info,"gene");
      }
      gene_name_raw=get_aux_field("gene_name",aux_data,afileds_info,"gene");
      this->gene_id=regex_replace(gene_id_raw,gene_regex,"");
      this->gene_name=regex_replace(gene_name_raw,gene_regex,"");
      this->gene_type=regex_replace(gene_type_raw,gene_regex,"");
    }else if(core_data[2]=="transcript"){
      string i_gene_id=regex_replace(aux_data[0],gene_regex,"");
      transcript_id_raw=get_aux_field("transcript_id",aux_data,afileds_info,"transcript");
      string transcript_id=regex_replace(transcript_id_raw,transcript_regex,"");
      transcript_feature i_transcript;
      if(this->gene_id==i_gene_id){
        i_transcript.start=stoi(core_data[3]);
        i_transcript.end=stoi(core_data[4]);
        transcript_name_raw=get_aux_field("transcript_name",aux_data,afileds_info,"transcript");
        i_transcript.transcript_name=regex_replace(transcript_name_raw,transcript_regex,"");
        i_transcript.transcript_id=transcript_id;
        this->transcripts[transcript_id]=i_transcript;
      }else{
        fprintf(stderr,"fill_in_feature_info: ERROR: fill_in_feature_info: type: transcript, gene_id: %s does not match to transcript_id: %s, check GTF\n", \
                this->gene_id.c_str(),transcript_id.c_str());
        exit(0);
      }
    }else if(core_data[2]=="exon"){
      pair<int,int> this_exon_pos={stoi(core_data[3]),stoi(core_data[4])};
      int exon_key=ceil((this_exon_pos.first+this_exon_pos.second)/2);
      
      bool loose_exon_recorded=false;
      int loose_recorded_id;
      //fprintf(stderr,"debug1\n");
      
      if(this->loose_exons.find(exon_key)==this->loose_exons.end()){
        auto exon=this->loose_exons.lower_bound(exon_key);
        //fprintf(stderr,"debug2\n");
        if(exon==this->loose_exons.begin()){
            //fprintf(stderr,"debug1: %i-%i\n",this_exon_pos.first,this_exon_pos.second);
          if(this->loose_exons.empty()){
            loose_exon_recorded=false;
            this->add_new_loose_exon(exon_key,this_exon_pos);
          }else if(!region_is_overlap(this_exon_pos,exon->second.pos)){
            loose_exon_recorded=false;
            this->add_new_loose_exon(exon_key,this_exon_pos);
          }else{
            loose_exon_recorded=true;
            this->update_loose_exon_info(exon->first,this_exon_pos);
            loose_recorded_id=exon->second.id;
          }
        }else if(exon==this->loose_exons.end()){
            //fprintf(stderr,"debug2: %i-%i\n",this_exon_pos.first,this_exon_pos.second);
          exon--;
          if(!region_is_overlap(this_exon_pos,exon->second.pos)){
            loose_exon_recorded=false;
            this->add_new_loose_exon(exon_key,this_exon_pos);
          }else{
            loose_exon_recorded=true;
            this->update_loose_exon_info(exon->first,this_exon_pos);
            loose_recorded_id=exon->second.id;
          }
        }else{
            //fprintf(stderr,"debug3: %i-%i->%i-%i\n",this_exon_pos.first,this_exon_pos.second,exon->second.pos.first,exon->second.pos.second);
          if(!region_is_overlap(this_exon_pos,exon->second.pos)){
            exon--;
              //fprintf(stderr,"debug3.1: %i-%i->%i-%i\n",this_exon_pos.first,this_exon_pos.second,exon->second.pos.first,exon->second.pos.second);
            if(!region_is_overlap(this_exon_pos,exon->second.pos)){
                //fprintf(stderr,"debug3.1.1: %i-%i->%i-%i\n",this_exon_pos.first,this_exon_pos.second,exon->second.pos.first,exon->second.pos.second);
              loose_exon_recorded=false;
              this->add_new_loose_exon(exon_key,this_exon_pos);
            }else{
                //fprintf(stderr,"debug3.1.2: %i-%i->%i-%i\n",this_exon_pos.first,this_exon_pos.second,exon->second.pos.first,exon->second.pos.second);
              loose_exon_recorded=true;
                this->update_loose_exon_info(exon->first,this_exon_pos);
              loose_recorded_id=exon->second.id;
            }
          }else{
              //has existing-exon match, attention, overlap larger than 0 of recorded exon will be considered a match, 
              //this is a very loose standard
              //better modification is needed.
              //fprintf(stderr,"debug3.1: %i-%i->%i-%i\n",this_exon_pos.first,this_exon_pos.second,exon->second.pos.first,exon->second.pos.second);
              //fprintf(stderr,"debug3.2: %i-%i->%i-%i\n",this_exon_pos.first,this_exon_pos.second,exon->second.pos.first,exon->second.pos.second);
            loose_exon_recorded=true;
            this->update_loose_exon_info(exon->first,this_exon_pos);
            loose_recorded_id=exon->second.id;
          }
        }
      }else{
        loose_exon_recorded=true;
        loose_recorded_id=this->loose_exons.find(exon_key)->second.id;
      }
        
      bool exon_recorded=false;
      int recorded_id;
      if(this->exons.find(exon_key)==this->exons.end()){
        exon_feature this_exon;
        this_exon.id=this->exons.size();
        this_exon.pos=this_exon_pos;
        this->exons[exon_key]=this_exon;
      }else{
        exon_recorded=true;
        recorded_id=this->exons.find(exon_key)->second.id;
      }
    
      string i_gene_id=regex_replace(aux_data[0],gene_regex,"");
      transcript_id_raw=get_aux_field("transcript_id",aux_data,afileds_info,"exon");
      string i_transcript_id=regex_replace(transcript_id_raw,transcript_regex,"");
      if(this->gene_id==i_gene_id&&this->transcripts.find(i_transcript_id)!=transcripts.end()){
        pair<int,int> i_pos={stoi(core_data[3]),stoi(core_data[4])};
        if(exon_recorded){
          this->transcripts[i_transcript_id].exons.push_back(recorded_id);
        }else{
          this->transcripts[i_transcript_id].exons.push_back(this->exons[exon_key].id);
        }
        if(loose_exon_recorded){
          this->transcripts[i_transcript_id].loose_exons.push_back(loose_recorded_id);
        }else{
          this->transcripts[i_transcript_id].loose_exons.push_back(this->loose_exons.size()-1);
        }
      }else{
        fprintf(stderr,"fill_in_feature_info: ERROR: fill_in_feature_info: type: exon, gene_id: %s does not match to transcript_id: %s, or transcript_id %s is not a key, check GTF\n",
                this->gene_id.c_str(),i_transcript_id.c_str(),i_transcript_id.c_str());
        exit(0);
      }
    }
  }
  
  void debugPrintInfo(bool stdout=false){
      string transcripts_string="";
      string exon_string="";
      if(transcripts.empty()){
        transcripts_string="\n";
      }else{
        for(auto feature:this->transcripts){
          string exon_str="";
          for(int exon:feature.second.exons){
            exon_str=exon_str+to_string(exon)+",";
          }
          transcripts_string=transcripts_string+"\n\\_______"+feature.second.transcript_id+'\t'+\
            feature.second.transcript_name+'\t'+to_string(feature.second.start)+"\t"+to_string(feature.second.end)+"\t"+exon_str;
        }
      }
      if(exons.empty()){
        exon_string="\n";
      }else{
        for(auto exon:exons){
          exon_string=exon_string+"\n"+gene_name+"\t"+to_string(exon.second.id)+"="+to_string(exon.second.pos.first)+"-"+to_string(exon.second.pos.second)+";";
        }
      }
      if(stdout){
        cout << gene_name.c_str() << "\t" << gene_id.c_str() << "\t" << gene_type.c_str() << "\t";
        cout << chromosome.c_str() << "\t" << start << "\t" << end << "\t" << exon_string << "\t";
        cout << transcripts_string.c_str() << endl;
      }else{
        fprintf(stderr,"%s\t%s\t%s\t%s\t%i\t%i%s\n\n",
                gene_name.c_str(),gene_id.c_str(),gene_type.c_str(),chromosome.c_str(),start,end,transcripts_string.c_str());
      }
  }
};

/*!
 * mapping_info
 * 
 * @abstract: parsed read mapping info from BAM 
 * 
 * @ID                    | read id
 * @read_length
 * @final_annotation_id   | gene id of the final annotation of read, string form of 'features', (a read may have multiple alignment to ref).
 * @final_annotation_name | gene name of the final annotation of read.
 * @read_orientation      | orientations of this read inferred from its mappings, "forward" or "reverse" or "*" 
 *                          (for read not mapped), separated by ',', this data member is the string form of 'read_orientations'.
 * @qcovs**               | query coverage of hits. 
 * @qhits**               | vector of pair<int,int>, indicate mapping start and end position on read, 1-based.
 * @hit_strands**         | orientation of hits corresponding to ref, "forward" or "reverse". 
 * @is_spliced**          | is hits spliced. 
 * @aligned_segments**    | some reads are aligned by segmentation (spliced), using vector of pair to represent them, if they are
 *                          aligned without segmentation, using its rhit as the only item in vector.   
 * @refs**                | chromosome id (in header of genome fasta). 
 * @features**            | gene ids of hits according to GTF file (a read may have multiple alignment to ref), 
 *                          notice that features.size() == qhits.size(), as some of the hits may not overlap with 
 *                          certain gene, i.e. not mapped to transcriptome, these hits would be annotated as 
 *                          "non-coding". 
 * @rcovs**               | coverage of gene (for non-coding hits, rcov=0). 
 * @rhits**               | vector of pair<int,int>, indicate mapping start and end positions on chromosome, 1-based.
 * @read_orientations**   | orientation of read inferred from its mappings, "forward" or "reverse" or "*".
 * @occupiedRegion        | regions that occupied by hit(s).
 * @freeRegion            | vector of pair<int,int>, indicate regions that not mapped to genome, 1-based.
 * @mergedRegions**       | vector of vector of int, indicate index (0-based) of qhits that were merged due to overlaping 
 *                          during calculation of occupiedRegion.
 * @mapped_to_genome
 * @mapped_to_transcriptome
 * 
 * 
 * note: the lengths of data members with suffix '**' are expected to be equal.
 */

struct mapping_info{
  string ID;
  int read_length;
  string final_annotation_id="default";
  string final_annotation_name="default";
  string final_transcript_name="defalut";
  string read_orientation="default";
  vector<float> qcovs;
  vector<pair<int, int>> qhits;
  vector<string> hit_strands;
  vector<bool> is_spliced;
  vector<vector<pair<int,int>>> aligned_segments;
  vector<string> refs;
  vector<string> features;
  vector<string> transcripts;
  vector<float> rcovs;
  vector<pair<int, int>> rhits;
  vector<int> annotation_scores;
  vector<string> read_orientations;
  vector<pair<int,int>> occupiedRegion;
  vector<pair<int,int>> freeRegion;
  vector<vector<int>> mergedRegions;
  bool mapped_to_genome;
  bool mapped_to_transcriptome;
  
  void debugPrintInfo(){
    string ref_str="";
    string strand_str="";
    for(string ref:this->refs){
      ref_str=ref_str+ref+",";
    }
    for(string strand:this->hit_strands){
      strand_str=strand_str+strand+",";
    }
    string hits_str=vector_pair_int2str(this->qhits);
    string free_region_str=vector_pair_int2str(this->freeRegion);
    string is_spliced_str=stringCat(this->is_spliced);
    fprintf(stderr,"%s\t%s\t%s\t%s\t%s\t%s\n",this->ID.c_str(),ref_str.c_str(),strand_str.c_str(),\
            hits_str.c_str(),free_region_str.c_str(),is_spliced_str.c_str());
  }//void debugPrintInfo()
  
  void writeAnnotation(ofstream& File){
    string free_region_str=vector_pair_int2str(this->freeRegion);
    string qhits_str=vector_pair_int2str(this->qhits);
    string qcov_str=stringCat(this->qcovs);
    string merged_hits_str="";
    string occupied_region_str=vector_pair_int2str(this->occupiedRegion);
    string rhits_str=vector_pair_int2str(this->rhits);
    string rcov_str=stringCat(this->rcovs);
    string is_spliced_str=stringCat(this->is_spliced);
    string hit_strand_str=stringCat(this->hit_strands);
    string annotation_score_str=stringCat(this->annotation_scores);
    if(this->mergedRegions.empty()){
      merged_hits_str="*";
    }else{
      for(vector<int> indices:this->mergedRegions){
        merged_hits_str=stringCat(indices)+";";
      }
      merged_hits_str.erase(merged_hits_str.end()-1);
    }
    string aligned_segments="";
    for(vector<pair<int,int>> hit:this->aligned_segments){
      aligned_segments=aligned_segments+regions2string(hit,',')+";";
    }
    aligned_segments.erase(aligned_segments.end()-1);
    File << this->ID << "\t" << this->final_annotation_id << "\t" << this->final_annotation_name << "\t";
    File << this->read_orientation << "\t" << free_region_str << "\t" << qhits_str << "\t" << qcov_str << "\t";
    File << rhits_str << "\t" << rcov_str << "\t" <<  merged_hits_str << "\t";
    File << occupied_region_str << "\t" << is_spliced_str << "\t" << aligned_segments << "\t" << final_transcript_name << "\t";
    File << hit_strand_str << "\t" << annotation_score_str << "\n";
  }
  
  
  /*!
   * @abstract some genes share same start or end position (a rare case in overlapped genes), this function is 
   *           used to choose the best one as the annotation with a specific read with its hit position on reference.
   * 
   */
  int choose_from_overlapping_feature(pair<int,int>& rhit, vector<seq_feature> overlapping_features, bool debug=false){
    if(overlapping_features.size()==0){
      return(0);
    }
    int min_overhang_index;
    float min_overhang;
    for(int i=0;i<overlapping_features.size();i++){
      pair<int,int> feat_pos={overlapping_features[i].start,overlapping_features[i].end};
      float overhang=abs(calculate_cov(feat_pos,rhit)-1);
      if(debug){
        fprintf(stderr,"choose_from_overlapping_feature: %i\t%s\t%f\n",i,overlapping_features[i].id.c_str(),overhang);
      }
      if(i==0){
        min_overhang_index=i;
        min_overhang=overhang;
      }else{
        if(overhang<min_overhang){
          min_overhang_index=i;
        }else{
          continue;
        }
      }
    }
    return(min_overhang_index);
  }
  
  
  /*!
   * @abstract annotate read with gene features
   * 
   * @param sorted_gene_features: .first chromosome
   *                              .second <map>, end position on its genome, sorted key
   *                              .second[i].first feature start
   *                              .second[i].second gene_id
   * @param gene_features: .first gene_id
   *                       .second mapping_info
   * @param debug
   * 
   */
  
  void annotateRead(unordered_map<string, map<int,vector<seq_feature>>>& sorted_gene_features, \
                    unordered_map<string, gene_feature>& gene_features, int tech, bool debug=false){
    if(debug){
      fprintf(stderr,"annotateRead: debug1: %s: start annotation\n",this->ID.c_str());
    }
    if(this->qhits.empty()){
      if(debug){
        fprintf(stderr,"annotateRead: debug1: %s: unmapped\n",this->ID.c_str());
      }
      this->mapped_to_genome=false;
      this->mapped_to_transcriptome=false;
      this->final_annotation_id="*";
      this->final_annotation_name="*";
      this->read_orientation="*";
      this->annotation_scores.push_back(9999);
      return;
    }else{
      if(debug){
        fprintf(stderr,"annotateRead: debug1: %s: mapped to genome\n",this->ID.c_str());
      }
      this->mapped_to_genome=true;
    }
    
    if(debug){
      fprintf(stderr,"annotateRead: debug1.1: %s: initiating supplementary hits\n",this->ID.c_str());
    }
    vector<float> sup_qcovs;
    vector<pair<int, int>> sup_qhits;
    vector<string> sup_hit_strands;
    vector<bool> sup_is_spliced;
    vector<vector<pair<int,int>>> sup_aligned_segments;
    vector<string> sup_refs;
    vector<string> sup_features;
    vector<float> sup_rcovs;
    vector<pair<int,int>> sup_rhits;
    vector<int> sup_annotation_scores;
    vector<string> sup_read_orientations;
    
    int sup_hit_index=this->rhits.size()+1;
    for(int hit_index=0;hit_index<this->rhits.size();hit_index++){
      int qstart=this->qhits[hit_index].first;
      int qend=this->qhits[hit_index].second;
      float qcov=static_cast<float>(qend-qstart+1)/this->read_length;
      int query_key1=this->rhits[hit_index].first;
      int query_key2=this->rhits[hit_index].second;
      auto candidate_feature1=sorted_gene_features[this->refs[hit_index]].lower_bound(query_key1);
      auto candidate_feature2=sorted_gene_features[this->refs[hit_index]].lower_bound(query_key2);
      pair<int,int> rhit={this->rhits[hit_index].first,this->rhits[hit_index].second};
      vector<seq_feature> candidate_feature_set;
      
      if(debug){
        fprintf(stderr,"annotateRead: debug2: hit: %i/%lu, %i-%i\n",hit_index+1,this->rhits.size(),rhit.first,rhit.second);
      }

      if(candidate_feature1==sorted_gene_features[this->refs[hit_index]].end()){
        if(debug){
          fprintf(stderr,"annotateRead: debug2.1: %s: non-coding\n",this->ID.c_str());
        }
        this->features.push_back("*");
        this->read_orientations.push_back("*");
        this->qcovs.push_back(qcov);
        this->rcovs.push_back(0);
        this->annotation_scores.push_back(9999);
        continue;
      }else{
        if(debug){
          fprintf(stderr,"annotateRead: debug2.2: %s: adding candidate1 features\n",this->ID.c_str());
        }
        int candidate1_shift_time=0;
        while(candidate1_shift_time<20&&candidate_feature1!=sorted_gene_features[this->refs[hit_index]].end()){
          //candidate_feature_set.push_back(candidate_feature1->second[choose_from_overlapping_feature(rhit,candidate_feature1->second,debug)]);
          candidate_feature_set.insert(candidate_feature_set.end(),candidate_feature1->second.begin(),candidate_feature1->second.end());
          candidate1_shift_time++;
          candidate_feature1++;
        }
      }
      
      if(debug){
        fprintf(stderr,"annotateRead: debug2.3: %s: adding candidate2 features\n",this->ID.c_str());
      }
      if(candidate_feature2!=sorted_gene_features[this->refs[hit_index]].end()){
        candidate_feature_set.push_back(candidate_feature2->second[choose_from_overlapping_feature(rhit,candidate_feature2->second,debug)]);
      }
      int candidate2_shift_time=0;
      while(candidate2_shift_time<20&&candidate_feature2!=sorted_gene_features[this->refs[hit_index]].begin()){
          candidate_feature2--;
          //candidate_feature_set.push_back(candidate_feature2->second[choose_from_overlapping_feature(rhit,candidate_feature2->second,debug)]);
          candidate_feature_set.insert(candidate_feature_set.end(),candidate_feature2->second.begin(),candidate_feature2->second.end());
          candidate2_shift_time++;
      }
      
      if(debug){
        fprintf(stderr,"annotateRead: debug2.4: %s: iterating candidate features\n",this->ID.c_str());
      }
      //int min_overhang=1000000;
      //float max_rcov=0;
      //string best_hit_id;
      //seq_feature best_feature;
      bool has_at_least_one_hit=false;
      vector<string> valid_hit_set;
      for(seq_feature candidate_feature:candidate_feature_set){
        if(debug){
          fprintf(stderr,"annotateRead: debug2.4.1: %s: %s strand:%i position: %i-%i\n",this->ID.c_str(),
                  candidate_feature.id.c_str(),candidate_feature.strand,candidate_feature.start,candidate_feature.end);
        }
        int feature_strand=candidate_feature.strand;
        pair<int,int> feature_hit={candidate_feature.start,candidate_feature.end};
        float feature_cov=calculate_cov(feature_hit,rhit);
        float rhit_cov=calculate_cov(rhit,feature_hit);
        int rhit_overhang=calculate_min_overhang(rhit,feature_hit,feature_strand,tech);
        if(debug){
          fprintf(stderr,"annotateRead: debug2.4.2: %s: %s: %f\t%i\n",this->ID.c_str(),
                  candidate_feature.id.c_str(),feature_cov,rhit_overhang);
        }
      
        if(feature_cov==0){
          continue;
        }
        
        if(std::count(valid_hit_set.begin(),valid_hit_set.end(),candidate_feature.id)){
          continue;
        }
        
        if(!has_at_least_one_hit){
          if(debug){
            fprintf(stderr,"annotateRead: debug2.4.3: found one hit\n");
          }
          valid_hit_set.push_back(candidate_feature.id);
          this->features.push_back(candidate_feature.id);
          string feature_strand=num2rna[gene_features[candidate_feature.id].strand];
          string hit_strand=this->hit_strands[hit_index];
          if(feature_strand==hit_strand){
            this->read_orientations.push_back("forward");
          }else{
            this->read_orientations.push_back("reverse");
          }
          this->qcovs.push_back(qcov);
          this->rcovs.push_back(feature_cov);
          this->annotation_scores.push_back(floor(rhit_overhang/10));
          has_at_least_one_hit=true;
        }else{
          if(debug){
            fprintf(stderr,"annotateRead: debug2.4.3: found one supplementary hit\n");
          }
          valid_hit_set.push_back(candidate_feature.id);
          sup_qcovs.push_back(qcov);
          sup_qhits.push_back(this->qhits[hit_index]);
          string feature_strand=num2rna[gene_features[candidate_feature.id].strand];
          string hit_strand=this->hit_strands[hit_index];
          sup_hit_strands.push_back(hit_strand);
          if(feature_strand==hit_strand){
            sup_read_orientations.push_back("forward");
          }else{
            sup_read_orientations.push_back("reverse");
          }
          sup_is_spliced.push_back(this->is_spliced[hit_index]);
          sup_refs.push_back(this->refs[hit_index]);
          sup_features.push_back(candidate_feature.id);
          sup_rcovs.push_back(feature_cov);
          sup_rhits.push_back(this->rhits[hit_index]);
          for(int merged_region_index=0;merged_region_index<this->mergedRegions.size();merged_region_index++){
            if(std::count(this->mergedRegions[merged_region_index].begin(),this->mergedRegions[merged_region_index].end(),hit_index)){
              this->mergedRegions[merged_region_index].push_back(sup_hit_index);
              sup_hit_index++;
            }else{
              continue;
            }
          }//for(int merged_region_index=0;merged_region_index<this->mergedRegions.size();merged_region_index++)
          sup_annotation_scores.push_back(floor(rhit_overhang/10));
          sup_aligned_segments.push_back(this->aligned_segments[hit_index]);
        }//if(!has_at_least_one_hit)
      }//for(seq_feature candidate_feature:candidate_feature_set)
      
      //if(max_rcov==0){
      if(!has_at_least_one_hit){
        if(debug){
          fprintf(stderr,"annotateRead: debug2.5: %s: no feature for this hit\n",this->ID.c_str());
        }
        this->features.push_back("*");
        this->read_orientations.push_back("*");
        this->qcovs.push_back(qcov);
        this->rcovs.push_back(0);
        this->annotation_scores.push_back(9999);
        continue;
      }else{
        continue;
      }
      
      //this->features.push_back(best_hit_id);
      //string feature_strand=num2rna[gene_features[best_hit_id].strand];
      //string hit_strand=this->hit_strands[hit_index];
      //this->qcovs.push_back(qcov);
      //this->rcovs.push_back(max_rcov);
      //this->annotation_scores.push_back(floor(min_overhang/10));
      //if(feature_strand==hit_strand){
      //  this->read_orientations.push_back("forward");
      //}else{
      //  this->read_orientations.push_back("reverse");
      //}
      //continue;
    } //for(int hit_index=0;hit_index<this->rhits.size();hit_index++)
    
    if(debug){
      fprintf(stderr,"annotateRead: debug3.0: %s: %li primary genes have hit\n",this->ID.c_str(),this->features.size());
    }
    this->qcovs.insert(this->qcovs.end(),sup_qcovs.begin(),sup_qcovs.end());
    this->qhits.insert(this->qhits.end(),sup_qhits.begin(),sup_qhits.end());
    this->hit_strands.insert(this->hit_strands.end(),sup_hit_strands.begin(),sup_hit_strands.end());
    this->is_spliced.insert(this->is_spliced.end(),sup_is_spliced.begin(),sup_is_spliced.end());
    this->refs.insert(this->refs.end(),sup_refs.begin(),sup_refs.end());
    this->features.insert(this->features.end(),sup_features.begin(),sup_features.end());
    this->rcovs.insert(this->rcovs.end(),sup_rcovs.begin(),sup_rcovs.end());
    this->rhits.insert(this->rhits.end(),sup_rhits.begin(),sup_rhits.end());
    this->read_orientations.insert(this->read_orientations.end(),sup_read_orientations.begin(),sup_read_orientations.end());
    this->aligned_segments.insert(this->aligned_segments.end(),sup_aligned_segments.begin(),sup_aligned_segments.end());
    this->annotation_scores.insert(this->annotation_scores.end(),sup_annotation_scores.begin(),sup_annotation_scores.end());
    if(debug){
      fprintf(stderr,"annotateRead: debug3.1: %s: %li genes have hit\n",this->ID.c_str(),this->features.size());
    }
    
    this->identifyIsoforms(gene_features,false,debug);
    if(debug){
      fprintf(stderr,"annotateRead: debug3.2: %s: isoforms identified\n",this->ID.c_str());
    }
    this->final_annotation_id=stringCat(this->features);
    this->final_annotation_name="";
    for(string id:this->features){
      if(gene_features[id].gene_name.empty()){
        this->final_annotation_name=this->final_annotation_name+"*"+",";
        continue;
      }
      this->final_annotation_name=this->final_annotation_name+gene_features[id].gene_name+",";
    }
    this->final_annotation_name.erase(this->final_annotation_name.end()-1);
    this->final_transcript_name=stringCat(this->transcripts);
    this->read_orientation=stringCat(this->read_orientations);
  } //annotateRead
  
  /*! 
   *  @abstract identify possible isoforms of each hit
   *  
   *  @param gene_features: gene feature map, created from function: readFeaturesFromGTF.
   *  
   *  @return modificate this->transcripts, format: gene/transcript(#TAG=exon1>exon2...exonN),
   *          the order of exons is 5' -> 3'. #TAG: see 'glossary'.
   *  
   *  @glossary unspliced: no splicing behavior was observed on read, and the hit region 
   *                       exceed known exonic regions or covers non-exonic regions of the gene.
   *                       
   *            mono_exonic: no splicing behavior was observed on read, but the hit region 
   *                         is in the range of certain known exon of a multi-exon gene.
   *                         
   *            mono-exon_match: no splicing behavior was observed on read, but the hit 
   *                             region is in the range of a mono-exon gene.
   *                             
   *            fsm/ism: full spliced match/incomplete spliced match, splicing behavior 
   *                     was observed on read, and the exons were uniquely matched to known 
   *                     transcript, though it may be truncated(ism).
   *                     
   *            ambiguous: splicing behavior was observed on read, but the exons were not
   *                       uniquely matched to known transcripts.
   *                       
   *            novel: splicing behavior was observed on read, but the arrange of exons were not 
   *                   matched to any known transcripts. notice that in some case the exon may be 
   *                   -1, which indicate the exon is not matched to any known exon (the exon is 
   *                   not recorded in GTF file).
   *                   
   */
  void identifyIsoforms(unordered_map<string, gene_feature>& gene_features, bool loose_exon=false, bool debug=false){
    for(int i=0;i<this->features.size();i++){
      this->transcripts.push_back("*");
    }
    int index=0;
    for(string feature:this->features){
      float exon_determine_threshold=0.9;
      if(loose_exon){
        exon_determine_threshold=0;
      }
      string gene_name=gene_features[feature].gene_name;
      map<int, exon_feature>& exon_table_used=gene_features[feature].exons;
      if(loose_exon){
        exon_table_used=gene_features[feature].loose_exons;
      }
      if(gene_name=="*"){
        gene_name=feature;
      }
      //not mapped or mapped to non-transcript region
      if(feature=="*"){
        //this->transcripts.push_back("*");
        index++;
        continue;
      }
      if(debug){
        fprintf(stderr,"identifyIsoforms: debug1: %s: inital trial using strict exon match\n",this->ID.c_str());
      }
      this->M_identifyIsoform(index,exon_table_used,exon_determine_threshold,gene_features,feature,gene_name,loose_exon,debug);
      if((!loose_exon && this->transcripts[index].find("novel")!=std::string::npos)||
         (this->transcripts[index].find(">")==std::string::npos && this->transcripts[index].find("unspliced")==std::string::npos)){
        if(debug){
          fprintf(stderr,"identifyIsoforms: debug2: %s: secondary trial using loose exon match\n",this->ID.c_str());
        }
        this->M_identifyIsoform(index,gene_features[feature].loose_exons,0,gene_features,feature,gene_name,true,debug);
        this->transcripts[index]=this->transcripts[index]+'*';
      }
      index++;
    } //for(string feature:this->features)
  } //identifyIsoforms
  
  /*!
   * @abstract module for isoform identification
   * 
   * 
   */
  void M_identifyIsoform(int index, map<int, exon_feature>& exon_table_used, float exon_determine_threshold,
                         unordered_map<string, gene_feature>& gene_features, string feature, string gene_name,
                         bool loose_exon, bool debug){
    //mapped to transcript region and spliced
    if(this->is_spliced[index]){
      vector<int> exons;
      bool has_novel_exon=false;
      for(pair<int,int> segment:this->aligned_segments[index]){
        int key=floor((segment.first+segment.second)/2);
        auto candidate=exon_table_used.lower_bound(key);
        if(candidate==exon_table_used.end()||calculate_cov(segment,candidate->second.pos)<=exon_determine_threshold){
          if(candidate!=exon_table_used.begin()){
            candidate--;
            //fprintf(stderr,"debug2: %i->%i-%i\n",key,candidate->second.pos.first,candidate->second.pos.second);
          }else{
            exons.push_back(-1);
            has_novel_exon=true;
            continue;
          }
        }
        if(calculate_cov(segment,candidate->second.pos)>exon_determine_threshold&&calculate_cov(candidate->second.pos,segment)>exon_determine_threshold){
          //fprintf(stderr,"debug3: %i->%i-%i\n",key,candidate->second.pos.first,candidate->second.pos.second);
          exons.push_back(candidate->second.id);
          continue;
        }else{
          exons.push_back(-1);
          has_novel_exon=true;
          continue;
        }
      } //for(pair<int,int> segment:this->aligned_segments[index])
      if(gene_features[feature].strand==3){
        reverse(exons.begin(),exons.end());
      }
      string arrange=stringCat(exons,'>',false);
      if(has_novel_exon){
        this->transcripts[index]=gene_name+"(novel="+arrange+")";
      }else{
        int match=0;
        string latest_transcript_name;
        if(debug){
          fprintf(stderr,"exon arrange: %s\n",arrange.c_str());
        }
        for(auto transcript:gene_features[feature].transcripts){
          vector<int>& transcript_exons_used=transcript.second.exons;
          if(debug){
            fprintf(stderr,"M_identifyIsoform: debug1: comparing %s: %s\n",
                    transcript.second.transcript_name.c_str(),
                    stringCat(transcript_exons_used,'>',false).c_str());
          }
          if(loose_exon){
            transcript_exons_used=transcript.second.loose_exons;
          }
          //if(exon_is_contained(transcript_exons_used,exons)){
          if(exon_is_contained(stringCat(transcript_exons_used,'>',false),arrange)){
            if(debug){
              fprintf(stderr,"M_identifyIsoform: debug2: hit found %s\n",
                      transcript.second.transcript_name.c_str());
            }
            match++;
            latest_transcript_name=transcript.second.transcript_name;
            if(latest_transcript_name.empty()){
              latest_transcript_name=transcript.second.transcript_id;
            }
            continue;
          }else{
            continue;
          }
        } //for(auto transcript:gene_features[feature].transcripts)
        if(match==0){
          this->transcripts[index]=gene_name+"(novel="+arrange+")";
        }else if(match==1){
          this->transcripts[index]=latest_transcript_name+"(fsm/ism="+arrange+")";
        }else{
          this->transcripts[index]=gene_name+"(ambiguous="+arrange+")";
        }
      }//if(has_novel_exon)
      
      ////mapped to transcript region but no spliced mark was observed, i.e., no 'N's was observed in cigar, 
      //mono exon match, mono_exonic, and unspliced transcript is included
    }else{
      if(exon_table_used.size()==1){
        auto the_only_transcript=gene_features[feature].transcripts.begin();
        string the_only_transcript_name=the_only_transcript->second.transcript_name;
        this->transcripts[index]=the_only_transcript_name+"(mono-exon_match)";
      }else{
        int key=floor((this->aligned_segments[index][0].first+this->aligned_segments[index][0].second)/2);
        auto candidate=exon_table_used.lower_bound(key);
        if(candidate==exon_table_used.end()||calculate_cov(this->aligned_segments[index][0],candidate->second.pos)<=exon_determine_threshold){
          if(candidate!=exon_table_used.begin()){
            candidate--;
          }else{
            this->transcripts[index]=gene_features[feature].gene_name+"(unspliced)";
          }
        }
        if(calculate_cov(this->aligned_segments[index][0],candidate->second.pos)>exon_determine_threshold){
          int match=0;
          string latest_transcript_name;
          for(auto transcript:gene_features[feature].transcripts){
            vector<int>& transcript_exons_used=transcript.second.exons;
            if(loose_exon){
              transcript_exons_used=transcript.second.loose_exons;
            }
            if(exon_is_contained(transcript_exons_used,candidate->second.id)){
              match++;
              latest_transcript_name=transcript.second.transcript_name;
              if(latest_transcript_name.empty()){
                latest_transcript_name=transcript.second.transcript_id;
              }
              continue;
            }else{
              continue;
            }
          } //for(auto transcript:gene_features[feature].transcripts)
          //fprintf(stderr,"%s\t%i\n",latest_transcript_name.c_str(),match);
          //in practice this condition not exist
          if(match==0){
            this->transcripts[index]=gene_name+"(novel="+to_string(candidate->second.id)+")";
            //mono_exonic
          }else if(match==1){
            this->transcripts[index]=latest_transcript_name+"(mono_exonic="+to_string(candidate->second.id)+")";
            //the only exon matched to multiple transcripts
          }else{
            this->transcripts[index]=gene_name+"(ambiguous="+to_string(candidate->second.id)+")";
          }
        }else{
          this->transcripts[index]=gene_features[feature].gene_name+"(unspliced)";
        }
      }
    }
  } //void M_identifyIsoform
  
  void findFreeRegion(){
    if(this->qhits.size()==0){
      this->freeRegion.push_back({1,this->read_length});
    }else if(this->qhits.size()==1){
      this->occupiedRegion=this->qhits;
      this->freeRegion=calculateFreeRegions(this->qhits);
    }else{
      vector<pair<int,int>> occupied_regions=mergeRegions(this->qhits);
      this->occupiedRegion=occupied_regions;
      this->freeRegion=calculateFreeRegions(occupied_regions);
    }
  }
  
  /*! 
   *  @abstract merge several regions that may overlap, if certain region does not overlap with any region,
   *            it will be returned as raw, else it will be merged with other regions, the merged region starts
   *            at the left-most position of overlapped regions and ends at the right-most position.
   *  
   *  @param raw_regions: vector of raw regions, .first=start, .second=end
   *  
   *  @return vector of merged regions, .first=merged start, .second=merged end
   */
  private: vector<pair<int, int>> mergeRegions(vector<pair<int,int>>& raw_regions){
    vector<pair<int, int>> res;
    map<int,vector<int>> ends_ordered;
    vector<int> regions_merged;
    int i=0;
    for(pair<int,int> raw_region:raw_regions){
      ends_ordered[raw_regions[i].first].push_back(i);
      ends_ordered[raw_regions[i].second].push_back(i);
      i++;
    }
    unordered_map<int,bool> recorded_region;
    pair<int,int> merged_region={-1,-1};;
    for(auto end:ends_ordered){
      if(recorded_region.empty()){
        for(int this_region:end.second){
          recorded_region[this_region]=true;
          regions_merged.push_back(this_region);
        }
        merged_region.first=end.first;
      }else{
        for(int this_region:end.second){
          if(recorded_region.find(this_region)==recorded_region.end()){
            recorded_region[this_region]=true;
            regions_merged.push_back(this_region);
          }else{
            recorded_region.erase(this_region);
          }
        }
        if(recorded_region.empty()){
          merged_region.second=end.first;
          res.push_back(merged_region);
          if(regions_merged.size()>1){
            this->mergedRegions.push_back(regions_merged);
          }
          regions_merged.clear();
          merged_region.first=-1;
          merged_region.second=-1;
        }
      } //else
    } //for(auto end:ends_ordered)
    return(res);
  }
    
    
  /*!
   * @abstract for a read that have non-overlapped hit(s), return region(s) that are not occupied by hit region, notice occupied_regions
   *           should be ordered by its start position.
   *           
   * @param occupied_regions: regions ordered by its start position on read.
   * 
   * @return vector of regions not occupied by hits (free region).
   * 
   */
  
  private: vector<pair<int, int>> calculateFreeRegions(vector<pair<int,int>>& occupied_regions){
    vector<pair<int, int>> res;
    if(occupied_regions.size()==1){
      if(occupied_regions[0].first==1&&occupied_regions[0].second==this->read_length){
        res.push_back({-1,-1});
        return(res);
      }else if(occupied_regions[0].first!=1&&occupied_regions[0].second==this->read_length){
        res.push_back({1,occupied_regions[0].first-1});
        return(res);
      }else if(occupied_regions[0].first==1&&occupied_regions[0].second!=this->read_length){
        res.push_back({occupied_regions[0].second+1,this->read_length});
        return(res);
      }else{
        res.push_back({1,occupied_regions[0].first-1});
        res.push_back({occupied_regions[0].second+1,this->read_length});
        return(res);
      }
    }else{
      for(int i=0;i+1<occupied_regions.size();i++){
        res.push_back({occupied_regions[i].second+1,occupied_regions[i+1].first-1});
      }
      if(occupied_regions[0].first==1&&occupied_regions[occupied_regions.size()-1].second==this->read_length){
        sleep(0);
      }else if(occupied_regions[0].first!=1&&occupied_regions[occupied_regions.size()-1].second==this->read_length){
        res.insert(res.begin(),{1,occupied_regions[0].first-1});
      }else if(occupied_regions[0].first==1&&occupied_regions[occupied_regions.size()-1].second!=this->read_length){
        res.push_back({occupied_regions[occupied_regions.size()-1].second+1,this->read_length});
      }else{
        res.insert(res.begin(),{1,occupied_regions[0].first-1});
        res.push_back({occupied_regions[occupied_regions.size()-1].second+1,this->read_length});
      }
      return(res);
    }
  }
    
  private: string vector_pair_int2str(vector<pair<int,int>>& vect){
    string res="";
    for(pair<int,int> i_pair:vect){
      res=res+to_string(i_pair.first)+"-"+to_string(i_pair.second)+";";
    }
    return(res);
  }
  
};

/*!
 * @abstract read gene, transcript, exon data from GTF file
 * 
 * @param GTFfilename
 * @param gene_features: empty map, gene_id-gene_feature, data will be filled in while reading and processing GTF file
 */
void readFeaturesFromGTF(const char * GTFfilename, unordered_map<string,
                         gene_feature>& gene_features, bool loose_transcript=false){
  fprintf(stderr,"Parsing GTF: %s\n",GTFfilename);
  ifstream gtffile;
  gtffile.open(GTFfilename);
  int index=0;
  vector<string> head_gene_ids;
  string last_gene_id;
  unordered_map<string,int> afileds_info;
  afileds_info["gene_gene_type_field"]=-2;
  afileds_info["gene_gene_biotype_field"]=-2;
  afileds_info["gene_gene_name_field"]=-2;
  afileds_info["transcript_transcript_id_field"]=-2;
  afileds_info["transcript_transcript_name_field"]=-2;
  afileds_info["exon_transcript_id_field"]=-2;

  while(!gtffile.eof()){
    string line;
    string aux_data;
    getline(gtffile,line,'\n');
    if(line==""||line.find("#")==0){
      continue;
    }
    vector<string> core_fields;
    vector<string> aux_fields;
    core_fields=tokenize(line,'\t');
    if(core_fields[2]=="Selenocysteine"||core_fields[2]=="start_codon"||core_fields[2]=="stop_codon"||core_fields[2]=="UTR"){
      continue;
    }
    aux_data=core_fields[8];
    aux_fields=tokenize(aux_data,';');
    if(head_gene_ids.empty()){
      head_gene_ids.push_back(regex_replace(aux_fields[0],gene_regex,""));
    }
    if(core_fields[2]=="gene"){
      string i_gene_id=regex_replace(aux_fields[0],gene_regex,"");
      if(gene_features.find(i_gene_id)==gene_features.end()){
        gene_feature feature;
        index++;
        feature.index=index;
        feature.fill_in_feature_info(core_fields,aux_fields,afileds_info);
        gene_features[feature.gene_id]=feature;
      }else{
        gene_features[i_gene_id].fill_in_feature_info(core_fields,aux_fields,afileds_info);
      }
    }else if(core_fields[2]=="transcript"){
      if(aux_data.find("transcript_support_level \"1\"")==string::npos&&!loose_transcript){
        continue;
      }
      string i_gene_id=regex_replace(aux_fields[0],gene_regex,"");
      if(gene_features.find(i_gene_id)==gene_features.end()){
        gene_feature feature;
        index++;
        feature.index=index;
        feature.gene_id=i_gene_id;
        feature.fill_in_feature_info(core_fields,aux_fields,afileds_info);
        gene_features[i_gene_id]=feature;
      }else{
        gene_features[i_gene_id].fill_in_feature_info(core_fields,aux_fields,afileds_info);
      }
    }else if(core_fields[2]=="exon"){
      if(aux_data.find("transcript_support_level \"1\"")==string::npos&&!loose_transcript){
        continue;
      }
      string i_transcript_id_raw=get_aux_field("transcript_id",aux_fields,afileds_info,"exon");
      string i_gene_id=regex_replace(aux_fields[0],gene_regex,"");
      string i_transcript_id=regex_replace(i_transcript_id_raw,transcript_regex,"");
      if(gene_features.find(i_gene_id)==gene_features.end()){
        gene_feature feature;
        transcript_feature t_feature;
        index++;
        feature.index=index;
        feature.gene_id=i_gene_id;
        t_feature.transcript_id=i_transcript_id;
        feature.transcripts[i_transcript_id]=t_feature;
        feature.fill_in_feature_info(core_fields,aux_fields,afileds_info);
        gene_features[i_gene_id]=feature;
      }else if(gene_features.find(i_gene_id)!=gene_features.end()&&gene_features[i_gene_id].transcripts.find(i_transcript_id)==gene_features[i_gene_id].transcripts.end()){
        transcript_feature t_feature;
        t_feature.transcript_id=i_transcript_id;
        gene_features[i_gene_id].transcripts[i_transcript_id]=t_feature;
        gene_features[i_gene_id].fill_in_feature_info(core_fields,aux_fields,afileds_info);
      }else{
        gene_features[i_gene_id].fill_in_feature_info(core_fields,aux_fields,afileds_info);
      }
    }else{
      continue;
    }
    if(index<=5&&head_gene_ids[head_gene_ids.size()-1]!=regex_replace(aux_fields[0],gene_regex,"")){
      head_gene_ids.push_back(regex_replace(aux_fields[0],gene_regex,""));
    }
    last_gene_id=regex_replace(aux_fields[0],gene_regex,"");
  }
  for(string gene_id:head_gene_ids){
    gene_features[gene_id].debugPrintInfo();
  }
  fprintf(stderr,"......\n");
  gene_features[last_gene_id].debugPrintInfo();
  fprintf(stderr,"------------------------------------\nAll features loaded, totally %i genes.\n",index);
}

/*!
 * @abstract sort features by end position in each chromosome
 * 
 * @param gene_features: reference (c++) of raw data from GTF filled in by readFeaturesFromGTF().
 * @param sorted_features: reference (c++) of sorted data, sorting are performed within 
 *        chromosome (unordered_map.first) by end position (map.first). map.second.first = start position,
 *        map.second.second = gene_name.
 */
void sortFeaturesByCoord(unordered_map<string, gene_feature>& gene_features, unordered_map<string, map<int, vector<seq_feature>>>& sorted_features){
  fprintf(stderr,"Start sorting features\n");
  for(auto feature:gene_features){
    int key=feature.second.end;
    seq_feature i_feature;
    i_feature.start=feature.second.start;
    i_feature.end=feature.second.end;
    i_feature.strand=feature.second.strand;
    i_feature.id=feature.second.gene_id;
    if(sorted_features[feature.second.chromosome].empty()){
      map<int, vector<seq_feature>> i_chr;
      sorted_features[feature.second.chromosome]=i_chr;
      sorted_features[feature.second.chromosome][key].push_back(i_feature);
    }else{
      sorted_features[feature.second.chromosome][key].push_back(i_feature);
    }
  }
  fprintf(stderr,"Finished sorting features\n");
}

vector<pair<int,int>> get_segments_to_ref(uint32_t* CIGAR, int ncigar, int rstart){
  vector<pair<int,int>> res;
  int num_seg=0;
  int seg_start=rstart;
  for(int i=1;i<=ncigar;i++){
    if(bam_cigar_opchr(CIGAR[i-1])=='N'){
      pair<int,int> this_seg;
      this_seg.first=seg_start;
      this_seg.second=rstart+bam_cigar2rlen(i-1,CIGAR)-1;
      seg_start=rstart+bam_cigar2rlen(i,CIGAR);
      res.push_back(this_seg);
      num_seg++;
    }
    if(i==ncigar&&num_seg!=0){
      res.push_back({seg_start,rstart+bam_cigar2rlen(ncigar,CIGAR)-1});
    }
  } //for(int i=0;i<ncigar;i++)
  return(res);
}

void readMappingFile(const char * Bamfilename, unordered_map<string, mapping_info>& ReadMapping, unordered_map<string, map<int,vector<seq_feature>>>& sorted_gene_features, \
                     unordered_map<string, gene_feature>& gene_features, int tech, bool debug, ofstream& File){
  fprintf(stderr,"start to read BAM file\n");
  int r=0;
  unsigned long int num_hits=0;
  unsigned long int num_read=0;
  string latest_read;
  htsFile *testfile=hts_open(Bamfilename,"r");
  bam1_t *b = NULL;
  bam_hdr_t *h = NULL;
  h = sam_hdr_read(testfile);
  if (h == NULL) {
    cout << "Cound not read header" << endl;
    exit(0);
  }
  b=bam_init1();
  while ((r = sam_read1(testfile, h, b)) >= 0) {
    num_hits++;
    mapping_info thismapping;
    string flag;
    if(bam_is_umapped(b)){
      num_hits--;
      num_read++;
      if(num_read%5000000==0){
        fprintf(stderr,"processed %li reads\n",num_read);
      }
      continue;
    }
    if(bam_is_rev(b)){
      flag="reverse";
    }else{
      flag="forward";
    }
    uint32_t* CIGAR=bam_get_cigar(b);
    int q_start;
    int q_end;
    vector<pair<int,int>> aligned_segments_to_ref;
    int ncigar=b->core.n_cigar;
    char left_copchr=bam_cigar_opchr(CIGAR[0]);
    int left_coplen=bam_cigar_oplen(CIGAR[0]);
    char right_copchr=bam_cigar_opchr(CIGAR[ncigar-1]);
    int right_coplen=bam_cigar_oplen(CIGAR[ncigar-1]);
    int qlen=bam_cigar2qlen(ncigar,CIGAR);
    int rlen=bam_cigar2rlen(ncigar,CIGAR);
    bool splicing_stat=isSpliced(CIGAR,ncigar);
    if(splicing_stat==true){
      aligned_segments_to_ref=get_segments_to_ref(CIGAR,ncigar,b->core.pos+1);
    }
    if(bam_is_rev(b)){
      std::swap(left_coplen,right_coplen);
      std::swap(left_copchr,right_copchr);
    }
    const char * ref_name=sam_hdr_tid2name(h,b->core.tid);
    if(left_copchr=='H'){
      qlen=qlen+left_coplen;
    }
    if(right_copchr=='H'){
      qlen=qlen+right_coplen;
    }
    if(left_copchr=='H'||left_copchr=='S'){
      q_start=left_coplen+1;
    }else{
      q_start=1;
    }
    if(right_copchr=='H'||right_copchr=='S'){
      q_end=qlen-right_coplen;
    }else{
      q_end=qlen;
    }
    //if read_id not exists, i.e. this read have not been recorded in "ReadMapping" structure
    if(ReadMapping.find(bam_get_qname(b))==ReadMapping.end()){
      thismapping.ID=bam_get_qname(b);
      thismapping.hit_strands.push_back(flag);
      pair<int, int> qhit={q_start,q_end};
      thismapping.qhits.push_back(qhit);
      thismapping.refs.push_back(ref_name);
      pair<int, int> rhit={b->core.pos+1,b->core.pos+rlen};
      if(splicing_stat==false){
        aligned_segments_to_ref.push_back(rhit);
      }
      thismapping.aligned_segments.push_back(aligned_segments_to_ref);
      thismapping.rhits.push_back(rhit);
      thismapping.is_spliced.push_back(splicing_stat);
      thismapping.read_length=qlen;
      ReadMapping[bam_get_qname(b)]=thismapping;
      num_read++;
      if(latest_read.length()!=0){
        ReadMapping[latest_read].findFreeRegion();
        ReadMapping[latest_read].annotateRead(sorted_gene_features,gene_features,tech,debug);
        ReadMapping[latest_read].writeAnnotation(File);
        ReadMapping.erase(latest_read);
      }
      latest_read=thismapping.ID;
      if(num_read%5000000==0){
        fprintf(stderr,"processed %li reads\n",num_read);
      }
    }else{
      ReadMapping[bam_get_qname(b)].hit_strands.push_back(flag);
      pair<int, int> qhit={q_start,q_end};
      ReadMapping[bam_get_qname(b)].qhits.push_back(qhit);
      ReadMapping[bam_get_qname(b)].refs.push_back(ref_name);
      pair<int, int> rhit={b->core.pos+1,b->core.pos+rlen};
      if(splicing_stat==false){
        aligned_segments_to_ref.push_back(rhit);
      }
      ReadMapping[bam_get_qname(b)].aligned_segments.push_back(aligned_segments_to_ref);
      ReadMapping[bam_get_qname(b)].rhits.push_back(rhit);
      ReadMapping[bam_get_qname(b)].is_spliced.push_back(splicing_stat);
    } //if(ReadMapping.find(bam_get_qname(b))==ReadMapping.end())
  } //while ((r = sam_read1(testfile, h, b)) >= 0)
  bam_destroy1(b);
  bam_hdr_destroy(h);
  ReadMapping[latest_read].findFreeRegion();
  ReadMapping[latest_read].annotateRead(sorted_gene_features,gene_features,tech,debug);
  ReadMapping[latest_read].writeAnnotation(File);
  ReadMapping.erase(latest_read);
  fprintf(stderr,"------------------------------------\nAll mappings read, totally %li reads, %li hits\n",num_read, num_hits);
}

  
  
  