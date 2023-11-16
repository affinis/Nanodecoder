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

#define bam_is_umapped(b) (((b)->core.flag&BAM_FUNMAP) != 0)

regex gene_regex("(gene_id \"| gene_type \"| gene_biotype \"| gene_name \"|\")");
regex transcript_regex("( transcript_id \"| transcript_type \"| transcript_name \"|\")");


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
 * @abstract transcript feature
 * 
 * @start 1-based
 * @end 1-based
 * @transcript_id
 * @transcript_name
 * @exons pairs of 1-based start&end sites of this gene indicating exon positions
 */
struct transcript_feature{
  int start;
  int end;
  string transcript_id;
  string transcript_name;
  vector<pair<int,int>> exons;
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
  unordered_map<string,transcript_feature> transcripts;
  
  /*!
   * @abstract fill in data from parsed GTF fields, notice when input data are from transcript/exon rather than gene, 
   *           related gene_id/transcript_id should already exists. 
   *            
   * @param core_data, tokenized (with '\t') fields, corresponding to first 8 column of GTF
   * @param aux_data, last (9th) column of GTF, tokenized with ';', the info in this fields depends on core_data[2] (feature type)
   */
  void fill_in_feature_info(vector<string>& core_data,vector<string>& aux_data, unordered_map<string, int>& afileds_info){
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
      string i_gene_id=regex_replace(aux_data[0],gene_regex,"");
      transcript_id_raw=get_aux_field("transcript_id",aux_data,afileds_info,"exon");
      string i_transcript_id=regex_replace(transcript_id_raw,transcript_regex,"");
      if(this->gene_id==i_gene_id&&this->transcripts.find(i_transcript_id)!=transcripts.end()){
        pair<int,int> i_pos={stoi(core_data[3]),stoi(core_data[4])};
        this->transcripts[i_transcript_id].exons.push_back(i_pos);
      }else{
        fprintf(stderr,"fill_in_feature_info: ERROR: fill_in_feature_info: type: exon, gene_id: %s does not match to transcript_id: %s, or transcript_id %s is not a key, check GTF\n",
                this->gene_id.c_str(),i_transcript_id.c_str(),i_transcript_id.c_str());
        exit(0);
      }
    }
  }
  
  void debugPrintInfo(){
      string transcripts_string="";
      if(transcripts.empty()){
        transcripts_string="\n";
      }else{
        for(auto feature:this->transcripts){
          string exon_str="";
          for(pair<int,int> exon:feature.second.exons){
            exon_str=exon_str+to_string(exon.first)+","+to_string(exon.second)+";";
          }
          transcripts_string=transcripts_string+"\n\\_______"+feature.second.transcript_id+'\t'+\
            feature.second.transcript_name+'\t'+to_string(feature.second.start)+"\t"+to_string(feature.second.end)+"\t"+exon_str;
        }
      }
      fprintf(stderr,"%s\t%s\t%s\t%s\t%i\t%i%s\n\n",gene_name.c_str(),gene_id.c_str(),gene_type.c_str(),chromosome.c_str(),start,end,transcripts_string.c_str());
  }
};

/*!
 * mapping_info
 * 
 * @abstract: parsed read mapping info from BAM 
 * 
 * @ID                    read id
 * @read_length
 * @final_annotation_id   gene id of the final annotation of read (a read may have multiple alignment to ref).
 * @final_annotation_name gene name of the final annotation of read.
 * @read_orientation      orientation of this read infered from its mapping, "forward" or "reverse" or "*" 
 *                        (for read not mapped).
 * @qcovs                 query coverage of hits.
 * @qhits                 vector of pair<int,int>, indicate mapping start and end position on read, 1-based.
 * @hit_strands           orientation of hits corresponding to ref, "forward" or "reverse".
 * @refs                  chromosome id (in header of genome fasta).
 * @features              gene ids of hits according to GTF file (a read may have multiple alignment to ref), 
 *                        notice that features.size() == qhits.size(), as some of the hits may not overlap with 
 *                        certain gene, i.e. not mapped to transcriptome, these hits would be annoatated as 
 *                        "non-coding".
 * @rcovs                 coverage of gene (for non-coding hits, rcov=0).
 * @rhits                 vector of pair<int,int>, indicate mapping start and end positions on chromosome, 1-based.
 * @read_orientations     orientation of read infered from its mappings, "forward" or "reverse" or "*".
 * @occupiedRegion        regions that occupied by hit(s).
 * @freeRegion            vector of pair<int,int>, indicate regions that not mapped to genome, 1-based.
 * @mergedRegions        vector of vector of int, indicate index (0-based) of qhits that were merged due to overlaping 
 *                        during calculation of occupiedRegion.
 * @mapped_to_genome
 * @mapped_to_transcriptome
 * 
 */

struct mapping_info{
  string ID;
  int read_length;
  string final_annotation_id="default";
  string final_annotation_name="default";
  string read_orientation="default";
  vector<float> qcovs;
  vector<pair<int, int>> qhits;
  vector<string> hit_strands;
  vector<string> refs;
  vector<string> features;
  vector<float> rcovs;
  vector<pair<int, int>> rhits;
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
    
    fprintf(stderr,"%s\t%s\t%s\t%s\t%s\n",this->ID.c_str(),ref_str.c_str(),strand_str.c_str(),\
            hits_str.c_str(),free_region_str.c_str());
  }//void debugPrintInfo()
  
  void writeAnnotation(ofstream& File){
    string free_region_str=vector_pair_int2str(this->freeRegion);
    string qhits_str=vector_pair_int2str(this->qhits);
    string qcov_str=stringCat(this->qcovs);
    string merged_hits_str="";
    string occupied_region_str=vector_pair_int2str(this->occupiedRegion);
    string rhits_str=vector_pair_int2str(this->rhits);
    string rcov_str=stringCat(this->rcovs);
    if(this->mergedRegions.empty()){
      merged_hits_str="*";
    }else{
      for(vector<int> indices:this->mergedRegions){
        merged_hits_str=stringCat(indices)+";";
      }
      merged_hits_str.erase(merged_hits_str.end()-1);
    }
    File << this->ID << "\t" << this->final_annotation_id << "\t" << this->final_annotation_name << "\t";
    File << this->read_orientation << "\t" << free_region_str << "\t" << qhits_str << "\t" << qcov_str << "\t";
    File << rhits_str << "\t" << rcov_str << "\t" <<  merged_hits_str << "\t";
    File << occupied_region_str << "\n";
  }
  
  
  /*!
   * @abstract annotate read with gene features
   * 
   * @param sorted_gene_features: .first chromosome
   *                              .second <map>, end position on its genome, sorted key
   *                              .second.first feature start
   *                              .second.second gene_id
   * @param gene_features: .first gene_id
   *                       .second mapping_info
   * @param debug
   * 
   */
  
  void annotateRead(unordered_map<string, map<int,pair<int,string>>>& sorted_gene_features, \
                    unordered_map<string, gene_feature>& gene_features, bool debug){
    if(false){
      debug=true;
    }
    vector<int> hits_transcript_idx;
    if(this->qhits.empty()){
      if(debug){
        cout << this->ID << ": unmapped" << endl;
      }
      this->mapped_to_genome=false;
      this->mapped_to_transcriptome=false;
      this->final_annotation_id="*";
      this->final_annotation_name="*";
      this->read_orientation="*";
      return;
    }else{
      this->mapped_to_genome=true;
    }
    for(int hit_index=0;hit_index<this->rhits.size();hit_index++){
      auto feature1=sorted_gene_features[this->refs[hit_index]].upper_bound(this->rhits[hit_index].first);
      auto feature2=sorted_gene_features[this->refs[hit_index]].upper_bound(this->rhits[hit_index].second);
      
      //no hits on genes throughout genome
      if(feature1==sorted_gene_features[this->refs[hit_index]].end()){
        if(debug){
          cout << this->ID << ": non-coding" << endl;
          continue;
        }else{
          this->features.push_back("*");
          this->read_orientations.push_back("*");
          int qstart=this->qhits[hit_index].first;
          int qend=this->qhits[hit_index].second;
          float qcov=static_cast<float>(qend-qstart+1)/this->read_length;
          this->qcovs.push_back(qcov);
          this->rcovs.push_back(0);
          continue;
        }
      }
      
      int feature_start=feature1->second.first;
      int feature_end=feature1->first;
      int hit_start=this->rhits[hit_index].first;
      int hit_end=this->rhits[hit_index].second;
      int qstart=this->qhits[hit_index].first;
      int qend=this->qhits[hit_index].second;
      string feature_id=feature1->second.second;
      string feature_strand=num2rna[gene_features[feature_id].strand];
      string hit_strand=this->hit_strands[hit_index];
      float qcov=static_cast<float>(qend-qstart+1)/this->read_length;
      this->qcovs.push_back(qcov);
      
      /*
       * feature->second.first: start position of feature
       *        feature->first: end position of feature
       *    map.rhits[x].first: start position of hit
       *   map.rhits[x].second: end position of hit
       * 
       * hit                +-----------------------+
       * feature    +-----------------------+
       */
      if(feature2==sorted_gene_features[this->refs[hit_index]].end()||feature_id!=feature2->second.second){
        float rcov=static_cast<float>(feature_end-hit_start+1)/(feature_end-feature_start+1);
        if(debug){
          cout << this->ID << ": hit: " << hit_start << "-" << hit_end << ", feature: ";
          cout << feature_id << ": " << feature_start << "-" << feature_end << " " << qcov << " " << rcov <<endl;
          continue;
        }else{
          this->features.push_back(feature_id);
          this->rcovs.push_back(rcov);
          hits_transcript_idx.push_back(hit_index);
          if(feature_strand==hit_strand){
            this->read_orientations.push_back("forward");
          }else{
            this->read_orientations.push_back("reverse");
          }
          continue;
        }
      }else{
        /*
         * hit     +-----------------------+  
         * feature                              +--------------------------+
         */
        if(feature_start>=hit_end){
          float rcov=0;
          if(debug){
            cout << this->ID << ": hit: " << hit_start << "-" << hit_end << ", feature: ";
            cout << feature_id << ": " << feature_start << "-" << feature_end;
            cout << " non-coding" << endl;
            continue;
          }else{
            this->features.push_back("*");
            this->rcovs.push_back(rcov);
            this->read_orientations.push_back("*");
            continue;
          }
          /*
           * hit            +--------------------+
           * feature     +--------------------------+
           * 
           * hit         +--------------------+
           * feature           +--------------------------+
           */
        }else{
          float rcov;
          if(hit_start>=feature_start){
            rcov=static_cast<float>(hit_end-hit_start+1)/(feature_end-feature_start+1);
          }else if(hit_start<feature_start){
            rcov=static_cast<float>(hit_end-feature_start+1)/(feature_end-feature_start+1);
          }else{
            fprintf(stderr,"annotateRead: ERROR, unexpected condition, hit_start-hit_end;\
            feature_start-feature_end: %i-%i;%i-%i",hit_start,hit_end,feature_start,feature_end);
            exit(0);
          }
          if(debug){
            cout << this->ID << ": hit: " << hit_start << "-" << hit_end << ", feature: ";
            cout << feature_id << ": " << feature_start << "-" << feature_end << " "  << qcov << " " << rcov <<endl;
            continue;
          }else{
            this->features.push_back(feature_id);
            this->rcovs.push_back(rcov);
            hits_transcript_idx.push_back(hit_index);
            if(feature_strand==hit_strand){
              this->read_orientations.push_back("forward");
            }else{
              this->read_orientations.push_back("reverse");
            }
            continue;
          }
        }
      }
    } //for(int hit_index=0;hit_index<this->rhits.size();hit_index++)
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
    this->read_orientation=stringCat(this->read_orientations);
    /*
    if(hits_transcript_idx.size()>0){
      this->mapped_to_transcriptome=true;
      if(hits_transcript_idx.size()==1){
        this->final_annotation_id=this->features[hits_transcript_idx[0]];
        this->final_annotation_name=gene_features[this->final_annotation_id].gene_name;
        this->read_orientation=this->read_orientations[hits_transcript_idx[0]];
      }else{
        this->final_annotation_id=stringCat(this->features);
        this->final_annotation_name="";
        for(string id:this->features){
          this->final_annotation_name=this->final_annotation_name+gene_features[id].gene_name+",";
        }
        this->final_annotation_name.erase(this->final_annotation_name.end()-1);
        this->read_orientation=stringCat(this->read_orientations);
      }
    }else{
      this->mapped_to_transcriptome=false;
      this->final_annotation_id="*";
      this->final_annotation_name="*";
      this->read_orientation="*";
    }
     */
  } //annotateRead
  
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
void readFeaturesFromGTF(const char * GTFfilename, unordered_map<string, gene_feature>& gene_features){
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
 * @abstract sort features by start site in each chromosome
 * 
 * @param gene_features: reference (c++) of raw data from GTF filled in by readFeaturesFromGTF().
 * @param sorted_features: reference (c++) of sorted data, sorting are performed within 
 *        chromosome (unordered_map.first) by end position (map.first). map.second.first = start position,
 *        map.second.second = gene_name.
 */
void sortFeaturesByCoord(unordered_map<string, gene_feature>& gene_features, unordered_map<string, map<int, pair<int,string>>>& sorted_features){
  fprintf(stderr,"Start sorting features\n");
  for(auto feature:gene_features){
    if(sorted_features[feature.second.chromosome].empty()){
      map<int, pair<int,string>> i_chr;
      sorted_features[feature.second.chromosome]=i_chr;
      sorted_features[feature.second.chromosome][feature.second.end]={feature.second.start,feature.second.gene_id};
    }else{
      sorted_features[feature.second.chromosome][feature.second.end]={feature.second.start,feature.second.gene_id};
    }
  }
  fprintf(stderr,"Finished sorting features\n");
}

void readMappingFile(const char * Bamfilename, unordered_map<string, mapping_info>& ReadMapping, unordered_map<string, map<int,pair<int,string>>>& sorted_gene_features, \
                     unordered_map<string, gene_feature>& gene_features, bool debug, ofstream& File){
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
      if(num_read%50000000==0){
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
    int ncigar=b->core.n_cigar;
    char left_copchr=bam_cigar_opchr(CIGAR[0]);
    int left_coplen=bam_cigar_oplen(CIGAR[0]);
    char right_copchr=bam_cigar_opchr(CIGAR[ncigar-1]);
    int right_coplen=bam_cigar_oplen(CIGAR[ncigar-1]);
    int qlen=bam_cigar2qlen(ncigar,CIGAR);
    int rlen=bam_cigar2rlen(ncigar,CIGAR);
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
    //if read_id not exists
    if(ReadMapping.find(bam_get_qname(b))==ReadMapping.end()){
      thismapping.ID=bam_get_qname(b);
      thismapping.hit_strands.push_back(flag);
      pair<int, int> qhit={q_start,q_end};
      thismapping.qhits.push_back(qhit);
      thismapping.refs.push_back(ref_name);
      pair<int, int> rhit={b->core.pos+1,b->core.pos+rlen};
      thismapping.rhits.push_back(rhit);
      thismapping.read_length=qlen;
      ReadMapping[bam_get_qname(b)]=thismapping;
      num_read++;
      if(latest_read.length()!=0){
        if(num_read<7){
          ReadMapping[latest_read].debugPrintInfo();
        }
        if(num_read==7){
          fprintf(stderr,"............\n");
        }
        ReadMapping[latest_read].findFreeRegion();
        ReadMapping[latest_read].annotateRead(sorted_gene_features,gene_features,false);
        ReadMapping[latest_read].writeAnnotation(File);
        ReadMapping.erase(latest_read);
      }
      latest_read=thismapping.ID;
      if(num_read%50000000==0){
        fprintf(stderr,"processed %li reads\n",num_read);
      }
    }else{
      ReadMapping[bam_get_qname(b)].hit_strands.push_back(flag);
      pair<int, int> qhit={q_start,q_end};
      ReadMapping[bam_get_qname(b)].qhits.push_back(qhit);
      ReadMapping[bam_get_qname(b)].refs.push_back(ref_name);
      pair<int, int> rhit={b->core.pos+1,b->core.pos+rlen};
      ReadMapping[bam_get_qname(b)].rhits.push_back(rhit);
    }
  }
  bam_destroy1(b);
  bam_hdr_destroy(h);
  ReadMapping[latest_read].debugPrintInfo();
  fprintf(stderr,"------------------------------------\nAll mappings read, totally %li reads, %li hits\n",num_read, num_hits);
}
  
  
  