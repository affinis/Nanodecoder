#include <unistd.h>
#include <stdlib.h>
#include <cstring>
#include <string>
#include <fstream>
#include <cstring>
#include <vector>
#include <unordered_map>
#include <iostream>

using namespace std;
/*
using std::string;
using std::unordered_map;
using std::ifstream;
using std::vector;
using std::fstream;
using std::endl;
using std::cout;
 */

ifstream samfile;
ifstream id2namefile;
fstream out_count_file;
fstream chimeraReadlist;
fstream chimeraReadMappingStat;
unordered_map<string,string> id2name;

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

pair<int,int> findClip(string const &cigar){
  size_t start;
  size_t end = 0;
  vector<string> splits;
  pair<int,int> res;
  res.first=-1;
  res.second=-1;
//  cout << "findClip.debug1" << endl;
  while ((start = cigar.find_first_not_of("SHDIMPX=N", end)) != string::npos)
  {
//    cout << "findClip.while.debug1" << endl;
    end = cigar.find_first_of("SHDIMPX=N", start);
//    cout << "findClip.while.debug2" << endl;
    splits.push_back(cigar.substr(start, end-start));
  }
//  cout << "findClip.debug2" << endl;
  int seg_num=splits.size();
  if(seg_num==1){
    return res;
  }
  size_t first_flag=cigar.find_first_of("SHDIMPX=N");
  size_t last_flag=cigar.find_last_of("SHDIMPX=N");
  if(cigar[first_flag]=='S'||cigar[first_flag]=='H'){
    res.first=stoi(splits[0]);
  }
  if(cigar[last_flag]=='S'||cigar[last_flag]=='H'){
    res.second=stoi(splits[splits.size()-1]);
  }
  return res;
}

int main(int argc, char *argv[]){
  int option;
  const char* Samfilename;
  const char* Id2namefilename;
  while((option=getopt(argc,argv,"f:m:"))!=-1){
    switch(option){
    case 'f':
      Samfilename=optarg;
      break;
    case 'm':
      Id2namefilename=optarg;
      break;
    } // end switch
  }// end while
  cout << "Loading ID to gene name mappings." << endl;
  id2namefile.open(Id2namefilename);
  while(!id2namefile.eof()){
    string line;
    getline(id2namefile,line,'\n');
    if(line==""){
      continue;
    }
    vector<string> fields=tokenize(line,'\t');
    id2name[fields[0]]=fields[1];
//    cout << fields[0] << "\t" << fields[1] << endl;
    }
  id2namefile.close();
  cout << "Mappings loaded" << endl;
  samfile.open(Samfilename);
  chimeraReadlist.open("chimeric_read_ids.lst",fstream::in|fstream::out|fstream::app);
  out_count_file.open("reads_properly_mapped.tsv",fstream::in|fstream::out|fstream::app);
  chimeraReadMappingStat.open("chimera_mapping_stat.tsv",fstream::in|fstream::out|fstream::app);
  unordered_map<string,string> read_bc2gene;
  unordered_map<string,string> chimera;
  while(!samfile.eof()){
    string status;
    getline(samfile,status,'\n');
    if(status==""){
      continue;
    }
    int isHeader=status.find('@');
    if(isHeader==0){
      continue;
    }else{
      string ref;
      vector<string> fields=tokenize(status,'\t');
      string read_bc=fields[0];
      int flag=stoi(fields[1]);
      if(chimera[read_bc].length()!=0){
        continue;
      }
      if(fields[2]!="*"){
        ref=id2name[fields[2]];
      }else{
        continue;
      }
      if(read_bc2gene[read_bc].length()==0){
        read_bc2gene[read_bc]=ref;
      }else{
        if(flag>=2048){
          int default_SA_fields=21;
//          cout << read_bc << endl;
          vector<string> SA_fields;
          pair<int,int> clip_pos;
          pair<int,int> SA_clip_pos;
          read_bc2gene.erase(read_bc);
          string readid=(tokenize(read_bc,'_'))[0];
//          cout << "debug1" << endl;
          string cigar=fields[5];
//          cout << "debug2" << endl;
          string SA_info=fields[default_SA_fields];
          while(SA_info.find("SA:Z:")!=0){
            SA_info=fields[default_SA_fields++];
          }
//          cout << "debug3" << endl;
          SA_fields=tokenize(SA_info,',');
//          cout << "debug4" << endl;
          clip_pos=findClip(cigar);
//         cout << "debug5" << endl;
          SA_clip_pos=findClip(SA_fields[3]);
//          cout << "debug6" << endl;
          if(clip_pos.first>=SA_clip_pos.first&&clip_pos.second>=SA_clip_pos.second){
            read_bc2gene[read_bc]=id2name[SA_fields[0].replace(0,5,"")];
            continue;
          }else if(clip_pos.first<=SA_clip_pos.first&&clip_pos.second<=SA_clip_pos.second){
            read_bc2gene[read_bc]=ref;
            continue;
          }else{
            chimeraReadMappingStat << read_bc << "\t" << clip_pos.first << "\t" << cigar << "\t" << clip_pos.second << "\t" << SA_clip_pos.first << "\t" << SA_fields[3] << "\t" << SA_clip_pos.second << endl;
            chimera[read_bc]=read_bc;
            chimeraReadlist << readid << endl;  
          }
        }else{
          continue;
        } //if(flag>=2048)
      } //if(read_bc2gene[read_bc].length()==0&&chimera[read_bc].length()==0)
    } //if(isHeader==0)
  } //while(!samfile.eof())
  for(auto it:read_bc2gene){
    out_count_file << it.first << "\t" << it.second << endl;
  }
  chimeraReadlist.close();
  out_count_file.close();
  chimeraReadMappingStat.close();
}
