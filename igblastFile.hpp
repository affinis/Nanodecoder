#include <fstream>
#include <string>
#include <vector>

using std::string;
using std::vector;

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

class igblastCloneRecord{
  string clone_id;
  string reppresent_read;
  int count;
  float frequency;
  string CDR3nucl;
  string CDR3aa;
  bool productive;
  string chain_type;
  string V_gene;
  string D_gene;
  string J_gene;
  string clone_def;
  
  igblastCloneRecord(string line){
    vector<string> fields;
    fields=tokenize(line,'\t');
    this->clone_id=fields[0];
    this->reppresent_read=fields[1];
    this->count=stoi(fields[2]);
    this->frequency=stof(fields[3]);
    this->clone_id=fields[0];
    this->clone_id=fields[0];
    this->clone_id=fields[0];
    this->clone_id=fields[0];
    this->clone_id=fields[0];
  }
};
