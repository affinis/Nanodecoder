#include <string.h>
#include <iostream>
#include <cstring>
#include <vector>
#include <bits/stdc++.h>
#include <thread>
#include <cmath>
#include <algorithm>

using std::string;
using std::vector;
//using std::max;
//using std::min;
using std::cout;
using std::endl;
using std::unordered_map;
using std::to_string;
using std::pair;

#define GAP_PENALTY 0
#define MATCH_SCORE 1
#define MISMATCH_SCORE 0

vector<string> genrateKmerFromSeq(string& seq, const int k);
vector<int> kmerDistances(string& queryseq, vector<string>& kmers);
void editDistance(int id, string seq1, string seq2, vector<int> results);
int editDistance(string& seq1, string& seq2);

char base_to_number[26] {1,-1,2,-1,-1,-1,3,-1,-1,-1,-1,-1,-1,0,-1,-1,-1,-1,-1,4,-1,-1,-1,-1,-1,-1};
char number_to_base[26] {'A','C','G','T'};

unordered_map<int, string> num_to_strand{
  {0,"original"},{1,"reverse"},{2,"complement"},{3,"reverse_complement"},{-1,"*"}
};

unordered_map<string, int> strand2num{
  {"+",0},{"-",3},{"forward",0},{"reverse",3},{"*",-1}
};

unordered_map<int, string> num2rna{
  {0,"forward"},{3,"reverse"},{-1,"*"}
};

unordered_map<int, int> reverse_strand{
  {0,3},{3,0},{-1,-1}
};

struct alignment {
  vector<char> aligned_seq1 {};
  vector<char> aligned_symbol {};
  vector<char> aligned_seq2 {};
  };

struct simple_hit{
  int hstart;
  int strand;
  int mismatch;
};

struct kmerDistancesThreadData{
  int threadID;
  string queryseq;
  vector<string> iblockkmers;
  };

struct kmerCandidate{
  int index;
  int strandness;  //0=original;1=reverse;2=complement;3=reverse_complement
  };

string reverse(string seq){
  string seq_rev=seq;
  for(int i=0;i<seq.length();++i){
    seq_rev[i]=seq[seq.length()-1-i];
    }
  return seq_rev;
  }

string complement(string seq){
  string seq_complement=seq;
  for(int i=0;i<seq.length();++i){
    seq_complement[i]=number_to_base[4-base_to_number[seq[i]-'A']];
    }
  return seq_complement;
  }

string reverse_complement(string seq){
  string seq_rc=seq;
  for(int i=0;i<seq.length();++i){
    seq_rc[i]=number_to_base[4-base_to_number[seq[seq.length()-1-i]-'A']];
    }
  return seq_rc;
  }

int w(const int base1,const int base2){
//  cout << base1 << endl;
//  cout << base2 << endl;
  if(base1==base2){
    return MATCH_SCORE;
  }else{
    return MISMATCH_SCORE;
    }
  }

int nuc_global_alignment(int* seq_1, int* seq_2,const int i,const int j){
  int scoreMat[(i+1)*(j+1)];
  scoreMat[0]=0;
//  cout << scoreMat[0] << "\t";
  for(int n=1;n<=j;n++){
    scoreMat[n]=n*GAP_PENALTY;
//    cout << scoreMat[n] << "\t";
    }
//  cout << endl;
  for(int n=1;n<=i;n++){
    scoreMat[n*(j+1)]=n*GAP_PENALTY;
    }
//  cout << scoreMat[1] << endl;
  
  for(int n=1;n<=i;n++){
//    cout << scoreMat[n*(j+1)+0] << "\t";
    for(int m=1;m<=j;m++){
      int left=scoreMat[n*(j+1)+m-1]+GAP_PENALTY;
      int upper=scoreMat[(n-1)*(j+1)+m]+GAP_PENALTY;
      int upperleft=scoreMat[(n-1)*(j+1)+m-1];
//      cout << scoreMat[n*(j+1)+m-i-1] << endl;
      int siteScore=w(seq_1[n-1],seq_2[m-1]);
      upperleft+=siteScore;
      scoreMat[n*(j+1)+m]=std::max(std::max(upper,left),upperleft);
//      cout << scoreMat[n*(j+1)+m] << "\t";
      }
//    cout << endl;
    }
  return scoreMat[(i+1)*(j+1)-1];
  }


vector<string> genrateKmerFromSeq(string& seq, const int k){
  vector<string> Kmers;
  const int len=seq.length();
  int i=0;
  while(i+k<=len){
    string kmer=seq.substr(i,k);
    Kmers.push_back(kmer);
    i++;
    }
  return Kmers;
  }

int basediff(string& seq1, string& seq2){
  int dist=0;
  if(seq1.length()!=seq2.length()){
    cout << "Length not equal while using function 'basediff()'";
    exit(1);
    }
  for(int i=0;i<seq1.length();i++){
      if(seq1[i]!=seq2[i]){
        dist++;
        }
    }
  return dist;
  }

bool seqEqual(string& seq1, string& seq2){
  if(seq1.length()!=seq2.length()){
    return false;
  }else{
    for(int i;i<=seq1.length();i++){
      if(seq1[i]!=seq2[i]){
        return false;
        }
      }
    return true;
    }
  }

int hamming_distance(string& A, string& B)
{
  int dist = 0;
  for (int i = 0; i < B.length(); ++i)
  {
    dist += (A[i] != B[i]);
  }
  return dist;
}

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

int editDistance(string& seq1, string& seq2, int indel_num){
  
  int n = seq1.length();
  int m = seq2.length();
  
  if(indel_num>=n||indel_num>=m){
    return 0;
  }
  vector<int> candidates;
  int dp[n + 1][m + 1];
  for (int i = 0; i <= n; i++)
  {
    for (int j = 0; j <= m; j++)
    {
      if (i == 0){
        dp[i][j] = j;
        }
      else if (j == 0){
        dp[i][j] = i;
        }
      else if (seq1[i-1] == seq2[j-1]){
        dp[i][j] = dp[i-1][j-1];
        }
      else{
        dp[i][j] = 1 + std::min(dp[i][j-1], std::min(dp[i-1][j], dp[i-1][j-1]));
        }
      if((i>=n-indel_num&&j==m)||(j>=m-indel_num&&i==n)){
        candidates.push_back(dp[i][j]);
      }
    }
  }
  int dist=*std::min_element(candidates.begin(),candidates.end());
  return dist;
}


int minDistance(string w1, string w2) {
  int n = w1.size();
  int m =w2.size();
  int** dp = new int*[n+1];
  for(int i =0;i<=n;i++){
    dp[i] = new int[m+1];
    for(int j=0;j<=m;j++){
      dp[i][j]=0;
      if(i==0)dp[i][j]=j;
      else if(j==0)dp[i][j] = i;
    }
  }
  w1 = " " + w1;
  w2 = " " + w2;
  for(int i =1;i<=n;i++){
    for(int j = 1;j<=m;j++){
      if(w1[i] !=w2[j]){
        dp[i][j] = 1+std::min({dp[i-1][j],dp[i][j-1],dp[i-1][j-1]});
      } else {
        dp[i][j] = dp[i-1][j-1];
      }
    }
  }
  return dp[n][m];
}


void editDistances(int id, string seq1, string seq2, vector<int> results){
  int n = seq1.length();
  int m = seq2.length();
  int dp[n + 1][m + 1];
  for (int i = 0; i <= n; i++)
  {
    for (int j = 0; j <= m; j++)
    {
      if (i == 0){
        dp[i][j] = j;
      }
      else if (j == 0){
        dp[i][j] = i;
      }
      else if (seq1[i-1] == seq2[j-1]){
        dp[i][j] = dp[i-1][j-1];
      }
      else{
        dp[i][j] = 1 + std::min(dp[i][j-1], std::min(dp[i-1][j], dp[i-1][j-1]));
      }
    }
  }
  results[id]=dp[n][m];
}

/*
vector<int> kmerDistances(string& queryseq, vector<string>& kmers){
  vector<int> distances;
  string revquery=reverse(queryseq);
  string comquery=complement(queryseq);
  string rcquery=reverse_complement(queryseq);
  for(vector<string>::iterator it = kmers.begin(); it != kmers.end(); ++it){
    int index=distance(kmers.begin(), it);
    int oriseqDist=editDistance(queryseq, *it);
    int revDist=editDistance(revquery, *it);
    int comDist=editDistance(comquery, *it);
    int rcDist=editDistance(rcquery, *it);
    int dist=min(min(oriseqDist,revDist),min(comDist,rcDist));
    distances.push_back(dist);
  }
  return distances;
}
 */

/*
kmerCandidate minKmerDistance(string& queryseq, vector<string>& kmers,
                          bool r=false, bool c=false, bool rc=false){
  string revquery=NULL;
  string comquery=NULL;
  string rcquery=NULL;
  if(r){
    revquery=reverse(queryseq);
    kmerCandidate revCandidate;
    
    }
  if(c){
    comquery=complement(queryseq);
    }
  if(rc){
    rcquery=reverse_complement(queryseq);
    }
  }
 */

vector<int> getMaxComplexitySegments(int total_length, int num_segments){
  int base_length=floor(total_length/num_segments);
//  cout << "base length: " << base_length << endl;
  int num_longer=total_length%num_segments;
//  cout << "num_longer: " << num_longer << endl;
  vector<int> segments_length;
  for(int i=0;i<num_segments;i++){
    segments_length.push_back(base_length);
    }
  for(int j=0;j<num_longer;j++){
    segments_length[j]++;
    }
//  cout << "num of segments: " << segments_length.size() << endl;
  return segments_length;
  }

vector<string> getMaxComplexitySegments(string& sequence, vector<int> segments){
  if(accumulate(segments.begin(),segments.end(),0)!=sequence.length()){
//    cout << sequence << endl;
//    for(vector<int>::iterator it=segments.begin();it!=segments.end();++it){
//      cout << *it << " ";
//      }
//    cout << endl;
    cout << "Segment lengths did not add up to sequence length" << endl;
    exit(0);
    }
  vector<string> segment_sequences;
  int start=0;
  for(int i=0;i<segments.size();i++){
    string segment=sequence.substr(start,segments[i]);
    segment_sequences.push_back(segment);
    start=start+segments[i];
    }
  return segment_sequences;
  }

vector<int> getMaxComplexityWordSizes(int total_length, int num_segments){
  vector<int> word_sizes;
  int base_length=floor(total_length/num_segments);
  int longer_length=base_length+1;
  word_sizes.push_back(longer_length);
  word_sizes.push_back(base_length);
  return word_sizes;
  }

vector<int> seq2CurrentLevels(string& seq, std::unordered_map<string, int>& current_table){
  vector<string> Sixmers=genrateKmerFromSeq(seq,6);
  vector<int> current_levels;
  for(string kmer:Sixmers){
    current_levels.push_back(current_table[kmer]);
  }
  return(current_levels);
}

string stringCat(vector<string>& vec,char sep=',',bool trim=true){
  string res="";
  if(!trim){
    res=res+sep;
  }
  for(string str:vec){
    res=res+str+sep;
  }
  if(!trim){
    return(res);
  }else{
    res.erase(res.end()-1);
    return(res);
  }
}

string stringCat(vector<float>& vec,char sep=',',bool trim=true){
  string res="";
  if(!trim){
    res=res+sep;
  }
  for(float flo:vec){
    res=res+to_string(flo)+sep;
  }
  if(!trim){
    return(res);
  }else{
    res.erase(res.end()-1);
    return(res);
  }
}
  
string stringCat(vector<int>& vec,char sep=',',bool trim=true){
  string res="";
  if(!trim){
    res=res+sep;
  }
  for(int integer:vec){
    res=res + to_string(integer) + sep;
  }
  if(!trim){
    return(res);
  }else{
    res.erase(res.end()-1);
    return(res);
  }
}

string stringCat(vector<bool>& vec, char sep=',',bool trim=true){
  string res="";
  if(!trim){
    res=res+sep;
  }
  for(bool boo:vec){
    res=res+to_string(boo)+sep;
  }
  if(!trim){
    return(res);
  }else{
    res.erase(res.end()-1);
    return(res);
  }
}

vector<pair<int,int>> string2regions(string& str){
  vector<pair<int,int>> res;
  vector<string> regions=tokenize(str,';');
  int i=0;
  for(string region:regions){
    if(region==""){
      continue;
    }else if(region=="-1--1"){
      pair<int,int> region_int={0,0};
      res.push_back(region_int);
    }else{
      pair<int,int> region_int;
      region_int.first=stoi(tokenize(region,'-')[0]);
      region_int.second=stoi(tokenize(region,'-')[1]);
      res.push_back(region_int);
    }
  }
  return(res);
}

vector<vector<int>> string2hit_overlaps(string& str){
  vector<vector<int>> res;
  if(str=="*"||str.empty()){
    return(res);
  }
  vector<int> hit_group;
  for(int i=0;i<str.length();i++){
    if(str[i]==','){
      continue;
    }
    if(str[i]==';'){
      res.push_back(hit_group);
      hit_group.clear();
    }
    hit_group.push_back(str[i]);
  }
  return(res);
}

string regions2string(vector<pair<int,int>>& regions,char sep=';'){
  if(regions.empty()){
    return("*");
  }
  string res="";
  for(pair<int,int> region:regions){
    res=res+to_string(region.first)+"-"+to_string(region.second)+sep;
  }
  res.erase(res.end()-1);
  return(res);
}

vector<int> string2strands(string& str){
  vector<int> res;
  if(str.empty()){
    return(res);
  }
  vector<string> fields=tokenize(str,',');
  for(string strand_str:fields){
    if(strand_str==""){
      continue;
    }
    res.push_back(strand2num[strand_str]);
  }
  return(res);
}

vector<float> string2covs(string& str){
  vector<float> res;
  vector<string> fields=tokenize(str,',');
  for(string cov_str:fields){
    if(cov_str==""){
      continue;
    }else if(cov_str.find("-")==0){
      res.push_back((float)0);
      continue;
    }
    res.push_back(stof(cov_str));
  }
  return(res);
}

vector<int> string2ints(string& str){
  vector<int> res;
  vector<string> fields=tokenize(str,',');
  for(string int_str:fields){
    if(int_str==""){
      continue;
    }
    res.push_back(stoi(int_str));
  }
  return(res);
}

vector<int> string2splice(string& str){
  vector<int> res;
  vector<string> fields=tokenize(str,',');
  for(string splice_str:fields){
    if(splice_str==""){
      continue;
    }
    res.push_back(stoi(splice_str));
  }
  return(res);
}

string strands2string(vector<int>& strands){
  string res="";
  if(strands.empty()){
    return("*");
  }
  for(int strand:strands){
    res=res+num2rna[strand]+",";
  }
  res.erase(res.end()-1);
  return(res);
}

bool isContaining(int query, pair<int,int>& target){
  if(query>=target.first&&query<=target.second){
    return(true);
  }else{
    return(false);
  }
}

bool isContaining(pair<int,int>& query, pair<int,int>& target){
  if(query.first>=target.first&&query.second<=target.second){
    return(true);
  }else{
    return(false);
  }
}


//calculate query coverage
float calculate_cov(pair<int,int>& query, pair<int,int>& subject){
  float cov;
  if(isContaining(query,subject)){
    cov=1;
  }else if(isContaining(query.first,subject)&&!isContaining(query.second,subject)){
    cov=static_cast<float>(subject.second-query.first+1)/(query.second-query.first+1);
  }else if(!isContaining(query.first,subject)&&isContaining(query.second,subject)){
    cov=static_cast<float>(query.second-subject.first+1)/(query.second-query.first+1);
  }else if(isContaining(subject,query)){
    cov=static_cast<float>(subject.second-subject.first+1)/(query.second-query.first+1);
  }else{
    cov=0;
  }
  return(cov);
}

int calculate_min_overhang(pair<int,int>& query, pair<int,int>& subject, int feature_strand, int tech){
  int query_length=query.second-query.first+1;
  int subject_length=subject.second-subject.first+1;
  int res;
  if(query.second<subject.first||query.first>subject.second){
    res=99999;
  }else{
    int case_key=tech+feature_strand;
    //fprintf(stderr,"%i\n",case_key);
    switch(case_key){
    //3' sequencing, gene is on positive strand
      case 3:
        res=abs(query.second-subject.second);
        break;
        //fprintf(stderr,"3' sequencing, gene is on positive strand\n");
    //5' sequencing, gene is on positive strand
      case 5:
        res=abs(query.first-subject.first);
        break;
        //fprintf(stderr,"5' sequencing, gene is on positive strand\n");
    //3' sequencing, gene is on negative strand
      case 6:
        res=abs(query.first-subject.first);
        break;
        //fprintf(stderr,"3' sequencing, gene is on negative strand\n");
    //5' sequencing, gene is on negative strand
      case 8:
        res=abs(query.second-subject.second);
        break;
        //fprintf(stderr,"5' sequencing, gene is on negative strand\n");
    }
    //return(std::min(abs(query.first-subject.first),abs(query.second-subject.second)));
  }
  return(res);
}

int pair_length(pair<int,int>& pair_in){
  if(pair_in.first>pair_in.second){
    return(pair_in.first-pair_in.second+1);
  }
  return(pair_in.second-pair_in.first+1);
}

bool exon_is_same(std::vector<int> &first, std::vector<int> &second){
  if (first.size()!=second.size()) {
    return false;
  }
  return std::is_permutation(first.begin(), first.end(), second.begin());
}

bool exon_is_contained(std::vector<int> first, std::vector<int> second){
  //we do not need to sort as the exon order is fixed in GTF
  
  // Sort first vector
  //std::sort(first.begin(), first.end());
  // Sort second vector
  //std::sort(second.begin(), second.end());
  // Check if  all elements of a second vector exists in first vector
  return std::includes(first.begin(), first.end(), second.begin(), second.end());
}

bool exon_is_contained(string first, string second){
  //we do not need to sort as the exon order is fixed in GTF
  if(first.find(second)!=string::npos){
    return(true);
  }else{
    return(false);
  }
}

bool exon_is_contained(std::vector<int> first, int single_element){
  vector<int> second;
  second.push_back(single_element);
  // Sort first vector
  std::sort(first.begin(), first.end());
  // Sort second vector
  std::sort(second.begin(), second.end());
  // Check if  all elements of a second vector exists in first vector
  return std::includes(first.begin(), first.end(), second.begin(), second.end());
}

bool region_is_overlap(pair<int,int> first, pair<int,int> second){
  if(first.first>second.second||second.first>first.second){
    return(false);
  }else{
    return(true);
  }
}
