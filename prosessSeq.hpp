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
using std::max;
using std::min;
using std::cout;
using std::endl;
using std::unordered_map;

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
  {"+",0},{"-",3}
};

unordered_map<int, string> num2rna{
  {0,"forward"},{3,"reverse"}
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
      scoreMat[n*(j+1)+m]=max(max(upper,left),upperleft);
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
        dp[i][j] = 1 + min(dp[i][j-1], min(dp[i-1][j], dp[i-1][j-1]));
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
        dp[i][j] = 1+min({dp[i-1][j],dp[i][j-1],dp[i-1][j-1]});
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
        dp[i][j] = 1 + min(dp[i][j-1], min(dp[i-1][j], dp[i-1][j-1]));
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

/*
vector<int> kmerDistancesMultithreads(string& queryseq, vector<string>& kmers){
  int numKmers=kmers.size();
  vector<int> distances[numKmers];
  for(int i=0;i<numKmers;++i){
    string kmer=kmers[i];
    thread threadApp {editDistance,queryseq,kmer};
    }
  }
*/

/*
int main(){
  alignment align_main;
  int seq1num[seq1.length()];
  int seq1num_r[seq1.length()];
  int seq1num_rc[seq1.length()];
  int seq2num[seq2.length()];
  for(int i=0;i<seq1.length();i++){
    seq1num[i]=base_to_number[seq1[i]-'A'];
//    cout << seq1num[i];
    }
//  cout << endl;
  for(int j=0;j<seq2.length();j++){
    seq2num[j]=base_to_number[seq2[j]-'A'];
//    cout << seq2num[j];
  }
//  cout << endl;
  ///get reverse and reverse complement sequence of seq1
  for(int i=0;i<seq1.length();i++){
    seq1num_r[i]=base_to_number[seq1[seq1.length()-1-i]-'A'];
    seq1num_rc[i]=3-base_to_number[seq1[seq1.length()-1-i]-'A'];
    }
  int score=nuc_global_alignment(seq1num,seq2num,seq1.length(),seq2.length());
  int score_r=nuc_global_alignment(seq1num_r,seq2num,seq1.length(),seq2.length());
  int score_rc=nuc_global_alignment(seq1num_rc,seq2num,seq1.length(),seq2.length());
  cout << "query score to subject: " << score << endl;
  cout << "reverse query score to subject: " << score_r << endl;
  cout << "reverse complement query score to subject: " << score_rc << endl;
  }
*/
