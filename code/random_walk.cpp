/*
Performs the random walk simulation on a given landscape, under a given genetic code.
Parameters:
	[1] Name of the file containing the genotype-phenotype map. Format of the file: First line header, each other line contains one sequence and its phenotype, tab-delimited.
	[2] Name of the file containing the genetic code. Format of the file: First line header, each other line tab-delimited aa (one-letter code) and a corresponding codon.
	[3] Name of the output file.
	[4] Population size (integer).

Compile as
g++ -std=c++14 -fmax-errors=1 -O3 -Wall -o random_walk random_walk.cpp
*/

#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <set>
#include <map> 
#include <unordered_map>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <string>

using namespace std;

// scores of the genotypes
vector<double> score;
// the genetic code				
map<string, char> code; //GCA -> A
// the input data
unordered_map<string, double> data; //ABCD -> 0.1
// map int -> nucleotide
map<int, char> nucs;
// the minimum score in data - will be used as score of sequences containing stop codons
// the exact value will be found later
double min_score = 1000;

// function to read genetic code from file
void read_code(string file_name, unordered_map<string, char>& code){
	cerr << "reading code" << endl;	
	const int lines_cnt = 64;
	ifstream ifile;
	ifile.open(file_name); 

	string s;
	ifile >> s >> s; // the header
	for(int i = 0; i < lines_cnt; ++i){
		char aa;
		string codon;
		ifile >> aa >> codon;
		code[codon] = aa;
	}

	ifile.close();
}

// translates a nucleotide sequence to protein, using the given genetic code
string translate(string s){
	assert(s.size() % 3 == 0);
	string ret;
	for(int i = 0; i < s.size(); i += 3){
		ret += code[ s.substr(i, 3) ];
	}
	return ret;
}

// generates random mRNA sequence of a given length
string gen_random_seq(int length, int seed) {
	srand(seed);
	string res = "";

	for (int i = 0; i<length; ++i) {
		char new_char = nucs[rand() % 4];
		res += new_char;
	}
	return res;
}

// performs the random walk
// seq = starting sequence, num_steps = length of the walk, pop_size = population size
pair< vector<double >, string > random_walk(string seq, int num_steps, int pop_size){ 
	vector<double > ret;
	string t = translate(seq);
	double score = (data.count(t)) ? data[t] : min_score;
	ret.push_back(score);

	for(int i = 0; i < num_steps; ++i){
		string new_seq = seq;
		while (seq == new_seq) {
			// generate random mutation
			int mutation = rand() % (4*seq.size());
			int pos = mutation / 4;
			char nuc = nucs[mutation % 4];
			new_seq[pos] = nuc;
		}
		t = translate(new_seq);
		// score of the mutated sequence
		double new_score = (data.count(t)) ? data[t] : min_score;
		// acceptance probability
		double p_accept;
		if (score==new_score) {
			p_accept = 1.0/pop_size;
		} else {	
			p_accept =  (1 - exp(score - new_score)) / (1 - exp( pop_size * (score - new_score)));
		}
		
		if ((float) rand()/RAND_MAX <= p_accept) {	// accept the mutation
			seq = new_seq;
			score = new_score;
		}
		ret.push_back(score);
		
	}
	return {ret, translate(seq)};
}

// run num_starts random walks of length walk_length
pair< vector<double>, map<string,int>> run_random_walks(int num_starts, int walk_length, int pop_size){
	int len = 0;
	vector<double> scores;
	map<string,int> reached_seqs;
	for(int i = 0; i < walk_length; i+=10) scores.push_back(0);
	// run num_starts random walks
	for (int j = 0; j < num_starts; ++j) {
		auto it = data.begin();
		string seq = gen_random_seq(3*(it->first.size()), j+1); // start in a randomly chosen sequence
		pair< vector< double >, string> walk = random_walk(seq, walk_length, pop_size);
		++len;
		// save the scores reached after each 10 steps
		for(int i = 0; i < walk_length/10; ++i){
			scores[i] += walk.first[10*i];			
		}
		// store which sequence was reached
		if (reached_seqs.count(walk.second)) ++reached_seqs[walk.second];
		else reached_seqs[walk.second] = 1;
	}

	for (int i =0; i<scores.size(); ++i) scores[i] /= len;

	return {scores, reached_seqs};
}

////////////////////////////////////////////
int main(int argc, char** argv){
	// command line parameters
	string data_file_str = string(argv[1]);
	string code_file_str = string(argv[2]);
	string output_file_str = string(argv[3]);
	int pop_size = atoi(argv[4]);

	// read the input data
	cerr << "Reading data" << endl;
	ifstream data_file;
	data_file.open(data_file_str); 
	string s;
	data_file >> s >> s;	// header	
	while(true){
		string sseq;
		double val;
		data_file >> sseq >> val;
		if(sseq.size() == 0){
			break;
		}
		data[sseq] = val;
		if (val < min_score) min_score=val;	// compute the minimal score found in the data
	}
	cerr << "Num of nodes: " << data.size() << endl;

	// output file
	ofstream out_file;
	out_file.open(output_file_str);

	srand(0); // seed
	nucs[0]='A'; nucs[1]='C'; nucs[2]='G'; nucs[3]='U';

	// the genetic code
	cerr << "Reading code" << endl;	
	read_code(code_file_str);

	// run the random walks
	int num_starts = 100000;
	int walk_length = 1000;
	pair< vector<double>, map<string,int>> res = run_random_walks(num_starts, walk_length, pop_size);
		
	// output
	for(int i = 0; i < (res.first).size(); ++i){
		out_file << res.first[i] << "\t";	 
	} 
	// uncomment this if you also want to output which sequences were reached (caution! easily generates very large files!)
	/*
	for(auto it=(res.second).begin(); it!=(res.second).end(); ++it) {
		out_file << it->first << ":" << it->second << ",\n"[it==(res.second).end()];
	}*/
	out_file << endl;
}