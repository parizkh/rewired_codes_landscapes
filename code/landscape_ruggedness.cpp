/*
Computes the landscape ruggedness characteristics (number of peaks; prevalence of epistasis; peak accessibility) of a GB1 landscape
assuming a given genetic code.
Parameters:
	[1] Name of the file containing the genotype-phenotype map. Format of the file: First line header, each other line contains one sequence and its phenotype, tab-delimited.
	[2] Name of the file containing the genetic code. Format of the file: First line header, each other line tab-delimited aa (one-letter code) and a corresponding codon.
	[3] Name of the output file.

Compile as
g++ -std=c++14 -fmax-errors=1 -O3 -Wall -o landscape_ruggedness landscape_ruggedness.cpp
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
#include <queue>
#include <climits>
#include <numeric>

using namespace std;

// This value is used for sequences containing stop codons.
const double EMPTY_VAL = -100.0;

struct Node{
	vector<int> neighbors;  
};


// Converts an integer to a DNA sequence of length 12.
string int_to_string(int u, int L){
	string ret = "";
	for(int i = 0; i < L*3; ++i){
		switch(u % 4){
			case 0: ret += "A"; break;
			case 1: ret += "C"; break;
			case 2: ret += "G"; break;
			case 3: ret += "U"; break;
		}
		u >>= 2;
	}
	return ret;
}

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
string translate(string nuc_seq, unordered_map<string, char>& code){
	string res = "";
	for (int i = 0; i<nuc_seq.size()/3; ++i) {
		string codon = nuc_seq.substr(3*i,3);
		res += string(1, code[codon]);
	}
	return res;
}

//read the phenotypic scores of the sequences
void read_scores(string file_name, vector<double>& scores, unordered_map<string, char>& code, int L){
	cerr << "reading scores" << endl;	
	// here we will store the scores of each nucleotide sequence
	scores.clear();
	scores.resize(1 << (6*L), EMPTY_VAL); // -100 for empty entries because of the stop codon

	// read the input data - scores for the protein sequences
	ifstream data_file;
	data_file.open(file_name); 
	string s;
	data_file >> s >> s;
	unordered_map<string, double> aa_scores;
	while(true){
		string sseq;
		double val;
		data_file >> sseq >> val;
		if(sseq.size() == 0){
			break;
		}
		aa_scores[sseq] = val;
	}

	cerr << "scores loaded, converting" << endl;
	// for each nucleotide sequence, find its translation and store the corresponding score
	for(int i = 0; i < (1 << (6*L)); ++i){
		string str = int_to_string(i, L);
		string aa_str = translate(str, code);
		if (aa_str.find('*') == string::npos) {
			scores[i] = aa_scores[aa_str];
		}

	}

	cerr << "num of nodes: " << aa_scores.size() << " " << scores.size() << endl;
}

// build the genotype network
void build_graph(vector<Node >& G, unordered_map<string, char>& code, unordered_map<string, int>& node_names, vector<double>& scores, int L){
	cerr << "building the graph" << endl;	
	G.clear();
	G.resize((1 << (6*L)), Node());

	for(int i = 0; i < G.size(); ++i){
		// generate all possible neighbors
		for(int j = 0; j < (3*L); ++j){
			char cur_val = (i >> (2 * j)) % 4;
			int base = i - (cur_val << (2*j));

			for(char new_val : {0, 1, 2, 3}){
				if(new_val == cur_val) continue;				
				int neighbor = base + (new_val << (2*j));
				G[i].neighbors.push_back(neighbor);
			}
		}
	}
}

// helper function: compute mean and variance of a vector of doubles
pair<double, double> mean_var(vector<double> v) {
	double sum = std::accumulate(v.begin(), v.end(), 0.0);
	double mean = sum / v.size();

	double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
	double var = sq_sum / v.size() - mean * mean;

	return {mean, var};
}

//////////////////////////////////// Landscape ruggedness functions

// Computes the number of peaks in the landscape and their mean score
pair<int, double> count_peaks(vector<Node>& G, vector<double>& scores){
	int cnt = 0;
	vector<double> peak_scores;
	
	vector<int> visited(scores.size(), -1); // have we seen this node before?
	// sort the sequences based on the scores, in descending order
	vector<pair<double, int> > sorted_seqs;
	for(int i = 0; i < scores.size(); ++i){
		sorted_seqs.push_back({scores[i], i});
	}
	sort(sorted_seqs.begin(), sorted_seqs.end(), greater<>());

	// iterate over the sequences in the order of descending score
	for(int i = 0; i < sorted_seqs.size(); ++i){
		// if we haven't visited this sequence before
		if(visited[sorted_seqs[i].second] == -1){			
			// this is a candidate peak	
			bool isPeak = true;		
			double peak_score = sorted_seqs[i].first;
			
			// explore the whole plateau (in case some neighbors have the same fitness, check that none of them was already visited)
			vector<int> peak_stack;
			peak_stack.push_back(sorted_seqs[i].second);			
			while(!peak_stack.empty()){
				int u = peak_stack[peak_stack.size()-1];
				peak_stack.pop_back();
				visited[u] = i;

				for(int v : G[u].neighbors){
					if(scores[v] == peak_score){
						if(visited[v] == -1){
							peak_stack.push_back(v);
						}	
						if(visited[v]!=-1 && visited[v]!=i) {
							isPeak = false;
						}				
					}
					visited[v] = i;		// this node was visited in iteration i
				}
			}

			if (isPeak==true) {	//this is a peak
				++cnt;
				peak_scores.push_back(peak_score);
			}
			
		}
		else {
			// not a peak
			// just visit all neighbors
			for(int v : G[sorted_seqs[i].second].neighbors){
				visited[v] = i;
			}
		}
	}

	// compute the mean peak score
	pair<double, double> scores_mean = mean_var(peak_scores);

	return {cnt, scores_mean.first};
}

// runs BFS on the landscape and computes the number of paths and the number of accessible paths from each node to the global peak
void BFS(vector<Node>& G, vector<double>& scores, vector<int>& dist, vector<long long>& num_paths, vector<long long>& num_accessible, vector<int>& global_max, bool oriented){
	// number of paths (regardless of accessibility) from each node to the global peak
	num_paths = vector<long long>(G.size(), 0);
	// number of accessible paths
	num_accessible = vector<long long>(G.size(), 0);
	// distances of the nodes from the global peak
	dist = vector<int>(G.size(), -1);

	// initialize the queue; add the global peak
	queue<int> q;
	for (int i = 0; i<global_max.size(); ++i) {
		q.push(global_max[i]);
		dist[global_max[i]] = 0;
		num_paths[global_max[i]] = 1;
		num_accessible[global_max[i]] = 1;
	}	

	// run BFS
	while (!q.empty()) {
		int u = q.front();
		q.pop();
		int d = dist[u];
		// add all neighbors to the queue
		for (int v : G[u].neighbors) {
			if (dist[v] == -1) {
				q.push(v);
				dist[v] = d+1;
			}
			// if this is the shortest possible way how to reach node v, update number of paths
			if (dist[v] == d+1) {
				num_paths[v] += num_paths[u];	// all paths
				if (scores[v] <= scores[u]) {	// accessible paths
					num_accessible[v] += num_accessible[u];
				}
				
			}
		}
	}
}

// Samples num_squares random squares of sequences (a wild-type, two single mutants and a corresponding double mutant) from the genotype network.
vector<vector<int> > sample_squares(vector<Node>& G, vector<double>& scores, int num_squares, int L){
	vector<vector<int>> squares;	// vector to store the results
	for (int i = 0; i<num_squares; ++i) {
		// the wildtype
		int seq0 = rand() % scores.size();
		int seq1; // first single mutant 
		int seq2; // second single mutant
		int seq12;	// double mutant
		
		// sample two random mutations
		// first mutation
		int pos1; int mut1; char cur_val1; 
		// position of the mutation
		pos1 = rand() % (L*3);
		while(true) {
			// mutates to
			mut1 = rand() % 4;
			cur_val1 = (seq0 >> (2 * pos1)) % 4;	
			int base = seq0 - (cur_val1 << (2*pos1));
			seq1 = base + (mut1 << (2*pos1));
			if (seq1!=seq0) break;
		}		
			
		// second mutation
		int pos2; int mut2; char cur_val2; 
		while(true) {
			// position of the mutation
			pos2 = rand() % (L*3);
			// the mutation is in a different position than the first one
			if(pos2 != pos1) {
				// mutates to
				mut2 = rand() % 4;
				cur_val2 = (seq0 >> (2 * pos2)) % 4;
				int base = seq0 - (cur_val2 << (2*pos2));
				seq2 = base + (mut2 << (2*pos2));
				if (seq2 != seq0) break;
			}			
		}

		// combination of the two mutations
		int base = seq0 - (cur_val1 << (2*pos1)) - (cur_val2 << (2*pos2));
		seq12 = base + (mut1 << (2*pos1)) + (mut2 << (2*pos2));			
		
		squares.push_back({seq0, seq1, seq2, seq12});
	}
	return squares;
}

// Computes the prevalence of epistasis type (no epistasis, magnitude, simple-sign or reciprocal-sign) among the squares in vectore squares.
vector<double> epistasis_types(vector<vector<int>> squares, vector<double>& scores) {
	int num_squares = squares.size();
	int mag = 0; int ss = 0; int rs = 0; int no_epi = 0;
	for (int i = 0; i < squares.size(); i++) {
		vector<int> square = squares[i];
		
		// scores of the sequences
		double score_AB = scores[square[0]];
		double score_aB = scores[square[1]];
		double score_ab = scores[square[3]];
		double score_Ab = scores[square[2]];

		// score differences
		double delta_ab_Ab = score_Ab - score_ab;
		double delta_ab_aB = score_aB - score_ab;
		double delta_Ab_AB = score_AB - score_Ab;
		double delta_aB_AB = score_AB - score_aB;

		// no epistasis
		if (score_AB + score_ab - score_Ab - score_aB == 0) {
			no_epi++;
		}
		// magnitude epistasis
		else if (delta_ab_Ab*delta_aB_AB >= 0 && delta_ab_aB*delta_Ab_AB>=0) {
			mag++;
		}
		// reciprocal sign epistasis
		else if (delta_ab_Ab*delta_aB_AB < 0 && delta_ab_aB*delta_Ab_AB<0) {
			rs++;
		} else ss++; // otherwise it simple-sign epistasis
	}
	vector<double> res = {((double)no_epi)/num_squares, ((double)mag)/num_squares, ((double)ss)/num_squares, ((double)rs)/num_squares};
	return res;
}

///////////////////////////////////////////////////////////
int main(int argc, char** argv){
	srand(1);
	// command line parameters
	string data_file_str = string(argv[1]);
	string code_file_str = string(argv[2]);
	string output_file_str = string(argv[3]);
	int L = atoi(argv[4]);
	
	//read the genetic code
	unordered_map<string, char> code; //GCA -> A
	read_code(code_file_str, code);
	
	// read the genotype-phenotype landscape
	vector<double> scores;
	read_scores(data_file_str, scores, code, L);

	//construct the genotype network
	vector<Node> G;
	unordered_map<string, int> node_names;
	build_graph(G, code, node_names, scores, L);

	// find global maxima
	vector<int> global_max;
	double max_score = -100;
	for (int i = 0; i<scores.size(); ++i ) {
		if (scores[i] >= max_score) {
			if (scores[i] == max_score) {
				global_max.push_back(i);
			}
			else {
				global_max = {i};
			}			
			max_score = scores[i];
		}
	}

	// output file
	ofstream out_file;
	out_file.open(output_file_str);
	
	// number of peaks
	cout << "Peaks" << endl;
	pair<int, double> res_peaks = count_peaks(G, scores);
	out_file << res_peaks.first << "\t" << res_peaks.second << "\t";	

	// epistasis
	cout << "Epistasis" << endl;
	vector<vector<int> > squares = sample_squares(G, scores, 1000000, L);
	vector<double> res = epistasis_types(squares, scores);
	out_file << res[0] << "\t" << res[1] << "\t" << res[2] << "\t" << res[3] << "\t";

	// proportion of accessible paths
	cout << "Accessible paths" << endl;
	vector<long long> num_paths;
	vector<long long> num_accessible;
	vector<int> dist;
	BFS(G, scores, dist, num_paths, num_accessible, global_max, false);
	double prob_acc = 0.0;
	int num_nonempty = 0;
	for (int i=0; i<scores.size(); ++i) {
		if (scores[i]!=EMPTY_VAL) {
			prob_acc += 1.0*num_accessible[i]/num_paths[i];
			num_nonempty++;
		}		
	}
	prob_acc /= num_nonempty;

	out_file << prob_acc << endl;
}