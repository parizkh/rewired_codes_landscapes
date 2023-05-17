/*
Performs the greedy walk simulation on a ParD3 landscape, under a given genetic code.
Parameters:
	[1] Name of the file containing the genotype-phenotype map. Format of the file: First line header, each other line contains one sequence and its phenotype, tab-delimited.
	[2] Name of the file containing the genetic code. Format of the file: First line header, each other line tab-delimited aa (one-letter code) and a corresponding codon.
	[3] Name of the output file.

Compile as
g++ -std=c++14 -fmax-errors=1 -O3 -Wall -o greedy_walk_ParD3 greedy_walk_ParD3.cpp
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

using namespace std;

const double EMPTY_VAL = -100.0;

struct Node{
	vector<int> neighbors;
	
	//computed in the algorithm
	vector<int> peaks; // indices of accessible peaks
	vector<double> path; // we use lexicographical order on these values in the alg.  
};

struct Event{
	vector<double> path; // path from a peak to node u
	int u; //node to relax
	int par_u; // parent that relaxed u
};

// how to compare two paths
class Compare
{
public:
    bool operator() (Event u, Event v)
    {
    	return lexicographical_compare(u.path.begin(), u.path.end(), v.path.begin(), v.path.end());
	}
};

// Converts an integer to a DNA sequence of length 9.
string int_to_string(int u){
	string ret = "";
	for(int i = 0; i < 9; ++i){
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

//read the phenotypic scores of the sequences
void read_scores(string file_name, vector<double>& scores, unordered_map<string, char>& code){
	cerr << "reading scores" << endl;	
	// here we will store the scores of each nucleotide sequence
	scores.clear();
	scores.resize(1 << 18, EMPTY_VAL); // -100 for empty entries because of the stop codon

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
	// for each nucleotide sequence of length 12, find its translation and store the corresponding score
	for(int i = 0; i < (1 << 18); ++i){
		string str = int_to_string(i);
		string s0 = str.substr(0,3);
		string s1 = str.substr(3,3);
		string s2 = str.substr(6,3);

		if(string(1,code[s0]) != "*" &&  string(1,code[s1]) != "*" && string(1,code[s2]) != "*" ){
			scores[i] = aa_scores[
				string(1, code[s0]) + string(1, code[s1]) + string(1, code[s2]) 
				];
		} 
	}

	cerr << "num of nodes: " << aa_scores.size() << " " << scores.size() << endl;
}

// build the genotype network
void build_graph(vector<Node >& G, unordered_map<string, char>& code, unordered_map<string, int>& node_names, vector<double>& scores){
	cerr << "building the graph" << endl;	
	G.clear();
	G.resize((1 << 18), Node());

	for(int i = 0; i < G.size(); ++i){
		// generate all possible neighbors
		for(int j = 0; j < 9; ++j){
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

// translates a nucleotide sequence to protein, using the given genetic code
string translate(string nuc_seq, unordered_map<string, char>& code){
	string res = "";
	for (int i = 0; i<nuc_seq.size()/3; ++i) {
		string codon = nuc_seq.substr(3*i,3);
		res += string(1, code[codon]);
	}
	return res;
}


// find the greedy paths
void find_paths(vector<Node>& G, vector<double>& scores){
	cerr << "finding paths" << endl;
	priority_queue<Event, vector<Event>, Compare> q;
	vector<double> parent_score(1<<18, -1);

	// initialize the queue: add all nodes; the queue is sorted based on the scores
	for(int i = 0; i < (1 << 18); ++i){
		if(scores[i] == EMPTY_VAL) continue;
		Event e = Event();
		e.path = {scores[i]};
		e.u = i;
		e.par_u = -1;
		q.push(e);

		parent_score[i] = scores[i];
	}

	vector<bool> relaxed(1<<18, false);

	while(!q.empty()){
		// the path in the queue beginning with the highest score
		Event e = q.top();
		q.pop();

		if(relaxed[e.u]){	
			// if there are two equally good paths starting from u
			if(G[e.u].path == e.path){ // by construction of the priority queue, e.path is always <= G[e.u].path
				// add the peaks that are reached by the current parent of u
				G[e.u].peaks.insert(G[e.u].peaks.end(), G[e.par_u].peaks.begin(), G[e.par_u].peaks.end());
				sort(G[e.u].peaks.begin(), G[e.u].peaks.end());
				auto it = unique(G[e.u].peaks.begin(), G[e.u].peaks.end());
				G[e.u].peaks.erase(it, G[e.u].peaks.end());
			}
		}
		else{
			// relax
			relaxed[e.u] = true;
			// update G
			G[e.u].path = e.path;
			if(e.par_u == -1){
				G[e.u].peaks = {e.u};
			}
			else{
				G[e.u].peaks = G[e.par_u].peaks;
			}

			// add neighbors to the queue, if 
			// 1) we are not just moving within a local peak
			// 2) neighbor v is not already in the queue with a better parent
			for(int v : G[e.u].neighbors){
				if(!(G[e.u].path.size() == 1 && scores[v] == G[e.u].path[0]) &&	parent_score[v] <= scores[e.u]){
					Event new_e = Event();
					
					new_e.u = v;
					new_e.par_u = e.u;

					vector<double> new_path = {scores[v]};
					new_path.insert(new_path.end(), G[e.u].path.begin(), G[e.u].path.end());
					new_e.path = new_path;

					// only add if v is not relaxed, or if the current path of v is the same as new_path 
					if (!relaxed[v] || new_path == G[v].path) {
						q.push(new_e);
						parent_score[v] = scores[e.u];
					}				
				} 
			}
		}
	}

}





int main(int argc, char** argv){
	// command line parameters
	string data_file_str = string(argv[1]);
	string code_file_str = string(argv[2]);
	string output_file_str = string(argv[3]);
	
	//read the genetic code
	unordered_map<string, char> code; //GCA -> A
	read_code(code_file_str, code);

	
	// genotype-phenotype map
	vector<double> scores;
	read_scores(data_file_str, scores, code);
	

	// the genotype network
	vector<Node> G;
	unordered_map<string, int> node_names;
	build_graph(G, code, node_names, scores);


	// the main computation
	find_paths(G, scores);

	// compute mean fitness reached and report which peaks are reached
	cout << "Computing mean number of steps, mean fitness and counting peaks." << endl;
	double sum = 0;
	double sum_len = 0;
	int num_vert = 0;
	map<string,double> reached_peaks;
	for (int i = 0; i < 1<<18; ++i) {
		if(scores[i] == EMPTY_VAL) continue;
		sum += G[i].path[G[i].path.size()-1];
		sum_len += G[i].path.size()-1;
		num_vert += 1;

		for (int j = 0; j<G[i].peaks.size(); ++j) {
			string translation = translate(int_to_string(G[i].peaks[j]), code);
			if (reached_peaks.count(translation)) reached_peaks[translation] += 1.0/G[i].peaks.size();
			else reached_peaks[translation] = 1.0/G[i].peaks.size();
		}

		assert(G[i].path[G[i].path.size()-1] == scores[G[i].peaks[0]]);


	}
	
	// output
	ofstream out_file;
	out_file.open(output_file_str);
	out_file << sum/num_vert << "\t" << sum_len/num_vert << "\t";
	for (auto it = reached_peaks.begin(); it!= reached_peaks.end(); ++it) {
		out_file << std::fixed << setprecision(2) <<  it->first << ":" << it->second << ",";
	}
	out_file << endl;


}




