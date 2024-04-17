#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<set>
#include <unordered_map>
#include <time.h>
#include<algorithm>

#include <iostream>
#include <windows.h>
#include <psapi.h>
#pragma comment(lib,"psapi.lib")

using namespace std;

string datatest_str = "../TestData/";

string path_predata = datatest_str + "dataset/875713/";// Data set
string path_result_root = datatest_str + "result/APSim/875713/";// Output result

string path_query = datatest_str + "query/875713/query100.txt";//query file



bool cmp(const pair<int, double>& pair_1, const pair<int, double>& pair_2)
{
	return pair_1.second > pair_2.second;
}


class GraphPush {
public:
	int n;    //number of nodes
	unsigned long long m;    //number of edges
	int** inAdjList; // = NULL;
	int** outAdjList; // = NULL;
	int* indegree; // = NULL;
	int* outdegree; // = NULL;
	double* indegRecip; // = NULL;

	GraphPush() {
	}


	~GraphPush() {
		if (indegree != nullptr)
			delete[]indegree;
		if (indegRecip != nullptr)
			delete[]indegRecip;
		if (outdegree != nullptr)
			delete[]outdegree;
		if (inAdjList != nullptr) {
			for (int i = 0; i < n; ++i) {
				delete[] inAdjList[i];
			}
			delete[] inAdjList;
		}
		if (outAdjList != nullptr) {
			for (int i = 0; i < n; ++i) {
				delete[] outAdjList[i];
			}
			delete[] outAdjList;
		}
	}

	void loadGraph(string file, string file_deg) {
		string datafilename = file;
		ifstream datafile(datafilename.c_str());
		datafile >> n;
		string degfilename = file_deg;
		ifstream degfile(degfilename.c_str());
		indegree = new int[n];
		outdegree = new int[n];
		int nid = 0;
		while (degfile >> indegree[nid] >> outdegree[nid]) {
			nid++;
		}
		degfile.close();

		cout << n;
		inAdjList = new int* [n];
		outAdjList = new int* [n];
		for (int vid = 0; vid < n; ++vid) {
			if (indegree[vid] > 0)
				inAdjList[vid] = new int[indegree[vid]];
			if (outdegree[vid] > 0)
				outAdjList[vid] = new int[outdegree[vid]];
		}

		int* inAdjIdxes = new int[n];
		int* outAdjIdxes = new int[n];
		for (int i = 0; i < n; i++) {
			inAdjIdxes[i] = 0;
			outAdjIdxes[i] = 0;
		}
		int from;
		int to;
		m = 0;
		while (datafile >> from >> to) {
			outAdjList[from][outAdjIdxes[from]] = to;
			outAdjIdxes[from]++;
			inAdjList[to][inAdjIdxes[to]] = from;
			inAdjIdxes[to]++;
			m++;
		}

		datafile.close();
		delete[]inAdjIdxes;
		delete[]outAdjIdxes;
		indegRecip = new double[n];
		for (int i = 0; i < n; ++i) {
			if (indegree[i] > 0) indegRecip[i] = 1.0 / (double)indegree[i];
		}

	}
	int getInSize(int vert) {
		return indegree[vert];
	}
	int getInVert(int vert, int pos) {
		return inAdjList[vert][pos];
	}
	int getOutSize(int vert) {
		return outdegree[vert];
	}
	int getOutVert(int vert, int pos) {
		return outAdjList[vert][pos];
	}
};

class GraphTest
{
public:
	int n;
	long m;
	int** inAdjList = nullptr;
	int** outAdjList = nullptr;
	int* indegree = nullptr;
	int* outdegree = nullptr;
	double* indegRecip;
	GraphTest() {
	}

	void InputGraph(string file)
	{
		m = 0;
		ifstream infile(file);
		infile >> n;
		indegree = new int[n];
		outdegree = new int[n];
		for (int i = 0; i < n; i++)
		{
			indegree[i] = 0;
			outdegree[i] = 0;
		}
		int from;
		int to;
		while (infile >> from && infile >> to)
		{
			outdegree[from]++;
			indegree[to]++;
		}
		inAdjList = new int* [n];
		outAdjList = new int* [n];
		for (int i = 0; i < n; i++)
		{
			inAdjList[i] = new int[indegree[i]];
			outAdjList[i] = new int[outdegree[i]];
		}
		int* pointer_in = new int[n];
		int* pointer_out = new int[n];
		for (int i = 0; i < n; i++)
		{
			pointer_in[i] = 0;
			pointer_out[i] = 0;
		}
		infile.clear();
		infile.seekg(0);
		infile >> n;
		while (infile >> from >> to) {
			outAdjList[from][pointer_out[from]] = to;
			pointer_out[from]++;
			inAdjList[to][pointer_in[to]] = from;
			pointer_in[to]++;
			m++;
		}
		indegRecip = new double[n];
		for (int i = 0; i < n; ++i) {
			if (indegree[i] > 0) indegRecip[i] = 1.0 / (double)indegree[i];
			else
				indegRecip[i] = 0;
		}
		infile.close();
		delete[] pointer_in;
		delete[] pointer_out;
	}

	~GraphTest() {

		if (indegree != nullptr)
			delete[]indegree;
		if (outdegree != nullptr)
			delete[] outdegree;
		if (indegRecip != nullptr)
			delete[] indegRecip;
		if (inAdjList != nullptr) {
			for (int i = 0; i < n; ++i) {
				delete[] inAdjList[i];
			}
			delete[] inAdjList;
		}
		if (outAdjList != nullptr) {
			for (int i = 0; i < n; ++i) {
				delete[] outAdjList[i];
			}
			delete[] outAdjList;
		}

	}

	int getInSize(int vert) {
		return indegree[vert];
	}
	int getInVert(int vert, int pos) {
		return inAdjList[vert][pos];
	}
	int getOutSize(int vert) {
		return outdegree[vert];
	}
	int getOutVert(int vert, int pos) {
		return outAdjList[vert][pos];
	}

};

class APSimStruct {
public:
	GraphTest g;
	APSimStruct(string file) {
		g.InputGraph(file);
	}
	~APSimStruct() {
	}

	void APSim(ofstream& ofs_Result, vector <int>& vec_query, int PathLength, int top_k, double DecayFactor, double delta) {


		stringstream str_PathLength;
		str_PathLength << PathLength;
		stringstream str_delta;
		str_delta << delta;
		string time_str = path_result_root + "L" + str_PathLength.str() + "_delta" + str_delta.str() + "/" + "TimeCost.txt";
		ofstream ofs_time(time_str);

		double time_sum = 0;
		for (vector<int>::iterator it_queryid = vec_query.begin(); it_queryid != vec_query.end(); ++it_queryid) {
			int queryid = *it_queryid;
			clock_t start_bg = clock();
			unordered_map<int, unordered_map<int, double>> map_SM_t_j, map_pst_SM_t_j;///xxxx

			map_SM_t_j.clear();
			map_pst_SM_t_j.clear();
			unordered_map<int, double> map_Sim;
			/************************query processing begin****************************************/
			unordered_map<int, set<int>> map_Len_NodeSet;
			map_Len_NodeSet[0].insert(queryid);
			for (int len = 1; len <= PathLength - 1; ++len)
				for (set<int>::iterator it_i = map_Len_NodeSet[len - 1].begin(); it_i != map_Len_NodeSet[len - 1].end(); ++it_i) {
					int indgr = g.indegree[*it_i];
					for (int in = 0; in < indgr; ++in) {
						int inid = g.inAdjList[*it_i][in];
						map_Len_NodeSet[len].insert(inid);
					}
				}
			for (int len = PathLength - 1; len >= 0; --len) {
				map_SM_t_j.clear();
				for (set<int>::iterator it_t = map_Len_NodeSet[len].begin(); it_t != map_Len_NodeSet[len].end(); ++it_t) {
					if (len == PathLength - 1) {
						int indgr = g.indegree[*it_t];
						for (int i = 0; i < indgr; ++i) {
							int iid = g.inAdjList[*it_t][i];
							int outdgr = g.outdegree[iid];
							for (int j = 0; j < outdgr; ++j) {
								int jid = g.outAdjList[iid][j];
								if (map_SM_t_j[*it_t].find(jid) != map_SM_t_j[*it_t].end())
									map_SM_t_j[*it_t][jid] += DecayFactor * g.indegRecip[*it_t] * g.indegRecip[jid];
								else
									map_SM_t_j[*it_t][jid] = DecayFactor * g.indegRecip[*it_t] * g.indegRecip[jid];
							}
						}
						map_SM_t_j[*it_t][*it_t] = 1;

					}
					else {
						unordered_map<int, double> Map_psm_key;
						int indgr = g.indegree[*it_t];
						for (int i = 0; i < indgr; ++i) {
							int iid = g.inAdjList[*it_t][i];
							for (unordered_map<int, double>::iterator it_key = map_pst_SM_t_j[iid].begin(); it_key != map_pst_SM_t_j[iid].end(); ++it_key)
								Map_psm_key[it_key->first] += it_key->second * g.indegRecip[*it_t];
						}

						for (unordered_map<int, double>::iterator it_key = Map_psm_key.begin(); it_key != Map_psm_key.end(); ++it_key)
							if (it_key->second >= delta) {
								int outdgr = g.outdegree[it_key->first];
								for (int j = 0; j < outdgr; ++j) {
									int jid = g.outAdjList[it_key->first][j];
									map_SM_t_j[*it_t][jid] += it_key->second * DecayFactor * g.indegRecip[jid];
								}

							}
						map_SM_t_j[*it_t][*it_t] = 1;
					}
				}
				for (unordered_map<int, unordered_map<int, double>>::iterator y = map_SM_t_j.begin(); y != map_SM_t_j.end(); ++y)
					for (unordered_map<int, double>::iterator z = y->second.begin(); z != y->second.end(); ++z)
						if (z->second < delta) {
							y->second.erase(z->first);
							if (y->second.size() == 0)
								map_SM_t_j.erase(y->first);
						}

				map_pst_SM_t_j.clear();
				map_pst_SM_t_j = map_SM_t_j;
			}
			clock_t end_bg = clock();
			map_Sim = map_SM_t_j[queryid];
			map_SM_t_j.clear();
			map_pst_SM_t_j.clear();
			map_Sim.erase(queryid);

			double time = (double)(end_bg - start_bg);
			time_sum += time;

			vector<pair<int, double>> map_QueryIDCenterIDSim;
			map_QueryIDCenterIDSim.clear();
			for (unordered_map<int, double>::iterator it = map_Sim.begin(); it != map_Sim.end(); ++it)
				map_QueryIDCenterIDSim.push_back(make_pair(it->first, it->second));
			map_Sim.clear();
			sort(map_QueryIDCenterIDSim.begin(), map_QueryIDCenterIDSim.end(), cmp);

			string stra = to_string(queryid);
			string str11 = path_result_root + "L" + str_PathLength.str() + "_delta" + str_delta.str() + "/" + stra + "klist.txt";
			ofstream ofs_klist(str11);


			ofs_time << queryid << "\tnode\tcosts\t" << (double)(end_bg - start_bg) / CLOCKS_PER_SEC << "s" << endl;

			ofs_Result << queryid << "\t";
			for (vector<pair<int, double>>::iterator it = map_QueryIDCenterIDSim.begin(); it != map_QueryIDCenterIDSim.end(); ++it) {
				ofs_Result << it->first << "\t" << it->second << "\t";
				ofs_klist << it->first << "\t" << it->second << endl;
			}
			ofs_Result << endl;
			ofs_klist.close();
			map_QueryIDCenterIDSim.clear();

		}

		ofs_time << "Total cost for all nodes " << (double)time_sum / CLOCKS_PER_SEC << "s" << endl;
		ofs_time << "Average cost per node " << (double)time_sum / vec_query.size() / CLOCKS_PER_SEC << "s" << endl;


		ofs_time.close();
		ofs_Result.close();
	}

};


int main() {

	ifstream ifs_queryNode(path_query);
	vector<int>vec_query;
	string buf;
	while (getline(ifs_queryNode, buf))
	{
		int a = atoi(buf.c_str());
		vec_query.push_back(a);
	}
	ifs_queryNode.close();

	int PathLength = 1;//1,2,3,4,5
	int top_k = 50;
	double DecayFactor = 0.6;
	double delta = 0.001;//  1-0.001  2-0.0005 3-0.0001  4-0.00005   5-0.00001 (Pruning_APSim)

	stringstream str_PathLength;
	str_PathLength << PathLength;
	stringstream str_top_k;
	str_top_k << top_k;
	stringstream str_DecayFactor;
	str_DecayFactor << DecayFactor;
	stringstream str_delta;
	str_delta << delta;

	ofstream ofs_Result(path_result_root + "Result_TopFast_C" + str_DecayFactor.str() + "_L" + str_PathLength.str() + "_delta" + str_delta.str() + "_top" + str_top_k.str() + ".txt");
	string file = path_predata + "relationshipmatrix.txt";

	APSimStruct TT(file);
	TT.APSim(ofs_Result, vec_query, PathLength, top_k, DecayFactor, delta);

	cout << "ok" << endl;
	system("pause");
	return 0;

}