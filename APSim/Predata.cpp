#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <unordered_map>
#include <vector>
#include <string>
#include <ctime>

using namespace std;
 
string dataset_str = "../TestData/dataset/1005/";// Data set

void WeightGraphOfCenterIDCenterID(ifstream& ifs_CenterCenter, ofstream& ofs_CenterIDCenterIDWeight, ofstream& ofs_CenterIDCenter, ofstream& ofs_deg, ofstream& ofs_stat, ofstream & ofs_G_rev) {
	set<string> set_value;//v
	unordered_map<string, int> map_vk;
	string a, b;
	while (ifs_CenterCenter >> a >> b) {
		set_value.insert(a);
		set_value.insert(b);
	}

	int count_wordkv = 0;
	for (set<string>::iterator it = set_value.begin(); it != set_value.end(); ++it) {
		ofs_CenterIDCenter << count_wordkv << "\t" << *it << endl;
		count_wordkv++;
	}
	long vexcount = set_value.size();

	int count_word = 0;
	for (set<string>::iterator it = set_value.begin(); it != set_value.end(); ++it) {
		map_vk.insert(unordered_map<string, int>::value_type(*it, count_word));
		count_word++;
	}

	ifs_CenterCenter.clear();
	ifs_CenterCenter.seekg(0);
	string str_line2;
	set<pair<int, int>> set_link;
	while (ifs_CenterCenter >> a >> b) {
		int r, c;
		r = map_vk[a];
		c = map_vk[b];
		set_link.insert(make_pair(r, c));
	}
	vector<int> indg, outdg;
	for (int i = 0; i < set_value.size(); ++i) {
		indg.push_back(0);
		outdg.push_back(0);
	}

	ofs_CenterIDCenterIDWeight << set_value.size() <<endl;
	for (set<pair<int, int>>::iterator link = set_link.begin(); link != set_link.end(); ++link) {
		ofs_G_rev << link->second << "\t" << link->first << endl;
		ofs_CenterIDCenterIDWeight << link->first << "\t" << link->second << endl;
		outdg[link->first]++;
		indg[link->second] ++;
	}


	ofs_stat << "node #:\t" << set_value.size() << "\tedge #:\t" << set_link.size() << endl;
	for (int i = 0; i < set_value.size(); ++i)
		ofs_deg << indg[i] << "\t" << outdg[i] << endl;

	set_value.clear();
	set_link.clear();

}

void CheckUDG(ifstream& ifs_CenterCenter, ofstream& ofs_CenterIDCenterIDWeight, ofstream& ofs_CenterIDCenter, ofstream& ofs_deg, ofstream& ofs_stat) {
	set<string> set_value;//v
	unordered_map<string, int> map_vk;
	string a, b;
	while (ifs_CenterCenter >> a >> b) {
		set_value.insert(a);
		set_value.insert(b);
	}

	int count_wordkv = 0;
	for (set<string>::iterator it = set_value.begin(); it != set_value.end(); ++it) {
		ofs_CenterIDCenter << count_wordkv << "\t" << *it << endl;
		count_wordkv++;
	}
	long vexcount = set_value.size();

	int count_word = 0;
	for (set<string>::iterator it = set_value.begin(); it != set_value.end(); ++it) {
		map_vk.insert(unordered_map<string, int>::value_type(*it, count_word));
		count_word++;
	}

	ifs_CenterCenter.clear();
	ifs_CenterCenter.seekg(0);
	string str_line2;
	set<pair<int, int>> set_link;
	int x = 0;
	while (ifs_CenterCenter >> a >> b) {
		int r, c;
		r = map_vk[a];
		c = map_vk[b];
		if (r == c) {
			set_link.insert(make_pair(r, c));
			x++;
		}
		else {
			set_link.insert(make_pair(c, r));
			set_link.insert(make_pair(r, c));
		}
	}
	cout <<"num of undirected edges:\t"<<(set_link.size()+x)/2;

	vector<int> indg, outdg;
	for (int i = 0; i < set_value.size(); ++i) {
		indg.push_back(0);
		outdg.push_back(0);
	}

	ofs_CenterIDCenterIDWeight << set_value.size() <</*"\t"<<set_link.size()<<*/endl;
	for (set<pair<int, int>>::iterator link = set_link.begin(); link != set_link.end(); ++link) {
		ofs_CenterIDCenterIDWeight << link->first << "\t" << link->second << endl;
		outdg[link->first]++;
		indg[link->second] ++;
	}

	ofs_stat << "node #:\t" << set_value.size() << "\tedge #:\t" << set_link.size() << endl;
	for (int i = 0; i < set_value.size(); ++i)
		ofs_deg << indg[i] << "\t" << outdg[i] << endl;

	set_value.clear();
	set_link.clear();

}


//digraph
void DataProcessing(string dataset_str) {
	ifstream ifs_CenterCenter(dataset_str + "G.txt");//dataset : "G.txt"	
	ofstream ofs_ObjectIDObjectIDMatrix(dataset_str + "relationshipmatrix.txt");
	ofstream ofs_CenterIDCenter(dataset_str + "id_node.txt");
	ofstream ofs_deg(dataset_str + "in_out_deg.txt");
	ofstream ofs_stat(dataset_str + "stat.txt");
	ofstream ofs_G_rev(dataset_str + "G_rev.txt");
	WeightGraphOfCenterIDCenterID(ifs_CenterCenter, ofs_ObjectIDObjectIDMatrix, ofs_CenterIDCenter, ofs_deg, ofs_stat, ofs_G_rev);
	ifs_CenterCenter.close();
	ofs_ObjectIDObjectIDMatrix.close();
	ofs_CenterIDCenter.close();
	ofs_deg.close();
	ofs_stat.close();
}

//undigraph
void CheckUndirectGraph(string dataset_str) {
	ifstream ifs_CenterCenter(dataset_str + "G.txt");//dataset : "G.txt"
	ofstream ofs_ObjectIDObjectIDMatrix(dataset_str + "relationshipmatrix.txt");
	ofstream ofs_CenterIDCenter(dataset_str + "id_node.txt");
	ofstream ofs_deg(dataset_str + "in_out_deg.txt");
	ofstream ofs_stat(dataset_str + "stat.txt");

	CheckUDG(ifs_CenterCenter, ofs_ObjectIDObjectIDMatrix, ofs_CenterIDCenter, ofs_deg, ofs_stat);

	ifs_CenterCenter.close();
	ofs_ObjectIDObjectIDMatrix.close();
	ofs_CenterIDCenter.close();
	ofs_deg.close();
	ofs_stat.close();

}



int main()
{
	clock_t start = clock();
	DataProcessing(dataset_str);
	//CheckUndirectGraph(string dataset_str);
	clock_t end = clock();
	cout << "Cost" << (double)(end - start) / CLOCKS_PER_SEC << "s" << endl;
	system("pause");
	return 0;

}
