#include <iostream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <string>
#include <cstring>

using namespace std;

#define MAX_LENGTH 1000000

class Prune {
private:
	string bamfile;
	string table;
	string refSeq;
	unordered_map< string, unordered_map <string, string>> pairdb;
	unordered_map<string, long > ctgdb;
	unordered_map<string, long > bamdb;
	unordered_map<string, long> removedb;

	bool Split(string source, string delim, vector<string>&target);
	void SortCtg(string ctg1, string ctg2, string &sctg1, string &sctg2);
	void output_nonBest(ofstream &fout2, vector<string>&lines);
public:
	Prune() {
		bamfile = "";
		table = "";
		refSeq = "";
	}
	Prune(string bamfile, string table, string refSeq) {
		this->bamfile = bamfile;
		this->table = table;
		this->refSeq = refSeq;
	}
	~Prune(){}
	void SetParameter(string bamfile, string table, string refSeq);
	bool GeneratePairsAndCtgs();
	bool CreateLogAndRemovedbFiles();
	bool CreatePrunedBam();
};

//Split string by delimiter
bool Prune::Split(string source, string delim, vector<string>&target) {
	target.clear();
	char *p;
	p = strtok(const_cast<char*>(source.c_str()), delim.c_str());
	if (!p) {
		return false;
	}
	while (p) {
		target.push_back(p);
		p = strtok(NULL, delim.c_str());
	}
	return true;
}

//Sort 2 ctgs by ascii, and assign to sctgs
void Prune::SortCtg(string ctg1, string ctg2, string &sctg1, string &sctg2) {
	if (ctg1.compare(ctg2) < 0) {
		sctg1 = ctg1;
		sctg2 = ctg2;
	}
	else {
		sctg1 = ctg2;
		sctg2 = ctg1;
	}
}

//Create contents for removedb_nonBest.txt, and add strings for remove to a map: removedb
void Prune::output_nonBest(ofstream &fout, vector<string>&lines) {
	vector<string> data;
	unordered_map<string, unordered_map<string, string>> hashdb;
	vector<string> rnamedb;
	hashdb.clear();
	for (long i = 0; i < lines.size(); i++) {
		Split(lines[i], "\t", data);
		if (hashdb.count(data[0]) == 0) {
			hashdb[data[0]]["retain"] = data[1];
			hashdb[data[0]]["num"] = data[2];
		}
		else if (hashdb.count(data[0]) > 0 && stol(data[2]) > stol(hashdb[data[0]]["num"])) {
			hashdb[data[0]]["retain"] = data[1];
			hashdb[data[0]]["num"] = data[2];
		}
	}
	for (long i = 0; i < lines.size(); i++) {
		Split(lines[i], "\t", data);
		if (hashdb[data[0]]["retain"] == data[1]) {
			continue;
		}
		else {
			fout << data[0] << "\t" << data[1] << "\t" << data[2] << "\tremove\t" << data[3] << "\n";
			Split(data[3], ",", rnamedb);
			for (long k = 0; k < rnamedb.size(); k++) {
				removedb[rnamedb[k]]++;
			}
		}
	}
	lines.clear();
}

void Prune::SetParameter(string bamfile, string table, string refSeq) {
	this->bamfile = bamfile;
	this->table = table;
	this->refSeq = refSeq;
}


//Read bamfiles and read them by samtools, then create pairdbs and ctgdbs;
bool Prune::GeneratePairsAndCtgs() {
	if (bamfile == "" || table == "" || refSeq == "") {
		return false;
	}
	else {
		string tempsam;
		string cmd;
		vector<string>data;
		string ctg1, ctg2;
		string sctg1, sctg2;

		cmd = "samtools view " + bamfile;
		FILE * fp = popen(cmd.c_str(), "r");
		if (fp) {
			char buffer[MAX_LENGTH];
			while (!feof(fp)) {
				if (fgets(buffer, sizeof(buffer), fp) != NULL) {
					tempsam = buffer;
					Split(tempsam, "\t", data);
					ctg1 = data[2];
					ctg2 = data[6];
					if (ctg2 == "=") {
						continue;
					}

					SortCtg(ctg1, ctg2, sctg1, sctg2);
					pairdb[sctg1][sctg2] += data[0] + ",";
					ctgdb[ctg1]++;
					ctgdb[ctg2]++;
				}
			}
			pclose(fp);
		}
		else{
			cout<<"Cannot open "<<bamfile<<"\n";
			return false;
		}
	}
}

//Create removedb_Allele.txt, removedb_nonBest.txt and log.txt;
bool Prune::CreateLogAndRemovedbFiles() {
	ifstream fin;
	ofstream fallele, fnonbest, flog;
	unordered_map<string, long>tempdb;
	vector<string>data;
	unordered_map<string, long >::iterator iter_ctgdb;
	vector<string>rnamedb;
	vector<string>lines;
	string key;
	string temp;
	string ctg1, ctg2;
	string sctg1, sctg2;
	long num_r;

	fallele.open("removedb_Allele.txt");
	fnonbest.open("removedb_nonBest.txt");
	flog.open("log.txt");
	if (fallele&&flog) {
		fin.open(table);
		if (fin) {
			while (getline(fin, temp)) {
				tempdb.clear();
				Split(temp, "\t", data);
				if (data.size() <= 3) {
					continue;
				}
				for (long i = 2; i < data.size() - 1; i++) {
					ctg1 = data[i];
					for (long j = i + 1; j < data.size(); j++) {
						ctg2 = data[j];
						SortCtg(ctg1, ctg2, sctg1, sctg2);
						key = sctg1 + "," + sctg2;
						tempdb[key]++;
						if (pairdb.count(sctg1) > 0 && pairdb[sctg1].count(sctg2) > 0) {
							fallele << sctg1 << "\t" << sctg2 << "\t" << pairdb[sctg1][sctg2] << "\n";
							rnamedb.clear();
							Split(pairdb[sctg1][sctg2], ",", rnamedb);
							for (long k = 0; k < rnamedb.size(); k++) {
								removedb[rnamedb[k]]++;
							}
						}
					}
				}
				flog << ">" << temp << "\n";
				lines.clear();
				for (long i = 2; i < data.size(); i++) {
					ctg1 = data[i];
					for (iter_ctgdb = ctgdb.begin(); iter_ctgdb != ctgdb.end(); iter_ctgdb++) {
						ctg2 = iter_ctgdb->first;
						SortCtg(ctg1, ctg2, sctg1, sctg2);
						key = sctg1 + "," + sctg2;
						if (tempdb.count(key) > 0) {
							continue;
						}
						if (pairdb.count(sctg1) == 0) {
							continue;
						}
						if (pairdb[sctg1].count(sctg2) == 0) {
							continue;
						}
						rnamedb.clear();
						Split(pairdb[sctg1][sctg2], ",", rnamedb);
						num_r = rnamedb.size();
						temp= ctg2 + "\t" + ctg1 + "\t" + to_string(num_r) + "\t" + pairdb[sctg1][sctg2];
						lines.push_back(temp);
						flog <<temp<<"\n";
					}
				}
				//Create contents for removedb_nonBest.txt
				output_nonBest(fnonbest, lines);
			}
			fin.close();
		}
		else {
			cout<<"Cannot open "<<table<<"\n";
			return false;
		}
		fallele.close();
		fnonbest.close();
		flog.close();
		return true;
	}
	else {
		return false;
	}
}

//Directly to create prunning.bam through pipe with samtools
bool Prune::CreatePrunedBam() {
	FILE *fp, *fout;
	string cmd;
	unordered_map<string, long >::iterator iter_bamdb;
	vector<string>data;
	long num_of_remove_reads;
	char buffer[MAX_LENGTH];
	string temp;
	string rname,ctg2;

	cmd = "samtools faidx " + refSeq;
	system(cmd.c_str());
	num_of_remove_reads = removedb.size();
	cout << "Removing " << num_of_remove_reads << " reads\n";
	
	//cmd = "samtools view -bt " + refSeq + ".fai - >> prunning.bam";
	//Create a pipe for bamfile
	//fout = popen(cmd.c_str(), "w");
	fout = fopen("prunning.sam","w");

	cmd = "samtools view " +  bamfile;
	fp = popen(cmd.c_str(), "r");
	while (!feof(fp)) {
		if (fgets(buffer, sizeof(buffer), fp)) {
			temp = buffer;
			memset(buffer,0,sizeof(buffer));
			data.clear();
			Split(temp, "\t", data);
			rname = data[0];
			ctg2 = data[6];
			if (ctg2 == "*") {
				continue;
			}
			if (removedb.count(rname) == 0) {
				//temp += "\n";
				fputs(temp.c_str(), fout);
			}
		}
	}
	pclose(fp);
	//pclose(fout);
	fclose(fout);
	cmd = "samtools view -bt " + refSeq + ".fai prunning.sam > prunning.bam";
	fp = popen(cmd.c_str(), "r");
	pclose(fp);
}

int main(int argc, char* argv[]) {
	if (argc != 7) {
		cout << "************************************************************************\n";
		cout << "    Usage: ./prune -i Allele.ctg.table -b bamfile -r draft.asm.fasta\n";
		cout << "      -h : help and usage.\n";
		cout << "      -i : Allele.ctg.table\n";
		cout << "      -b : input bam file\n";
		cout << "      -r : draft.asm.fasta\n";
		cout << "************************************************************************\n";
	}
	else {
		string bamfile;
		string table;
		string refSeq;
		for (long i = 1; i < 7; i += 2) {
			if (strcmp(argv[i], "-i") == 0) {
				table = argv[i + 1];
				continue;
			}
			if (strcmp(argv[i], "-b") == 0) {
				bamfile = argv[i + 1];
				continue;
			}
			if (strcmp(argv[i], "-r") == 0) {
				refSeq = argv[i + 1];
				continue;
			}
		}
		Prune prune;
		prune.SetParameter(bamfile, table, refSeq);
		prune.GeneratePairsAndCtgs();
		prune.CreateLogAndRemovedbFiles();
		prune.CreatePrunedBam();
	}
	return 0;
}
