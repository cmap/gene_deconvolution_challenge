#include<bits/stdc++.h>
#include<dirent.h>
#include <thread>

using namespace std;

#define MAX 1000
//namespace fs = std::filesystem;

vector<string> fnames; // ok for global
vector<vector<int> > ret[MAX]; // each have own space to save output
map<int, double> fi[399]; // own space for each
map<int, int> mp0, mp1; // barcode --> gene pair

vector<string> readfile(string fname) {
	vector<string> ret;
	string line;
	ifstream myfile (fname);
	if (myfile.is_open()) {
		while(getline(myfile, line)) {
			ret.push_back(line);
		}
		myfile.close();
	}
	else{
		cout << " Unable to open file\n";
	}
	return ret;
}

void f(int i) {
	int FN = i;
	string fname = fnames[i];
	bool is_dpk = false;
	if (fname[7] == 'D') {
		is_dpk = true;
	}
//	cout << fname << " -- here\n" ;
	vector<string> lines = readfile(fname);
	vector<pair<int, int> > barcode_fi;
	int n = lines.size();
	for(int i =1; i < n; i++) {
		string line = lines[i];
		int sz = line.size();

		int num[2];
		num[0] = 0; num[1] = 0;
		int k = 0;
		int j = 0;
		while(j < sz) {
			if(line[j] != '\t') {
				num[k]*=10;
				num[k] += (line[j] - '0');
			}else{
				k=1;
			}
			j++;
		}
		if(num[1]==0) {	}
		else barcode_fi.push_back(make_pair(num[0], num[1]));
	}
	n = barcode_fi.size();	
	sort(barcode_fi.begin(), barcode_fi.end());
	for(int i = 0; i < barcode_fi.size(); i++) {
	//	cout << barcode_fi[i].first << "," << barcode_fi[i].second << endl;
	}
	int start = 0;
	int tot_one_peaks = 0;
	for(int i = 1; i < n; i++){
		if((barcode_fi[start].first != barcode_fi[i].first) || (i == n-1)) {
			int end = i - 1;
			int tot = end - start + 1;
			int median = barcode_fi[(start+end)/2].second;
			//cout << tot << "," << start << "..." << end << endl;
			if(tot == 0) cout << "MMMMMMM\n";


/*
			// Single Peak Identify
			int no_of_vals_in_1_bar = 8;
			int st_val = barcode_fi[start].second;
			int end_val = barcode_fi[end].second;
			int per_val = ((end_val - st_val) / tot);
			int bar_width = per_val * no_of_vals_in_1_bar;
			map<int,int> bar_height;
			for(int j = start; j <= end; j++) {
				int bar_no = (barcode_fi[j].second - st_val) / bar_width;
				bar_height[bar_no] +=1;
			}
			bool one_peak = true;
			int k = 0;
			while(bar_height[k] <= bar_height[k+1]) k++;	
			int max_bar_no = (end_val - st_val) / bar_width;
			while((k < max_bar_no+1) && (bar_height[k] >= bar_height[k+1])) k++;
			if(k != max_bar_no+1) one_peak = false;

			if(one_peak) {
				tot_one_peaks++;
				
				
				cout << " -----------------  \n";
				for(int j = start; j <= end; j++){
					cout << barcode_fi[j].second <<",";
				}
				cout << endl;
				for(int j = 0; j <= max_bar_no; j++){
					cout << bar_height[j] <<",";
				}
				cout << " ------------------- \n";
				
				int gene1 = mp0[barcode_fi[start].first];
        	    int gene2 = mp1[barcode_fi[start].first];
           		start = i;
				int median = barcode_fi[(start+end)/2].second;
            	fi[FN][gene1] = 1.3*median; 
            	fi[FN][gene2] = 1.3*median;
				continue;
			}
*/
	 
		    double pc1[] = { 0.14,  0.25, 0.3};
            double pc2[] = {  0.76, 0.86};
			int range = barcode_fi[start + (int)(pc1[0]*tot) ].second;
			vector<pair<double, pair<int, int> > > candidates;
			for(int t = 0; t < 3; t++) {
				for(int t2=0; t2 < 2; t2++) { // 1-->2 [D:673.4k, lit:784.5k]
					// x, 2x, tot 3x, x/2/3x is 16.6%, 66.66%
					int v1 = barcode_fi[start + (int)(pc1[t]*tot) ].second;
					int v2 = barcode_fi[start + (int)(pc2[t2]*tot) ].second;
					double sz1, sz2;
					double mv1, mv2;
					double sum1 = 0;
					double sum2 = 0;
					for(int z = 0; z < 2; z++){// 2-->3 525.9k [D:672.5k, lit: 782.7k]
												// -->4 525.7k (due to time) [D: 672.8k, lit:783.2k]
						// v1 and v2 are guesses
						// lets refine our guesses (also known as kmeans algo)
						vector<int> gr1, gr2;
						for(int j = start; j <= end; j++){
							int fi = barcode_fi[j].second;
							double d1 = abs(v1-fi);
							double d2 = abs(v2-fi);
						/*	
							if(d1 <= range) gr1.push_back(fi);
							if(d2 <= range) gr2.push_back(fi);
						*/	
					/*	
							if(d1/(v1+1) < d2/(v2+1)) gr1.push_back(fi);
							else gr2.push_back(fi);
					*/	
							
							if(d1 < 0.82*d2) {gr1.push_back(fi); sum1 += fi;} // 1.0->0.9, 520.2k, 522.9k [DPK:670k, lit:780k],
																			// ->0.8, 524.7k [D:673k, lit:780k]
																			// -> 0.82 525.5k [D:672.7k, lit:781.1k]
		
                            else {gr2.push_back(fi); sum2+=fi;}
							
							
						}
						sz1 = gr1.size();
						
						if (sz1 == 0) { 
							cout << tot << ","<<start <<"," << end<< " ====== \n ";
							for(int j = start; j <= end; j++) cout << barcode_fi[j].second<< endl;
						cout << "v1="<<v1<<endl;
						}

						sz2 = gr2.size();
/*
						//cout << "SZ " << sz1 << ", " << sz2 << endl; 	
						double cnt1 = 0;
						double cnt2 = 0; 
						map<int, int> mp1, mp2;
						map<int, int> cnts;
						set<int> val_set;
						int pc10 = gr1[(int)(sz1/(2))]/10;
						if(pc10 == 0) pc10 = 1;
						for(int x = 0; x < sz1; x++) {
							int v = (gr1[x] / pc10)*pc10;
							val_set.insert(v);
							cnts[v]+=1;
						}
						vector<int> vals(val_set.begin(), val_set.end());
						sort(vals.begin(), vals.end());
						v1 =  gr1[(int)(sz1/(2))]; // default
						int max_cnt = 0;
						// Find first PEAK
						for(int x = 0; x < vals.size(); x++){
							//int prev = vals[x-1];
							int now = vals[x];
							//int next = vals[x+1];

//							cout << "fi:" << now << ", cnt:"<< cnts[now] << endl;
							
							if(cnts[now] > 2) // atleast 3 height peaks
							if ((cnts[prev] < cnts[now]) && (cnts[now] > cnts[next])) {
								v1 = now;
								break;
							}
							if (cnts[now] > max_cnt){
								max_cnt = cnts[now];
								v1 = now;
							}
						} 
						max_cnt = 0;
						map<int, int> cnts2;
						set<int> val_set2;
						pc10 = gr2[(int)(sz2/(2))]/10;
						if(pc10 == 0) pc10 = 1;
						for(int x = 0; x < sz2; x++) {
							int v = (gr2[x] / pc10)*pc10;
							val_set2.insert(v);
							cnts2[v]+=1;
						}
						vector<int> vals2(val_set2.begin(), val_set2.end());
						sort(vals2.begin(), vals2.end());
						v2 =  gr2[(int)(sz2/(2))]; // default

						// Find first PEAK
						for(int x = 0; x < vals2.size(); x++){
							//int prev = vals2[x-1];
							int now = vals2[x];
							//int next = vals2[x+1];
							
							if(cnts2[now] > 2) // atleast 3 height peaks
							if ((cnts2[prev] < cnts2[now]) && (cnts2[now] > cnts2[next])) {
								v2 = now ;
								break;
							}
							if (cnts2[now] > max_cnt){
                                max_cnt = cnts2[now];
                                v2 = now;
                            }
						} 

*/
						v1 = gr1[(int)(sz1/2)];
						v2 = gr2[(int)(sz2/2)];
						mv1 = sum1/sz1; //gr1[(int)(sz1/2)];
						mv2 = sum2/sz2; //gr2[(int)(sz2/2)];
						/*
						// peak height
						for(int x = 0; x < sz1; x++) {
							if(v1 == gr1[x]) cnt1++;
							mp1[gr1[x]]++;
						}
						v2 = gr2[(int)(sz2/2)];
						//peak height
						for(int x = 0; x < sz2; x++) {
							if(v2 == gr2[x]) cnt2++;
							mp2[gr2[x]]++;
						}
						*/


						/*
						// count vals
						for(map<int, int>::iterator it = mp1.begin(); it != mp1.end(); it++){
							//cout << it->first << "->" << it->second << endl;
							for(int w = 0; w < it->second; w++) {
								cout <<"||";
							}
							cout << endl;
						}
						cout <<" ---------------------------- \n";
						 for(map<int, int>::iterator it = mp2.begin(); it != mp2.end(); it++){
                            //cout << it->first << "->" << it->second << endl;
                            for(int w = 0; w < it->second; w++) {
                                cout <<"||";
                            }
                            cout << endl;
                        }*/
						//cout << cnt1 << " " << cnt2 << endl;	
						
						/* Experiment 
						double r1 = sz1/sz2;
                    	if (r1 < 1.5) r1 -= 0.5;
                    	double r2 = sz2/sz1;
                    	if(r2 < 1.5) r2 -= 0.5;
                    	if(r1 > 2.0) r1 -= 0.6;
                    	if(r2 > 2.0) r2 -= 0.6;
                    	candidates.push_back(make_pair(abs(2.0 - r1), make_pair(v1,v2)));
                    	candidates.push_back(make_pair(abs(2.0 - r2), make_pair(v2,v1)));	
						*/
					}
					double r1 = sz1/sz2;
					if (r1 < 1.7) r1 -= 0.5; // 1.5-->1.7, 519.5k to 520.2k
					double r2 = sz2/sz1;
					if(r2 < 1.5) r2 -= 0.5;
					if(r1 > 2.0) r1 -= 0.3;
					if(r2 > 2.0) r2 -= 0.3; 
					candidates.push_back(make_pair(abs(2.3 - r1), make_pair(v1,v2))); // 2.0 --> 2.3, 505k to 518.5k
					candidates.push_back(make_pair(abs(2.1 - r2), make_pair(v2,v1))); // 2.0 --> 2.1, 518.5k to 519.5k 
						//cout <<"Z = " <<z << " -- " << v1 << "," << v2 << " ["<<sz1<<","<<sz2<<"]\n";  
				
					
				}
				
			}
			sort(candidates.begin(), candidates.end());
			int higher = 0.5*(candidates[0].second.first + candidates[0].second.first); // higher membership
			int lower = 0.5*(candidates[0].second.second + candidates[0].second.second) ; // lower membership
			
			/************************************************************************
			reduce higher proportioned one if its value is less (see 199 barcode)
			like: 165,168,170,177,222,223,224,227,300,302,302,1201,1256,1376,1454
			[but reduce didnt work, somehow increasing gives better answer(?)]

			increase "lower" if its like:
			122,124,166,176,1122,1123,1234,1345,1348,1455,1456,1477,1488,1655,1766

			So, essesntially low FIs (~122..600) are increased by 30-40%
			*************************************************************************/
			if (is_dpk) {
				if (higher < lower) { higher = (int)(1.3*higher); lower = (int)(1.02*lower);}
				else {lower = (int)(1.4*lower);  higher =(int) (1.1*higher);}
			}else{
				if (higher < lower) { higher = (int)(1.3*higher); lower = (int)(.95*lower);}
                else {lower = (int)(1.4*lower);  }	
			
			}
			//if (candidates[0].first > 1) {cout << "[>0.9]\n" ; higher =  lower = median;} 
			//cout << candidates[0].first << "  -- " << higher << " " << lower << endl;
			int gene1 = mp0[barcode_fi[start].first];
			int gene2 = mp1[barcode_fi[start].first];
			start = i;
			fi[FN][gene1] = 0.5*lower ; // 1->0.5 .. 527k, [D:677k, lit:785.2k]
			fi[FN][gene2] = 0.5*higher ;	

/*
			if(one_peak) {
				fi[FN][gene1] = (lower+higher) * 0.5;
           	    fi[FN][gene2] = fi[FN][gene1];
			}
*/

		}	
	}
//	cout << "TOT ONE PEAKS: " << tot_one_peaks << endl;

}

int main(int argc, char* argv[]) {
	vector<string> args;
	for(int i = 0; i < argc; i++){
		stringstream ss;
        ss << argv[i];
        string s = ss.str();
		args.push_back(s);
	}
	ifstream myfile("barcode_to_gene_map.txt");
	vector<string> vec;
	vector<string> colname;
	vector<int> gene_list;
	string line;
	if (myfile.is_open()) {
		while(getline(myfile, line)) {
				vec.push_back(line);
		}
	}else{
		cout << " Not open map file\n";
	}
    
	for(int i = 1; i < vec.size(); i++) {
		string line = vec[i];
		int sz = line.size();
		int k = 0;
		int num[3];
		num[0] = num[1] = num[2] = 0;
		for (int j = 0; j < sz; j++) {
			if(line[j] != '\t') {
                num[k]*=10;
                num[k] += (line[j] - '0');
            }else{
                	k+=1;
            }
            
		}
		if ((num[0] != 11) && (num[0] != 499))
			gene_list.push_back(num[1]);
		cout << num[0] << "," << num[1] << "," <<  num[2] << endl;
		if(num[2] == 1)
			mp1[num[0]] = num[1]; // map<int, int> mp[2];
		else mp0[num[0]] = num[1];

	}
	string s = args[2]; //"/Users/manish/Work/topcoder-solutions/cmap/competitor_pack_v2/input/DPK.CP001_A549_24H_X1_B42";
	stringstream ss;
	DIR *dir;
	struct dirent *ent;
	cout << s << " [opening dir] " <<  endl;
	if ((dir = opendir (s.c_str())) != NULL) {
	  while ((ent = readdir (dir)) != NULL) {
		stringstream ss;
		ss << ent->d_name;
		if(ss.str().size() > 5){
			fnames.push_back(s+"/"+ss.str());
		}
	  }
	  closedir (dir);
	} else {
	  /* could not open directory */
	  perror ("");
	  return EXIT_FAILURE;
	}
	sort(fnames.begin(), fnames.end());
	for(int i = 0; i < fnames.size(); i++) {
	
		//cout << "calling " << i << endl;
		int len = fnames[i].size();
        colname.push_back(fnames[i].substr(len - 7, 3));
	//	f(i);
	}
	
	for(int i = 0; i < fnames.size(); i+=6+6){
		thread th1(f, i);
		thread th2(f, i+1);
		thread th3(f, i+2);
		thread th4(f, i+3);
		thread th5(f, i+4);
		thread th6(f, i+5);
		 thread th11(f, i+6);
        thread th21(f, i+7);
        thread th31(f, i+8);
        thread th41(f, i+9);
        thread th51(f, i+10);
        thread th61(f, i+11);
		th1.join();
		th2.join();
		th3.join();
		th4.join();
		th5.join();
		th6.join();
		   th11.join();
        th21.join();
        th31.join();
        th41.join();
        th51.join();
        th61.join();
	}
	
	cout << gene_list.size() << endl;	
	cout <<"fnames szi" << fnames.size();
	ofstream myfile2;
	string outfile = args[4]+"/"+args[8]+".gct";
    cout <<"[Wriiten] " << outfile << endl;
	myfile2.open (outfile.c_str());
	myfile2 << "#1.3\n";
	myfile2 << "976\t"<<fnames.size()<<"\t0\t0\n";
	myfile2 << "id\t";
	for(int i = 0; i < colname.size(); i++) {
		myfile2 << colname[i];
		if(i != colname.size()-1) myfile2 << "\t";
	}
	myfile2 <<"\n";

	for(int j = 0; j < gene_list.size(); j++) {
		  myfile2 << gene_list[j] <<"\t";
		   for(int i = 0; i < fnames.size(); i++){
				
				myfile2 << fi[i][gene_list[j]];
				if(i != fnames.size() - 1) myfile2 << "\t";
			}myfile2 << endl;
	}
	myfile2.close();

	return 0;
}
