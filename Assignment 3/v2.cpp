// Function Names: FunctionName
// Variable Name: variable_name

#include<bits/stdc++.h>
#include <mpi.h>
#include <chrono>
#include <omp.h>

using namespace std;

float CosineDistance(vector<float> v1, vector<float> v2){
	float mod_x = 0;
	float mod_y = 0;
	float num = 0;
	for(int i = 0; i < v1.size(); i++){
		mod_x = mod_x + v1[i]*v1[i];
		mod_y = mod_y + v2[i]*v2[i];
		num = num + v1[i]*v2[i];
	}
	float sim = num / (sqrt(mod_x) * sqrt(mod_y));
	return 1 - sim;
}

priority_queue<pair<float, int>> SearchLayer(int k, vector<float> &query, priority_queue<pair<float, int>> &candidates, vector<int> &indptr, vector<int> &index, vector<int> &level_offset, int level, unordered_map<int, bool>& visited, vector<vector<float>>& vect) {
	priority_queue<pair<float, int>> topk(candidates);
	while(candidates.size() > 0) {
		pair<float, int> pair = candidates.top(); candidates.pop();
		int ep = pair.second;
		int start = indptr[ep] + level_offset[level];
		int end = indptr[ep] + level_offset[level + 1];
		for(int i = start; i < end; i++) {
			int px = index[i];
			if(visited[px] || px == -1) {
				continue;
			}
			
			visited[px] = true;
			float dist = CosineDistance(query, vect[px]);
			if(dist > topk.top().first && topk.size() == k){
				continue;
			}
			topk.push({dist, px});
			if(topk.size() > k){
				topk.pop();
			}
			candidates.push({dist, px});
		}
	}
	return topk;
}

priority_queue<pair<float, int>> QueryHNSW(int k, vector<float> &query, int ep, vector<int> &indptr, vector<int> &index, vector<int> &level_offset, int max_level, vector<vector<float>> &vect){
	priority_queue<pair<float, int>> topk;
	topk.push({CosineDistance(query, vect[ep]), ep});
	unordered_map<int, bool> visited;
	visited[ep] = true;
	int L = max_level;
	for(int level = max_level; level >= 0; level--){
		topk = SearchLayer(k, query, topk, indptr, index, level_offset, level, visited, vect);
	}
	return topk;
}

int main(int argc, char *argv[]) {
	auto begin = std::chrono::high_resolution_clock::now();
	string out_path = argv[1];
	int k = stoi(argv[2]);
	string user_file = argv[3];
	string user_output_file = argv[4];

	int max_level, ep, lines, dim, index_size, indptr_size, level_offset_size;

	ifstream in;
	in.open(out_path + "/max_level.bin", ios::in | ios::binary);
	in.read((char *)&max_level, sizeof(int));
	in.close();

	in.open(out_path + "/ep.bin", ios::in | ios::binary);
	in.read((char *)&ep, sizeof(int));
	in.close();

	in.open(out_path + "/extra.bin", ios::in | ios::binary);
	in.read((char *)&lines, sizeof(int));
	in.read((char *)&dim, sizeof(int));
	in.read((char *)&index_size, sizeof(int));
	in.read((char *)&indptr_size, sizeof(int));
	in.read((char *)&level_offset_size, sizeof(int));
	in.close();

	int id;
	float f;
	in.open(out_path + "/index.bin", ios::in | ios::binary);
	vector<int> index;
	for(int i = 0; i < index_size; i++){
		
		in.read((char *)&id, sizeof(int));
		index.push_back(id);
	}
	in.close();

	in.open(out_path + "/indptr.bin", ios::in | ios::binary);
	vector<int> indptr;
	for(int i = 0; i < indptr_size; i++){
		in.read((char *)&id, sizeof(int));
		indptr.push_back(id);
	}
	in.close();

	in.open(out_path + "/level.bin", ios::in | ios::binary);
	vector<int> level;
	for(int i = 0; i < lines; i++){
		in.read((char *)&id, sizeof(int));
		level.push_back(id);
	}
	in.close();

	in.open(out_path + "/level_offset.bin", ios::in | ios::binary);
	vector<int> level_offset;
	for(int i = 0; i < level_offset_size; i++){
		in.read((char *)&id, sizeof(int));
		level_offset.push_back(id);
	}
	in.close();

	in.open(out_path + "/vect.bin", ios::in | ios::binary);
	vector<vector<float>> vect;
	for(int i = 0; i < lines; i++){
		vector<float> v;
		for(int j = 0; j < dim; j++){
			in.read((char *)&f, sizeof(int));
			v.push_back(f);
		}
		vect.push_back(v);
	}
	in.close();
	
	vector<vector<float>> user;
	in.open(user_file);
	while(in >> f) {
		vector<float> v;
		v.push_back(f);
		for(int iter = 0; iter < dim - 1; iter++){
			in >> f;
			v.push_back(f);
		}
		user.push_back(v);
	}
	in.close();

	int rank, size;
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int chunk = user.size()/size;
    int start = rank*chunk;
    int end;

    if(rank == size - 1){
        end = user.size();
    }
    else{
        end = start + chunk;
    }

    int ans[(end - start)*k];

    ofstream fout;
    int send = 1;
    int receive = 1;
    if(rank == 0){
        fout.open("output.dat", ios::binary | ios::out);
        fout.seekp(user.size()*k*4-1);
        fout.write("",1);
        fout.close();
        for(int q = 1; q < size; q++){
            MPI_Send(&send, 1, MPI_INT, q, 0, MPI_COMM_WORLD);
        }
    }
    else{
        MPI_Recv(&receive, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    fout.open("output.dat", ios::binary | ios::out | ios::in);
    #pragma omp parallel
    {
    	#pragma omp for
    	for(int i = start; i < end; i++){
	    	vector<float> query = user[i];
	    	priority_queue<pair<float, int>> topk = QueryHNSW(k, query, ep, indptr, index, level_offset, max_level, vect);
			
			// cout << "Rank: " << rank << " User Id: " << i << " Thread Id: " << omp_get_thread_num() << endl;
			// for(int i = 0; i < k; i++){
			// 	pair<float, int> pair = topk.top();
			// 	topk.pop();
			// 	for(int j = 0; j < vect[pair.second].size(); j++){
			// 		cout << vect[pair.second][j] << " ";
			// 	}
			// 	cout << endl;
			// }
			
			int cur = k - 1;
			for(int j = 0; j < k; j++){
				pair<float, int> p = topk.top();
				topk.pop();
				ans[(i - start)*k + cur] = p.second;

				cur--;
			}
	    }
	}

	fout.seekp(4*start*k, ios::beg);
	int temp = 0;
	for(int i = start; i < end; i++){
		for(int j = 0; j < k; j++){
			int peer = ans[temp];
			fout.write((char*) &peer, sizeof(int));
			temp++;
		}
	}
	fout.close();
	
    
    if(rank == 0){
    	for(int q = 1; q < size; q++){
    		MPI_Recv(&receive, 1, MPI_INT, q, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    	}
    	in.open("output.dat", ios::in | ios::binary);
    	ofstream fout1;
    	fout1.open(argv[4]);
    	for(int i = 0; i < user.size(); i++){
    		for(int j = 0; j < k; j++){
    			int peer;
    			in.read((char *)&peer, sizeof(int));
    			fout1 << peer;
    			if(j != k - 1){
    				fout1 << " ";
    			}
    		}
    		fout1 << endl;
    	}
    	in.close();
    	remove("output.dat");
    }
    else{
    	MPI_Send(&send, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }

	MPI_Finalize();
	auto end_time = std::chrono::high_resolution_clock::now();
	cout << "Time taken: " << (1e-6 * (std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - begin)).count()) << "ms" << endl;
}