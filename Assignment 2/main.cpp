#include <string>
#include <mpi.h>
#include <assert.h>
#include "randomizer.hpp"
#include <bits/stdc++.h>
#include <chrono>

using namespace std; 

int SWAP_INT32(int x){
    return ((x) >> 24) | (((x) & 0x00FF0000) >> 8) | (((x) & 0x0000FF00) << 8) | ((x) << 24);
}

void wtf(int tid, int size, int num_nodes, int num_steps, int num_walks, int num_rec, vector<vector<int>> adj, Randomizer r){
    ofstream fout;
    int send = 1;
    int receive = 1;
    if(tid == 0){
        fout.open("output.dat", ios::binary | ios::out);
        fout.seekp(num_nodes*(2*num_rec+1)*4-1);
        fout.write("",1);
        fout.close();
        for(int q = 1; q < size; q++){
            MPI_Send(&send, 1, MPI_INT, q, 0, MPI_COMM_WORLD);
        }
    }
    if(tid != 0){
        MPI_Recv(&receive, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    fout.open("output.dat", ios::binary | ios::out | ios::in);
    int chunk = num_nodes/size;
    int start = tid*chunk;
    int end;

    if(tid == size - 1){
        end = num_nodes;
    }
    else{
        end = start + chunk;
    }
    int s = 0;
    string s1 = "NULL";
    for(int i = start; i < end; i++){
        vector<int> m(num_nodes, 0);
        for(int j = 0; j < adj[i].size(); j++){
            for(int k = 0; k < num_walks; k++){
                int start_node = adj[i][j];
                int cur_node = adj[i][j];
                for(int l = 0; l < num_steps; l++){
                    if(adj[cur_node].size() == 0){
                        cur_node = start_node;
                    }
                    else{
                        int next_step = r.get_random_value(i);
                        if(next_step < 0){
                            cur_node = start_node;
                        }
                        else{
                            int node_id = next_step%(adj[cur_node].size());
                            cur_node = adj[cur_node][node_id];
                        }
                    }
                    m[cur_node]++;
                }
            }
        } 
        for(int j = 0; j < adj[i].size(); j++){
            m[adj[i][j]] = 0;
        }
        m[i] = 0;
        vector<pair<int, int>> q;
        for(int j = 0; j < m.size(); j++){
            if(m[j] != 0){
                q.push_back({m[j], -1*j});
            }
        }
        sort(q.begin(), q.end(), greater<pair<int,int>>());
        fout.seekp(4*i*(2*num_rec + 1), ios::beg);
        int j = 0;
        int outdegree = SWAP_INT32(adj[i].size());
        fout.write((char*) &outdegree, 4);
        while(j < num_rec && j < q.size()){
            int node = -1*q[j].second;
            int big_node = SWAP_INT32(node);
            int big_score = SWAP_INT32(q[j].first);
            fout.write((char*) &big_node, 4);
            fout.write((char*) &big_score, 4);
            j++;
        }
        while(j < num_rec){
            fout.write(s1.c_str(), 4);
            fout.write(s1.c_str(), 4);
            j++;
        }
    }
    fout.close();
}

int main(int argc, char* argv[]){
    assert(argc > 8);
    std::string graph_file = argv[1];
    int num_nodes = std::stoi(argv[2]);
    int num_edges = std::stoi(argv[3]);
    float restart_prob = std::stof(argv[4]);
    int num_steps = std::stoi(argv[5]);
    int num_walks = std::stoi(argv[6]);
    int num_rec = std::stoi(argv[7]);
    int seed = std::stoi(argv[8]);

    vector<vector<int>> adj(num_nodes);
    unsigned char buf[4];
    FILE *fp = fopen(graph_file.c_str(),"rb");
    int i = 0, a;
    while (fread(buf, 4, 1, fp) != 0){
        int num = buf[0]<<24 | buf[1]<<16 | buf[2]<<8 | buf[3];
        if(i%2 == 0){
            a = num;
        }
        else{
            adj[a].push_back(num);
        }
        i+=1;
    }

    int rank, size;
    Randomizer random_generator(seed, num_nodes, restart_prob);

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    wtf(rank, size, num_nodes, num_steps, num_walks, num_rec, adj, random_generator);
    
    MPI_Finalize();
}