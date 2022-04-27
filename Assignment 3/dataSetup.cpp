#include<bits/stdc++.h>

using namespace std;

int main(int argc, char *argv[]) {
	string in_path = argv[1];
	string out_path = argv[2];

	ifstream in;
	int d;
	float f;
	int lines, dim, index_size, indptr_size, level_offset_size;

	ofstream max_level;
	max_level.open(out_path + "/max_level.bin", ios::out | ios::binary);

	in.open(in_path + "/max_level.txt");
	while(in >> d) {
		max_level.write((char*)&d, sizeof(int));
	}
	in.close();

	ofstream ep;
	ep.open(out_path + "/ep.bin", ios::out | ios::binary);

	in.open(in_path + "/ep.txt");
	while(in >> d) {
	  ep.write((char*)&d, sizeof(int));
	}
	in.close();

	ofstream level;
	level.open(out_path + "/level.bin", ios::out | ios::binary);

	in.open(in_path + "/level.txt");
	while(in >> d) {
	  level.write((char*)&d, sizeof(int));
	}
	in.close();

	ofstream index;
	index.open(out_path + "/index.bin", ios::out | ios::binary);

	in.open(in_path + "/index.txt");
	while(in >> d) {
	  index.write((char*)&d, sizeof(int));
	  index_size++;
	}
	in.close();

	ofstream indptr;
	indptr.open(out_path + "/indptr.bin", ios::out | ios::binary);

	in.open(in_path + "/indptr.txt");
	while(in >> d) {
	  indptr.write((char*)&d, sizeof(int));
	  indptr_size++;
	}
	in.close();

	ofstream level_offset;
	level_offset.open(out_path + "/level_offset.bin", ios::out | ios::binary);

	in.open(in_path + "/level_offset.txt");
	while(in >> d) {
	  level_offset.write((char*)&d, sizeof(int));
	  level_offset_size++;
	}
	in.close();

	lines;
	string line;
	in.open(in_path + "/vect.txt");
	for(lines = 0; std::getline(in, line); lines++);
	in.close();
	
	int n = 0;
	ofstream vect;
	vect.open(out_path + "/vect.bin", ios::out | ios::binary);

	in.open(in_path + "/vect.txt");
	while(in >> f) {
	  vect.write((char*)&f, sizeof(float));
	  n++;
	}
	in.close();

	dim = n/lines;
	ofstream extra;
	extra.open(out_path + "/extra.bin", ios::out | ios::binary);
	extra.write((char*)&lines, sizeof(int));
	extra.write((char*)&dim, sizeof(int));
	extra.write((char*)&index_size, sizeof(int));
	extra.write((char*)&indptr_size, sizeof(int));
	extra.write((char*)&level_offset_size, sizeof(int));
}