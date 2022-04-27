// TODO: Efficient Checking in Interpolation
// Function Names: functionName
// Variable Name: variable_name
// Filtering values not matching

#include<bits/stdc++.h>
#include <chrono>

using namespace std;

int pad = 255;

// x --> x-coordinate, y --> y-coordinate
__device__
bool inLimits(int x, int y, int x_max, int y_max){
	if(x >= 0 && x < x_max && y >= 0 && y < y_max){
		return true;
	}
	return false;
}

__device__
float sine(int angle) {
	if(angle == 0){
		return 0;
	}
	else if(angle == 45){
		return 0.70710678118;
	}
	return -0.70710678118;
}

__device__
float cosine(int angle) {
	if(angle == 0){
		return 1;
	}
	return 0.70710678118;
}

// point[0] --> x-coordinate, point[1] --> y-coordinate
// pivot[0] --> x-coordinate, pivot[1] --> y-coordinate
// res[0] --> x-coordinate, res[1] --> y-coordinate
__device__
void getRotatedCoordinate(float *res, int *point, int *pivot, int angle) {
	float s = sine(angle);
	float c = cosine(angle);

	float x_new = (point[0] - pivot[0])*c - (point[1] - pivot[1])*s;
	float y_new = (point[0] - pivot[0])*s + (point[1] - pivot[1])*c;

	res[0] = x_new + pivot[0];
	res[1] = y_new + pivot[1];
}

float sineHost(int angle) {
	if(angle == 0){
		return 0;
	}
	else if(angle == 45){
		return 1/sqrt(2);
	}
	return -1*(1/sqrt(2));
}

float cosineHost(int angle) {
	if(angle == 0){
		return 1;
	}
	return 1/sqrt(2);
}

// point[0] --> x-coordinate, point[1] --> y-coordinate
// pivot[0] --> x-coordinate, pivot[1] --> y-coordinate
// res[0] --> x-coordinate, res[1] --> y-coordinate

void getRotatedCoordinateHost(float *res, int point[], int pivot[], int angle) {
	float s = sineHost(angle);
	float c = cosineHost(angle);

	float x_new = (point[0] - pivot[0])*c - (point[1] - pivot[1])*s;
	float y_new = (point[0] - pivot[0])*s + (point[1] - pivot[1])*c;

	res[0] = x_new + pivot[0];
	res[1] = y_new + pivot[1];
}

__device__
void interpolate(int *X, int point[], int pivot[], int angle, float *interpolated, int x_m, int x_n){
	float v[2];
	getRotatedCoordinate(v, point, pivot, angle);
	float x = v[0], y = v[1];
	int x1 = floor(v[0]), y1 = floor(v[1]); // corresponds to (0, 0)
	int x2 = floor(v[0]), y2 = floor(v[1] + 1); // corresponds to (0, 1)
	int x3 = floor(v[0] + 1), y3 = floor(v[1]); // corresponds to (1, 0)
	int x4 = floor(v[0] + 1), y4 = floor(v[1] + 1); // corresponds to (1, 1)
	if(!inLimits(x1, y1, x_n, x_m) || !inLimits(x4, y4, x_n, x_m)){
		interpolated[0] = -1;
		return;
	}
	for(int i = 0; i < 3; i++){
		float z00 = X[y1*x_n*3 + x1*3 + i];
		float z01 = X[y2*x_n*3 + x2*3 + i];
		float z10 = X[y3*x_n*3 + x3*3 + i];
		float z11 = X[y4*x_n*3 + x4*3 + i];
		interpolated[i] = z00*(x4-x)*(y4-y) + z10*(x-x2)*(y2-y) + z01*(x3-x)*(y-y3) + z11*(x-x1)*(y-y1);
	}
}

// start[0] --> x-coordinate
// start[1] --> y-coordinate
// end[0] --> x-coordinate
// end[1] --> y-coordinate
__device__
float rmsd(int *X, int *Y, int start[], int end[], int angle, int q_m, int q_n, int x_m, int x_n) {
	int s_x = start[0];
	int s_y = start[1];
	int e_x = end[0];
	int e_y = end[1];
	int m = e_y - s_y;
	int n = e_x - s_x;
	float rmsd = 0;
	for(int i = s_y; i < e_y; i++){
		for(int j = s_x; j < e_x; j++){
			float interpolated[3];
			int point[2] = {j, i};
			interpolate(X, point, start, angle, interpolated, x_m, x_n);
			if(interpolated[0] == -1){
				return FLT_MAX;
			}
			for(int k = 0; k < 3; k++){
				rmsd += (interpolated[k] - Y[(i - s_y)*q_n*3 + (j - s_x)*3 + k])*(interpolated[k] - Y[(i - s_y)*q_n*3 + (j - s_x)*3 + k]);
			}	
		}
	}
	rmsd /= m*n*3;
	return sqrt(rmsd);
}

void get_m_n(string file, int &m, int &n){
	ifstream in;
	in.open(file);
	in >> m;
	in >> n;
	in.close();
}

void readImage(string file, int *A) {
	int m, n, pixel;

	ifstream in;
	in.open(file);
	in >> m;
	in >> n;

	for(int i = 0; i < m; i++){
		for(int j = 0; j < n; j++){
			for(int k = 0; k < 3; k++){
				in >> pixel;
				A[(m - 1 - i)*n*3 + j*3 + k] = pixel;
			}
		}
	}
	in.close();
}


void filterBox(int *res, int angle, int q_m, int q_n){
	int *pivot, *br, *tl, *tr;
	pivot = new int[2];
	br = new int[2];
	tl = new int[2];
	tr = new int[2];

	pivot[0] = 0; pivot[1] = 0;
	br[0] = q_n - 1; br[1] = 0;
	tl[0] = 0; tl[1] = q_m - 1;
	tr[0] = q_n - 1; tr[1] = q_m - 1;

	float *tl_n, *br_n, *tr_n;
	tl_n = new float[2];
	br_n = new float[2];
	tr_n = new float[2];
	getRotatedCoordinateHost(tl_n, tl, pivot, angle);
	getRotatedCoordinateHost(br_n, br, pivot, angle);
	getRotatedCoordinateHost(tr_n, tr, pivot, angle);

	int x1, x2, y1, y2;

	if(angle == 45){
		x1 = ceil(tl_n[0]);
		x2 = floor(br_n[0]);
		y1 = pivot[1];
		y2 = floor(tr_n[1]);
	}
	else if(angle == -45){
		x1 = pivot[0];
		x2 = floor(tr_n[0]);
		y1 = ceil(br_n[1]);
		y2 = floor(tl_n[1]);
	}
	else{
		x1 = pivot[0];
		x2 = br[0];
		y1 = pivot[1];
		y2 = tl[1];
	}
	res[0] = x1;
	res[1] = x2;
	res[2] = y1;
	res[3] = y2;
	delete pivot;
	delete br;
	delete tl;
	delete tr;
	delete tl_n;
	delete br_n;
	delete tr_n;
}

// box[0] --> left-x, box[1] --> right-x, box[2] --> bottom-y, box[3] --> top-y
// pivot[0] --> x-coordinate, pivot[1] --> y-coordinate
__device__
float filter(float *grey_data_img, float query_grey, int box[], int pivot[], float t2, int x_m, int x_n){
	// return true;

	int x1 = pivot[0] + box[0];
	int x2 = pivot[0] + box[1];
	int y1 = pivot[1] + box[2];
	int y2 = pivot[1] + box[3];
	if(!inLimits(x1, y1, x_n, x_m) || !inLimits(x2, y2, x_n, x_m)){
		return FLT_MAX; // Actually false
	}

	float data_grey = 0;
	int count = 0;
	for(int i = y1; i <= y2; i++){
		for(int j = x1; j <= x2; j++){
			if(inLimits(j, i, x_n, x_m)){
				data_grey += (float)grey_data_img[i*x_n + j]; // Check this?
				count++;
			}
		}
	}
	data_grey /= (x2 - x1 + 1)*(y2 - y1 + 1);
	return abs(data_grey - query_grey);
}

void getGrey(int *X, float *grey_data_img, int x_m, int x_n){
	for(int i = 0; i < x_m; i++){
		for(int j = 0; j < x_n; j++){
			float x = 0;
			for(int k = 0; k < 3; k++){
				x += X[i*x_n*3 + j*3 + k];
			}
			grey_data_img[i*x_n + j] = (float)x/3;
		}
	}
}

__global__
void kernel(float *res, int *X, int *Q, float *grey_data_img, int *box0, int *box1, int *box2, int x_m, int x_n, int q_m, int q_n, float t1, float t2, float query_grey) {
	int angles[3] = {0, 45, -45};
	int offset = blockDim.x*blockIdx.x + threadIdx.x;
	if(offset >= x_m*x_n){
		return;
	}
	int i = offset / x_n; // y-coordinate
	int j = offset % x_n; // x-coordinate
	int pivot[2] = {j, i};
	for(int iter = 0; iter < 3; ++iter){
		int angle = angles[iter];
		res[i*x_n*3 + j*3 + iter] = -1;
		int box[4];
		if(iter == 0){
			box[0] = box0[0]; box[1] = box0[1]; box[2] = box0[2]; box[3] = box0[3];
		}
		else if(iter == 1){
			box[0] = box1[0]; box[1] = box1[1]; box[2] = box1[2]; box[3] = box1[3];
		}
		else{
			box[0] = box2[0]; box[1] = box2[1]; box[2] = box2[2]; box[3] = box2[3];
		}
		float diff = filter(grey_data_img, query_grey, box, pivot, t2, x_m, x_n);
		if(diff <= t2){
			int start[2] = {j, i};
			int end[2] = {j + q_n, i + q_m};
			float f = rmsd(X, Q, start, end, angle, q_m, q_n, x_m, x_n);
			if(f <= t1){
				res[i*x_n*3 + j*3 + iter] = (float)f;
			}
		}
	}
}

int main(int argc, char *argv[]) {
	// auto begin = std::chrono::high_resolution_clock::now();
	string data_image_path = argv[1];
	string query_image_path = argv[2];
	float t1 = stof(argv[3]);
	float t2 = stof(argv[4]);
	int topn = stoi(argv[5]);

	int x_m, x_n, q_m, q_n; // q_m is y-max, q_n is x-max, x_m is y-max, x_n is x-max

	int *X, *Q, *dX, *dQ;
	float *grey_data_img, *res, *dgrey_data_img, *dres;

	get_m_n(data_image_path, x_m, x_n);
	X = new int[x_m*x_n*3];
	res = new float[x_m*x_n*3];
	readImage(data_image_path, X);

	get_m_n(query_image_path, q_m, q_n);
	Q = new int[q_m*q_n*3];
	readImage(query_image_path, Q);

	grey_data_img = new float[x_m*x_n];
	getGrey(X, grey_data_img, x_m, x_n);

	// Finding Average Greyscale value of Query Image
	float query_grey = 0;
	for(int i = 0; i < q_m; i++){
		for(int j = 0; j < q_n; j++){
			float x = 0;
			for(int k = 0; k < 3; k++){
				x += Q[i*q_n*3 + j*3 + k];
			}
			query_grey += (float)x/3;
		}
	}
	query_grey /= q_m*q_n;

	// Finding axis-aligned bounding box
	int *boxes[3];
	int *dbox0, *dbox1, *dbox2;
	for(int i = 0; i < 3; i++){
		boxes[i] = new int[4];
	}
	filterBox(boxes[0], 0, q_m, q_n);
	filterBox(boxes[1], 45, q_m, q_n);
	filterBox(boxes[2], -45, q_m, q_n);

	// Allocate on GPU
	cudaMalloc(&dX, x_m*x_n*3*sizeof(int));
	cudaMalloc(&dQ, q_m*q_n*3*sizeof(int));
	cudaMalloc(&dgrey_data_img, x_m*x_n*sizeof(float));
	cudaMalloc(&dres, x_m*x_n*3*sizeof(float));
	cudaMalloc(&dbox0, 4*sizeof(int));
	cudaMalloc(&dbox1, 4*sizeof(int));
	cudaMalloc(&dbox2, 4*sizeof(int));

	// Copy on GPU
	cudaMemcpy(dX, X, x_m*x_n*3*sizeof(int), cudaMemcpyDefault);
	cudaMemcpy(dQ, Q, q_m*q_n*3*sizeof(int), cudaMemcpyDefault);
	cudaMemcpy(dgrey_data_img, grey_data_img, x_m*x_n*sizeof(float), cudaMemcpyDefault);
	cudaMemcpy(dbox0, boxes[0], 4*sizeof(int), cudaMemcpyDefault);
	cudaMemcpy(dbox1, boxes[1], 4*sizeof(int), cudaMemcpyDefault);
	cudaMemcpy(dbox2, boxes[2], 4*sizeof(int), cudaMemcpyDefault);

	int num_blocks = ceil((float)(x_n*x_m)/1024);
	kernel<<<num_blocks, 1024>>>(dres, dX, dQ, dgrey_data_img, dbox0, dbox1, dbox2, x_m, x_n, q_m, q_n, t1, t2, query_grey);

	cudaMemcpy(res, dres, x_m*x_n*3*sizeof(float), cudaMemcpyDeviceToHost);
	int angles[3] = {0, 45, -45};

	priority_queue<pair<float, tuple<int, int, int>>> q; 

	for(int i = 0; i < x_m; i++){ // i correcposds to y-coordinate
		for(int j = 0; j < x_n; j++){ // j corresponds to x-coordinate
			for(int k = 0; k < 3; k++){
				float f = res[i*x_n*3 + j*3 + k];
				if(f != -1){
					if(q.size() >= topn && f > q.top().first){
						continue;
					}
					q.push({f, {i, j, angles[k]}});
					if(q.size() > topn){
						q.pop();
					}
				}
			}
		}
	}
	
	ofstream outdata;
	outdata.open("output.txt");
	priority_queue<pair<float, tuple<int, int, int>>> output; 
	while(!q.empty()){
		pair<float, tuple<int, int, int>> p = q.top();
		tuple<int, int, int> tup = p.second;
		q.pop();
		output.push({-p.first, tup});
	}
	while(!output.empty()){
		pair<float, tuple<int, int, int>> p = output.top();
		output.pop();
		tuple<int, int, int> tup = p.second;
		//cout << get<0>(tup) << " " << get<1>(tup) << " " << get<2>(tup) << " " << -1*p.first << endl; // Remove this
		outdata << get<0>(tup) << " " << get<1>(tup) << " " << get<2>(tup) << endl;
	}

	delete X;
	delete Q;
	delete grey_data_img;

	cudaFree(dX);
	cudaFree(dQ);
	cudaFree(dgrey_data_img);
	cudaFree(dres);
	cudaFree(dbox0);
	cudaFree(dbox1);
	cudaFree(dbox2);

	outdata.close();
	//auto end1 = std::chrono::high_resolution_clock::now();
	//cout << "Time taken for completion: " << (1e-6 * (std::chrono::duration_cast<std::chrono::nanoseconds>(end1 - begin)).count()) << "ms" << endl;
	// Remove this
}