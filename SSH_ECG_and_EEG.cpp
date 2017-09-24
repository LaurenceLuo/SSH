// DTW.cpp : Defines the entry point for the console application.
//

//#include "stdafx.h"
// **************************************************************
// *
// *  Preproccessing.cpp
// *
// *  Based on LB_Keogh and LSH algorithm, perform fast Dynamic Time Wrapping
// *
// **************************************************************
#include <stdio.h>
#include <cstdlib>
#include <math.h>
#include <time.h>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <list>
#include <algorithm>
#include <cstdlib>
#include <utility>
#include <random>
//#include "RTree_HD.h"
#include <queue>
#include <ctime>
#include <iterator>


using namespace std;


#define M 5000//The number of time series
#define T 256 //The length of a time serie 256
#define D 64 //The dimension of a time point
#define bandwidth  0.12*T//Used for Sakoe-Chiba Band restriction, 0<=bandwidth<=T
//#define slope_variance 1 //Used for Itakura parallelogram restriction
#define constraint 4 //Itakura parallelogram window, must be smaller than T
#define PAAReductionFactor 1//the equal amount for each cell, must be a factor of T
#define L 3//The number of LSH function groups
#define K 6//The number of LSH functions in each group
#define W 1//Numerical constant
#define threshold 120//410
#define BlockNum T/PAAReductionFactor//The number of blocks left after using PAA
#define Epsilon 0.03*sqrt(D)*sqrt(T)//Threshold used in RTree RangeSearch
#define R 0.5*sqrt(D)
#define KNN 50
#define FLAG 50
#define LSHfunc 256
#define NUM 0.1*M
#define ratio 0.99
#define SD 0.01
#define ClusterNum 2

//Parameters for SSH contrast experiment
#define Delta 3 //Step Size
#define Size_W 80 //Vector Size
#define NB (T-Size_W)/Delta //Stream Size
#define N_grams 13 //The number of grams
#define Length_N pow(2,N_grams) // #type of shingles
#define HashtableNum 20 //The number of hash tables





//compute Euclidean distance between two datasets with length N
float distance(float p[], float q[], int N)
{
	float dist = 0;
	for (int i = 0; i<N; i++)
	{
		dist += (p[i] - q[i])*(p[i] - q[i]);
	}
	return sqrt(dist);
}

//compute Euclidean distance between two points with dimension d
float distance_HD(float p[], float q[], int d)
{
	float dist = 0;
	for (int i = 0; i<d; i++) {
		dist += (p[i] - q[i])*(p[i] - q[i]);
	}
	return dist;
}


//compute Euclidean distance between two points with dimension d
float Euclidean_disatance(float p[], float q[], int d)
{
	float dist = 0;
	for (int i = 0; i<d; i++) {
		dist += (p[i] - q[i])*(p[i] - q[i]);
	}
	return sqrt(dist);
}


//compute Euclidean distance between two series with dimension d
float distance_HD(float** p, float** q)
{
	float dist = 0;
	for (int i = 0; i<T; i++) {
		for (int j = 0; j<D; j++) {
			dist += (p[i][j] - q[i][j])*(p[i][j] - q[i][j]);
		}
	}
	return sqrt(dist);
}


float*** z_normalization(float***data, int row, int col)
{
    float sum=0;
    float square_root_sum=0;
    for (int i = 0; i < M; i++){
        for (int j = 0; j<row; j++) {
            for (int k = 0; k<col; k++) {
                sum+=data[i][j][k];
                square_root_sum+=data[i][j][k]*data[i][j][k];
            }
        }
    }
    float mean=sum/(M*row*col);
    float square_root_mean=square_root_sum/(M*row*col);
    square_root_mean=sqrt(square_root_mean-mean*mean);
    for (int i = 0; i < M; i++){
        for (int j = 0; j<row; j++) {
            for (int k = 0; k<col; k++) {
                data[i][j][k]=(data[i][j][k]-mean)/square_root_mean;
            }
        }
    }
    return data;
}

//load dataset from files
float*** load_and_z_normalization(const char* filename, int row, int col)
{
    ifstream file(filename); // declare file stream:
    string value;
    string num;
    int i, j;
    float*** data = new float**[M];
    for (int i = 0; i < M; i++){
        data[i] = new float*[row];
        for (int j = 0; j < row; j++){
            data[i][j] = new float[col];
        }
    }
    float sum=0;
    float square_root_sum=0;
    for (int i = 0; i < M; i++){
        for (int j = 0; j<row; j++) {
            for (int k = 0; k<col; k++) {
                getline(file, value);
                num = string(value);
                //cout<<num<<endl;
                data[i][j][k] = ::atof(num.c_str());
                sum+=data[i][j][k];
                square_root_sum+=data[i][j][k]*data[i][j][k];
            }
        }
    }
    float mean=sum/(M*row*col);
    float square_root_mean=square_root_sum/(M*row*col);
    square_root_mean=sqrt(square_root_mean-mean*mean);
    for (int i = 0; i < M; i++){
        for (int j = 0; j<row; j++) {
            for (int k = 0; k<col; k++) {
                data[i][j][k]=(data[i][j][k]-mean)/square_root_mean;
            }
        }
    }
    file.close();
    return data;
}

/*float** load_query(const char* filename, int row, int col)
{
    ifstream file(filename); // declare file stream:
    string value;
    string num;
    int i, j;
    float** data = new float*[row];
    for (int i = 0; i < row; i++){
        data[i] = new float[col];
    }
    
    for (int j = 0; j<row; j++) {
        for (int k = 0; k<col; k++) {
            getline(file, value);
            num = string(value);
            //cout<<num<<endl;
            data[j][k] = ::atof(num.c_str());
        }
    }
    file.close();
    return data;
}



float** load_query_and_z_normalization(const char* filename, int row, int col)
{
    ifstream file(filename); // declare file stream:
    string value;
    string num;
    int i, j;
    float** data = new float*[row];
    for (int i = 0; i < row; i++){
        data[i] = new float[col];
    }
    
    float sum=0;
    float square_root_sum=0;
    for (int j = 0; j<row; j++) {
        for (int k = 0; k<col; k++) {
            getline(file, value);
            num = string(value);
            //cout<<num<<endl;
            data[j][k] = ::atof(num.c_str());
            sum+=data[j][k];
            square_root_sum+=data[j][k]*data[j][k];
        }
    }
    float mean=sum/(row*col);
    float square_root_mean=square_root_sum/(row*col);
    square_root_mean=sqrt(square_root_mean-mean*mean);
    for (int j = 0; j<row; j++) {
        for (int k = 0; k<col; k++) {
            data[j][k]=(data[j][k]-mean)/square_root_mean;
        }
    }
    
    file.close();
    return data;
}*/

//load dataset from files (multivariates)
float**  load_data(const char* filename, int row, int col)
{
	ifstream file(filename); // declare file stream:
	string value;
	string num;
	int i, j;
	int count = -1;
	float **data = new float*[row];
	for (int i = 0; i < row; i++)
		data[i] = new float[col];

	for (int i = 0; i<row; i++) {
		for (int j = 0; j<col; j++) {
			getline(file, value, ' ');
			num = string(value, 0, value.length());
			data[i][j] = ::atof(num.c_str());
		}
		getline(file, value, '\n');
	}
	file.close();
	return data;
}

//normalize input datasets to the range of [0,1]
void normalization_HD(float***&p) {
	float max[D] = { -INFINITY };
	float min[D] = { INFINITY };

	for (int d = 0; d<D; d++) {
		for (int i = 0; i<M; i++) {
			for (int j = 0; j<T; j++) {
				if (p[i][j][d] >= max[d])
					max[d] = p[i][j][d];
				if (p[i][j][d]<min[d])
					min[d] = p[i][j][d];
			}
		}
	}
	for (int i = 0; i<M; i++) {
		for (int j = 0; j<T; j++) {
			for (int d = 0; d<D; d++) {
				p[i][j][d] = (p[i][j][d] - min[d]) / (max[d] - min[d]);
			}
		}
	}

}


//Basic one dimensional DTW
float DTW_Basic(float* p, float* q)
{
	float gamma[T][T];
	float dist[T][T];
	for (int i = 0; i<T; i++) {
		for (int j = 0; j<T; j++) {
			dist[i][j] = (p[i] - q[j])*(p[i] - q[j]);//distance(p[i],q[j]);
		}
	}
	gamma[0][0] = dist[0][0];
	for (int i = 1; i<T; i++) {
		gamma[0][i] = dist[0][i] + gamma[0][i - 1];
		//gamma[0][i]=INFINITY;
	}
	for (int i = 1; i<T; i++) {
		gamma[i][0] = dist[i][0] + gamma[i - 1][0];
		//gamma[i][0]=INFINITY;
	}
	for (int i = 1; i<T; i++) {
		for (int j = 1; j<T; j++) {
			if ((j - i < bandwidth) && (j - i > -bandwidth))//Rectangle restriction
				gamma[i][j] = min(gamma[i - 1][j - 1], min(gamma[i - 1][j], gamma[i][j - 1])) + dist[i][j];
			else gamma[i][j] = dist[i][j];
		}
	}

	//cout<<gamma[95][95]<<endl;
	vector<pair<int, int> > pair_vector;
	int i = 0;
	int j = 0;

	while (i<T - 1 && j<T - 1) {
		if (i == T - 2 && j != T - 2)
			j += 1;
		else if (j == T - 2 && i != T - 2)
			i += 1;
		else if (i == T - 2 && i == T - 2) {
			i += 1;
			j += 1;
		}
		else {
			if (gamma[i + 1][j + 1] - dist[i + 1][j + 1] == gamma[i + 1][j])
				i += 1;
			else if (gamma[i + 1][j + 1] - dist[i + 1][j + 1] == gamma[i][j + 1])
				j += 1;
			else {
				i += 1;
				j += 1;
			}
		}
		pair_vector.push_back(make_pair(i, j));
	}
	//cout<<"(p, q)"<<endl;
	float cost = 0;
	for (int i = 0; i<pair_vector.size(); i++) {
		//cout << "Pair: "<<pair_vector[i].first << ", " << pair_vector[i].second << endl;
		cost = cost + (p[pair_vector[i].first] - q[pair_vector[i].second])*(p[pair_vector[i].first] - q[pair_vector[i].second]);
		//cout<<cost<<endl;
	}
	return sqrt(cost);
}




//Multi-dimension DTW by calculating DTW in every dimension and sum them up, using DTW_Basic function
float DTW_1D(float** p, float** q) {
	float gamma[D][D];
	float dist[D][D];
	float** p_new = new float*[D];
	float** q_new = new float*[D];
	for (int i = 0; i < D; i++) {
		p_new[i] = new float[T];
		q_new[i] = new float[T];
	}
	for (int i = 0; i<D; i++) {
		for (int j = 0; j<T; j++) {
			p_new[i][j] = p[j][i];
			q_new[i][j] = q[j][i];
		}
	}
	float cost = 0;
	for (int i = 0; i<D; i++) {
		cost += DTW_Basic(p_new[i], q_new[i]);
	}


	for (int i = 0; i < D; i++) {
		delete[] p_new[i];
		delete[] q_new[i];
	}
	delete[] p_new;
	delete[] q_new;

	return cost;//square root? Already did in DTW_Basic
}


//Basic multiple dimensional DTW
float DTW_HD(float** p, float** q)
{
	float gamma[T][T];
	float dist[T][T];
	for (int i = 0; i<T; i++) {
		for (int j = 0; j<T; j++) {
			dist[i][j] = distance_HD(p[i], q[j], D); //no sqrt
		}
	}
	gamma[0][0] = dist[0][0];
	for (int i = 1; i<T; i++) {
		gamma[0][i] = dist[0][i] + gamma[0][i - 1];
		gamma[i][0] = dist[i][0] + gamma[i - 1][0];
		//gamma[i][0]=INFINITY;
	}
	float temp = 0;
	for (int i = 1; i<T; i++) {
		for (int j = 1; j<T; j++) {
			if ((j - i < bandwidth) && (j - i > -bandwidth))//Rectangle restriction
				gamma[i][j] = min(gamma[i - 1][j - 1], min(gamma[i - 1][j], gamma[i][j - 1])) + dist[i][j];
			else gamma[i][j] = dist[i][j];
			if (gamma[i][j] >= temp) {
				temp = gamma[i][j];
			}
		}
	}
	//cout<<gamma[95][95]<<endl;
	vector<pair<int, int> > pair_vector;
	int i = 0;
	int j = 0;
	while (i<T - 1 && j<T - 1) {
		if (i == T - 2 && j != T - 2)
			j += 1;
		else if (j == T - 2 && i != T - 2)
			i += 1;
		else if (i == T - 2 && i == T - 2) {
			i += 1;
			j += 1;
		}
		else {
			if (gamma[i + 1][j + 1] - dist[i + 1][j + 1] == gamma[i + 1][j])
				i += 1;
			else if (gamma[i + 1][j + 1] - dist[i + 1][j + 1] == gamma[i][j + 1])
				j += 1;
			else {
				i += 1;
				j += 1;
			}
		}
		//pair_vector.push_back(make_pair(i,j));
		pair_vector.push_back(make_pair(i, j));
	}
	float cost = 0;
	for (int i = 0; i<pair_vector.size(); i++) {
		//cout << "Pair: "<<pair_vector[i].first << ", " << pair_vector[i].second << endl;
		cost = cost + distance_HD(p[pair_vector[i].first], q[pair_vector[i].second], D);
	}
	return sqrt(cost);
}

vector<int> DTW_GroundTruth_Range(float**query, float*** datasets) {
// find series which DTW<value and fall in the r-envelop of query series
	vector<int >candidate;
	float dtw_dist = 0;
	for (int i = 0; i<M; i++) {
		dtw_dist = DTW_HD(query, datasets[i]);
		if (dtw_dist <= Epsilon) 
		{
			bool isTrue = true;
			for (int j = 0; j<T; j++)
			{
				if (Euclidean_disatance(query[j], datasets[i][j], D)> R)
				{
					isTrue = false;
					break;
				}			
			}
			if (isTrue)			
			candidate.push_back(i);
		}
	}
	return candidate;
}

vector<int> DTW_GroundTruth_KNN(float** query, float*** datasets){
    struct sort_int_StoL {
        bool operator()(int left, int right) {
            return left < right;
        }
    };//Int from small to large
    vector<pair<int,float> > candidate_KNN;
    int count=0;
    struct sort_pred {
        bool operator()(const std::pair<int,float> &left, const std::pair<int,float> &right) {
            return left.second < right.second;
        }
    };
    for(int m=0;m<M;m++){
        if(count<KNN){
            candidate_KNN.push_back(make_pair(m,DTW_HD(query,datasets[m])));
        }
        else{
            sort(candidate_KNN.begin(),candidate_KNN.end(),sort_pred());
            float temp=DTW_HD(query,datasets[m]);
            if(temp<candidate_KNN.back().second){
                candidate_KNN.pop_back();
                candidate_KNN.push_back(make_pair(m,temp));
            }
        }
        count++;
    }
    vector<int> KNN_output;
    for(vector<pair<int,float> >::iterator it=candidate_KNN.begin();it!=candidate_KNN.end();++it){
        KNN_output.push_back((*it).first);
    }
    //sort(KNN_output.begin(), KNN_output.end(), sort_int_StoL());
    return KNN_output;
}



vector<int> REnvelope_GroundTruth(float** query, float*** datasets) {
	vector<int> candidate;
	for (int m = 0; m<M; m++)
	{
		bool isTrue = true;
		for (int i = 0; i<T; i++)
		{
			if (Euclidean_disatance(query[i], datasets[m][i], D)> R)
			{
				isTrue = false;
				break;
			}			
		}
		if (isTrue)
			candidate.push_back(m);
	}
	return candidate;
}



float accuracy(vector<int> truth, vector<int> results)
{
	float accuracy;
	float count,acc_count;
	acc_count = 0;
	vector<int>::iterator it;
	vector<int>::iterator it2;

	for (it = truth.begin(); it != truth.end(); ++it)
	{   
		count = 0;
		for (it2 = results.begin(); it2 != results.end(); ++it2)
		{
			if ((*it) == (*it2))
			{
				count += 1;
			}
				
		}
		if (count != 0)
			acc_count += 1;
	}
	//cout << acc_count << endl;
	int size=truth.size();
	//cout << truth.size() << endl;
	accuracy = acc_count / size;
	return accuracy;
}

float accuracy_KNN(vector<int> truth, vector<int> results, int query_id, float*** datasets)
{
    float sum_truth,sum_results=0;
    
    for(int i=0;i<KNN;i++){
        sum_truth+=DTW_HD(datasets[truth[i]],datasets[query_id]);
        sum_results+=DTW_HD(datasets[results[i]],datasets[query_id]);
    }
    float accuracy=sum_results/sum_truth;
    return accuracy;
}

float accuracy_KNN_classification(vector<int> truth, vector<int> results, int query_id, float*** datasets)
{
    int count_truth=0;
    int count_results=0;
    cout<<"Query class type number: "<<datasets[query_id][0][D]<<endl;
    for(int i=0;i<KNN;i++){
        if(datasets[query_id][0][D]==datasets[truth[i]][0][D]){
            count_truth++;
        }
        else {
            cout<<"Groud Truth Wrong Index: "<<results[i]<<", Label: "<<datasets[results[i]][0][D]<<endl;
        }
        if(datasets[query_id][0][D]==datasets[results[i]][0][D]){
            count_results++;
        }
        else {
            cout<<"Result Wrong index: "<<results[i]<<", Label: "<<datasets[results[i]][0][D]<<endl;
        }
    }
    cout<<"Ground Truth KNN classification number: "<<count_truth<<endl;
    cout<<"LSH KNN classification number: "<<count_results<<endl;
    
    float accuracy=float(count_results)/float(count_truth);
    return accuracy;
}


// all sliding window use a same r
float* generateRandom_vector_r(){
    //Generate r vector
    default_random_engine generator(time(NULL));
    normal_distribution<float> distribution(0.0, 1.0);
    float* random=new float[Size_W];
    for(int i=0;i<Size_W;i++){
        float temp = distribution(generator);
        //while (temp<0 || temp>1)
            //temp = distribution(generator);
        random[i]=temp;
    }
    return random;
}


// each window uses a different r vector
float** generateRandom_vector_r_w(){
    //Generate r vector
    default_random_engine generator(time(NULL));
    normal_distribution<float> distribution(0.0, 1.0);
    float** random_r=new float*[NB];
    for(int m=0;m<NB;m++){
        random_r[m]=new float[Size_W];
        for(int i=0;i<Size_W;i++){
            float temp = distribution(generator);
            //while (temp<0 || temp>1)
                //temp = distribution(generator);
            random_r[m][i]=temp;
        }
    }
    //cout<<"RANDOM: "<<endl;
    for(int m=0;m<NB;m++){
        for(int i=0;i<Size_W;i++){
            cout<<random_r[m][i]<<" ";
        }
        //cout<<endl;
    }
    return random_r;
}

/*
vector<int> bit_profile_extraction_onevalue(float* vector_r, float** q){
    vector<int> bit_profile(NB);
    for(int b=0;b<NB;b++){
        float sum=0;
        for(int d=0;d<D;d++){
            for(int w=0;w<Size_W;w++){
                sum+=q[b*Delta+w][d]*vector_r[w];
            }
        }
        if(sum>=0){
            bit_profile[b]=1;
        }
        else{
            bit_profile[b]=-1;
        }
    }
    cout<<"Bit Profile 1D: "<<endl;
    for(vector<int>::iterator it = bit_profile.begin(); it != bit_profile.end(); ++it){
        cout<<(*it)<<" ";
    }
    cout<<endl;
    return bit_profile;
}*/

vector<int> bit_profile_extraction_onevalue_hd(float** vector_r, float** q){
    vector<int> bit_profile(NB);
    for(int b=0;b<NB;b++){
        float sum=0;
        for(int d=0;d<D;d++){
            for(int w=0;w<Size_W;w++){
                sum+=q[b*Delta+w][d]*vector_r[b][w];
            }
        }
        if(sum>=0){
            bit_profile[b]=1;
        }
        else{
            bit_profile[b]=-1;
        }
    }
    //cout<<"Bit Profile 1D: "<<endl;
    for(vector<int>::iterator it = bit_profile.begin(); it != bit_profile.end(); ++it){
        //cout<<(*it)<<" ";
    }
    //cout<<endl;
    return bit_profile;
}

vector< vector<int> > bit_profile_extraction_dvalue(float** vector_r, float** q){
    vector<vector<int > > bit_profile(NB,vector<int >(D));
    for(int b=0;b<NB;b++){
        for(int d=0;d<D;d++){
            float sum=0;
            for(int w=0;w<Size_W;w++){
                sum+=q[b*Delta+w][d]*vector_r[b][w];
            }
            if(sum>=0){
                bit_profile[b][d]=1;
            }
            else{
                bit_profile[b][d]=-1;
            }
        }
    }
    /*cout<<"Bit Profile HD: "<<endl;
    for(vector< vector<int> >::iterator it = bit_profile.begin(); it != bit_profile.end(); ++it){
        for(vector<int>::iterator it2 = (*it).begin(); it2 != (*it).end(); ++it2){
            cout<<(*it2)<<" ";
        }
        cout<<endl;
    }*/
    return bit_profile;
}


vector<int> generate_n_grams_onevalue(vector<int> bit_vector){
    vector<int> n_grams(Length_N,0);
    for(int n=0;n<NB-N_grams;n++){
        int index=0;
        int count=N_grams-1;
        for(int i=0;i<N_grams;i++){
            if(bit_vector[n+i]==1){
                index+=pow(2,count);
            }
            count--;
        }
        //cout<<"index: "<<index<<endl;
        n_grams[index]++;
    }
    //cout<<"N Grams 1D: "<<endl;
    for(vector<int>::iterator it = n_grams.begin(); it != n_grams.end(); ++it){
        //cout<<(*it)<<" ";
    }
    //cout<<endl;
    return n_grams;
}

vector<int> generate_n_grams_dvalue(vector< vector<int> > bit_vector){
    vector<int> n_grams(Length_N,0);
    for(int n=0;n<NB-N_grams;n++){
        for(int d=0;d<D;d++){
            int index=0;
            int count=N_grams-1;
            for(int i=0;i<N_grams;i++){
                if(bit_vector[n+i][d]==1){
                    index+=pow(2,count);
                }
                count--;
            }
            //cout<<"index: "<<index<<endl;
            n_grams[index]++;
        }
    }
    //cout<<"N Grams 1D: "<<endl;
    //for(vector<int>::iterator it = n_grams.begin(); it != n_grams.end(); ++it){
        //cout<<(*it)<<" ";
    //}
    //cout<<endl;
    return n_grams;
}



int* generate_seed()
{
    int* seed= new int[HashtableNum];
    for(int i=0; i<HashtableNum;i++){
	seed[i]=rand() % 100000;
	//cout<<seed[i]<<endl;
	}
    return seed;
}


vector<int> gen_minhash_value(vector<int> n_grams, int * seed){
    vector<int> hash_value(HashtableNum);
    //initialise hashes to 0
    for(int i=0; i<HashtableNum; i++)
	hash_value[i]=0;

    for(int h=0;h<HashtableNum;h++){
        int randomseed=seed[h];
        //cout<<"????"<<endl;
        int count=0;
        while(1){
            //cout<<"seed: "<<randomseed<<endl;
            if(randomseed<0){
                cout<<"seed: "<<randomseed<<endl;
                continue;
            }
            default_random_engine generator(randomseed);
            uniform_real_distribution<float> distribution(0.0, 1.0);
            
            float r = (Length_N)*distribution(generator);
            //check whether r falls in the green zone
            float integer;
            float decimal=modf(r,&integer);
            //cout<<r<<": "<<decimal<<endl;
            //cout<<n_grams[integer]<<"/"<<NB-N_grams<<"="<<float(n_grams[integer])/float(NB-N_grams)<<endl;
            if(decimal<float(n_grams[integer])/float(NB-N_grams)){
                //cout<<r<<": "<<decimal<<endl;
                //cout<<n_grams[integer]<<"/"<<NB-N_grams<<"="<<float(n_grams[integer])/float(NB-N_grams)<<endl;
                break;
            }
            if(count>100000){
                //cout<<"random"<<": "<<1000000<<endl;
                //cout<<n_grams[integer]<<"/"<<NB-N_grams<<"="<<float(n_grams[integer])/float(NB-N_grams)<<endl;
                //default_random_engine generator(time(NULL));
                //uniform_real_distribution<float> distribution(0.0, 1.0);
                hash_value[h]=0;
                break;
            }
            count++;
            hash_value[h]++;
            randomseed=r*1000;
        }
    }
    cout<<"Hash value: "<<endl;
    for(vector<int>::iterator it = hash_value.begin(); it != hash_value.end(); ++it){
        cout<<(*it)<<" ";
    }
    cout<<endl;
    return hash_value;
}





//one window-> one value
vector<vector<int> > full_minhash_values_type1(float***datasets, float** vector_r, int* seed)
{
    vector<vector<int > > full_hash_values;
    vector<int> bit_profile;
    vector<int> grams;
    for(int m=0;m<M; m++)
    {
        cout<<"Series "<<m<<": "<<endl;
        bit_profile=bit_profile_extraction_onevalue_hd(vector_r,datasets[m]);
        grams=generate_n_grams_onevalue(bit_profile);
        full_hash_values.push_back(gen_minhash_value(grams, seed));
    }
    return full_hash_values;
}



// one window--> d values
vector<vector<int> > full_minhash_values_type2(float***datasets, float** vector_r, int* seed)
{
    vector<vector<int > > full_hash_values;
    vector<vector<int> > bit_profile;
    vector<int> grams;
    for(int m=0;m<M; m++)
    {
        bit_profile=bit_profile_extraction_dvalue(vector_r,datasets[m]);
        grams=generate_n_grams_dvalue(bit_profile);
        full_hash_values.push_back(gen_minhash_value(grams, seed));
    } 
    return full_hash_values;
}



vector<int> SSH_generate_candidate(float** query, float*** datasets,float** vector_r, int* seed, int type)
{   
    // type1: one window->one value;     type2: one window-> d values
    vector<vector<int> >full_minhash_values;
    vector<int> query_hashvalue;
    if(type==1)
    {  //calculate minhash value for query time series and datasets
        vector<int> query_bit_profile=bit_profile_extraction_onevalue_hd(vector_r, query);
        vector<int> grams=generate_n_grams_onevalue(query_bit_profile);
        query_hashvalue=gen_minhash_value(grams, seed);
        full_minhash_values=full_minhash_values_type1(datasets, vector_r, seed);
    }
    else
    {
        vector<vector<int> > bit_profile=bit_profile_extraction_dvalue(vector_r,query);
        vector<int> grams2=generate_n_grams_dvalue(bit_profile);
        query_hashvalue=gen_minhash_value(grams2, seed);
        full_minhash_values=full_minhash_values_type2(datasets, vector_r, seed);
    }
    //generate candidates
    vector<int> candidate;
    for(int m=0; m<M; m++)
    {
        for(int k=0; k<HashtableNum; k++){
            if((query_hashvalue[k]==0)||(full_minhash_values[m][k]==0)){
                continue;
            }
            if(query_hashvalue[k]==full_minhash_values[m][k]){
                candidate.push_back(m);
                break;
            }
        }
    }
    cout<<"the size of candiate is:"<<candidate.size()<<endl;
    return candidate;
}


vector<int>  SSH_KNN(vector<int> candidate, float** query, float*** datasets) {
    int count = 0;
    struct sort_pred {
        bool operator()(const std::pair<int, float> &left, const std::pair<int, float> &right) {
            return left.second < right.second;
        }
    };
    vector<pair<int, float> > candidate_KNN;
    for (vector<int>::iterator it = candidate.begin(); it != candidate.end(); ++it) {
        if (count<KNN) {
            candidate_KNN.push_back(make_pair((*it), DTW_HD(query, datasets[*it])));
            sort(candidate_KNN.begin(), candidate_KNN.end(), sort_pred());
        }
        else {
            sort(candidate_KNN.begin(), candidate_KNN.end(), sort_pred());
            float temp = DTW_HD(query, datasets[*it]);
            if (temp<candidate_KNN.back().second) {
                candidate_KNN.pop_back();
                candidate_KNN.push_back(make_pair((*it), temp));
            }
        }
        count++;
    }
    vector<int> KNN_output;
    for (vector<pair<int, float> >::iterator it = candidate_KNN.begin(); it != candidate_KNN.end(); ++it) {
        KNN_output.push_back((*it).first);
        //cout<<(*it).second<<endl;
    }
    
    return  KNN_output;
    
}

int main()
{
    //load data
    int num = 1;
    int i, j;
    clock_t beginRQ, endRQ;
    double time;
    int query_id=10;

    beginRQ = clock();

    // load ECG dataset and Z-normalization
    //string input="./ECG.txt";
    //float*** datasets = load_and_z_normalization(input.c_str(), T, D);
    
    // load EEG dataste and Z-normaliazation
    float*** dataset = new float**[10000];
    for (i = 0; i<6000; i++) {
		string count = to_string(num);
		count = "../EEG/class_a/" + count + ".txt";
		num++;
		dataset[i] = load_data(count.c_str(), T, D);
	}
	num = 1;
	for (i = M/2; i<10000; i++) {
		string count = to_string(num);
		count = "../EEG/class_c/" + count + ".txt";
		num++;
		dataset[i] = load_data(count.c_str(), T, D);
	}
    float*** data=z_normalization(dataset,T, D);	
    cout << "finish normalization" << endl;
    float*** datasets = new float**[M];
    for(i=0;i<M;i++)
	datasets[i]=data[i];

   ////////////////////////////////////////////////////////////////////////////
	
    cout << "DTW_HD Test: " << DTW_HD(datasets[query_id], datasets[2]) << endl;	
    cout << "Euclidean distance Test: " << distance_HD(datasets[query_id], datasets[2])<< endl;
    endRQ = clock();
    time = double(endRQ - beginRQ) / CLOCKS_PER_SEC;
    cout << "The time for load data: " << time<< " seconds." << endl;
    
    cout << "finish z-normalization" << endl;
    

    

    float** vectorR=generateRandom_vector_r_w();
    int* seeds=generate_seed();

    //calculate hashvalue for query set
    vector<int> test=bit_profile_extraction_onevalue_hd(vectorR,datasets[query_id]);
    cout<<"Query Bit Profile Stream Size: "<<test.size()<<endl;
	vector<int> grams=generate_n_grams_onevalue(test);
	for(vector<int>::iterator it = grams.begin(); it != grams.end(); ++it){
        cout<<(*it)<<" ";
    }
    cout<<endl;
    cout<<"Query N Grams Size2 and Hash table: "<<grams.size()<<endl;
	
    vector<int> query_hashvalue=gen_minhash_value(grams, seeds);

    cout<<endl;
    
    cout<<"Query Hash value Size3: "<<query_hashvalue.size()<<endl;
    
    vector<int> grams1=generate_n_grams_onevalue(test);
	for(vector<int>::iterator it1 = grams1.begin(); it1 != grams1.end(); ++it1){
       // cout<<(*it1)<<" ";
    }

   //full_minhash_values_type1(datasets, vectorR, seeds);
    clock_t SSH_begin = clock();
    vector<int> SSH_candidate=SSH_generate_candidate(datasets[query_id], datasets,vectorR, seeds, 1);
    clock_t SSH_end = clock();
    double elapsed_secsSSH = double(SSH_end - SSH_begin) / CLOCKS_PER_SEC;
    cout<<"The time spent for SSH indexing HD: "<< elapsed_secsSSH<<" seconds."<<endl;
    cout<<"****************************************************************************"<<endl;


    
/*	//cout << "LB_PAA Test: " << compute_LB_PAA_HD(datasets[0], datasets[2]) << endl;
	clock_t DTW_onceBegin=clock();
	cout << "DTW_HD Test: " << DTW_HD(datasets[query_id], datasets[query_id]) << endl;
	clock_t DTW_onceEnd=clock();
	double DTW_once = double(DTW_onceEnd - DTW_onceBegin) / CLOCKS_PER_SEC;
	cout << "The time for DTW calculation once: " << DTW_once<< " seconds." << endl;
	//cout << "DTW_1D Test: " << DTW_1D(datasets[0], datasets[2]) << endl;
	cout << "Euclidean distance Test: " << distance_HD(datasets[query_id], datasets[query_id])<< endl;
	cout<<"****************************************************************************"<<endl;
    
 
	//DTW Range Query Ground Truth
	cout << "DTW Range Query Ground Truth: " << endl;
        beginRQ = clock();
	int count = 0;
	vector<int> DTW_groundtruth_Range = DTW_GroundTruth_Range(datasets[query_id], datasets);
	for (vector<int >::iterator it = DTW_groundtruth_Range.begin(); it != DTW_groundtruth_Range.end(); ++it) {
		cout << "Candidate series number for DTW range query ground truth: " << (*it)<< endl;
		count++;
	}
	cout << "The total number of series for DTW range query ground truth: " << count << endl;
	endRQ = clock();
	double elapsed_secsRQ = double(endRQ - beginRQ) / CLOCKS_PER_SEC;
	cout << "The time spent for DTW range query ground truth: " << elapsed_secsRQ << " seconds." << endl;
	//cout << "***************************************************************************" << endl;
*/

    //DTW KNN Query Ground Truth
    cout<<"DTW KNN Query Ground Truth: "<<endl;
    clock_t beginDTWKNN = clock();
    int countDTWKNN=0;
    vector<int> DTW_groundtruth_KNN=DTW_GroundTruth_KNN(datasets[query_id],datasets);
    for(vector<int>::iterator it=DTW_groundtruth_KNN.begin();it!=DTW_groundtruth_KNN.end();++it){
        cout<<"series number for DTW KNN Query ground truth: "<<(*it)<<endl;
        countDTWKNN++;
    }
    cout<<"The total number of candidate series for DTW KNN Query: "<<countDTWKNN<<endl;
    clock_t endDTWKNN = clock();
    double elapsed_secsDTWKNN = double(endDTWKNN - beginDTWKNN) / CLOCKS_PER_SEC;
    cout<<"The time spent for DTW KNN Query ground truth: "<< elapsed_secsDTWKNN<<" seconds."<<endl;
    cout<<"****************************************************************************"<<endl;
   
 
    //SSH KNN
    cout<<"SSH KNN: "<<endl;
    clock_t SSHKNN_begin = clock();
    vector<int> ssh_KNN=SSH_KNN(SSH_candidate,datasets[query_id],datasets);
    cout << "The number of SSH KNN: " << ssh_KNN.size()<< "." << endl;
    for (vector<int>::iterator it = ssh_KNN.begin(); it != ssh_KNN.end(); ++it) {
        //cout << "Candidate series number for SSH KNN: " << (*it) << endl;
        cout << (*it) << endl;
    }
    clock_t SSHKNN_end = clock();
    double elapsed_secs_SSHKNN = double(SSHKNN_end - SSHKNN_begin) / CLOCKS_PER_SEC;
    cout << "the accuracy of direct SSH KNN: " << accuracy(DTW_groundtruth_KNN, ssh_KNN) << endl;
    //cout << "the accuracy of direct SSH KNN Label classification: " << accuracy_KNN_classification(DTW_groundtruth_KNN, ssh_KNN, query_id, datasets) << endl;
    cout<<"The time spent for SSH KNN: "<< elapsed_secs_SSHKNN<<" seconds."<<endl;

  
	for (int i = 0; i < M; i++) {
		for (int j = 0; j<T; j++) {
			delete[]datasets[i][j];
		}
		delete[] datasets[i];
	}
	delete[] datasets;


	return 0;
}

	
	








