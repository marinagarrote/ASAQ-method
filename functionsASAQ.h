// functionsSAQ.h V.1.0
// Date: Dec/2020

// This file is part of the "SAQ" program
/*
This code has been written by  M. Garrote-López. 
Please quote the paper: M. Casanellas, J. Fernández-Sánchez, M. Garrote-López
"SAQ: semi-algebraic quartet reconstruction method", https://arxiv.org/abs/2011.13968
*/

#if ! defined (FUNCTIONS_SAQ_H) 
#define FUNCTIONS_SAQ_H

#include "typedefs.h"
#include <cmath>
#include <cstdio>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <armadillo>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace arma;

struct Alignment {
	unsigned int num_taxa;
	unsigned int seq_len;
	vector<string> taxa; // the taxa names
  map<string,string> seqs;  // the sequences
};
vector<int> get_species();
Alignment readFASTA (string fname);
map<string, int> Getcolumns (Alignment &align);
T4F convert_tensor_4 (map <string, int> tensor, int n);
mat marginalization_3 (cube Mtensor, int m);
cube marginalization_4 (T4F Mtensor, int m);
MM double_marginalizations(T4F tensor);
vec SAQ_score(T4F P, vector<int> top, int k, double filter);
vec SAQ(T4F Mtensor, vector<int> top, MM Ns, int k, double filter);
mat Flattening_M (T4F P, vector<int> split, vector <int> othersp);
double projeccio_rankdist_sl(mat A, int k, bool normalize);
int numb (int a, int b, int c);
vec_tensor tensor_transformations(T4F tensor, vector<int> top, MM N);
T4F operador_ast_M (T4F P, int x, mat M);
string convertInt(int number);
mat Flattening(map<string,int> tensor, vector<int> split, vector <int> othersp);
vector<long int> Decomp(string it, vector <int> split, vector <int> othersp);
float Erik2(mat flat, int n_partitions, vector<int> split, vector<int> count_split);
VM logDetNs(MM N);
mat condNs(MM N);
bool useErik(MM Ns, mat logDet, mat cond, mat signs, vector<int> top, double tol);
double estGamma(mat logDet, vector<int> top);
#endif
