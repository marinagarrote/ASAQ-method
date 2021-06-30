// SAQ.cpp V.1.0
// Date: Dec/2020

// This file is part of the "SAQ" program
/*
This code has been written by  M. Garrote-López. 
Please quote the paper: M. Casanellas, J. Fernández-Sánchez, M. Garrote-López
"SAQ: semi-algebraic quartet reconstruction method", https://arxiv.org/abs/2011.13968
*/

#include "functionsASAQ.h"
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <armadillo>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

using namespace std;	
using namespace arma; // for 'armadillo';

const string s[] = {"A", "C", "G", "T"};
const double epsilon = 1e-16;

int main(int argc, char **argv)
{
  // INPUT PATH
  if (argc < 2) {
  	cerr << "Missing filename" << endl;
  	return 1;
  }
  string filename = argv[1];
  double filter, tol;
  
  if (argv[2] == NULL){ filter = -1;
  }else{filter = stod(argv[2]);}
  if (argv[3] == NULL){ tol = 5000;
  }else{tol = stod(argv[3]);}
  

	int i,j,k,q;
	int numb_taxa=4;
	vector <int> power4=get_species(); 
	map<string,string>::const_iterator pos;

  
  // READ DATA
  Alignment align = readFASTA(filename);
  map<string, int> tensor = Getcolumns(align);
  
  map <string,string> dna=align.seqs; 
	map<string,string> myseqs=align.seqs;
	tensor=Getcolumns(align);
	int length=align.seq_len;	
	string st_length=convertInt(length);
  
  // Erik+2 method 
		
	mat flat_01 = Flattening(tensor,{0,1},{2,3});
	mat flat_02 = Flattening(tensor,{0,2},{1,3});
	mat flat_03 = Flattening(tensor,{0,3},{1,2});
 
	double erik2_01 = Erik2(flat_01, 2, {0,1},{2,3});
	double erik2_02 = Erik2(flat_02, 2, {0,2},{1,3});
	double erik2_03 = Erik2(flat_03, 2, {0,3},{1,2});
  
  //cerr << erik2_01 << " " << erik2_02 << " " << erik2_03 << endl;

	// ASAQ method 
	T4F Mtensor = convert_tensor_4(tensor, length);
	T4F zeroTens = T4F(4, T3F(4, MF(4, VF(4,0))));
	int rank_matrix = 4; double eps = 1e-7;
  
  MM Ns = double_marginalizations(Mtensor);
  mat cond = condNs(Ns);
  VM LS = logDetNs(Ns);
  mat logDet = LS[0];
  mat signs = LS[1];
  
  
  bool useErik2;
   
  if(erik2_01 < min(erik2_02, erik2_03)){
    useErik2 = useErik(Ns, logDet, cond, signs, {0,1,2,3}, tol);
  }else if(erik2_02 < min(erik2_01, erik2_03)){
    useErik2 = useErik(Ns, logDet, cond, signs, {0,2,1,3}, tol);
  }else{
    useErik2 = useErik(Ns, logDet, cond, signs, {0,3,1,2}, tol);
  }
  
  string sp0,sp1,sp2,sp3;
  sp0 = align.taxa[0];
  sp1 = align.taxa[1];
  sp2 = align.taxa[2];
  sp3 = align.taxa[3];
   
  cout << endl; 
  cout << "Taxa 1: " << sp0 << endl; 
  cout << "Taxa 2: " << sp1 << endl; 
  cout << "Taxa 3: " << sp2 << endl; 
  cout << "Taxa 4: " << sp3 << endl; 
  cout << "Topologies : \t" <<  "12|34 \t \t" << " 13|24 \t \t" << " 14|23" << endl;
    
  if(useErik2){
  
    double erik2_01_inv, erik2_02_inv, erik2_03_inv;
    
    if(erik2_01==0){
      erik2_01_inv = 1; erik2_01_inv = 1; erik2_01_inv = 1;
    }else if(erik2_02==0){
      erik2_02_inv = 1; erik2_02_inv = 1; erik2_02_inv = 1;
    }else if(erik2_03==0){
      erik2_03_inv = 1; erik2_03_inv = 1; erik2_03_inv = 1;
    }else{
      erik2_01_inv = 1./erik2_01;
      erik2_02_inv = 1./erik2_02;
      erik2_03_inv = 1./erik2_03;
    }   
    
	  double sumERIK2 = erik2_01_inv + erik2_02_inv + erik2_03_inv;
     
    cout << "ASAQ weights: \t " << erik2_01_inv/sumERIK2 << "\t \t" <<  erik2_02_inv/sumERIK2 << "\t \t" <<  erik2_03_inv/sumERIK2 << endl;  
    cout << "Method  used: Erik+2" << endl;
  }else{
  
    vec SAQ_01 = SAQ(Mtensor, {0,1,2,3}, Ns, rank_matrix, filter);
  	vec SAQ_02 = SAQ(Mtensor, {0,2,1,3}, Ns, rank_matrix, filter);
    vec SAQ_03 = SAQ(Mtensor, {0,3,1,2}, Ns, rank_matrix, filter);
   
    double score_01, score_02, score_03;
    
    score_01 = SAQ_01[0];
    score_02 = SAQ_02[0];
    score_03 = SAQ_03[0];
    
  
    double sumSAQ = score_01 + score_02 + score_03; 
  	cout << "ASAQ weights: \t" << score_01/sumSAQ << "\t \t" << score_02/sumSAQ << "\t \t" <<  score_03/sumSAQ << endl;   
    cout << "Method  used: SAQ" << endl;              
  }
  return(0);
}












