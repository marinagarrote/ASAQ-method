// functionsSAQ.cpp V.1.0
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

using namespace std;	
using namespace arma; // for 'armadillo';

const string s[4] = {"A", "C", "G", "T"};
const double epsilon = 1e-16;
vector <int> power4; 

vector<int> get_species()
{
	int aux=1; 
	int numb_sp=4; 
	for (int i=0;i<=numb_sp;i++) {	
		power4.push_back(aux); 
		aux=4*aux;
	}
	return(power4); 
}

Alignment readFASTA (string fname) { 	
	Alignment align;
	string line;
	// vector<string> species;
  fstream file;
  string speciesname;
  string dna="";
  string dnaaux;
  int numchar;
  int i;
  vector<int> al_length;
  //opening file
  file.open(fname.c_str(), fstream::in); //open(pointer to filename, type i.e. in/out)
 	if (file == NULL) {
 		cout <<"cannot open file " << fname.c_str() << "\n";
 		exit(1);
	}
  // get the first line from the FASTA file
  getline(file, line, '\n');//reads the first line
  if (line[0]!='>'){
	 	cout <<"error: not a FASTA file \n";
  	exit(0);
  }
  numchar=line.length()-1;
  speciesname=line.substr(1,numchar);
  align.taxa.push_back(speciesname);
  if (align.seqs[speciesname].length() > 0) {
		cerr << "Warning - found 2+ sequences with same name " << speciesname << endl;		
	}
	while(! file.eof()){
		getline(file, line, '\n');
	  
		if (line[0]!='>'){
   		dna += line;
   	}
		else{ 
   		align.seqs[speciesname]=dna;
   		al_length.push_back(dna.length());
   		numchar=line.length()-1;
   		speciesname=line.substr(1,numchar);
   		align.taxa.push_back(speciesname);
   		dna.clear();
   		if (align.seqs[speciesname].length() > 0) {
				cerr << "Warning - found 2+ sequences with same name " << speciesname << endl;		
			}
		}
	}
	align.seqs[speciesname]=dna;
	al_length.push_back(dna.length());
	file.close();
	align.num_taxa=align.taxa.size();
	if(al_length.size() != align.num_taxa)	cout << "error in getting the alignment!!!";
	for (i=0; i< al_length.size()-1; i++){
		if (al_length[i] != al_length[i+1]){
	   	cout << "sequence " <<align.taxa[i+1] <<" doesn't have same length \n";
	   	exit(0);
 		}
	}
	align.seq_len=al_length.back();
	return align;
}

// from an alignment, it obtains the relative frequencies of the columns
map<string, int> Getcolumns (Alignment &align){
	int i,numb_taxa,al_length=align.seq_len; 
	long long int numb_obs;
	float aux=1./al_length; 
	long int numb_col;
	string col;

	numb_taxa=align.num_taxa;

	map<string,string> myseqs=align.seqs;
	map<string,string>::const_iterator pos;
	map<long int,int>::const_iterator pattern;
	map<string, int> columnscount;
	string mystring;

	numb_obs=power4[numb_taxa];

	for (i=0; i < al_length; i++) {
		col="";

	for (int j=0;j<align.num_taxa;j++) {	
		mystring=myseqs[align.taxa[j]];
		col+=mystring[i];
	}

	//	for (pos = myseqs.begin(); pos != myseqs.end(); ++pos) {
  //   			col.push_back(pos->second[i]);
	//}
	columnscount[col]++;
	}
	return (columnscount);
}

mat Flattening(map<string,int> tensor, vector<int> split, vector <int> othersp)
{

 	int i,j,k,a,b; 
	long int coord_row,coord_col;
	map<long int,map<long int,int>> matrix;
	vector<long int> coord_it; 
	map<long int,int> mat_rows, mat_cols;

	int i_row=0;
	int i_col=0;

	for (map<string,int>::iterator it=tensor.begin(); it!=tensor.end(); ++it) {
		coord_it=Decomp(it->first,split,othersp);
		coord_row=coord_it[0];
		coord_col=coord_it[1];
		matrix[coord_row][coord_col]=it->second;

		if (mat_rows.find(coord_row)==mat_rows.end()) {
			i_row++; 
			mat_rows[coord_row]=i_row-1;
		}			
		if (mat_cols.find(coord_col)==mat_cols.end()) {
			i_col++; 
			mat_cols[coord_col]=i_col-1;
		}
	}
	map<long int,int>::iterator index_i;
	map<long int,int>::iterator index_j;

	long int n_rows=mat_rows.size(); 
	long int n_cols=mat_cols.size();
	mat flatmatrix=zeros<mat>(n_rows,n_cols);

	for (index_i=mat_rows.begin();index_i!=mat_rows.end();++index_i) {
		for (index_j=mat_cols.begin();index_j!=mat_cols.end();++index_j) {		
			flatmatrix(index_i->second,index_j->second)=matrix[index_i->first][index_j->first];
		}
	}
	return flatmatrix;
}

vector<long int> Decomp(string it, vector <int> split, vector <int> othersp) 
{
	int i; 
	int aa=split.size();
	int bb=othersp.size();;
	int numb=aa+bb;

	long int coord_row=0;
	long int coord_col=0;

	vector <long int> coords; 
	for(i=0;i<aa;i++) {
		coord_row+=(it[split[i]]-'0')*power4[numb-1-split[i]];
	}	
	for(i=0;i<bb;i++) {
		coord_col+=(it[othersp[i]]-'0')*power4[numb-1-othersp[i]];
	}	

	coords={coord_row,coord_col};
	return(coords);	
}

float Erik2(mat flat, int n_partitions, vector<int> split, vector<int> count_split)
{
	int i,k,j,a,b;
	mat flat_colcorrected,flat_rowcorrected;
	mat colsum,rowsum;
	int obs_row=flat.n_rows; 
	int obs_col=flat.n_cols; 
	vector<int> cols2remove; 
	vector<int> rows2remove; 

	a=split.size();
	b=count_split.size();
	mat unit_split = ones<mat>(1,obs_row);
	mat unit_othersp = ones<mat>(obs_col,1);
	vec singval_col, singval_row;
	colsum=unit_split*flat; 
	rowsum=flat*unit_othersp; 

	for (i=0;i<colsum.size();i++) {
		if (colsum[i]<=2) cols2remove.push_back(i); 
		else {colsum[i]=1./colsum[i];}
	}
	for (i=0;i<rowsum.size();i++) {
		if (rowsum[i]<=2) rows2remove.push_back(i); 
		else {rowsum[i]=1./rowsum[i];}
	}
	int n_cols2remove=cols2remove.size();
	int n_rows2remove=rows2remove.size();
	flat_colcorrected=flat*diagmat(colsum);
	flat_rowcorrected=diagmat(rowsum)*flat;
	for(j=n_cols2remove-1;j>=0;j--) {
		flat_colcorrected.shed_col(cols2remove[j]);
	}
	// cout << "col corrected " << endl; 
	// cout << flat_colcorrected << endl; 
	for(j=n_rows2remove-1;j>=0;j--) {
		flat_rowcorrected.shed_row(rows2remove[j]);
	}
	// flat_rowcorrected.shed_col(1);
	// cout << "row corrected " << endl; 
	// cout << flat_rowcorrected << endl; 

	singval_col=svd(flat_colcorrected);	
	// cout << "singular values columns: " << singval_col << endl; 
	singval_row=svd(flat_rowcorrected); 
	// cout << "singular values rows: " << singval_row << endl; 

	float tail_row=0;
	float tail_col=0;
	for(k=4*n_partitions;k<singval_col.n_elem;k++) {
		tail_col+=singval_col[k]*singval_col[k];
	}
	for(k=4*n_partitions;k<singval_row.n_elem;k++) {
		tail_row+=singval_row[k]*singval_row[k];
	}
	return((sqrt(tail_col)+sqrt(tail_row))/2);
}


T4F convert_tensor_4 (map <string, int> tensor, int n){
		T4F Mtensor(4, T3F(4, MF(4, VF(4,0))));
		for (int i = 0; i < 4; ++i){
			for (int j = 0; j < 4; ++j){
				for (int k = 0; k < 4; ++k){
					for (int l = 0; l < 4; ++l){
						string ss= s[i]+s[j]+s[k]+s[l];
						if ( tensor.find(ss) != tensor.end() ) Mtensor[i][j][k][l]= tensor[ss]/(double)n;
						//cerr << ss << ": "  << Mtensor[i][j][k][l] << " ";
					}
				}
			}
		}
		return Mtensor;
}

mat marginalization_3 (cube Mtensor, int m){ // m in[1,3] ----------P*m1-----------
	mat M = zeros<mat>(4,4);
	for (int i = 0; i < 4; ++i){
		for (int j = 0; j < 4; ++j){
			for (int k = 0; k < 4; ++k){
				if (m == 1) M(i,j)+=Mtensor(k,i,j);
				else if (m == 2)M(i,j)+=Mtensor(i,k,j);
				else if (m==3) M(i,j)+=Mtensor(i,j,k);
				else cout << "WARNING" << endl;
			}
		}
	}
	for (int i = 0; i < 4; ++i){
		for (int j = 0; j < 4; ++j){
			if( abs(M(i,j)) < epsilon) M(i,j) = 0;
		}
	}
	return M;
}

cube marginalization_4 (T4F Mtensor, int m){ // m in[1,4] ----------P*m1-----------
	//T3F T(4, MF(4, VF(4,0)));
	cube T = zeros<cube>(4,4,4);
	for (int i = 0; i < 4; ++i){
		for (int j = 0; j < 4; ++j){
			for (int k = 0; k < 4; ++k){
				for (int l = 0; l < 4; ++l){
					if (m == 1) T(i,j,k)+=Mtensor[l][i][j][k];
					else if (m == 2) T(i,j,k)+=Mtensor[i][l][j][k];
					else if (m == 3) T(i,j,k)+=Mtensor[i][j][l][k];
					else if (m == 4) T(i,j,k)+=Mtensor[i][j][k][l];
					else cout << "WARNING" << endl;
				}
			}
		}
	}
	
	for (int i = 0; i < 4; ++i){
		for (int j = 0; j < 4; ++j){
			for (int k = 0; k < 4; ++k){
				if(abs(T(i,j,k)) < epsilon) T(i,j,k) = 0;
			}
		}
	}
	return T;
}

MM double_marginalizations(T4F tensor){
	MM N(4,VM(4));
	for (int i = 0; i < 4; ++i){
		for (int j = 0; j < 4; ++j){
			if (i < j){ N[i][j] = marginalization_3(marginalization_4(tensor,numb(i+1,j+1,5)),numb(i+1,j+1,numb(i+1,j+1,5)));
				//~ cout << i + 1 << " " << j + 1  << " " << det(N[i][j]) << endl;
				}
			else if (i > j) N[i][j] = (N[j][i]).t();
			else N[i][j]= zeros<mat>(4,4);
		}
	}
	return(N);
}




vec SAQ_score(T4F P, vector<int> top, int k, double filter){

	mat Flat1 = Flattening_M(P, {top[0], top[1]}, {top[2], top[3]}); mat Flat1_2 = Flattening_M(P, {top[1], top[0]}, {top[2], top[3]}); 
	mat Flat2 = Flattening_M(P, {top[0], top[2]}, {top[1], top[3]});
	mat Flat3 = Flattening_M(P, {top[0], top[3]}, {top[1], top[2]});
		
	double normF = norm(Flat1, "fro");

// Projectem a les matrius SPSD i calculem la distancia a les matrius de rank k	dels flatenings
	double scores1_proj_rank = (projeccio_rankdist_sl(Flat1, k, false)/normF + projeccio_rankdist_sl(Flat1_2, k, false)/normF)/2;
	double scores2_proj_rank = projeccio_rankdist_sl(Flat2, k, false)/normF;
	double scores3_proj_rank = projeccio_rankdist_sl(Flat3, k, false)/normF;

  //double score_dif = scores1_proj_rank - min(scores2_proj_rank, scores3_proj_rank);
  double score_quoInv = min(scores2_proj_rank, scores3_proj_rank)/scores1_proj_rank;
  double score_quo = scores1_proj_rank/min(scores2_proj_rank, scores3_proj_rank);
  vec scores = {NAN, NAN};
	//if(Flat1.min() > -1) scores = {score_dif, score_quo};
  if(Flat1.min() > filter) scores = {score_quoInv, score_quo};
  return(scores);
}


vec SAQ(T4F Mtensor, vector<int> top, MM Ns, int k, double filter){
  vec_tensor P = tensor_transformations(Mtensor, top, Ns);
  double meanScoreDif = 0; 
  double meanWeightDif = 0;
  double cont = 0;
  vec scoreTrans_i;
  for(int i = 0; i < 16; ++i){
    scoreTrans_i = SAQ_score(P[i], top, k, filter);
    if(!std::isnan(scoreTrans_i[0])){
      ++cont,
      meanScoreDif = meanScoreDif + scoreTrans_i[0];
      meanWeightDif = meanWeightDif + scoreTrans_i[1];
    }
  }
  vec scores_SAQ = {meanScoreDif/cont, meanWeightDif/cont};
  return(scores_SAQ);
}
    
mat Flattening_M (T4F P, vector<int> split, vector <int> othersp){//Flattening Arbres de 4 fulles ---> 16x16
	mat flatt = zeros<mat>(16,16);
	vector <int> v(4,0);
	for (int i = 0; i < 16; ++i){
		for (int j = 0; j < 16; ++j){
			v[split[0]] = (int)i/4;
			v[split[1]] = (int)i%4;
			v[othersp[0]] = (int)j/4;
			v[othersp[1]] = (int)j%4;
			//~ if (abs(P[v[0]][v[1]][v[2]][v[3]]) < epsilon) flatt(i,j) = 0;
			//~ else flatt(i,j) = P[v[0]][v[1]][v[2]][v[3]]; 
			flatt(i,j) = P[v[0]][v[1]][v[2]][v[3]]; 
		}
	}
	return flatt;
}


double projeccio_rankdist_sl(mat A, int k, bool normalize){
	mat B = (A + trans(A))/2;
	vec eigval;
	mat eigvec;
	eig_sym(eigval, eigvec, B);
  
  //escriu_vec(eigval);
	double norm;
	if(normalize){
		norm = 0;
		for(int i = 0; i < eigval.n_elem; ++i){
			if(eigval(i) < epsilon) eigval(i) = 0;
			norm = norm + eigval(i)*eigval(i);
		}
			
		//for(int i = 0; i < eigval.n_elem; ++i) eigval(i) = eigval(i)/sum;
	}else{
		norm = 1;
		for(int i = 0; i < eigval.n_elem; ++i){
			if(eigval(i) < epsilon) eigval(i) = 0;
		}
	}
  //escriu_vec(eigval);
	double aux = 0;
	for(int i = 0; i < (16-k); ++i) aux += eigval(i)*eigval(i);
	
	aux = sqrt(aux)/sqrt(norm);
  //cout << aux << endl;
	return(aux);
}



int numb (int a, int b, int c){
	if (a!=4 and b!=4 and c!=4) return 4;
	else if (a!=3 and b!=3 and c!=3) return 3;
	else if (a!=2 and b!=2 and c!=2) return 2;
	else if (a!=1 and b!=1 and c!=1) return 1;
}


vec_tensor tensor_transformations(T4F tensor, vector<int> top, MM N){

	vector<int> top_1 = {top[0], top[1]};
	vector<int> top_2 = {top[2], top[3]};
	
	vec_tensor P (16, T4F(4, T3F(4, MF(4, VF(4,0)))));
	T4F zeroTens = T4F(4, T3F(4, MF(4, VF(4,0))));
 
	int cont = 0;
	for(int i = 0; i < 2; ++i){
		for (int j = 2; j < 4; ++j){	
			for(int k = 0; k < 2; ++k){
				for (int l = 2; l < 4; ++l){
										
          //cerr << "1  --  " << det(N[top[l]][top[i]]) << endl;
					if(abs(det(N[top[l]][top[i]])) < epsilon){
            P[cont] = zeroTens;
            ++ cont;
            continue;
          }
          mat Nli = solve(N[top[l]][top[i]], eye<mat>(4,4));
					T4F P1 = operador_ast_M(tensor, top[i]+1, Nli*N[top[l]][top[1-i]]);
   	      
					if(abs(det(N[top[k]][top[j]])) < epsilon){
            P[cont] = zeroTens;
            ++ cont;
            continue;
          }
          //cerr << "2  --  " << det(N[top[k]][top[j]]) << endl;
					mat Nkj = solve(N[top[k]][top[j]], eye<mat>(4,4));
					T4F P2 = operador_ast_M(P1, top[j]+1, Nkj*N[top[k]][top[5-j]]);
					
					P[cont] = P2;
					++cont;
				}
			}
		}
	}
	return(P);
}


T4F operador_ast_M (T4F P, int x, mat M){

	T4F T(4, T3F(4, MF(4, VF(4,0))));
	if (x == 1){
		for (int i = 0; i < 4; ++i){
			for (int j = 0; j < 4; ++j){
				for (int k = 0; k < 4; ++k){
					for (int l = 0; l < 4; ++l){
						for (int m = 0; m < 4; ++m){
							//~ if(abs(P[m][j][k][l]*M(m,i)) >= epsilon)
							 T[i][j][k][l] += P[m][j][k][l]*M(m,i); 
							//~ T[i][j][k][l] += P[m][j][k][l]*M(m,i); 
						}
					}
				}
			}
		}
	}
	
	else if (x == 2){
		for (int i = 0; i < 4; ++i){
			for (int j = 0; j < 4; ++j){
				for (int k = 0; k < 4; ++k){
					for (int l = 0; l < 4; ++l){
						for (int m = 0; m < 4; ++m){
							//~ if(abs(P[i][m][k][l]*M(m,j)) >= epsilon) 
							T[i][j][k][l] += P[i][m][k][l]*M(m,j); 
							//~ T[i][j][k][l] += P[i][m][k][l]*M(m,j);
						}
					}
				}
			}
		}
	}
	
	else if (x ==3){
		for (int i = 0; i < 4; ++i){
			for (int j = 0; j < 4; ++j){
				for (int k = 0; k < 4; ++k){
					for (int l = 0; l < 4; ++l){
						for (int m = 0; m < 4; ++m){
							//~ if(abs(P[i][j][m][l]*M(m,k)) >= epsilon) 
							T[i][j][k][l]+= P[i][j][m][l]*M(m,k); 
							//~ T[i][j][k][l]+= P[i][j][m][l]*M(m,k); 
						}
					}
				}
			}
		}
	}

	else if ( x== 4 ) {
		for (int i = 0; i < 4; ++i){
			for (int j = 0; j < 4; ++j){
				for (int k = 0; k < 4; ++k){
					for (int l = 0; l < 4; ++l){
						for (int m = 0; m < 4; ++m){
							//~ if(abs(P[i][j][k][m]*M(m,l)) >= epsilon)
							 T[i][j][k][l] += P[i][j][k][m]*M(m,l); 
							 //~ T[i][j][k][l] += P[i][j][k][m]*M(m,l); 
						}
					}
				}
			}
		}
	}
	else cerr << "WARNING!!!" << endl;
	
	//~ for (int i = 0; i < 4; ++i){
		//~ for (int j = 0; j < 4; ++j){
			//~ for (int k = 0; k < 4; ++k){
				//~ for (int l = 0; l < 4; ++l){
					//~ if(abs(T[i][j][k][l]) < epsilon) T[i][j][k][l] = 0;
				//~ }
			//~ }
		//~ }
	//~ }
	return T;
}

string convertInt(int number) {
   std::stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}


VM logDetNs(MM N){
   mat logDet(4,4, fill::zeros); 
   mat signs(4,4, fill::zeros);
   double sign12; log_det(logDet(0,1), signs(0,1), N[0][1]);
   double sign13; log_det(logDet(0,2), signs(0,2), N[0][2]);
   double sign14; log_det(logDet(0,3), signs(0,3), N[0][3]);
   double sign23; log_det(logDet(1,2), signs(1,2), N[1][2]);
   double sign24; log_det(logDet(1,3), signs(1,3), N[1][3]);
   double sign34; log_det(logDet(2,3), signs(2,3), N[2][3]);
   logDet(1,0) = logDet(0,1); logDet(2,0) = logDet(0,2); logDet(3,0) = logDet(0,3);
   logDet(2,1) = logDet(1,2); logDet(3,1) = logDet(1,3); logDet(3,2) = logDet(2,3);
   signs(1,0) = signs(0,1); signs(2,0) = signs(0,2); signs(3,0) = signs(0,3);
   signs(2,1) = signs(1,2); signs(3,1) = signs(1,3); signs(3,2) = signs(2,3);
   
   //cerr << logDet(0,1) << " " << logDet(0,2) << " " << logDet(0,3) << " " << logDet(1,2) << " " << logDet(1,3) << " " << logDet(2,3) << endl;
   VM output = {logDet, signs};
   return(output);
}

mat condNs(MM N){
  mat condMat(4,4, fill::zeros);  
  condMat(0,1) = cond(N[0][1]); condMat(0,2) = cond(N[0][2]);  condMat(0,3) = cond(N[0][3]);
  condMat(1,2) = cond(N[1][2]); condMat(1,3) = cond(N[1][3]);  condMat(2,3) = cond(N[2][3]);
  //cerr << condMat(0,1) << " " << condMat(0,2) << " " << condMat(0,3) << " " << condMat(1,2) << " " << condMat(1,3) << " " << condMat(2,3) << endl;
  return(condMat);
}

bool useErik(MM Ns, mat logDet, mat cond, mat signs, vector<int> top, double tol){
  if(cond(0,1) > tol || cond(0,2) > tol || cond(0,3) > tol || cond(1,2) > tol || cond(1,3) > tol || cond(2,3) > tol) return(false);
  
  double s1, sign1;
  s1 = logDet(top[0],top[2]) + logDet(top[1],top[3]) - logDet(top[0],top[1]) - logDet(top[2],top[3]); 
  sign1 = signs(top[0],top[2])*signs(top[1],top[3])*signs(top[0],top[1])*signs(top[2],top[3]);
  //cerr << s1 << " " << sign1 << endl;
  if(s1 > 0) return(false); // || sign1 != 1
  
  double s2, sign2;
  s2 = logDet(top[0],top[3]) + logDet(top[1],top[2]) - logDet(top[0],top[1]) - logDet(top[2],top[3]);
  sign2 = signs(top[0],top[3])*signs(top[1],top[2])*signs(top[0],top[1])*signs(top[2],top[3]);
  //cerr << s2 << " " << sign2 << endl;
  if(s2 > 0) return(false); //|| sign2 != 1
  
  return(true);
}


double estGamma(mat logDet, vector<int> top){
  double s1, s2;
  s1 = logDet(top[0],top[2]) + logDet(top[1],top[3]) - logDet(top[0],top[1]) - logDet(top[2],top[3]); 
  s2 = logDet(top[0],top[3]) + logDet(top[1],top[2]) - logDet(top[0],top[1]) - logDet(top[2],top[3]);
  return(max(s1,s2));
}
  
  


  
  
  





















