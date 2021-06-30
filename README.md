# ASAQ: Algebraic and semi-algebraic quartet reconstruction method

ASAQ V.1.0
Update: June/202a

This code has been written by M. Garrote-López and designed by M. Casanellas, J. Fernández-Sánchez and M. Garrote-López.

## Description
The ASAQ method has been implemented in C++ to be run in a UNIX terminal. Given an alignment of four sequences, ASAQ computes three weights, one for each possible tree topology. 


### Input

The input of the method is a FASTA file containing the alignment of four species and the parameter _filter_ which is optional.


| Input | Description | Example |
|:------------- |:------------- | :-----: |
| fastaFile | Path to the fasta file with the alignment | /home/user/data/alignment.fa |
| filter (optional)| Minimum value accepted on the entries of the vectors obtained by a _T_-leaf transformation on the distribution obtained from the alignment. <br /> The default value is _filter = -1_ | -0.5 |
| threshold (optional)| Maximum condition number of N_ij accepted in the computation of I_AB. <br /> The default value is _threshold = 5000_ | 1000 |


### Output

The program outputs three weights provided by the method SAQ applied to the alignment. Each weight is computed according to one possible topology relating the sequences of the alignment _seq1, seq2, seq3, seq4_:

* Weight of tree _seq1,seq2 | seq3,seq4_
* Weight of tree _seq1,seq3 | seq2,seq4_
* Weight of tree _seq1,seq4 | seq2,seq3_ 


## Usage 

The executable files ASAQ have been written in C++ code and compiled using the c++11 compiler. To run ASAQ it is necessary to have installed the C++ linear algebra library Armadillo (http://arma.sourceforge.net/). 

### Compiling SAQ
```
  g++ -std=c++11 -c functionsASAQ.cpp -O1 -larmadillo
  g++ -std=c++11 ASAQ.cpp functionsASAQ.o -o ASAQ -O1 -larmadillo
```
    
### Running an experiment
```
  ./ASAQ fastaFile filter
```
The parameter _filter_ can be omitted. 

### Examples:
```
  ./ASAQ /home/user/data/alignment.fa -0.5 1000
  ./ASAQ /home/user/data/alignment.fa
```
