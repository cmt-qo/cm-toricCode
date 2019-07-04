/* 
Author: Agnes Valenti
this file contains useful functions for vector operations
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <string>
#include <sstream>
#include <random>
#include <iomanip>
#include <limits>
#include <ctime>
using namespace std;


//useful functions for vector operations

vector<double> operator* (const double a,vector<double> state_){        
  #pragma omp parallel for num_threads(72)
  for (int i=0; i<state_.size(); i++){
     state_[i]=a*state_[i];}
  return state_;
 }

void pluseq(vector<double> & state_,vector<double> state2){
  #pragma omp parallel for num_threads(72)
  for (int i=0; i<state_.size(); i++){
    state_[i]=state_[i]+state2[i];}
 }

void mineq(vector<double> & state_,vector<double> state2){
  #pragma omp parallel for num_threads(72)
  for (int i=0; i<state_.size(); i++){
    state_[i]=state_[i]-state2[i];}
 }
