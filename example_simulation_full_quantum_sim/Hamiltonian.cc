/*
Author: Agnes Valenti
This file contains the class 'Hamiltonian'
*/

#include "vectoroperations.cc"
using namespace std;


class Hamiltonian{
   //class implementing Hamiltonian and other operators and operations on a full quantum state.
   //The full quantum state is not stored as member variable, but passed as argument to functions.
   //The class loads field parameters from the files filenameA, filenameB and stores the fields in
   //the member variables lambdaA_, lambdaB


   //system length
   const int L;

   //number of spinconfigurations, n_=2**(2*L*L-1)
   unsigned long long int n_;

   //field configurations: 2*L*L elements
   //beta_z - fieldconfiguration
   vector<double> lambdaA_;
   //beta_x - fieldconfiguration
   vector<double> lambdaB_;
   
   //useful quantities for computation
   vector<unsigned long long int> pow2;


public:
   Hamiltonian(string filenameA, string filenameB,int length=2):L(length){     
     InitializePow2();
     LoadParameters(filenameA,filenameB);
     n_=pow2[2*L*L];
   }

 void InitializePow2(){
   cout<<"Initializing pow2..."<<endl;
   pow2.resize(2*L*L+1);
   for (int i=0; i<2*L*L+1; i++){
     if (i==0)
       pow2[0]=1;
     else
       pow2[i]= pow2[i-1]*2;
   }
   cout<<"Finished initializing pow2"<<endl;
 }

 //swap_spins: function used to implement sigma_x
 //spinconf: store full quantum state as vector. 
 //At index ''spinconf'', the basis vector corresponds to 
 //the spinconfiguration equal to the binary representation of spinconf
 //swap_spins: returns (as spinconf) the configuration 
 //(=index of the configuration) resulting from flipping spin spin_i in configuration spinconf
 inline void swap_spins(unsigned long long int & spinconf, int spin_i){
   unsigned long long int spinconf1=spinconf;
   spinconf1=spinconf1%(pow2[2*L*L-spin_i]);
   spinconf1=spinconf1/pow2[2*L*L-(spin_i+1)];  
   spinconf1=1-2*spinconf1;  
   spinconf=spinconf+spinconf1*pow2[2*L*L-(spin_i+1)];
 }

 //sign: function used to implement sigma_z
 //returns the eigenvalue of sigma_z at spin_i of the spinconfiguration spinconf: 
 //If spin is down at spin_i (coefficient 1 in binary repr.), it returns -1. 
 //If spin is up at spin_i (coefficient 0 in binary repr.), it returns 1
 inline int sign(unsigned long long int & spinconf, int spin_i){
   unsigned long long int spinconf1=spinconf;
   spinconf1=spinconf1%(pow2[2*L*L-spin_i]);
   spinconf1=spinconf1/pow2[2*L*L-(spin_i+1)];  //a_i
   spinconf1=1-2*spinconf1;  //sign
   return spinconf1;
 }

 //implements sigma_x on state ''state_''
 //void version
 inline void vsigma_x(vector<double> & state_,int spin){ 
   vector<double> initialstate_=state_;            
   #pragma omp parallel for num_threads(72)
   for(unsigned long long int i=0; i<n_; i++){                
      unsigned long long int i2=i;
      swap_spins(i2, spin);
      state_[i]=initialstate_[i2]; 
      }
 }

 //implements sigma_z on state ''state_''
 //void version
 inline void vsigma_z(vector<double> & state_,int spin){ 
   #pragma omp parallel for num_threads(72)
   for(unsigned long long int i=0; i<n_; i++){                
      int signc;
      signc=sign(i, spin);
      state_[i]=state_[i]*signc;}
 }
 
 //implements sigma_x on state ''state_'', without altering initialstate
 //returns state
 inline vector<double> sigma_x(vector<double> & initialstate_,int spin){ 
   vector<double> state_=initialstate_;            
   #pragma omp parallel for num_threads(72)
   for(unsigned long long int i=0; i<n_; i++){         
      unsigned long long int i2;
      i2=i;
      swap_spins(i2, spin);
      state_[i]=initialstate_[i2]; 
      }
 return state_;
 }
 
 //implements sigma_x on state ''state_'',
 //returns state
 inline vector<double> sigma_z(vector<double> state_,int spin){
   vsigma_z(state_,spin);
 return state_;
 }

 //implements tensor product of sigma_z on spins spin1 and spin2, returns state
 inline vector<double> sigma_z(vector<double> state_,int spin1,int spin2){             
   vsigma_z(state_,spin1);
   vsigma_z(state_,spin2);
 return state_;
 }

 //implements tensor product of sigma_z on spins spin1, spin2 and spin3, returns state
 inline vector<double> sigma_z(vector<double> state_,int spin1,int spin2, int spin3){           
   vsigma_z(state_,spin1);
   vsigma_z(state_,spin2);
   vsigma_z(state_,spin3);
 return state_;
 }

 //implements tensor product of sigma_z on spins spin1, spin2, spin3 and spin4, returns state
 inline vector<double> sigma_z(vector<double> state_,int spin1,int spin2, int spin3, int spin4){            
   vsigma_z(state_,spin1);
   vsigma_z(state_,spin2);
   vsigma_z(state_,spin3);
   vsigma_z(state_,spin4);
 return state_;
 }

 //implements tensor product of sigma_x on spins spin1 and spin2, returns state
 inline vector<double> sigma_x(vector<double> state_,int spin1,int spin2){             
   vsigma_x(state_,spin1);
   vsigma_x(state_,spin2);
 return state_;
 }

 //implements tensor product of sigma_x on spins spin1, spin2 and spin3, returns state
 inline vector<double> sigma_x(vector<double> state_,int spin1,int spin2, int spin3){           
   vsigma_x(state_,spin1);
   vsigma_x(state_,spin2);
   vsigma_x(state_,spin3);
 return state_;
 }

 //implements tensor product of sigma_x on spins spin1, spin2, spin3 and spin4, returns state
 inline vector<double> sigma_x(vector<double> state_,int spin1,int spin2, int spin3, int spin4){             
   vsigma_x(state_,spin1);
   vsigma_x(state_,spin2);
   vsigma_x(state_,spin3);
   vsigma_x(state_,spin4);
 return state_;
 }

 //implements vertex stabilizer operator on state (at vertex), does not alter state, returns state
 vector<double> A_s(vector<double> state_,int vertex){ 
   int l,o,r,u;
   get_indices_vertex(vertex,l,o,r,u);
   vsigma_x(state_,l);
   vsigma_x(state_,o);
   vsigma_x(state_,r);
   vsigma_x(state_,u);
 return state_;
 }

 //implements plaquette stabilizer operator on state (at plaquette), does not alter state, returns state
 vector<double> B_p(vector<double> state_,int plaquette){ 
   int l,o,r,u;
   get_indices_plaquette(plaquette,l,o,r,u);
   vsigma_z(state_,l);
   vsigma_z(state_,o);
   vsigma_z(state_,r);
   vsigma_z(state_,u);
 return state_;
 }

 //implements sigma_z-exponential on state, does not alter state
 vector<double> exp_vertex(const vector<double> & state_,int vertex){
   int l,o,r,u;
   get_indices_vertex(vertex,l,o,r,u);
   double lbda1=lambdaA_[l];
   double lbda2=lambdaA_[o];
   double lbda3=lambdaA_[r];
   double lbda4=lambdaA_[u];
   double alpha0=+cosh(lbda1)*cosh(lbda2)*cosh(lbda3)*cosh(lbda4);  
   double alpha1=-sinh(lbda1)*cosh(lbda2)*cosh(lbda3)*cosh(lbda4); 
   double alpha2=-cosh(lbda1)*sinh(lbda2)*cosh(lbda3)*cosh(lbda4); 
   double alpha3=-cosh(lbda1)*cosh(lbda2)*sinh(lbda3)*cosh(lbda4); 
   double alpha4=-cosh(lbda1)*cosh(lbda2)*cosh(lbda3)*sinh(lbda4); 
   double alpha5=+sinh(lbda1)*sinh(lbda2)*cosh(lbda3)*cosh(lbda4); 
   double alpha6=+sinh(lbda1)*cosh(lbda2)*sinh(lbda3)*cosh(lbda4); 
   double alpha7=+sinh(lbda1)*cosh(lbda2)*cosh(lbda3)*sinh(lbda4); 
   double alpha8=+cosh(lbda1)*sinh(lbda2)*sinh(lbda3)*cosh(lbda4); 
   double alpha9=+cosh(lbda1)*sinh(lbda2)*cosh(lbda3)*sinh(lbda4); 
   double alpha10=+cosh(lbda1)*cosh(lbda2)*sinh(lbda3)*sinh(lbda4);
   double alpha11=-sinh(lbda1)*sinh(lbda2)*sinh(lbda3)*cosh(lbda4); 
   double alpha12=-sinh(lbda1)*sinh(lbda2)*cosh(lbda3)*sinh(lbda4); 
   double alpha13=-sinh(lbda1)*cosh(lbda2)*sinh(lbda3)*sinh(lbda4); 
   double alpha14=-cosh(lbda1)*sinh(lbda2)*sinh(lbda3)*sinh(lbda4);
   double alpha15=+sinh(lbda1)*sinh(lbda2)*sinh(lbda3)*sinh(lbda4);
   vector <double> statesum_=alpha0*state_;
   pluseq(statesum_,alpha1*sigma_z(state_,l));       
   pluseq(statesum_,alpha2*sigma_z(state_,o));
   pluseq(statesum_,alpha3*sigma_z(state_,r));
   pluseq(statesum_,alpha4*sigma_z(state_,u));
   pluseq(statesum_,alpha5*sigma_z(state_,l,o));
   pluseq(statesum_,alpha6*sigma_z(state_,l,r));
   pluseq(statesum_,alpha7*sigma_z(state_,l,u));
   pluseq(statesum_,alpha8*sigma_z(state_,o,r));
   pluseq(statesum_,alpha9*sigma_z(state_,o,u));
   pluseq(statesum_,alpha10*sigma_z(state_,r,u));
   pluseq(statesum_,alpha11*sigma_z(state_,l,o,r));
   pluseq(statesum_,alpha12*sigma_z(state_,l,o,u));
   pluseq(statesum_,alpha13*sigma_z(state_,l,r,u));
   pluseq(statesum_,alpha14*sigma_z(state_,o,r,u));
   pluseq(statesum_,alpha15*sigma_z(state_,l,o,r,u));
   return statesum_;
 }

 //implements sigma_x-exponential on state, does not alter state
 vector<double> exp_plaquette(vector<double> & state_,int plaquette){
   int l,o,r,u;
   get_indices_plaquette(plaquette,l,o,r,u);
   double lbda1=lambdaB_[l];
   double lbda2=lambdaB_[o];
   double lbda3=lambdaB_[r];
   double lbda4=lambdaB_[u];
   double alpha0=+cosh(lbda1)*cosh(lbda2)*cosh(lbda3)*cosh(lbda4);  
   double alpha1=-sinh(lbda1)*cosh(lbda2)*cosh(lbda3)*cosh(lbda4); 
   double alpha2=-cosh(lbda1)*sinh(lbda2)*cosh(lbda3)*cosh(lbda4); 
   double alpha3=-cosh(lbda1)*cosh(lbda2)*sinh(lbda3)*cosh(lbda4); 
   double alpha4=-cosh(lbda1)*cosh(lbda2)*cosh(lbda3)*sinh(lbda4); 
   double alpha5=+sinh(lbda1)*sinh(lbda2)*cosh(lbda3)*cosh(lbda4); 
   double alpha6=+sinh(lbda1)*cosh(lbda2)*sinh(lbda3)*cosh(lbda4); 
   double alpha7=+sinh(lbda1)*cosh(lbda2)*cosh(lbda3)*sinh(lbda4); 
   double alpha8=+cosh(lbda1)*sinh(lbda2)*sinh(lbda3)*cosh(lbda4); 
   double alpha9=+cosh(lbda1)*sinh(lbda2)*cosh(lbda3)*sinh(lbda4); 
   double alpha10=+cosh(lbda1)*cosh(lbda2)*sinh(lbda3)*sinh(lbda4);
   double alpha11=-sinh(lbda1)*sinh(lbda2)*sinh(lbda3)*cosh(lbda4); 
   double alpha12=-sinh(lbda1)*sinh(lbda2)*cosh(lbda3)*sinh(lbda4); 
   double alpha13=-sinh(lbda1)*cosh(lbda2)*sinh(lbda3)*sinh(lbda4); 
   double alpha14=-cosh(lbda1)*sinh(lbda2)*sinh(lbda3)*sinh(lbda4);
   double alpha15=+sinh(lbda1)*sinh(lbda2)*sinh(lbda3)*sinh(lbda4);
   vector <double> statesum_=alpha0*state_;
   pluseq(statesum_,alpha1*sigma_x(state_,l));        
   pluseq(statesum_,alpha2*sigma_x(state_,o));
   pluseq(statesum_,alpha3*sigma_x(state_,r));
   pluseq(statesum_,alpha4*sigma_x(state_,u));
   pluseq(statesum_,alpha5*sigma_x(state_,l,o));
   pluseq(statesum_,alpha6*sigma_x(state_,l,r));
   pluseq(statesum_,alpha7*sigma_x(state_,l,u));
   pluseq(statesum_,alpha8*sigma_x(state_,o,r));
   pluseq(statesum_,alpha9*sigma_x(state_,o,u));
   pluseq(statesum_,alpha10*sigma_x(state_,r,u));
   pluseq(statesum_,alpha11*sigma_x(state_,l,o,r));
   pluseq(statesum_,alpha12*sigma_x(state_,l,o,u));
   pluseq(statesum_,alpha13*sigma_x(state_,l,r,u));
   pluseq(statesum_,alpha14*sigma_x(state_,o,r,u));
   pluseq(statesum_,alpha15*sigma_x(state_,l,o,r,u));
   return statesum_;
 }

 //applies sigma_z-exponential as in exactly solvable ground states (not Hamiltonian form)
 void apply_expA(vector<double> & state_){
   cout<<"Apply expA"<<endl;
   for (int i=0; i<2*L*L; i++){
     vector<double> statesum_=state_;
     double lbda1=cosh(0.5*lambdaA_[i]);
     double lbda2=sinh(0.5*lambdaA_[i]);
     statesum_=lbda1*state_;       
     pluseq(statesum_,lbda2*sigma_z(state_,i));  
     state_=statesum_;
     //cout<<"i: "<<i<<endl;
   }
 normalize(state_);
 }

 //applies sigma_x-exponential as in exactly solvable ground states (not Hamiltonian form)
 void apply_expB(vector<double> & state_){
   cout<<"Apply expB"<<endl;
   for (int i=0; i<2*L*L; i++){
     vector<double> statesum_=state_;
     double lbda1=cosh(0.5*lambdaB_[i]);
     double lbda2=sinh(0.5*lambdaB_[i]);
     statesum_=lbda1*state_;       
     pluseq(statesum_,lbda2*sigma_x(state_,i));  
     state_=statesum_;
   }
 normalize(state_);
 }

 //applies Hamiltonian to state, with offset -8*L*L and factor of (-1) to ensure that the ground state has the maximal energy
 //caution regarding numerical stability: if fields are high, offset might have to set to a larger value (if not, method might converge to an excited state with larger absolute value of energy)
 void apply_H(vector<double> & state_){
   vector<double> initialstate_=state_;
   for(int i=0; i<L*L; i++){
     if (i==0)
       state_=double(-1)*A_s(initialstate_,i);
     else
       mineq(state_, A_s(initialstate_,i));
     pluseq(state_, exp_vertex(initialstate_,i));
     mineq(state_, B_p(initialstate_,i));
     pluseq(state_, exp_plaquette(initialstate_,i));
     pluseq(state_, double(-8*L*L)*initialstate_);
   }
   state_=double(-1)*state_;
 }

 //applies Hamiltonian to state, with offset -8*L*L and factor of (-1) to ensure that the ground state has the maximal energy
   //caution regarding numerical stability: if fields are high, offset might have to set to a larger value (if not, method might converge to an excited state with larger absolute value of energy)
   //allows magnetic fields in x-direction of stength h
 void apply_H_noisy(vector<double> & state_,double h=0.0){
   vector<double> initialstate_=state_;
   for(int i=0; i<L*L; i++){
     if (i==0)
       state_=double(-1)*A_s(initialstate_,i);
     else
       mineq(state_, A_s(initialstate_,i));
     pluseq(state_, exp_vertex(initialstate_,i));
     mineq(state_, B_p(initialstate_,i));
     pluseq(state_, exp_plaquette(initialstate_,i));
     pluseq(state_,h*sigma_x(initialstate_,i));
     pluseq(state_, double(-8*L*L)*initialstate_);
   }
   state_=double(-1)*state_;
 }

 //normalizes the state ''state_''
 inline void normalize(vector<double> & state_){ 
   double norm=0;
   for (int i=0; i<n_; i++){
      norm+=state_[i]*state_[i];
   }
   norm=sqrt(norm);
   if (norm<1e-4)
      norm=1;
   for (int i=0; i<n_; i++){
      state_[i]=state_[i]/norm;
   }
 }

 //returns overlap between state_ and ground state
 double calculate_overlap(vector<double> & state_){
    double fid=0;
    vector<double> initialstate_;
    get_TCGS(initialstate_);
    for (int i=0; i<n_; i++){
      fid+=initialstate_[i]*state_[i];
    }
 return fid;
 }

 //returns overlap between state_ and state1
 double overlap(vector<double> & state_,vector<double> & state1){
    double fid=0;
    for (int i=0; i<n_; i++){
      fid+=state_[i]*state1[i];
    }
 return fid;
 }

 //returns toric ground state
 void get_TCGS(vector<double> & state_){
   state_.resize(n_,0);
   for (int N=0; N<int(pow2[L*L]);N++){
     int n=N;
     vector<int> spinconf(2*L*L,1);
     for (int i=0; n>0; i++){
       int a=n%2;
       if(a==1)
         flip_vertex(spinconf,i);
       n=n/2;
     }
     double weight=1;
     unsigned long long index=0;
     for (int s_i=0; s_i<2*L*L; s_i++){
       if (spinconf[2*L*L-1-s_i]==(-1)){
          index+=pow2[s_i];
       }
     } 
     state_[index]=weight;
   
   }
 normalize(state_);
 }

 //initializes state in exactly solvable ground state form, only z-fields
 void Initializestate(vector<double> & state_){
   cout<<"starting initializing: "<<endl;
   state_.resize(n_,0);
   for (int N=0; N<int(pow2[L*L]);N++){
     int n=N;
     vector<int> spinconf(2*L*L,1);
     for (int i=0; n>0; i++){
       int a=n%2;
       if(a==1)
         flip_vertex(spinconf,i);
       n=n/2;
     }
     double weight=0;
     for (int spin_i=0; spin_i<2*L*L; spin_i++){
      weight+=0.5*lambdaA_[spin_i]*spinconf[spin_i];  
    }
    weight=exp(weight);
    
    unsigned long long index=0;
    for (int s_i=0; s_i<2*L*L; s_i++){
      if (spinconf[2*L*L-1-s_i]==(-1)){
        index+=pow2[s_i];
        }
    } 
   state_[index]=weight;
   
 }
 normalize(state_);
 }

 //initializes state_ with all spins up
 vector<double> Initializestate_allspinsup(vector<double> & state_){
   cout<<"Starting initializing all spins up..."<<endl;
   state_.resize(n_,0);
   state_[0]=1;
   cout<<"Done"<<endl;
 return state_;  
 }

 //loads field condigurations form files filename_A and filename_B
 void LoadParameters(string filename_A, string filename_B){
   cout<<"Starting loading Parameters..."<<endl;
   lambdaA_.resize(2*L*L);
   ifstream fin(filename_A.c_str());

   if(!fin.good()){
      cerr<<"# Error : Cannot load from file "<<filename_A<<" : file not found."<<endl;
      abort();
    }
   for (int l_i=0; l_i<2*L*L; l_i++){
      fin>>lambdaA_[l_i];
   }
  
   lambdaB_.resize(2*L*L);
   ifstream fin_B(filename_B.c_str());

   if(!fin_B.good()){
      cerr<<"# Error : Cannot load from file "<<filename_B<<" : file not found."<<endl;
      abort();
    }
   for (int l_i=0; l_i<2*L*L; l_i++){
      fin_B>>lambdaB_[l_i];
   }
  cout<<"Finished loading parameters"<<endl;
 }

 //returns index of kronecker product
 inline int get_index1D(int x,int y,int AB){
    int z;
    z=L*2*x+y*2+AB;
    return z;
    }

 //returns (in call-by reference arguments) 1D indices of spins around plaquette p
 inline void get_indices_plaquette(int p, int & l1D, int & o1D, int & r1D, int & u1D){
    int x=int(p)/int(L);
    int y=p-L*x;
    int yp1=y+1;
    int xm1=x-1;
    if (y==(L-1))
       	yp1=0;
    if (x==0)
       xm1=L-1;
    l1D=get_index1D(x,y,0);
    o1D=get_index1D(xm1,y,1);
    r1D=get_index1D(x,yp1,0);
    u1D=get_index1D(x,y,1);
 }

 //returns (in call-by reference arguments) 1D indices of spins around vertex s
 inline void get_indices_vertex(int s, int & l1D, int & o1D, int & r1D, int & u1D){
    int x=int(s)/int(L);
    int y=s-L*x;
    int ym1=y-1;
    int xp1=x+1;
    if (y==0)
       ym1=L-1;
    if (x==(L-1))
       xp1=0;
    l1D=get_index1D(x,ym1,1);
    o1D=get_index1D(x,y,0);
    r1D=get_index1D(x,y,1);
    u1D=get_index1D(xp1,y,0);
 }

 //flips spins in configuration spinlattice (not state, but vector with -1 for spin down, 1 for spin up) around vertex s
 inline void flip_vertex(vector<int> & spinlattice, int s){
    int l1D, o1D, r1D, u1D;
    get_indices_vertex(s,l1D,o1D,r1D,u1D);
    spinlattice[l1D]=-spinlattice[l1D];
    spinlattice[o1D]=-spinlattice[o1D];
    spinlattice[r1D]=-spinlattice[r1D];
    spinlattice[u1D]=-spinlattice[u1D];
 }
 unsigned long long int get_n(){
    return n_;}
};

