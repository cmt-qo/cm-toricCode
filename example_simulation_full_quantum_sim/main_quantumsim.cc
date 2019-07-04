/*
Author: Agnes Valenti
Input: field configurations 'rbl_x.txt' and 'rbl_z.txt'
Output: expectation values of stabilizers 'expv_x.txt' and 'expv_z.txt' as well as of sigma_z and sigma_x: 'sign_x.txt' and 'sign_z.txt'
The expectation values are calculated via full quantum simulation. The ground state is found numerically.
*/

#include "Hamiltonian.cc"
#include "get_vp.cc"
using namespace std;

//lattice length
const int L_=3;

//number of iterations
int iterations=6;



int main(){

//create object Hamiltonian H with loaded field configurations: beta_z in the file rbl_z.txt and beta_x in the file rbl_x.txt
Hamiltonian H("rbl_z.txt","rbl_x.txt", L_);

//state, state2: full quantum states. The states are ordered as follows:
//The basis vector with 1 at index i and 0 otherwise corresponds to 
//the spinconfiguration equal to the binary representation of i (in the binary representation: 1 spin down, 0 spin up)
vector<double> state;
vector<double> state2;


//The state is initialized to be approximately equal to the ground state, such that faster convergence is obtained when finding the ground state. In particular: As it is intialized, it would correspond to the ground state if on every spin only one type of field would be present
H.get_TCGS(state);
H.apply_expA(state);
H.apply_expB(state);


//find ground state by applying the Hamiltonian often to the state. The Hamiltonian is modified, such that the ground state has the largest energy
cout<<"Finding the ground state numerically..."<<endl;
double h=0.0;
for (int i=0; i<1000; i++){   
    H.apply_H_noisy(state,h);
    H.normalize(state);
    }
cout<<"Done"<<endl;


//---------------------------------------
//calculate expectation values

cout<<"Calculate expectation values..."<<endl;
vector<double> statecopy=state;

ofstream file_a;
ostringstream fileNameStream_a("");
fileNameStream_a<<"expv_z.txt";
string fileName_a=fileNameStream_a.str();
file_a.open(fileName_a.c_str());

for (int spin_atlattice=0; spin_atlattice<2*L_*L_; spin_atlattice++){
  //calculate exp.values As, As+1, AsAs+1, write in file expv_x.txt
  int vertex1,vertex2;
  get_2vertices(spin_atlattice,vertex1,vertex2,L_);
  cout<<"vertex 1: "<<vertex1<<"vertex 2: "<<vertex2<<endl;
  vector<double> state1=H.A_s(state,vertex2);
  double expv2=H.overlap(state,state1);
  state1=H.A_s(state,vertex1);
  double expv1=H.overlap(state,state1);
  state1=H.A_s(state1,vertex2);
  double expv12=H.overlap(state,state1);
  file_a<<expv1<<" " <<expv2<<" "<<expv12<<" ";
  file_a<<endl;
  
 }

file_a.close();


ofstream file_a2;
ostringstream fileNameStream_a2("");
fileNameStream_a2<<"sign_z.txt";
string fileName_a2=fileNameStream_a2.str();
file_a2.open(fileName_a2.c_str());

for (int spin_atlattice=0; spin_atlattice<2*L_*L_; spin_atlattice++){
  //calculate exp.values of sigma_z to determine sign of beta_z. Write in file sign_x.txt
  vector<double> state1=H.sigma_z(state,spin_atlattice);
  double expv=H.overlap(state,state1);
  file_a2<<expv;
  file_a2<<endl;
  
 }

file_a2.close();



ofstream file_b;
ostringstream fileNameStream_b("");
fileNameStream_b<<"expv_x.txt";
string fileName_b=fileNameStream_b.str();
file_b.open(fileName_b.c_str());

for (int spin_atlattice=0; spin_atlattice<2*L_*L_; spin_atlattice++){
  //calculate exp.values B_p, Bp+1, BpBp+1, write in file expv_z.txt
  int plaquette1,plaquette2;
  get_2plaquettes(spin_atlattice,plaquette1,plaquette2,L_);
  vector<double> state1=H.B_p(state,plaquette2);
  double expv2=H.overlap(state,state1);
  state1=H.B_p(state,plaquette1);
  double expv1=H.overlap(state,state1);
  state1=H.B_p(state1,plaquette2);
  double expv12=H.overlap(state,state1);
  file_b<<expv1<<" " <<expv2<<" "<<expv12<<" ";
  file_b<<endl;
 }

file_b.close();


ofstream file_b2;
ostringstream fileNameStream_b2("");
fileNameStream_b2<<"sign_x.txt";
string fileName_b2=fileNameStream_b2.str();
file_b2.open(fileName_b2.c_str());

for (int spin_atlattice=0; spin_atlattice<2*L_*L_; spin_atlattice++){
  //calculate exp.values of sigma_x to determine sign of beta_x. Write in file sign_z.txt (z because z-loops)
  vector<double> state1=H.sigma_x(state,spin_atlattice);
  double expv=H.overlap(state,state1);
  file_b2<<expv;
  file_b2<<endl;
  
 }

file_b2.close();
cout<<"Done"<<endl;


 return 0;
 }
