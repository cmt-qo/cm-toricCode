//Author: Agnes Valenti
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


void Seed(mt19937 & gen){
    random_device rd;
    //gen_.seed(std::time(nullptr));
    gen.seed(rd());
    //choose seed: deterministic (seed with deterministic number) or random (processor time or random device), here random
    }

//returns index of kronecker product
inline int get_index1D(int x,int y,int AB, int L_){
    int z;
    z=L_*2*x+y*2+AB;
    return z;
    }

//returns (as call by reference variables s1, s2) the indices of the two vertices adjacent to spin
inline void get_2vertices(int spin, int & s1, int & s2, int L_){
   int spin_x,spin_y,spin_AB;
   spin_x=spin/(2*L_);
   spin_y=(spin-spin_x*2*L_)/2;
   spin_AB=(spin-spin_x*2*L_-spin_y*2);
   int xs1,xs2,ys1,ys2;
   if (spin_AB==0){
      ys1=spin_y;
      ys2=spin_y;
      if (spin_x==0)
        xs1=L_-1;
      else
        xs1=(spin_x-1);
      xs2=spin_x;
      }
   if (spin_AB==1){
      xs1=spin_x;
      xs2=spin_x;
      ys1=spin_y;
      ys2=(spin_y+1)%L_;
      }
   s1=L_*xs1+ys1;
   s2=L_*xs2+ys2;
   }

//returns (as call by reference variables p1, p2) the indices of the two plaquettes adjacent to spin
inline void get_2plaquettes(int spin, int & p1, int & p2, int L_){
   int spin_x,spin_y,spin_AB;
   spin_x=spin/(2*L_);
   spin_y=(spin-spin_x*2*L_)/2;
   spin_AB=(spin-spin_x*2*L_-spin_y*2);
   int xp1,xp2,yp1,yp2;
   if (spin_AB==1){
      yp1=spin_y;
      yp2=spin_y;
      xp1=spin_x;
      xp2=(spin_x+1)%L_;
      }
   if (spin_AB==0){
      xp1=spin_x;
      xp2=spin_x;
      yp1=spin_y;
      if (spin_y==0)
        yp1=L_-1;
      else
        yp1=(spin_y-1);
      yp2=spin_y%L_;
      }
   p1=L_*xp1+yp1;
   p2=L_*xp2+yp2;
   }


//calculates mean, statistical error and autocorrelation time for vector "value_"
//IN: value_: vector of f(x), where x are samples obtained from a Markov Chain, f the function whose integral we aim to estimate via Monte Carlo sampling
//wert_: variable to store mean
//werterror_: variable to store statistical error
void OutputValue(vector<double> value_,double & wert_, double & werterror_){
    int nbins=100;

    int blocksize=floor(double(value_.size())/double(nbins));
    //cout<<"sizeenergy"<<value_.size()<<endl;
    double value_average=0;
    double value_averagesq=0;

    double value_average_unblocked=0;
    double value_averagesq_unblocked=0;

    for(int i=0;i<nbins;i++){
      double eblock=0;
      for(int j=i*blocksize;j<(i+1)*blocksize;j++){
        eblock+=value_[j];
        assert(j<value_.size());

        double delta=value_[j]-value_average_unblocked;
        value_average_unblocked+=delta/double(j+1);
        double delta2=value_[j]-value_average_unblocked;
        value_averagesq_unblocked+=delta*delta2;
      }
      eblock/=double(blocksize);
      double delta=eblock-value_average;
      value_average+=delta/double(i+1);
      double delta2=eblock-value_average;
      value_averagesq+=delta*delta2;
    }

    value_averagesq/=(double(nbins-1));
    value_averagesq_unblocked/=(double((nbins*blocksize-1)));

    double estav=value_average; ///double(nspins_);
    wert_=estav;
    double esterror=sqrt(value_averagesq/double(nbins)); ///double(nspins_);
    werterror_=esterror;

  //cout<<"# Estimated autocorrelation time is ";
  //cout<<0.5*double(blocksize)*value_averagesq/value_averagesq_unblocked<<endl;
  }


//functions to define X-loops (to calculate expectation values of)

//IN: s1, s2: indices of two (for our case, adjacent) vertices
//OUT: loop: vector of length 2, contains s1 and s2
vector<int> make_As12(int s1, int s2){
    vector<int> loop(2);
    loop[0]=s1;
    loop[1]=s2;
    return loop;
    }

//IN: s: index of vertex
//OUT: loop: vector of length 1, contains s
vector<int> make_As(int s){
    vector<int> loop(1);
    loop[0]=s;
    return loop;
    }

//functions to define Z-loops (to calculate expectation values of)

//IN: p1, p2: indices of two (for our case, adjacent) vertices
//OUT: loop: vector of length 2, contains p1 and p2
vector<int> make_Bp12(int p1, int p2){
    vector<int> loop(2);
    loop[0]=p1;
    loop[1]=p2;
    return loop;
    }

//IN: p: index of vertex
//OUT: loop: vector of length 1, contains p
vector<int> make_Bp(int p){
    vector<int> loop(1);
    loop[0]=p;
    return loop;
    }


//flips spins around vertex s
//IN: spinconfiguration: vector of length 2*L*L, each entry can take the values +1 (spin up) or -1 (spin down)
//s: index of vertex
void flip_vertex(vector<int> & spinlattice, int s, int L_){
    int x=int(s)/int(L_);
    int y=s-L_*x;
    int ym1,xp1;
    ym1=y-1;
    xp1=x+1;
        
    if (y==0){
       ym1=L_-1;}
    if (x==(L_-1)){
       xp1=0;}
    int l1D, o1D, r1D, u1D;
    l1D=get_index1D(x,ym1,1,L_);
    o1D=get_index1D(x,y,0,L_);
    r1D=get_index1D(x,y,1,L_);
    u1D=get_index1D(xp1,y,0,L_);
    spinlattice[l1D]=-spinlattice[l1D];
    spinlattice[o1D]=-spinlattice[o1D];
    spinlattice[r1D]=-spinlattice[r1D];
    spinlattice[u1D]=-spinlattice[u1D];
    }

//flips spins around plaquette p
//IN: spinconfiguration: vector of length 2*L*L, each entry can take the values +1 (spin up) or -1 (spin down)
//p: index of plaquette
void flip_plaquette(vector<int> & spinlattice, int p, int L_){    
    int x=int(p)/int(L_);
    int y=p-L_*x;
    int yp1=y+1;
    int xm1=x-1;
    if (y==(L_-1))
       	yp1=0;
    if (x==0)
       xm1=L_-1;
    int l1D=get_index1D(x,y,0,L_);
    int o1D=get_index1D(xm1,y,1,L_);
    int r1D=get_index1D(x,yp1,0,L_);
    int u1D=get_index1D(x,y,1,L_);
    spinlattice[l1D]=-spinlattice[l1D];
    spinlattice[o1D]=-spinlattice[o1D];
    spinlattice[r1D]=-spinlattice[r1D];
    spinlattice[u1D]=-spinlattice[u1D];
    }


//reads field configuration from file 'rbl_x.txt' into variable passed by reference (logits)
//IN: logits: vector for field configuration, of size (2*L^2, num_ex)
//num_ex: number of examples (number of states), here 1
void read_rbl_x(vector<vector<double> > & logits_, int num_ex, int num_categories){
    vector<double> labels_(num_ex);
    string filename="rbl_x.txt";
    ifstream fin(filename.c_str());
    if(!fin.good()){
       cerr<<"# Error : Cannot load from file "<<filename<<" : file not found."<<endl;
       abort();
       }
    for (int fin_j=0; fin_j<num_categories; fin_j++){
       for (int fin_i=0; fin_i<num_ex; fin_i++){
          fin>>logits_[fin_i][fin_j];}}
    }

//reads field configuration from file 'rbl_z.txt' into variable passed by reference (logits)
//IN: logits: vector for field configuration, of size (2*L^2, num_ex)
//num_ex: number of examples (number of states), here 1
void read_rbl_z(vector<vector<double> > & logits_, int num_ex, int num_categories){
    vector<double> labels_(num_ex);
    string filename="rbl_z.txt";
    ifstream fin(filename.c_str());
    if(!fin.good()){
       cerr<<"# Error : Cannot load from file "<<filename<<" : file not found."<<endl;
       abort();
       }
    for (int fin_j=0; fin_j<num_categories; fin_j++){
       for (int fin_i=0; fin_i<num_ex; fin_i++){
          fin>>logits_[fin_i][fin_j];}}
    }





