//Author: Agnes Valenti
#include "useful_functions.cc"
using namespace std;


class MonteCarlo{
 int L_;
 uniform_real_distribution<> distu_;
 uniform_int_distribution<> distnx_;
 uniform_int_distribution<> distny_;
 bool vertex_;
public:
  MonteCarlo(int length, bool vertex):L_(length),distnx_(0,L_-1),distny_(0,L_-1),vertex_(vertex){}

 inline void get_indices_vertex(int x, int y, int & l1D, int & o1D, int & r1D, int & u1D){
  int ym1,xp1;
  ym1=y-1;
  xp1=x+1;
  if (y==0){
    ym1=L_-1;}
  if (x==(L_-1)){
    xp1=0;}
        
  l1D=get_index1D(x,ym1,1,L_);
  o1D=get_index1D(x,y,0,L_);
  r1D=get_index1D(x,y,1,L_);
  u1D=get_index1D(xp1,y,0,L_);
  }

 inline void get_indices_plaquette(int x, int y, int & l1D, int & o1D, int & r1D, int & u1D){
  int yp1=y+1;
  int xm1=x-1;
  if (y==(L_-1))
    yp1=0;
  if (x==0)
    xm1=L_-1;
  l1D=get_index1D(x,y,0,L_);
  o1D=get_index1D(xm1,y,1,L_);
  r1D=get_index1D(x,yp1,0,L_);
  u1D=get_index1D(x,y,1,L_);
  }

//MC: returns expectation value of X-loop specified by 'weights' using the Metropolis-Hastings algorithm 
//IN: 'weights: vector of length 2L^2, weights[i]=1 if i is included in the X-loop that is measured, 0 otherwise
//spinweights: it is measured on ground state with field configuration 'spinweights' 
//measured every messung*sweeps after thermalization time 'thermalize'
//beta: field strength
//index: helper index
//OUT: vector 'Werte', returning measured values along the Markov Chain to compute expectation value of X-loop
 vector<double> MC(double beta,unsigned int sweeps, int messung, int thermalize, vector<int> weights, int & index, vector<double> spinweights, mt19937 & gen_){ 
    //spin lattice, Markov chain flips spins on the lattice
    vector<int> spins(2*L_*L_); 
    
    //magnetization
    double M=0;
    for (int i=0; i<2*L_*L_; i++){   //initializing
         spins[i]=1;  //has to be initialized in cold state (eigenstate of all B_p, could also be a sum over contractible loops)
         M+=spins[i];
          }
           
    
    unsigned int maxit=sweeps*L_*L_;
    
    //values needed to calculate expectation value (to be returned)
    vector<double> wert1((sweeps)/messung-(thermalize)/messung);
    vector<double> wert2((sweeps)/messung-(thermalize)/messung);
    
    double Mloop1=0;
    for (int l1=0;l1<2*L_*L_;l1++){
        if (weights[l1]==1){
            Mloop1+=spins[l1];}
        }
    
    index=0; 
    
    for (int n=0; n<maxit; n++){
        
        //choose random vertex to be flipped
        int x;
        x=distnx_(gen_);//random vertex x-coordinate
        int y;
        y=distnx_(gen_);//random vertex y-Koordinate
        
        //calculate 'energy' difference when flipping the spins around the vertex specified by (x,y)
        double dE;
        int l1D, o1D, r1D, u1D;
        if (vertex_){
           get_indices_vertex(x,y,l1D,o1D,r1D,u1D);}
        else{
           get_indices_plaquette(x,y,l1D,o1D,r1D,u1D);}
        
        dE=2.0*beta*(spins[l1D]*spinweights[l1D]+spins[o1D]*spinweights[o1D]+spins[r1D]*spinweights[r1D]+spins[u1D]*spinweights[u1D]);  

        double rnumber;
        rnumber=distu_(gen_);
        //acceptance of vertex flip (Metropolis-Hastings test)
        if (rnumber<=exp(-dE)){
           spins[l1D]=(-1)*spins[l1D];
           spins[o1D]=(-1)*spins[o1D];
           spins[r1D]=(-1)*spins[r1D];
           spins[u1D]=(-1)*spins[u1D];
           M=M+2*(spins[l1D]+spins[o1D]+spins[r1D]+spins[u1D]);
           }
          
        //measure every L*L*messung steps
        if ((n%(L_*L_*messung)==0)&&(n>=L_*L_*thermalize)){
           double Mloop=0;
           for (int l1=0;l1<2*L_*L_;l1++){
               if (weights[l1]==1){
                  Mloop+=spins[l1]*spinweights[l1];}
               }
           wert1.at(index)=exp(-beta*Mloop);
           wert2.at(index)=1.0;
           index+=1;
           }
    }
    
    //write wert1 and wert2 into one vector
    vector<double> Werte(index+index);
    for (int i=0; i<index; i++){
        Werte.at(i)=wert1.at(i);
        Werte.at(i+index)=wert2.at(i);
    }
    
    return Werte;
    }


//MC: returns expectation value of sigma_z(position) via Metropolis-Hastings algorithm 
//IN: 'position': spin, on which sigma_z acts
//spinweights: it is measured on ground state with field configuration 'spinweights' 
//measured every messung*sweeps after thermalization time 'thermalize'
//beta: field strength
//index: helper index
//OUT: vector 'Werte', returning measured values along the Markov Chain to compute expectation value of sigma_z(position)
 vector<double> MC_sigma(double beta,unsigned int sweeps, int messung, int thermalize, int & index, vector<double> spinweights, mt19937 & gen_, int position){ 
    //spin lattice, Markov chain flips spins on the lattice
    vector<int> spins(2*L_*L_); 
    
    //magnetization
    double M=0;
    for (int i=0; i<2*L_*L_; i++){   //initializing
         spins[i]=1;  //has to be initialized in cold state (eigenstate of all B_p, could also be a sum over contractible loops)
         M+=spins[i];
          }
           
    
    unsigned int maxit=sweeps*L_*L_;
    
    //values needed to calculate expectation value (to be returned)
    vector<double> wert1((sweeps)/messung-(thermalize)/messung);
    vector<double> wert2((sweeps)/messung-(thermalize)/messung);
    
    index=0; 
    
    for (int n=0; n<maxit; n++){
        
        //choose random vertex to be flipped
        int x;
        x=distnx_(gen_);//random vertex x-coordinate
        int y;
        y=distnx_(gen_);//random vertex y-coordinate

        //calculate 'energy' difference when flipping the spins around the vertex specified by (x,y)
        double dE;
        int l1D, o1D, r1D, u1D;
        if (vertex_){
           get_indices_vertex(x,y,l1D,o1D,r1D,u1D);}
        else{
           get_indices_plaquette(x,y,l1D,o1D,r1D,u1D);}
        
        dE=2.0*beta*(spins[l1D]*spinweights[l1D]+spins[o1D]*spinweights[o1D]+spins[r1D]*spinweights[r1D]+spins[u1D]*spinweights[u1D]);  
        
        double rnumber;
        rnumber=distu_(gen_);
        //acceptance of vertex flip (Metropolis-Hastings test)
        if (rnumber<=exp(-dE)){
           spins[l1D]=(-1)*spins[l1D];
           spins[o1D]=(-1)*spins[o1D];
           spins[r1D]=(-1)*spins[r1D];
           spins[u1D]=(-1)*spins[u1D];
           M=M+2*(spins[l1D]+spins[o1D]+spins[r1D]+spins[u1D]);
           }
          
        //measure every L*L*messung steps
        if ((n%(L_*L_*messung)==0)&&(n>=L_*L_*thermalize)){
           wert1.at(index)=spins[position];
           wert2.at(index)=1.0;
           index+=1;
           }
    }
    
    //write wert1 and wert2 into one vector
    vector<double> Werte(index+index);
    for (int i=0; i<index; i++){
        Werte.at(i)=wert1.at(i);
        Werte.at(i+index)=wert2.at(i);
    }
    
    return Werte;
    }


};

