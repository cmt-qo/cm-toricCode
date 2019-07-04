//Author: Agnes Valenti
using namespace std;

//expvalues_xloops: calculates expectation values of A_s, A_s' and A_s A_s'.
//IN: reads in field configuration 'rbl_x.txt'
//OUT: outputs expectation values 'expv_x.txt'
void expvalues_z(){

//lattice length
const int L=4;

//number of expectation values to be calculated per spin
const int N_loops=3;

//number of spins
const int num_categories=2*L*L;
const int num_ex0=1;
double nspins_= 2*L*L;



//parameters for Monte Carlo sampling
unsigned int sweeps=100000;
int messung=20;
int thermalize=0.1*sweeps;
         
  
//if given, read in num_ex (number of examples)
int num_ex=num_ex0;
//field configuration (beta_z-fields)
vector<vector<double> > Logits_(num_ex,vector<double>(num_categories));
//read in field configuration
read_rbl_z(Logits_,num_ex,nspins_);

//open file to write expectation values
ofstream file_a;
ostringstream fileNameStream_a("");
fileNameStream_a<<"expv_z.txt";
string fileName_a=fileNameStream_a.str();
file_a.open(fileName_a.c_str());


for (int ne_i=0; ne_i<num_ex; ne_i++){
   //new variable to store field configuration in
   vector<double> spinweights(2*L*L,0);
   //only useful when num_ex>1
   vector<double> logits_(num_categories);
   for (int i=0; i<num_categories; i++){
       logits_[i]=Logits_[ne_i][i];}

   for(int w_i=0; w_i<2*L*L; w_i++){
      spinweights[w_i]=logits_[w_i];
      }
  
   //store expectation values of stabilizers in 'expvaluesloops'
   vector<vector<double>> expvaluesloops(2*L*L,vector<double>(N_loops));
   

   #pragma omp parallel for num_threads(4)
   for (int loop_i=0; loop_i<N_loops*2*L*L; loop_i++){   
      //loop over all expectation values to calculate: number of expectation values per spin times number of spins
      //notice that we calculate some expectation values twice, we can fasten the algorithm by 1/3, if we re-use already calculated expectation values

      //make loop------- weights[i]=1 if spin i is in loop that is measured, 0 else. Loop contains indices s of vertices, for which A_s is part of the X-loop
      vector<int> loop;
      //initialize weights 
      vector<int> weights(2*L*L);
      for (int i=0; i<2*L*L; i++){
          weights[i]=1;
          }

      //loop: which As in loop (make loop) (X-loop=product of vertex operators)
      //spin, for which expectation values are measured
      int spin_atlattice=loop_i/3;

      //adjacent vertices to spin_atlattice
      int vertex1,vertex2;
      get_2vertices(spin_atlattice,vertex1,vertex2,L);
      
      if (loop_i%3==0)
        loop=make_As(vertex1);
      if (loop_i%3==1)
        loop=make_As(vertex2);
      if (loop_i%3==2)
        loop=make_As12(vertex1,vertex2);
      
      
      //use the vector 'loop' to obtain 'weights': Flip the spins of a vertex s, if s is an element of 'loop'
      for(int i=0; i<loop.size(); i++){
         flip_vertex(weights,loop[i],L);
         }
      
      //change weights from a vector containing 1 and -1 to a vector containing 1 and 0
      for (int i=0; i<2*L*L; i++){
         if (weights[i]==1){
            weights[i]=0;}
         else if (weights[i]==-1){
            weights[i]=1;}
         }

      
      //------------------------------------------------------------------------------------------------------------
      //Monte Carlo: sampling

      //random number generator (inside parallelized loop, such that there is one generator for each thread)
      mt19937 gen_;
      Seed(gen_);   
     
      //field strength
      double beta=1;

      //helper variable
      int index=0;
      
      //bool value determines, if expectation values of X-loops of Z-loops is calculated. True for X-loops
      MonteCarlo sample(L,true);
      //run sampling process via Metropolis-Hastings algorithm
      vector<double> werte=sample.MC(beta,sweeps,messung,thermalize, weights, index, spinweights, gen_);
      int anzahl=index;
      vector<double> wert1(anzahl);
      vector<double> wert2(anzahl);
      for (int i_werte=0; i_werte<anzahl; i_werte++){
          wert1.at(i_werte)=werte.at(i_werte);
          wert2.at(i_werte)=werte.at(i_werte+anzahl);
          }
      double wert1average_;
      double wert1error_;

      //calculate mean and error from Markov chain samples
      OutputValue(wert1,wert1average_,wert1error_);
      expvaluesloops[spin_atlattice][loop_i%3]=wert1average_;
         
      //MC sampling up to here  
      //----------------------------------------------------------------------------------
      }//closes Nloops loop


   //write expectation values into file
   for(int s_i=0; s_i<2*L*L; s_i++){
     for(int l_i=0; l_i<N_loops; l_i++){
          file_a<<expvaluesloops[s_i][l_i]<<" ";}
     file_a<<endl;}
   file_a<<endl;
   
   }  //closes ne_i loop here

file_a.close();

 }

