//Author: Agnes Valenti
using namespace std;

//sign_x: calculates expectation values of sigma_z.
//IN: reads in field configuration 'rbl_x.txt'
//OUT: outputs expectation values 'sign_x.txt'
void sign_x(){

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
fileNameStream_a<<"sign_z.txt";
string fileName_a=fileNameStream_a.str();
file_a.open(fileName_a.c_str());

//new variable to store field configuration in
vector<double> spinweights(2*L*L,0);
//only useful when num_ex>1
vector<double> logits_(num_categories);
for (int i=0; i<num_categories; i++){
   logits_[i]=Logits_[0][i];}

for(int w_i=0; w_i<2*L*L; w_i++){
   spinweights[w_i]=logits_[w_i];
   }
   
//store expectation values of sigma_z for each spin in 'sigmas'
vector<double> sigmas(2*L*L);
   

#pragma omp parallel for num_threads(4)
for (int spin_atlattice=0; spin_atlattice<2*L*L; spin_atlattice++){   
   //loop over all expectation values to calculate = number of spins
    
   //------------------------------------------------------------------------------------------------------------
   //MC: sampling

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
   vector<double> werte=sample.MC_sigma(beta,sweeps,messung,thermalize, index, spinweights, gen_, spin_atlattice);
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
   sigmas[spin_atlattice]=wert1average_;
         
   //MC sampling up to here     
   //----------------------------------------------------------------------------------  
   }//closes spin_atlattice loop


//write expectation values into file  
for(int s_i=0; s_i<2*L*L; s_i++){
   file_a<<sigmas[s_i]<<" ";
   file_a<<endl;}
file_a<<endl;

file_a.close();

 }

