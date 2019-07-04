//Author: Agnes Valenti
#include "MC/sampler.cc"
#include "MC/expvalues_z.cc"
#include "MC/expvalues_x.cc"
#include "MC/sign_z.cc"
#include "MC/sign_x.cc"
using namespace std;


int main(){

//calculate expectation values of Xloops (A_s, A_s', A_s A_s'),
//input: 'rbl_z.txt', output: 'expv_z.txt'
expvalues_z();

//calculate expectation values of Zloops (B_p, B_p', B_p B_p')
//input: 'rbl_x.txt', output: 'expv_x.txt'
expvalues_x();

//calculate expectation values of sigma_z,
//input: 'rbl_z.txt', output: 'sign_z.txt'
sign_x();

//calculate expectation values of sigma_x,
//input: 'rbl_x.txt', output: 'sign_x.txt'
sign_z();


 return 0;}

