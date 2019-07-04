"""
Author: Agnes Valenti
This file calculates the Hamiltonian error as explained in appendix  . Input are the field configurations 'rbl_x.txt' and 'rbl_z.txt' and the initial configurations 'rbl_x0.txt' and 'rbl_z0.txt'.
"""
import numpy as np

#numerically found best polynomial fit for the error rate, given the vertex (plaquette) flip probability p
def error_rate(p):
  return (0.2187*p+0.72419*p**2-2.5398*p**3+4.90118*p**4)

#calculate the single qubit phase flip probability
for i in xrange(1):
  expv=np.loadtxt("expv_x.txt")
  expv2=np.sum(expv, axis=0)/np.shape(expv)[0]
  vertexprob= (1-expv2[0])/2.0
  a= error_rate(vertexprob)
  print "Probability of a single qubit phase flip: ", a
  
#calculate the single qubit bit flip probability
for i in xrange(1):
  expv=np.loadtxt("expv_z.txt".format(i))
  expv2=np.sum(expv, axis=0)/np.shape(expv)[0]
  vertexprob= (1-expv2[0])/2.0
  a= error_rate(vertexprob)
  print "Probability of a single qubit bit flip: ", a





  
