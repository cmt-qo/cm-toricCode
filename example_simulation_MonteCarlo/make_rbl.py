"""
This file outputs random field configurations 'rbl_x.txt' and 'rbl_z.txt' as well as copies 'rbl_x0.txt' and 'rbl_z0.txt' meeting the restriction of solubility, i.e. on every spin only a field in one direction can be present.
"""

import numpy as np
import random as rd


L=4
num_ex0=1
labels=np.arange(num_ex0)
run=np.zeros(num_ex0)
beta=np.random.rand(num_ex0)*0.4+0.1
logits=np.random.rand(num_ex0,2*L*L)+1.0

#for i in xrange(num_ex0):
#   logits[i,:]=softmax(logits[i,:])

A=np.random.randint(0,2,size=(num_ex0,2*L*L))*2-1

logits=logits*A

f1=open('rbl_x.txt', "w")
f2=open('rbl_z.txt', "w")
f3=open('rbl_x0.txt', "w")
f4=open('rbl_z0.txt', "w")
for j in range(2*L*L):
   for i in range(num_ex0):
      r=np.random.rand()
      #print r
      if r>0.5:
        f1.write(str(logits[i,j])+" ")
        f2.write(str(0.0)+" ")
        f3.write(str(logits[i,j])+" ")
        f4.write(str(0.0)+" ")
      else:
        f1.write(str(0.0)+" ")
        f2.write(str(logits[i,j])+" ")
        f3.write(str(0.0)+" ")
        f4.write(str(logits[i,j])+" ")
   f1.write("\n")
   f2.write("\n")
   f3.write("\n")
   f4.write("\n")
f1.close()
f2.close()
f3.close()
f4.close()
