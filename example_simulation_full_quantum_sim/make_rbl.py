import numpy as np
import random as rd

def softmax(x):
    """Compute softmax values for each sets of scores in x."""
    e_x = np.exp(x - np.max(x))
    return e_x / e_x.sum()


L=3
num_ex0=1
labels=np.arange(num_ex0)
run=np.zeros(num_ex0)
beta=np.random.rand(num_ex0)*0.4+0.1
logits=np.random.rand(num_ex0,2*L*L)+0.5
#for i in xrange(num_ex0):
#   logits[i,:]=softmax(logits[i,:])

A=np.random.randint(0,2,size=(num_ex0,2*L*L))*2-1
print(A)

logits=logits*A
np.savetxt("rbl_x.txt",logits)
np.savetxt("rbl_x0.txt", logits)

logits=np.random.rand(num_ex0,2*L*L)+0.5
#for i in xrange(num_ex0):
#   logits[i,:]=softmax(logits[i,:])

A=np.random.randint(0,2,size=(num_ex0,2*L*L))*2-1
print(A)

logits=logits*A

np.savetxt("rbl_z.txt",logits)
np.savetxt("rbl_z0.txt",logits)
