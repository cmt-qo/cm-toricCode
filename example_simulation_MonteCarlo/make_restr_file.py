"""
This code outputs a file 'correctx.txt' utilized to keep the model soluble during the correction process.
"""
import numpy as np

rblz=np.loadtxt("rbl_z.txt")
rblx=np.loadtxt("rbl_x.txt")

correctx=np.zeros_like(rblz)

for i in range(np.shape(rblz)[0]):
  if (abs(rblz[i])>0 and abs(rblx[i])>0):
     print("The given field configurations do not meet the restriction of solubility. On spin ", i, " fields in both x- and z-direction are present")
  if abs(rblz[i])>0:
     correctx[i]=1


np.savetxt("correctx.txt",correctx)
