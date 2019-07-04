"""
Author: Agnes Valenti
This file calculates the Hamiltonian error as explained in appendix  . Input are the field configurations 'rbl_x.txt' and 'rbl_z.txt' and the initial configurations 'rbl_x0.txt' and 'rbl_z0.txt'.
"""
import numpy as np

L=4

#returns index of kronecker product
def get_index1D(x,y,AB):  
    #returns index of kronecker product
    z=L*2*x+y*2+AB
    return z

#returns indices of the four spins around the plaquette specified by (x,y)
def get_indices_plaquette(x,y):  #x:0..k-1, y:0..k-1
    yp1=y+1
    xm1=x-1
    if y==(L-1):
       	yp1=0
    if x==0:
       xm1=L-1
    l=[x,y,0]
    o=[xm1,y,1]
    r=[x,yp1,0]
    u=[x,y,1]
    l1D=get_index1D(l[0],l[1],l[2])
    o1D=get_index1D(o[0],o[1],o[2])
    r1D=get_index1D(r[0],r[1],r[2])
    u1D=get_index1D(u[0],u[1],u[2])
    return l1D,o1D,r1D,u1D

#returns indices of the four spins around the vertex (star) specified by (x,y)
def get_indices_star(x,y):  #x:0..k-1, y:0..k-1
    ym1=y-1
    xp1=x+1
    if y==0:
       ym1=L-1
    if x==(L-1):
       xp1=0
    l=[x,ym1,1]
    o=[x,y,0]
    r=[x,y,1]
    u=[xp1,y,0]
    l1D=get_index1D(l[0],l[1],l[2])
    o1D=get_index1D(o[0],o[1],o[2])
    r1D=get_index1D(r[0],r[1],r[2])
    u1D=get_index1D(u[0],u[1],u[2])
    return l1D,o1D,r1D,u1D


#useful functions to calculate normalization

def add_error1(i1,i2,i3,i4):
    return -np.sinh(i1)*np.cosh(i2)*np.cosh(i3)*np.cosh(i4)

def add_error2(i1,i2,i3,i4):
    return np.sinh(i1)*np.sinh(i2)*np.cosh(i3)*np.cosh(i4)
    
def add_error3(i1,i2,i3,i4):
    return -np.sinh(i1)*np.sinh(i2)*np.sinh(i3)*np.cosh(i4)
    
def add_error4(i1,i2,i3,i4):
    return np.sinh(i1)*np.sinh(i2)*np.sinh(i3)*np.sinh(i4)
   

#returns normalization for Hamiltonian error
def calculatenorm(lambdas, lambdasz):
   norm1=0
   normsinglespinsr=np.zeros(2*L*L)
   normsinglespinsrz=np.zeros(2*L*L)
   for x in range(L):
       for y in range(L):
           l,o,r,u=get_indices_star(x,y)
           norm1+=(add_error2(lambdas[l],lambdas[o],lambdas[r],lambdas[u]))**2
           norm1+=(add_error2(lambdas[l],lambdas[r],lambdas[o],lambdas[u]))**2
           norm1+=(add_error2(lambdas[l],lambdas[u],lambdas[r],lambdas[o]))**2
           norm1+=(add_error2(lambdas[o],lambdas[r],lambdas[l],lambdas[u]))**2
           norm1+=(add_error2(lambdas[o],lambdas[u],lambdas[l],lambdas[r]))**2
           norm1+=(add_error2(lambdas[r],lambdas[u],lambdas[l],lambdas[o]))**2
        
           norm1+=(add_error3(lambdas[o],lambdas[r],lambdas[u],lambdas[l]))**2
           norm1+=(add_error3(lambdas[l],lambdas[r],lambdas[u],lambdas[o]))**2
           norm1+=(add_error3(lambdas[l],lambdas[o],lambdas[u],lambdas[r]))**2
           norm1+=(add_error3(lambdas[l],lambdas[o],lambdas[r],lambdas[u]))**2
        
           norm1+=(add_error4(lambdas[l],lambdas[o],lambdas[r],lambdas[u]))**2
        
           
           normsinglespinsr[l]+=add_error1(lambdas[l],lambdas[o],lambdas[r],lambdas[u])
           normsinglespinsr[o]+=add_error1(lambdas[o],lambdas[l],lambdas[r],lambdas[u])
           normsinglespinsr[r]+=add_error1(lambdas[r],lambdas[o],lambdas[l],lambdas[u])
           normsinglespinsr[u]+=add_error1(lambdas[u],lambdas[l],lambdas[r],lambdas[o])
           
           l,o,r,u=get_indices_plaquette(x,y)
           norm1+=(add_error2(lambdasz[l],lambdasz[o],lambdasz[r],lambdasz[u]))**2
           norm1+=(add_error2(lambdasz[l],lambdasz[r],lambdasz[o],lambdasz[u]))**2
           norm1+=(add_error2(lambdasz[l],lambdasz[u],lambdasz[r],lambdasz[o]))**2
           norm1+=(add_error2(lambdasz[o],lambdasz[r],lambdasz[l],lambdasz[u]))**2
           norm1+=(add_error2(lambdasz[o],lambdasz[u],lambdasz[l],lambdasz[r]))**2
           norm1+=(add_error2(lambdasz[r],lambdasz[u],lambdasz[l],lambdasz[o]))**2
       
           norm1+=(add_error3(lambdasz[o],lambdasz[r],lambdasz[u],lambdasz[l]))**2
           norm1+=(add_error3(lambdasz[l],lambdasz[r],lambdasz[u],lambdasz[o]))**2
           norm1+=(add_error3(lambdasz[l],lambdasz[o],lambdasz[u],lambdasz[r]))**2
           norm1+=(add_error3(lambdasz[l],lambdasz[o],lambdasz[r],lambdasz[u]))**2
        
           norm1+=(add_error4(lambdasz[l],lambdasz[o],lambdasz[r],lambdasz[u])-add_error4(lambdalz[l],lambdalz[o],lambdalz[r],lambdalz[u]))**2
        
           
           normsinglespinsrz[l]+=add_error1(lambdasz[l],lambdasz[o],lambdasz[r],lambdasz[u])
           normsinglespinsrz[o]+=add_error1(lambdasz[o],lambdasz[l],lambdasz[r],lambdasz[u])
           normsinglespinsrz[r]+=add_error1(lambdasz[r],lambdasz[o],lambdasz[l],lambdasz[u])
           normsinglespinsrz[u]+=add_error1(lambdasz[u],lambdasz[l],lambdasz[r],lambdasz[o])
   norm1+=np.dot((normsinglespinsr).T,  (normsinglespinsr)) 
   norm1+=np.dot((normsinglespinsrz).T,  (normsinglespinsrz))       
   return np.sqrt(norm1)
   




#initial field configuration in z-direction
lambdar=np.loadtxt("rbl_x0.txt")

#initial field configuration in x-direction
lambdarz=np.loadtxt("rbl_z0.txt")


for i in xrange(1):
   lambdal=lambdar-np.loadtxt("rbl_x.txt")
   lambdalz=lambdarz-np.loadtxt("rbl_z.txt")
   normr=calculatenorm(lambdar,lambdarz)
   if normr<1e-5:
     normr=1
   norml=calculatenorm(lambdal,lambdalz)
   if norml<1e-5:
     norml=1
   error=0
   
   #useful quantities for short-term storage 
   singlespinsr=np.zeros(2*L*L)
   singlespinsl=np.zeros(2*L*L)
   singlespinsrz=np.zeros(2*L*L)
   singlespinslz=np.zeros(2*L*L)

   #calculate Hamiltonian error
   for x in range(L):
       for y in range(L):
           l,o,r,u=get_indices_star(x,y)
           #add up terms resulting from the expansion of the exponential
           error+=(add_error2(lambdar[l],lambdar[o],lambdar[r],lambdar[u])/normr-add_error2(lambdal[l],lambdal[o],lambdal[r],lambdal[u])/norml)**2
           error+=(add_error2(lambdar[l],lambdar[r],lambdar[o],lambdar[u])/normr-add_error2(lambdal[l],lambdal[r],lambdal[o],lambdal[u])/norml)**2
           error+=(add_error2(lambdar[l],lambdar[u],lambdar[r],lambdar[o])/normr-add_error2(lambdal[l],lambdal[u],lambdal[r],lambdal[o])/norml)**2
           error+=(add_error2(lambdar[o],lambdar[r],lambdar[l],lambdar[u])/normr-add_error2(lambdal[o],lambdal[r],lambdal[l],lambdal[u])/norml)**2
           error+=(add_error2(lambdar[o],lambdar[u],lambdar[l],lambdar[r])/normr-add_error2(lambdal[o],lambdal[u],lambdal[l],lambdal[r])/norml)**2
           error+=(add_error2(lambdar[r],lambdar[u],lambdar[l],lambdar[o])/normr-add_error2(lambdal[r],lambdal[u],lambdal[l],lambdal[o])/norml)**2
        
           error+=(add_error3(lambdar[o],lambdar[r],lambdar[u],lambdar[l])/normr-add_error3(lambdal[o],lambdal[r],lambdal[u],lambdal[l])/norml)**2
           error+=(add_error3(lambdar[l],lambdar[r],lambdar[u],lambdar[o])/normr-add_error3(lambdal[l],lambdal[r],lambdal[u],lambdal[o])/norml)**2
           error+=(add_error3(lambdar[l],lambdar[o],lambdar[u],lambdar[r])/normr-add_error3(lambdal[l],lambdal[o],lambdal[u],lambdal[r])/norml)**2
           error+=(add_error3(lambdar[l],lambdar[o],lambdar[r],lambdar[u])/normr-add_error3(lambdal[l],lambdal[o],lambdal[r],lambdal[u])/norml)**2
        
           error+=(add_error4(lambdar[l],lambdar[o],lambdar[r],lambdar[u])/normr-add_error4(lambdal[l],lambdal[o],lambdal[r],lambdal[u])/norml)**2
        
           
           singlespinsr[l]+=add_error1(lambdar[l],lambdar[o],lambdar[r],lambdar[u])
           singlespinsl[l]+=add_error1(lambdal[l],lambdal[o],lambdal[r],lambdal[u])
           singlespinsr[o]+=add_error1(lambdar[o],lambdar[l],lambdar[r],lambdar[u])
           singlespinsl[o]+=add_error1(lambdal[o],lambdal[l],lambdal[r],lambdal[u])
           singlespinsr[r]+=add_error1(lambdar[r],lambdar[o],lambdar[l],lambdar[u])
           singlespinsl[r]+=add_error1(lambdal[r],lambdal[o],lambdal[l],lambdal[u])
           singlespinsr[u]+=add_error1(lambdar[u],lambdar[l],lambdar[r],lambdar[o])
           singlespinsl[u]+=add_error1(lambdal[u],lambdal[l],lambdal[r],lambdal[o])
           
           l,o,r,u=get_indices_plaquette(x,y)
           error+=(add_error2(lambdarz[l],lambdarz[o],lambdarz[r],lambdarz[u])/normr-add_error2(lambdalz[l],lambdalz[o],lambdalz[r],lambdalz[u])/norml)**2
           error+=(add_error2(lambdarz[l],lambdarz[r],lambdarz[o],lambdarz[u])/normr-add_error2(lambdalz[l],lambdalz[r],lambdalz[o],lambdalz[u])/norml)**2
           error+=(add_error2(lambdarz[l],lambdarz[u],lambdarz[r],lambdarz[o])/normr-add_error2(lambdalz[l],lambdalz[u],lambdalz[r],lambdalz[o])/norml)**2
           error+=(add_error2(lambdarz[o],lambdarz[r],lambdarz[l],lambdarz[u])/normr-add_error2(lambdalz[o],lambdalz[r],lambdalz[l],lambdalz[u])/norml)**2
           error+=(add_error2(lambdarz[o],lambdarz[u],lambdarz[l],lambdarz[r])/normr-add_error2(lambdalz[o],lambdalz[u],lambdalz[l],lambdalz[r])/norml)**2
           error+=(add_error2(lambdarz[r],lambdarz[u],lambdarz[l],lambdarz[o])/normr-add_error2(lambdalz[r],lambdalz[u],lambdalz[l],lambdalz[o])/norml)**2
        
           error+=(add_error3(lambdarz[o],lambdarz[r],lambdarz[u],lambdarz[l])/normr-add_error3(lambdalz[o],lambdalz[r],lambdalz[u],lambdalz[l])/norml)**2
           error+=(add_error3(lambdarz[l],lambdarz[r],lambdarz[u],lambdarz[o])/normr-add_error3(lambdalz[l],lambdalz[r],lambdalz[u],lambdalz[o])/norml)**2
           error+=(add_error3(lambdarz[l],lambdarz[o],lambdarz[u],lambdarz[r])/normr-add_error3(lambdalz[l],lambdalz[o],lambdalz[u],lambdalz[r])/norml)**2
           error+=(add_error3(lambdarz[l],lambdarz[o],lambdarz[r],lambdarz[u])/normr-add_error3(lambdalz[l],lambdalz[o],lambdalz[r],lambdalz[u])/norml)**2
        
           error+=(add_error4(lambdarz[l],lambdarz[o],lambdarz[r],lambdarz[u])/normr-add_error4(lambdalz[l],lambdalz[o],lambdalz[r],lambdalz[u])/norml)**2
        
           
           singlespinsrz[l]+=add_error1(lambdarz[l],lambdarz[o],lambdarz[r],lambdarz[u])
           singlespinslz[l]+=add_error1(lambdalz[l],lambdalz[o],lambdalz[r],lambdalz[u])
           singlespinsrz[o]+=add_error1(lambdarz[o],lambdarz[l],lambdarz[r],lambdarz[u])
           singlespinslz[o]+=add_error1(lambdalz[o],lambdalz[l],lambdalz[r],lambdalz[u])
           singlespinsrz[r]+=add_error1(lambdarz[r],lambdarz[o],lambdarz[l],lambdarz[u])
           singlespinslz[r]+=add_error1(lambdalz[r],lambdalz[o],lambdalz[l],lambdalz[u])
           singlespinsrz[u]+=add_error1(lambdarz[u],lambdarz[l],lambdarz[r],lambdarz[o])
           singlespinslz[u]+=add_error1(lambdalz[u],lambdalz[l],lambdalz[r],lambdalz[o])
  
   error+=np.dot((singlespinsr/normr-singlespinsl/norml).T,  (singlespinsr/normr-singlespinsl/norml)) 
   error+=np.dot((singlespinsrz/normr-singlespinslz/norml).T,  (singlespinsrz/normr-singlespinslz/norml))
   

   print "Hamiltonian error: ",error
   

