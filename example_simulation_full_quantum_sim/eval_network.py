"""
Author: Agnes Valenti
This file evaluates the Artificial Neural Network trained to predict field strengths from expectation values.
Evaluation: input: 'expv_x.txt', 'expv_z.txt' (expectation values for x-fields and z-fields), 'rbl_x.txt', 'rbl_z.txt' (current field configurations).
Faulty measurements can be simulated by setting the parameter 'faulty' to 'True'. The standard deviation 
(parameter 'coeff') can be chosen, default value: 0.05.
Output: corrected field configurations 'rbl_x.txt' and 'rbl_z.txt'
"""
from __future__ import division
from __future__ import absolute_import
from __future__ import print_function


import sys
sys.path.append('../ANN/')

import tensorflow as tf
import matplotlib.pyplot as plt
import numpy as np
import network as net

#tf.set_random_seed(5)
tf.set_random_seed(6)
np.random.seed(3)

#lattice size
L=3 

#set faulty=True to simulate faulty measurements. The parameter 'coeff' corresponds to the standard deviation
faulty=False
#choose standard deviation to simulate faulty measurements. Default value: 0.05
coeff=0.05
 
saver=tf.train.Saver()
sess = tf.Session()
sess.run(tf.global_variables_initializer())

   
if True:
    #evaluate network
    
    #load expectation values measured on state, fields in z-direction
    state=np.loadtxt('expv_z.txt')

    #add Gaussian noise to measurement outcomes to simulate faulty measurements
    if faulty:
       r=0
       rn=0
       sign=0
       signn=0
       for i in range(np.shape(state)[0]):
         r1=np.random.normal(0,coeff)

         r12=np.random.normal(0,coeff)

         if i==0:
           r=np.random.normal(0,coeff)
           rn=r
         if i==np.shape(state)[0]-1:
           r1=rn
         state[i][0]=state[i][0]+r
         state[i][1]=state[i][1]+r1
         state[i][2]=state[i][2]+r12
         r=r1


    #load field strengths for correction
    betas=np.loadtxt('rbl_z.txt')

    #load trained network (change if lattice length is changed)
    tf.reset_default_graph()
    saver.restore(sess, '../ANN/slstates_dense_k3/-49900')
    
    #network predictions of absolute values of field strengths
    beta_output_ = sess.run([net.beta_output], feed_dict={net.tf_x: state})
    beta_output_=np.reshape(beta_output_,(np.shape(betas)[0],))
    
    beta_output_z=beta_output_

    #determine signs 
    sign_z=np.loadtxt("sign_z.txt")
    betas_z=np.loadtxt("rbl_z.txt")
    for i in range(2*L*L):
       if sign_z[i]>0:
          betas_z[i]=betas_z[i]-abs(beta_output_z[i])
       else:
          betas_z[i]=betas_z[i]+abs(beta_output_z[i])
    
    #save corrected configuration
    np.savetxt("rbl_z.txt", betas_z)



    #load expectation values measured on state, fields in x-direction
    state=np.loadtxt('expv_x.txt')

    #add Gaussian noise to measurement outcomes to simulate faulty measurements
    if faulty:
       r=0
       rn=0
       sign=0
       signn=0
       for i in range(np.shape(state)[0]):
         r1=np.random.normal(0,coeff)

         r12=np.random.normal(0,coeff)

         if i==0:
           r=np.random.normal(0,coeff)
           rn=r
         if i==np.shape(state)[0]-1:
           r1=rn
         state[i][0]=state[i][0]+r
         state[i][1]=state[i][1]+r1
         state[i][2]=state[i][2]+r12
         r=r1

    #load field strengths for correction
    betas=np.loadtxt('rbl_x.txt')

    #network predictions of absolute values of field strengths
    beta_output_ = sess.run([net.beta_output], feed_dict={net.tf_x: state})
    beta_output_=np.reshape(beta_output_,(np.shape(betas)[0],))
    
    #determine signs 
    beta_output_x=beta_output_
    sign_x=np.loadtxt("sign_x.txt")
    betas_x=np.loadtxt("rbl_x.txt")
    for i in range(2*L*L):
       if sign_x[i]>0:
         betas_x[i]=betas_x[i]-abs(beta_output_x[i])
       else:
         betas_x[i]=betas_x[i]+abs(beta_output_x[i])
    
    #save corrected configuration
    np.savetxt("rbl_x.txt", betas_x)
