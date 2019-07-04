import tensorflow as tf

#learning rate used for training: here not important
LR=0.00001

#----------------------------------------------------------------------------------------------
#network

#stabilizer expectation value input
tf_x = tf.placeholder(tf.float32, [None, 3])    # value in the range of (0, 1)

#label: field strength
tf_y_beta=tf.placeholder(tf.float32, [None,1])

#hidden layers
en0 = tf.layers.dense(tf_x, 128, tf.nn.relu,name='en0_')  
en10 = tf.layers.dense(en0, 150, tf.nn.relu)
en1 = tf.layers.dense(en10, 128, tf.nn.relu)      

#output neuron
beta_output=tf.layers.dense(en1, 1)

tf_y_beta2=tf.cast(tf_y_beta,tf.float64)

#calculate loss of continuous parameter as mean squared error
loss=tf.losses.mean_squared_error(tf_y_beta2,beta_output)

predictions=beta_output

train = tf.train.AdamOptimizer(LR).minimize(loss)

#-----------------------------------------------------------------------------------------------
