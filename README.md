# cm-toricCode

To run a simulation, the user might choose between a correction process with the restriction of solubility
(contained in the folder 'example_simulation_MonteCarlo' or a full quantum simulation for the general
non-solvable model (contained in the folder 'example_simulation_full_quantum_sim').

The correction process with the restriction of solubility runs on a standard laptop in less than 5 minutes,
while full quantum simulation has to be run on a computation node with 72 threads and finishes in about
1.5 hours. We therefore suggest the solvable Monte Carlo simulation for testing purposes.

The user may choose one of the two simulations and proceed as follows:



--------------------------- MONTE CARLO SIMULATION ---------------------------------------------------

The complete description of the simulation including additional information and changes that can be 
made can be found in the file 'example_simulation_MonteCarlo/README_MonteCarlo'.

In the following, we only state the necessary steps to achieve minimal usage.


Change directory to example_simulation_MonteCarlo.

####################### COMPILING THE CODE ###########################################################

The provided codes ending with ".cc" use the standard c++ library and no additional packages. They can
be compiled using the command "make" in the folder 'example_simulation_MonteCarlo'.
The executable "run_MC" is produced.

A c++11 compliant compiler is assumed. Parallel programming using OpenMP is used.
These codes are used to calculate expectation values via Monte Carlo sampling.


The provided codes ending with ".py" use Python 2.7 and the TensorFlow library.
These codes are used to evaluate a pre-trained neural network.

######################################################################################################


########################## RUNNING A SIMULATION ######################################################
Minimal usage:
Once compiled, the necessary steps to run a simulation with 4 correction iterations are contained in 
the script "run.sh". 
To run, the user may give permission with "chmod +x run.sh" and execute the script with the command 
"./run.sh".

The simulation should run in less than 5 minutes on a standard laptop.


Outcome of the simulation are the following files:
'rbl_x0.txt', 'rbl_z0.txt': these are the initial field configurations
'rbl_x.txt', 'rbl_z.txt': these are the final field configurations after 4 iterations
The field configurations can be compared.
In addition, in every iteration step the probability of a single qubit error and the Hamiltonian error
are calculated and printed.


######################################################################################################

The complete description of the simulation including additional information and changes that can be 
made can be found in the file 'example_simulation_MonteCarlo/README_MonteCarlo'.

------------------------------------------------------------------------------------------------------




-------------------------- FULL QUANTUM SIMULATION ---------------------------------------------------

The complete description of the simulation including additional information and changes that can be 
made can be found in the file 'example_simulation_full_quantum_sim/README_quantum'.

In the following, we only state the necessary steps to achieve minimal usage.


Change directory to 'example_simulation_full_quantum_sim'.

####################### COMPILING THE CODE ###########################################################

The provided codes ending with ".cc" use the standard c++ library and no additional packages. They can
be compiled using the command "make" in the folder 'example_simulation_full_quantum_sim'.
The executable "run_quantumsim" is produced.

A c++11 compliant compiler is assumed. Parallel programming using OpenMP is used.
These codes are used to calculate expectation values via full quantum simulation.


The provided codes ending with ".py" use Python 2.7 and the TensorFlow library.
These codes are used to evaluate a pre-trained neural network.

######################################################################################################


########################## RUNNING A SIMULATION ######################################################

Minimal usage:
Once compiled, the necessary steps to run a simulation with 4 correction iterations are contained in 
the script "run.sh". 
To run, the user may give permission with "chmod +x run.sh" and execute the script with the command 
"./run.sh".
Alternatively, the user might run the steps detailed in the file 
'example_simulation_full_quantum_sim/README_quantum'.

IMPORTANT:
The simulation is parallelized and takes about 1.5 hours on a computation node with 72 threads and is 
not to be run on a laptop.

Outcome of the simulation are the following files:
'rbl_x0.txt', 'rbl_z0.txt': these are the initial field configurations
'rbl_x.txt', 'rbl_z.txt': these are the final field configurations after 4 iterations
The field configurations can be compared.
In addition, in every iteration step the probability of a single qubit error and the Hamiltonian error
are calculated and printed.

######################################################################################################


------------------------------------------------------------------------------------------------------


