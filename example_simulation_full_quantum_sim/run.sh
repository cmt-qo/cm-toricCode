printf "\n"
echo "Create random field configurations..."
python make_rbl.py
python make_restr_file.py
echo "Done"
printf "\n"
printf "\n"

echo "Start first iteration: calculate required expectation values..."
./run_quantumsim
echo "Done"
printf "\n"
echo "Calculate error measures with respect to the current field configuration:"
cp *.txt ../error_measures
cd ../error_measures
python Hamiltonianerror_L3.py  
python error_rate.py
rm *.txt
cd ../example_simulation_full_quantum_sim
printf "\n"
echo "Evaluate trained neural network and correct field configurations..."
python eval_network.py
echo "Done"
printf "\n"
printf "\n"

echo "Start second iteration: calculate required expectation values..."
./run_quantumsim
echo "Done"
printf "\n"
echo "Calculate error measures with respect to the current field configuration:"
cp *.txt ../error_measures
cd ../error_measures
python Hamiltonianerror_L3.py  
python error_rate.py
rm *.txt
cd ../example_simulation_full_quantum_sim
printf "\n"
echo "Evaluate trained neural network and correct field configurations..."
python eval_network.py
echo "Done"
printf "\n"
printf "\n"

echo "Start third iteration: calculate required expectation values..."
./run_quantumsim
echo "Done"
printf "\n"
echo "Calculate error measures with respect to the current field configuration:"
cp *.txt ../error_measures
cd ../error_measures
python Hamiltonianerror_L3.py  
python error_rate.py
rm *.txt
cd ../example_simulation_full_quantum_sim
printf "\n"
echo "Evaluate trained neural network and correct field configurations..."
python eval_network.py
echo "Done"
printf "\n"
printf "\n"

echo "Start fourth iteration: calculate required expectation values..."
./run_quantumsim
echo "Done"
printf "\n"
echo "Calculate error measures with respect to the current field configuration:"
cp *.txt ../error_measures
cd ../error_measures
python Hamiltonianerror_L3.py  
python error_rate.py
rm *.txt
cd ../example_simulation_full_quantum_sim
printf "\n"
echo "Evaluate trained neural network and correct field configurations..."
python eval_network.py
echo "Done"
printf "\n"
printf "\n"


rm expv_x.txt
rm expv_z.txt
rm sign_x.txt
rm sign_z.txt
