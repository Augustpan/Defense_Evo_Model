mkdir output
cd output

../od_worker od_output1.csv a 0 5&
../od_worker od_output2.csv a 5 5&
../od_worker od_output3.csv a 10 5&
../od_worker od_output4.csv a 15 5&

../mc_worker mc_output1.csv a 0 5&
../mc_worker mc_output2.csv a 5 5&
../mc_worker mc_output3.csv a 10 5&
../mc_worker mc_output4.csv a 15 5&

cd..