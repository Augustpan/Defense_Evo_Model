cd src

gcc od_worker.c utils.c model.c -o ../od_worker
gcc mc_worker.c utils.c model.c -o ../mc_worker

cd ..