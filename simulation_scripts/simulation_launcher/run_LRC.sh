for alpha in 0.0 0.3 0.7 0.9; do
	for delta in 0.0 0.2 0.4 0.6; do
		for migration in 0.01 0.05 0.1 0.5; do
		    nohup python3 run_LRC.py --model /home/samuel/selfing_sex_alloc/slim_models/LRC.slim --alpha $alpha --delta $delta --migration $migration --nrep 5 &
	done
done
