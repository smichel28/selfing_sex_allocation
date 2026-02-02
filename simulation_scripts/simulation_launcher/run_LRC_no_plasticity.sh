for alpha in 0.0 0.3 0.7; do
	for migration in 0.01 0.1 1.0; do
		nohup python3 run_LRC_no_plasticity.py --model /home/samuel/selfing_sex_alloc/slim_models/LRC_no_plasticity.slim --alpha $alpha --delta 0.2 --migration $migration --nrep 3 &
	done
done
