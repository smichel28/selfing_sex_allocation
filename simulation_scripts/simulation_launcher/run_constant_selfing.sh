for self_rate in 0.0 0.3 0.7 0.9; do
	nohup python3 run_constant_selfing.py --model /home/samuel/selfing_sex_alloc/slim_models/constant_selfing.slim --selfing_rate $self_rate --delta 0.2 --nrep 5 &
done
