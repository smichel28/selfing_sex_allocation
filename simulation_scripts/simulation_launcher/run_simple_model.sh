for alpha in 0.0 0.3 0.7 0.9; do
	for delta in 0.0 0.2 0.4 0.6; do
		nohup python3 run_simple_model.py --model simple_model.slim --alpha $alpha --delta $delta --nrep 1 &
	done
done
