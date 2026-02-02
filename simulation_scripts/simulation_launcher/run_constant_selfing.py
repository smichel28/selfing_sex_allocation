
import argparse
import subprocess

p = argparse.ArgumentParser()
p.add_argument("--model")
p.add_argument("--selfing_rate")
p.add_argument("--delta")
p.add_argument("--nrep", type=int)
args = p.parse_args()

for rep in range(args.nrep):
    subprocess.run(f"slim -d N_G=100000 -d POP_S=10000 -d SELF_R={args.selfing_rate} -d DEL={args.delta} -d C=0.01 -d MUT=0.005 -d SI=0.01 -d N_S=200 {args.model}", shell=True)

