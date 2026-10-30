#!/bin/bash
# Train a small MACE model for Ti on CPU (run inside mace/). Resumable: re-run the same command.
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
mace_run_train --name=ti_mace --train_file=train.xyz --valid_file=valid.xyz \
  --energy_key=REF_energy --forces_key=REF_forces --E0s=average \
  --model=MACE --r_max=5.0 --num_interactions=2 --num_channels=32 --max_L=0 --correlation=3 \
  --batch_size=4 --valid_batch_size=10 --max_num_epochs=100 --lr=0.01 --amsgrad \
  --ema --ema_decay=0.99 --swa --start_swa=75 --eval_interval=2 \
  --device=cpu --default_dtype=float32 --seed=1 --restart_latest \
  --results_dir=res --checkpoints_dir=ck --model_dir=mod --log_dir=log
