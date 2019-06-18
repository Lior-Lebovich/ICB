#!/bin/bash
job='job_'
output='output_'
dat='motor_both9_'
numrun=0
qsub -N "$job$dat$numrun" -o "$output$dat$numrun" -b y "source ~/envs/yon/bin/activate; python motor_both9.py "$numrun""
numrun=1
qsub -N "$job$dat$numrun" -o "$output$dat$numrun" -b y "source ~/envs/yon/bin/activate; python motor_both9.py "$numrun""
numrun=2
qsub -N "$job$dat$numrun" -o "$output$dat$numrun" -b y "source ~/envs/yon/bin/activate; python motor_both9.py "$numrun""
numrun=3
qsub -N "$job$dat$numrun" -o "$output$dat$numrun" -b y "source ~/envs/yon/bin/activate; python motor_both9.py "$numrun""
numrun=4
qsub -N "$job$dat$numrun" -o "$output$dat$numrun" -b y "source ~/envs/yon/bin/activate; python motor_both9.py "$numrun""
numrun=5
qsub -N "$job$dat$numrun" -o "$output$dat$numrun" -b y "source ~/envs/yon/bin/activate; python motor_both9.py "$numrun""
numrun=6
qsub -N "$job$dat$numrun" -o "$output$dat$numrun" -b y "source ~/envs/yon/bin/activate; python motor_both9.py "$numrun""
numrun=7
qsub -N "$job$dat$numrun" -o "$output$dat$numrun" -b y "source ~/envs/yon/bin/activate; python motor_both9.py "$numrun""
numrun=8
qsub -N "$job$dat$numrun" -o "$output$dat$numrun" -b y "source ~/envs/yon/bin/activate; python motor_both9.py "$numrun""
numrun=9
qsub -N "$job$dat$numrun" -o "$output$dat$numrun" -b y "source ~/envs/yon/bin/activate; python motor_both9.py "$numrun""
numrun=10
qsub -N "$job$dat$numrun" -o "$output$dat$numrun" -b y "source ~/envs/yon/bin/activate; python motor_both9.py "$numrun""
numrun=11
qsub -N "$job$dat$numrun" -o "$output$dat$numrun" -b y "source ~/envs/yon/bin/activate; python motor_both9.py "$numrun""
