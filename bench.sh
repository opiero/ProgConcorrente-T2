#!/bin/bash

for tam in 1000 5000 10000
do
	cp matriz$(echo $tam).txt matriz.txt
	cp vetor$(echo $tam).txt vetor.txt

	for it in 1 2 3 4 5
	do
		mpirun --hostfile hosts -np 1 ./gaussjordan.out 1
	done

	for no in 2 4 8
	do
		for it in 1 2 3 4 5
		do
			mpirun --hostfile hosts -np $no ./gaussjordan.out 1
		done
	done

	for thr in 4 8
	do
		for it in 1 2 3 4 5
		do
			mpirun --hostfile hosts -np 11 ./gaussjordan.out $thr
		done
	done
	
	mv time.txt time$(echo $tam).txt
done