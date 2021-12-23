#! /bin/sh
#just call di script with the first 2 arguments ($1, $2) as name of nodes involved
# remember to do a make before the run!

module load openmpi-4.1.1+gnu-9.3.0

for type in core socket node
do 
	for iteration in {1..5}
	do
		for btl in tcp vader
		do
			filename=$type-ob1-$btl
			mkdir -p "csv/openmpi/$filename"

			echo "mpirun --map-by $type --mca pml ob1 --mca btl self,$btl -np 2 ./IMB-MPI1 PingPong -msglog 29" > "csv/openmpi/$filename/$filename-$iteration".csv

			if [ "$type" == "node" ];
			then
				echo "$1, $2" >> "csv/openmpi/$filename/$filename-$iteration".csv
			else
				echo "$1" >> "csv/openmpi/$filename/$filename-$iteration".csv
			fi

			echo "lamba, bandwith computed by fitting data" >> "csv/openmpi/$filename/$filename-$iteration".csv
			echo "$filename"

			mpirun --map-by $type --mca pml ob1 --mca btl self,$btl -np 2 ~/mpi-bench/mpi-benchmarks/src_c/IMB-MPI1 PingPong -msglog 29 2>/dev/null | grep -v ^\# | grep -v '^$' | tr -s ' ' | sed 's/  */,/g' | cut -c 2- >> "csv/openmpi/$filename/$filename-$iteration".csv
		done

		filename=$type-ucx
		mkdir -p "csv/openmpi/$filename"

		echo "mpirun --map-by $type --mca pml ucx -np 2 ./IMB-MPI1 PingPong -msglog 29 " > "csv/openmpi/$filename/$filename-$iteration".csv

		if [ "$type" == "node" ];
		then
			echo "$1, $2" >> "csv/openmpi/$filename/$filename-$iteration".csv
		else
			echo "$1" >> "csv/openmpi/$filename/$filename-$iteration".csv
		fi

		echo "lamba, bandwith computed by fitting data" >> "csv/openmpi/$filename/$filename-$iteration".csv 
		echo "$filename"

		mpirun --map-by $type --mca pml ucx -np 2 ~/mpi-bench/mpi-benchmarks/src_c/IMB-MPI1 PingPong -msglog 29 2>/dev/null | grep -v ^\# | grep -v '^$' | tr -s ' ' | sed 's/  */,/g' | cut -c 2- >> "csv/openmpi/$filename/$filename-$iteration".csv
	done
done

