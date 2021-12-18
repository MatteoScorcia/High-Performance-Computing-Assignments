#! /bin/sh
#just call di script with the first 2 arguments ($1, $2) as name of nodes involved

module load openmpi-4.1.1+gnu-9.3.0
module load intel

for j in tcp shm
	do 
		do for iteration in {1..50}
		filename="core-$j"
		mkdir -p "csv/intel/$filename"
			do
				echo "mpiexec -n 2 -genv I_MPI_PIN_PROCESSOR_LIST=0,2 -genv I_MPI_OFI_PROVIDER $j ./IMB-MPI1 PingPong -msglog 29" > "csv/intel/$filename/$filename-$iteration".csv
				echo "$1" >> "csv/intel/$filename/$filename-$iteration".csv 
				echo "lamba, bandwith computed by fitting data" >> "csv/intel/$filename/$filename-$iteration".csv  
				mpiexec -n 2 -genv I_MPI_PIN_PROCESSOR_LIST=0,2 -genv I_MPI_OFI_PROVIDER $j ~/mpi-bench/mpi-benchmarks/src_c/IMB-MPI1 PingPong -msglog 29 2>/dev/null | grep -v ^\# | grep -v '^$' | tr -s ' ' | sed 's/  */,/g' | cut -c 2- >> "csv/intel/$filename/$filename-$iteration".csv
				echo "$filename-$iteration" 
			done
		done

for j in tcp shm
	do 
		do for iteration in {1..50}
		filename="socket-$j"
		mkdir -p "csv/intel/$filename"
			do
				echo "mpiexec -n 2 -genv I_MPI_PIN_PROCESSOR_LIST=0,1 -genv I_MPI_OFI_PROVIDER $j ./IMB-MPI1 PingPong -msglog 29" > "csv/intel/$filename/$filename-$iteration".csv 
				echo "$1" >> "csv/intel/$filename/$filename-$iteration".csv 
				echo "lamba, bandwith computed by fitting data" >> "csv/intel/$filename/$filename-$iteration".csv 
				mpiexec -n 2 -genv I_MPI_PIN_PROCESSOR_LIST=0,1 -genv I_MPI_OFI_PROVIDER $j ~/mpi-bench/mpi-benchmarks/src_c/IMB-MPI1 PingPong -msglog 29 2>/dev/null | grep -v ^\# | grep -v '^$' | tr -s ' ' | sed 's/  */,/g' | cut -c 2- >> "csv/intel/$filename/$filename-$iteration".csv 
				echo "$filename-$iteration"
			done
		done

for j in tcp shm
	do 
		do for iteration in {1..50}
			filename="node-$j"
			mkdir -p "csv/intel/$filename"
			do
				echo "mpiexec -n 2 -ppn 1 -genv I_MPI_OFI_PROVIDER $j ./IMB-MPI1 PingPong -msglog 29" > "csv/intel/$filename/$filename-$iteration".csv 
				echo "$1, $2" >> "csv/intel/$filename/$filename-$iteration".csv 
				echo "lamba, bandwith computed by fitting data" >> "csv/intel/$filename/$filename-$iteration".csv 
				mpiexec -n 2 -ppn 1 -genv I_MPI_OFI_PROVIDER $j ~/mpi-bench/mpi-benchmarks/src_c/IMB-MPI1 PingPong -msglog 29 2>/dev/null | grep -v ^\# | grep -v '^$' | tr -s ' ' | sed 's/  */,/g' | cut -c 2- >> "csv/intel/$filename/$filename-$iteration".csv 
				echo "$filename-$iteration"
			done
		done

#by core (shm -> shared memory transport)
# mpiexec -n 2 -genv I_MPI_PIN_PROCESSOR_LIST=0,2 ~/mpi-bench/mpi-benchmarks/src_c/IMB-MPI1 PingPong -msglog 29 2>/dev/null | grep -v ^\# | grep -v '^$' | tr -s ' ' | sed 's/  */,/g' | cut -c 2- >> "csv/intel/$filename/$filename-$iteration".csv 

# mpiexec -n 2 -genv I_MPI_PIN_PROCESSOR_LIST=0,2 -genv I_MPI_OFI_PROVIDER tcp ~/mpi-bench/mpi-benchmarks/src_c/IMB-MPI1 PingPong -msglog 29 2>/dev/null | grep -v ^\# | grep -v '^$' | tr -s ' ' | sed 's/  */,/g' | cut -c 2- >> "csv/intel/$filename/$filename-$iteration".csv 

# mpiexec -n 2 -genv I_MPI_PIN_PROCESSOR_LIST=0,2 -genv I_MPI_OFI_PROVIDER shm ~/mpi-bench/mpi-benchmarks/src_c/IMB-MPI1 PingPong -msglog 29 2>/dev/null | grep -v ^\# | grep -v '^$' | tr -s ' ' | sed 's/  */,/g' | cut -c 2- >> "csv/intel/$filename/$filename-$iteration".csv 

#by socket
# mpiexec -n 2 -genv I_MPI_PIN_PROCESSOR_LIST=0,1 ~/mpi-bench/mpi-benchmarks/src_c/IMB-MPI1 PingPong -msglog 29 2>/dev/null | grep -v ^\# | grep -v '^$' | tr -s ' ' | sed 's/  */,/g' | cut -c 2- >> "csv/intel/$filename/$filename-$iteration".csv 

# mpiexec -n 2 -genv I_MPI_PIN_PROCESSOR_LIST=0,1 -genv I_MPI_OFI_PROVIDER tcp ~/mpi-bench/mpi-benchmarks/src_c/IMB-MPI1 PingPong -msglog 29 2>/dev/null | grep -v ^\# | grep -v '^$' | tr -s ' ' | sed 's/  */,/g' | cut -c 2- >> "csv/intel/$filename/$filename-$iteration".csv 

# mpiexec -n 2 -genv I_MPI_PIN_PROCESSOR_LIST=0,1 -genv I_MPI_OFI_PROVIDER shm ~/mpi-bench/mpi-benchmarks/src_c/IMB-MPI1 PingPong -msglog 29 2>/dev/null | grep -v ^\# | grep -v '^$' | tr -s ' ' | sed 's/  */,/g' | cut -c 2- >> "csv/intel/$filename/$filename-$iteration".csv 

#by node
# mpiexec -n 2 -ppn 1 ~/mpi-bench/mpi-benchmarks/src_c/IMB-MPI1 PingPong -msglog 29 2>/dev/null | grep -v ^\# | grep -v '^$' | tr -s ' ' | sed 's/  */,/g' | cut -c 2- >> "csv/intel/$filename/$filename-$iteration".csv 

# mpiexec -n 2 -ppn 1 -genv I_MPI_OFI_PROVIDER tcp ~/mpi-bench/mpi-benchmarks/src_c/IMB-MPI1 PingPong -msglog 29 2>/dev/null | grep -v ^\# | grep -v '^$' | tr -s ' ' | sed 's/  */,/g' | cut -c 2- >> "csv/intel/$filename/$filename-$iteration".csv 

# mpiexec -n 2 -ppn 1 -genv I_MPI_OFI_PROVIDER shm ~/mpi-bench/mpi-benchmarks/src_c/IMB-MPI1 PingPong -msglog 29 2>/dev/null | grep -v ^\# | grep -v '^$' | tr -s ' ' | sed 's/  */,/g' | cut -c 2- >> "csv/intel/$filename/$filename-$iteration".csv 