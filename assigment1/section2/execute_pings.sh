#! /bin/sh
module load openmpi-4.1.1+gnu-9.3.0
for type in core socket node
	do for i in ob1 ucx
		 do for j in tcp vader
			do
				filename="$type-$i-$j"
				if ["$filename" != "node-ob1-vader"];
				then
					echo "mpirun  --map-by $type --mca pml $i --mca btl self,$j -np 2 ./IMB-MPI1 PingPong -msglog 29" > "csv/openmpi/$filename".csv
					if ["$type" == "node"];
					then
						echo "$1, $2" >> "csv/openmpi/$filename".csv
					else
						echo "$1" >> "csv/openmpi/$filename".csv
					fi
					echo "lamba, bandwith computed by fitting data" >> "csv/openmpi/$filename".csv
					mpirun  --map-by $type --mca pml $i --mca btl self,$j -np 2 ~/mpi-bench/mpi-benchmarks/src_c/IMB-MPI1 PingPong -msglog 29 2>/dev/null | grep -v ^\# | grep -v '^$' | tr -s ' ' | sed 's/  */,/g' | cut -c 2- >> "csv/openmpi/$filename".csv
					echo "$filename"
				fi
			done
		done
	done

