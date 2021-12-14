#! /bin/sh
#module load openmpi-4.1.1+gnu-9.3.0
for type in core socket node
	do for i in ob1 ucx
		 do for j in tcp vader
			do
				namefile="$type-$i-$j"
				echo "$type"
			done
		done
	done

