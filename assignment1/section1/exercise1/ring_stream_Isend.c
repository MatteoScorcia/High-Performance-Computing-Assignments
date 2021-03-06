#include <stdio.h>
#include <mpi.h>

void execute_mpi_ring(int numprocs, FILE *fptr);

int main(int argc, char *argv[])
{

	FILE *fptr;
	int numprocs;
	fptr = fopen("ring_results.txt", "a");
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	double start_time = MPI_Wtime();
	execute_mpi_ring(numprocs, fptr);
	double elapsed_time = MPI_Wtime();

	double total_time;
	MPI_Reduce(&elapsed_time, &total_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if (my_rank == 0)
		printf("# total time spent for the execution with %d processes: %10.8f \n", numprocs, total_time);

	fclose(fptr);
	MPI_Finalize();
	return 0;
}

void execute_mpi_ring(int numprocs, FILE *fptr)
{
	//1D ring array (periodic boundary)
	int ndims = 1;
	int dims[1] = {0};
	int periods[1] = {1};
	int reorder = 1;

	MPI_Comm ring_communicator;
	MPI_Dims_create(numprocs, ndims, dims);
	MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &ring_communicator);

	int my_rank;
	MPI_Comm_rank(ring_communicator, &my_rank);

	MPI_Request req_left[2];
	MPI_Request req_right[2];

	int direction = 0, disposition = 1;
	int left_neighbour, right_neighbour;

	//find left and right neighbours of the mpi process
	MPI_Cart_shift(ring_communicator, direction, disposition, &left_neighbour, &right_neighbour);

	MPI_Status msg_left_status;
	MPI_Status msg_right_status;
	int msg_left_buf, msg_left_tag;
	int msg_right_buf, msg_right_tag;
	int send_left_buf, send_right_buf;

	int counter_msg = 0;

	for (int iteration = 1; iteration <= numprocs; iteration++)
	{
		// find right tags of messages based on current iteration
		MPI_Cart_shift(ring_communicator, direction, iteration, &msg_left_tag, &msg_right_tag);
		msg_left_tag *= 10;
		msg_right_tag *= 10;

		if (iteration == 1)
		{
			//send message to left neighbour
			send_left_buf = my_rank;
			MPI_Isend(&send_left_buf, 1, MPI_INT, left_neighbour, my_rank * 10, ring_communicator, &req_left[0]);

			//send message to right neighbour
			send_right_buf = -1 * my_rank;
			MPI_Isend(&send_right_buf, 1, MPI_INT, right_neighbour, my_rank * 10, ring_communicator, &req_right[0]);
		}
		else
		{
			//send message to left neighbour with right tag
			send_left_buf = msg_right_buf - my_rank;
			MPI_Isend(&send_left_buf, 1, MPI_INT, left_neighbour, msg_right_status.MPI_TAG, ring_communicator, &req_left[0]);

			//send message to right neighbour with right tag
			send_right_buf = msg_left_buf + my_rank;
			MPI_Isend(&send_right_buf, 1, MPI_INT, right_neighbour, msg_left_status.MPI_TAG, ring_communicator, &req_right[0]);
		}

		//receive message from right neighbour
		MPI_Irecv(&msg_right_buf, 1, MPI_INT, right_neighbour, msg_right_tag, ring_communicator, &req_right[1]);

		//receive message from left neighbour
		MPI_Irecv(&msg_left_buf, 1, MPI_INT, left_neighbour, msg_left_tag, ring_communicator, &req_left[1]);

		//wait for eventual Isend to finish sending before starting a new iteration
		MPI_Wait(&req_left[0], MPI_STATUS_IGNORE);
		MPI_Wait(&req_right[0], MPI_STATUS_IGNORE);

		//wait for Irecv to finish before starting a new iteration
		MPI_Wait(&req_right[1], &msg_right_status);
		MPI_Wait(&req_left[1], &msg_left_status);
		counter_msg += 2;
	}

	fprintf(fptr, "I am process %d and i have received %d messages. My final messages have tag %d and value '%d'-left, '%d'-right\n", my_rank, counter_msg, msg_right_tag, msg_left_buf, msg_right_buf);
}