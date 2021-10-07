// Compile the code: mpic++ helloMPI.cpp
// Run the executable: mpirun -np 2 a.out

#include <iostream>
using namespace std;

#include <mpi.h>

int main()
{
    int err;
	int rank,size;

    MPI_Init(NULL, NULL);
    
	err = MPI_Comm_size(MPI_COMM_WORLD, &size);
	err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    cout << "Size of communicator: " << size << endl;
    cout << "Rank of communicator: " << rank << endl;

    MPI_Finalize();
    
}

