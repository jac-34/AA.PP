// Compile the code: mpic++ helloMPI.cpp
// Run the executable: mpirun -np 2 a.out

#include <iostream>
using namespace std;

#include <mpi.h>

int main()
{
    int err;
	char procname[MPI_MAX_PROCESSOR_NAME];
	int name_len;

    MPI_Init(NULL, NULL);
    
	err = MPI_Get_processor_name(procname, &name_len);
	cout << "Processor name: " << procname << endl;

    MPI_Finalize();
    
}

