// Compile the code: mpic++ helloMPI.cpp
// Run the executable: mpirun -np 2 a.out

#include <iostream>
using namespace std;

#include <mpi.h>

int main()
{
    int err;
    int version, subversion;

    MPI_Init(NULL, NULL);
    
    err = MPI_Get_version(&version, &subversion);

    cout << "The MPI version is: " << version << "." << subversion << endl;

    MPI_Finalize();
    
}

