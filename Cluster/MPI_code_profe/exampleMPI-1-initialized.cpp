// Compile the code: mpic++ helloMPI.cpp
// Run the executable: mpirun -np 2 a.out

#include <iostream>
using namespace std;

#include <mpi.h>

int main()
{
    int flagInit = -1;
    int flagFinalize = -1;
    int err = 0;

    err = MPI_Initialized(&flagInit);
    err = MPI_Finalized(&flagFinalize);

    cout << "MPI is initialized (before Init): " << flagInit << endl;
    cout << "MPI is finalized (before Init): " << flagFinalize << endl;

    MPI_Init(NULL, NULL);
    
    err = MPI_Initialized(&flagInit);
    err = MPI_Finalized(&flagFinalize);

    cout << "MPI is initialized (in MPI block): " << flagInit << endl;
    cout << "MPI is finalized (in MPI block): " << flagFinalize << endl;
    
    MPI_Finalize();
    
    err = MPI_Initialized(&flagInit);
    err = MPI_Finalized(&flagFinalize);

    cout << "MPI is initialized (after Finalize): " << flagInit << endl;
    cout << "MPI is finalized (after Finalize): " << flagFinalize << endl;
}

