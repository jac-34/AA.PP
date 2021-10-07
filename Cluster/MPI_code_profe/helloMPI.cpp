// Compile the code: mpic++ helloMPI.cpp
// Run the executable: mpirun -np 2 a.out

#include <iostream>
using namespace std;

#include <mpi.h>

int main()
{
    cout << "Hello before MPI initialisation" << endl;
    
    MPI_Init(NULL, NULL);
    
    cout << "Hello from MPI block" << endl;
    
    MPI_Finalize();

    cout << "Hello after MPI finalisation" << endl;
}

