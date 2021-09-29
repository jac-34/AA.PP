#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <random>
#include <omp.h>

int main() {
    //std::default_random_engine generator;
    //std::uniform_real_distribution<double> distribution(0, 1);


    int* arr = (int*) malloc(10 * sizeof(int));
    #pragma omp parallel for num_threads(10)
    for (int i = 0; i < 10; i++)
    {
        arr[i] = rand();
    }
    for (int i = 0; i < 10; i++)
    {
        printf("arr[%i] = %i\n", i, arr[i]);
    }
    
    free(arr);
    return 0;
}