#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main() {
    int* arr = (int*) malloc(5 * sizeof(int));
    for (int i = 0; i < 5; i++)
    {
        arr[i] = pow(i, 2);
        printf("arr[%i] = %i\n", i, arr[i]);
    }
    return 0;
}