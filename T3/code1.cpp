#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <math.h>
#include <random>
#include <iostream>
#include <mpi.h>
#include <cstdlib>
using namespace std;

float alpha(float x, float y){
    return x * (x - 1) * y * (y - 1) + 1;
}

float** create_matrix(int m, int n){
    float** matrix = (float**) calloc(m, sizeof(float*));
    for (int i = 0; i < m; i++)
    {
        matrix[i] = (float*) calloc(n, sizeof(float));
    }
    return matrix;
}

void free_matrix(float** matrix, int m){
    for (int i = 0; i < m; i++)
    {
        free(matrix[i]);
    }
    free(matrix);
}

float** north(int rows, int n, int first_row){
    float** matrix = create_matrix(rows, n);
    float h = (float) (n+1)/2;
    for (int i = 0; i < rows; i++)
    {   
        for (int j = 0; j < n; j++)
        {
            if (j == n - 1)
            {
                matrix[i][j] = 0;
            }
            else
            {
                float k = i + first_row + 1;
                float l = j + 1;
                matrix[i][j] = - alpha(k * h, (l + 0.5) * h) / pow(h, 2.0);
            }   
        }
    }
    return matrix;    
}

float** south(int rows, int n, int first_row){

    float** matrix = create_matrix(rows, n);
    float h = (float) (n+1)/2;
    for (int i = 0; i < rows; i++)
    {   
        for (int j = 0; j < n; j++)
        {
            if (j == 0)
            {
                matrix[i][j] = 0;
            }
            else
            {
                float k = i + first_row + 1;
                float l = j + 1;
                matrix[i][j] = - alpha(k * h, (l - 0.5) * h) / pow(h, 2.0);
            }   
        }
    }   
    return matrix; 
}

float** east(int rows, int n, int first_row, int world_rank, int world_size){
    float** matrix = create_matrix(rows, n);
    float h = (float) (n+1)/2;
    for (int i = 0; i < rows; i++)
    {   
        for (int j = 0; j < n; j++)
        {
            if (world_rank == world_size - 1 && i == rows - 1)
            {
                matrix[i][j] = 0;
            }
            else
            {
                float k = i + first_row + 1;
                float l = j + 1;
                matrix[i][j] = - alpha((k + 0.5) * h, l * h) / pow(h, 2.0);
            }   
        }
    }   
    return matrix; 
}

float** west(int rows, int n, int first_row){
    float** matrix = create_matrix(rows, n);
    float h = (float) (n+1)/2;
    for (int i = 0; i < rows; i++)
    {   
        for (int j = 0; j < n; j++)
        {
            if (j == 0)
            {
                matrix[i][j] = 0;
            }
            else
            {
                float k = i + first_row + 1;
                float l = j + 1;
                matrix[i][j] = - alpha((k - 0.5) * h, l * h) / pow(h, 2.0);
            }   
        }
    }   
    return matrix; 
}

float** center(int rows, int n, int first_row){
    float** matrix = create_matrix(rows, n);
    float h = (float) (n+1)/2;
    for (int i = 0; i < rows; i++)
    {   
        for (int j = 0; j < n; j++)
        {
            float k = i + first_row + 1;
            float l = j + 1;
            matrix[i][j] = (alpha((k - 0.5) * h, l * h) + alpha((k + 0.5) * h, l * h) +
                            alpha(k * h, (j - 0.5) * h) + alpha(k * h, (j + 0.5) * h)) / 
                            pow(h, 2.0);  
        }
    }   
    return matrix; 
}

float** create_x(int rows, int n, int world_rank, int world_size){
    x = create_matrix(rows, n + 2);
    std::default_random_engine generator;
    std::uniform_real_distribution<float> distribution(0.0, 1.0);
    if (world_rank == 0)
    {
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < n + 2; j++)
            {
                if (j == 0 || j == n + 1 || i == 0)
                {
                    x[i][j] = 0;
                }
                else
                {
                    x[i][j] = distribution(generator);
                }
            }
        }
    }
    else if (world_rank == world_size - 1)
    {
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < n + 2; j++)
            {
                if (j == 0 || j == n + 1 || i == rows - 1)
                {
                    x[i][j] = 0;
                }
                else
                {
                    x[i][j] = distribution(generator);
                }
            }
        }
    }
    else
    {
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < n + 2; j++)
            {
                if (j == 0 || j == n + 1)
                {
                    x[i][j] = 0;
                }
                else
                {
                    x[i][j] = distribution(generator);
                }
            }
        }
    }
    return x;
}

float** center_matvec(float** C, float** N, float** S, float** E, float** W, 
                float** x, float* upper_x, float* lower_x, int rows, int n){
    // Matvec para bloques al centro. 
    float** b = create_matrix(rows, n);
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < n; j++)
        {
            int k = i + 1;
            int l = j + 1;
            if (i == 0)
            {
                b[i][j] = C[i][j] * x[k][l] + N[i][j] * x[k][l + 1] + 
                          S[i][j] * x[k][l - 1] + E[i][j] * x[k + 1][l] + 
                          W[i][j] * upper_x[l];
            }
            else if (i == rows - 1)
            {
                b[i][j] = C[i][j] * x[k][j] + N[i][j] * x[k][l + 1] + 
                          S[i][j] * x[k][l - 1] + W[i][j] * x[k - 1][l] + 
                          E[i][j] * lower_x[l];
            }
            else
            {
                b[i][j] = C[i][j] * x[k][l] + N[i][j] * x[k][l + 1] + 
                          S[i][j] * x[k][l - 1] + W[i][j] * x[k - 1][l] + 
                          E[i][j] * x[k + 1][l];
            }   
        }
    }
    return b;
}

float** upper_matvec(float** C, float** N, float** S, float** E, float** W, 
                float** x, float* lower_x, int rows, int n){
    // Matvec para el primer bloque. 
    float** b = create_matrix(rows, n);
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < n; j++)
        {
            int k = i + 1;
            int l = j + 1;
            if (i == rows - 1)
            {
                b[i][j] = C[i][j] * x[k][l] + N[i][j] * x[k][l + 1] + 
                          S[i][j] * x[k][l - 1] + W[i][j] * x[k - 1][l] + 
                          E[i][j] * lower_x[l];
            }
            else
            {
                b[i][j] = C[i][j] * x[k][l] + N[i][j] * x[k][l + 1] + 
                          S[i][j] * x[k][l - 1] + W[i][j] * x[k - 1][l] + 
                          E[i][j] * x[k + 1][l];
            }
        }
    }
    return b;
}

float** lower_matvec(float** C, float** N, float** S, float** E, float** W, 
                float** x, float* lower_x, int rows, int n){
    // Matvec para el ultimo bloque. 
    float** b = create_matrix(rows, n);
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < n; j++)
        {
            int k = i + 1;
            int l = j + 1;
            if (i == 0)
            {
                b[i][j] = C[i][j] * x[k][l] + N[i][j] * x[k][l + 1] + 
                          S[i][j] * x[k][l - 1] + E[i][j] * x[k + 1][l] + 
                          W[i][j] * lower_x[l];
            }
            else
            {
                b[i][j] = C[i][j] * x[k][l] + N[i][j] * x[k][l + 1] + 
                          S[i][j] * x[k][l - 1] + W[i][j] * x[k - 1][l] + 
                          E[i][j] * x[k + 1][l];
            }
        }
    }
    return b;
}

float vecvec(float* x, float* y, int n){
    // Calcula el producto vector vector x^T y
    float z = 0;
    for (int i = 0; i < n; i++)
    {
        z += x[i] * y[i];
    }
    return z;
}

float** vecsum(float** x, float** y, int rows, int n){
    // Calcula la suma de dos vectores x + y
    float** z = create_matrix(rows, n);
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < n; j++)
        {
            z[i][j] = x[i][j] + y[i][j];
        }
    }
    return z;
}

float** vecdif(float** x, float** y, int rows, int n){
    // Calcula la resta de dos vectores x - y
    float** z = create_matrix(rows, n);
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < n; j++)
        {
            z[i][j] = x[i][j] - y[i][j];
        }
    }
    return z;
}

float* scalevec(float* v, float l, int n){
    // Pondera el vector v por l
    for (int i = 0; i < n; i++)
    {
        v[i] = v[i] * l;
    }
    return v;
}

float** create_b(int rows, int n, int world_rank){
    float** b = create_matrix(rows, n);
    if (world_rank == 0)
    {
        b[rows - 1][n - 1] = 1;
    }
    return b;
}   

void print_matrix(float** matrix, int filas, int columnas) {
  printf("\n");
  for (int i = 0; i < filas; i++) {
      for (int j = 0; j < columnas; j++) {
        printf("%f ", matrix[i][j]); 
      }
    printf("\n");
  }
}

int main(){
    MPI_Init(NULL,NULL);
	int world_size, world_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (world_rank == 0) {
        printf("Size: %i\n\n", world_size);
    }

    int first_row, local_rows, n, err;

    n = 10;

    local_rows = n / world_size;
    first_row = world_rank * local_rows;

    if (world_rank == world_size - 1) {
        local_rows += n % world_size;
    }

    // Creamos los arreglos correspondientes
    float** C = center(local_rows, n, first_row);
    float** N = north(local_rows, n, first_row);
    float** S = south(local_rows, n, first_row);
    float** E = east(local_rows, n, first_row, world_rank, world_size);
    float** W = west(local_rows, n, first_row);
    float** x = create_x(local_rows, n, world_rank, world_size);
    float** b = create_b(local_rows, n, world_rank);

    bool no_conv = true
    while (no_conv)
    {

    if (world_rank != 0 && world_rank != world_size - 1)
    {
        float* local_upper_x = x[0];
        float* local_lower_x = x[local_rows - 1];
        float* lower_x = (float*) calloc(n, sizeof(float));
        float* upper_x = (float*) calloc(n, sizeof(float));

        // Comunicarse con el vecino de arriba
        MPI_Sendrecv(local_upper_x, n, MPI_FLOAT, world_rank - 1, world_rank,
                    upper_x, n, MPI_FLOAT, world_rank - 1, world_rank - 1, 
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Comunicarse con el vecino de abajo
        MPI_Sendrecv(local_lower_x, n, MPI_FLOAT, world_rank + 1, world_rank,
                    lower_x, n, MPI_FLOAT, world_rank + 1, world_rank + 1, 
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Calcular el residuo
        float** aprox_b = center_matvec(C, N, S, E, W, x, upper_x, lower_x, rows, n)
        float** r = vecdif(aprox_b, b)
    }
        else if (world_rank==0)
    {
        float* local_lower_x = x[local_rows - 1];
        float* lower_x = (float*) calloc(n, sizeof(float));

        // Comunicarse con el vecino de abajo
        MPI_Sendrecv(local_lower_x, n, MPI_FLOAT, world_rank + 1, world_rank,
                    lower_x, n, MPI_FLOAT, world_rank + 1, world_rank + 1, 
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Calcular el residuo
        float** aprox_b = upper_matvec(C, N, S, E, W, x, lower_x, rows, n)
        float** r = vecdif(aprox_b, b)
    }
        else
    {
        float* local_upper_x = x[0];
        float* upper_x = (float*) calloc(n, sizeof(float));
        
        // Comunicarse con el vecino de arriba
        MPI_Sendrecv(local_upper_x, n, MPI_FLOAT, world_rank - 1, world_rank,
                    upper_x, n, MPI_FLOAT, world_rank - 1, world_rank - 1, 
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Calcular el residuo
        float** aprox_b = lower_matvec(C, N, S, E, W, x, upper_x, rows, n)
        float** r = vecdif(aprox_b, b);
    }    
}