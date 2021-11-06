#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
//#include <chrono>
#include <math.h>
//#include <random>
#include <iostream>
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
    float h = (float) 1 / (n + 1);
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
    float** x;
    srand(world_rank * 10);
    //std::default_random_engine generator;
    //std::uniform_real_distribution<float> distribution(0.0, 1.0);
    if (world_rank == 0)
    {
        x = create_matrix(rows + 1, n + 2);
        for (int i = 0; i < rows + 1; i++)
        {
            for (int j = 0; j < n + 2; j++)
            {
                if (j == 0 || j == n + 1 || i == 0)
                {
                    x[i][j] = 0;
                }
                else
                {
                    x[i][j] = (float) (rand()%100000) / 100000;
                }
            }
        }
    }
    else if (world_rank == world_size - 1)
    {
        x = create_matrix(rows + 1, n + 2);
        for (int i = 0; i < rows + 1; i++)
        {
            for (int j = 0; j < n + 2; j++)
            {
                if (j == 0 || j == n + 1 || i == rows)
                {
                    x[i][j] = 0;
                }
                else
                {
                    x[i][j] = (float) (rand()%100000) / 100000;
                }
            }
        }
    }
    else
    {
        x = create_matrix(rows, n + 2);
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
                    x[i][j] = (float) (rand()%100000) / 100000;
                }
            }
        }
    }
    return x;
}

float** center_matvec(float** C, float** N, float** S, float** E, float** W, 
                float** x, float* upper_x, float* lower_x, int rows, int n){
    // Matvec para bloques al centro. 
    float** b = create_matrix(rows, n + 2);
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < n; j++)
        {
            int k = i;
            int l = j + 1;
            if (i == 0)
            {
                b[k][l] = C[i][j] * x[k][l] + N[i][j] * x[k][l + 1] + 
                          S[i][j] * x[k][l - 1] + E[i][j] * x[k + 1][l] + 
                          W[i][j] * upper_x[l];
            }
            else if (i == rows - 1)
            {
                b[k][l] = C[i][j] * x[k][j] + N[i][j] * x[k][l + 1] + 
                          S[i][j] * x[k][l - 1] + W[i][j] * x[k - 1][l] + 
                          E[i][j] * lower_x[l];
            }
            else
            {
                b[k][l] = C[i][j] * x[k][l] + N[i][j] * x[k][l + 1] + 
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
    float** b = create_matrix(rows + 1, n + 2);
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < n; j++)
        {
            int k = i + 1;
            int l = j + 1;
            if (i == rows - 1)
            {
                b[k][l] = C[i][j] * x[k][l] + N[i][j] * x[k][l + 1] + 
                          S[i][j] * x[k][l - 1] + W[i][j] * x[k - 1][l] + 
                          E[i][j] * lower_x[l];
            }
            else
            {
                b[k][l] = C[i][j] * x[k][l] + N[i][j] * x[k][l + 1] + 
                          S[i][j] * x[k][l - 1] + W[i][j] * x[k - 1][l] + 
                          E[i][j] * x[k + 1][l];
            }
        }
    }
    return b;
}

float** lower_matvec(float** C, float** N, float** S, float** E, float** W, 
                float** x, float* upper_x, int rows, int n){
    // Matvec para el ultimo bloque. 
    float** b = create_matrix(rows + 1, n + 2);
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < n; j++)
        {
            int k = i;
            int l = j + 1;
            if (i == 0)
            {
                b[k][l] = C[i][j] * x[k][l] + N[i][j] * x[k][l + 1] + 
                          S[i][j] * x[k][l - 1] + E[i][j] * x[k + 1][l] + 
                          W[i][j] * upper_x[l];
            }
            else
            {
                b[k][l] = C[i][j] * x[k][l] + N[i][j] * x[k][l + 1] + 
                          S[i][j] * x[k][l - 1] + W[i][j] * x[k - 1][l] + 
                          E[i][j] * x[k + 1][l];
            }
        }
    }
    return b;
}

float** comm_arrays(float** x, int xrows, int n, int world_rank, int world_size){
    // Comunica las fronteras entre los bloques

    float* local_upper_x = x[0];
    float* local_lower_x = x[xrows - 1];
    float* boundaries[2];
    boundaries[0] = (float*) calloc(n, sizeof(float));
    boundaries[1] = (float*) calloc(n, sizeof(float));
    
    if (world_rank != 0)
    {
        // Comunicarse con el vecino de arriba
        MPI_Sendrecv(local_upper_x, n, MPI_FLOAT, world_rank - 1, world_rank,
                     boundaries[1], n, MPI_FLOAT, world_rank - 1, world_rank - 1, 
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    if (world_rank != world_size - 1)
    {
        // Comunicarse con el vecino de abajo
        MPI_Sendrecv(local_lower_x, n, MPI_FLOAT, world_rank + 1, world_rank,
                     boundaries[0], n, MPI_FLOAT, world_rank + 1, world_rank + 1, 
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    return boundaries;
}

float find_local_max(float** r, int rows, int n){
    int max;
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (r[i][j] > max)
            {
                max = r[i][j];
            }
        }
    }
    return max;
}

float find_max(float** r, int rows, int n){
    // Encuentra el elemento maximo de un array 2D distribuido en varios CPUs
    int local_max = find_local_max(r, rows, n);
    int global_max;
    MPI_Allreduce(&local_max, &global_max, 1, MPI_FLOAT, MPI_MAX,
              MPI_COMM_WORLD);
    return global_max;
}

bool is_ready(float** r, int rows, int n, float tol, int world_rank){
    float max = find_max(r, rows, n);
    if (world_rank==0){printf("residuo: %f\n", max);}
    if (max < tol)
    {
        return true;
    }
    return false;
}

float** mat_vec(float** C, float** N, float** S, float** E, float** W, float** x, 
                int rows, int xrows, int n, int world_rank, int world_size){
    
    // Implementa matriz vector de manera compacta
    float** boundaries = comm_arrays(x, xrows, n + 2, world_rank, world_size);
    float* lower_x = boundaries[0];
    float* upper_x = boundaries[1];
    float** b;
    if (world_rank == 0)
    {
        b = upper_matvec(C, N, S, E, W, x, lower_x, rows, n);
    }
    else if (world_rank == world_size - 1)
    {
        b = lower_matvec(C, N, S, E, W, x, lower_x, rows, n);
    }
    else
    {
        b = center_matvec(C, N, S, E, W, x, upper_x, lower_x, rows, n);        
    }
    free(upper_x);
    free(lower_x);
    return b;
}


float vecvec(float** x, float** y, int rows, int n){
    // Calcula el producto vector vector x^T y
    float z = 0;
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < n; j++)
        {
            z += x[i][j] * y[i][j];
        }
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

float** scalevec(float** v, float l, int rows, int n){
    // Pondera el vector v por l
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < n; j++)
        {
            v[i][j] = v[i][j] * l;
        }
    }
    return v;
}

float** create_b(int rows, int n, int world_rank, int world_size){
    float** b;
    if (world_rank == 0 || world_rank == world_size - 1)
    {
        b = create_matrix(rows + 1, n + 2);
        b[rows - 1][n / 2] = 1;
    }
    else{
        b = create_matrix(rows, n + 2);
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

float** copy(float** x, int rows, int n){
    float** copied = create_matrix(rows, n);
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < n; j++)
        {
            copied[i][j] = x[i][j];
        }
    }
    return copied;
}

void print_vector(float* vector, int n){
    for (int j = 0; j < n; j++) {
        printf("%f ", vector[j]); 
    }
    printf("\n");
}

int main(int argc, char** argv) {
    MPI_Init(NULL,NULL);
	int world_size, world_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int err;
	char procname[MPI_MAX_PROCESSOR_NAME];
	int name_len;

    err = MPI_Get_processor_name(procname, &name_len);
	cout << "Processor name: " << procname << endl;

    // Se particionan las filas
    int first_row, local_rows, n, xrows;

    n = 6;

    local_rows = n / world_size;
    xrows = local_rows;             // xrows indica la canditad de filas de la matriz
                                    // x, que es distinta para el primer y ultimo procesador
    if (world_rank == 0 || world_rank == world_size - 1)
    {
        xrows += 1;
    }
    
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
    float** new_x;
    float** b = create_b(local_rows, n, world_rank, world_size);
    float tol = 0.001;

    float** aprox_b = mat_vec(C, N, S, E, W, x, local_rows, xrows, n, world_rank, world_size);
    float** r = vecdif(aprox_b, b, xrows, n + 2);
    float** new_r;
    bool ready = false;
    int iter = 1;
    float beta;
    float rho = 0;
    float old_rho;
    float** p;
    float** new_p;
    float** q = create_matrix(xrows, n + 2);
    float** new_q;
    float delta;

    while (!ready && iter < 100)
    {
        // Implementacion del metodo GC
        old_rho = rho;
        rho = vecvec(r, r, xrows, n + 2);

        if (iter == 1){
            p = copy(r, xrows, n + 2);
        }
        else {
            beta = rho / old_rho;

            new_p = vecdif(r, scalevec(p, beta, xrows, n + 2), xrows, n + 2);
            free_matrix(p, xrows);
            p = new_p;
        }
        new_q = mat_vec(C, N, S, E, W, p, local_rows, xrows, n, world_rank, world_size);
        free_matrix(q, xrows);
        q = new_q;

        delta = rho / vecvec(q, p, xrows, n + 2);

        new_x = vecdif(x, scalevec(p, delta, xrows, n + 2), xrows, n + 2);
        free_matrix(x, xrows);
        x = new_x;

        new_r = vecdif(r, scalevec(q, delta, xrows, n + 2), xrows, n + 2);
        free_matrix(r, xrows);
        r = new_r;

        ready = is_ready(r, xrows, n + 2, tol, world_rank);
        iter++;
    }
    
    free_matrix(C, local_rows);
    free_matrix(N, local_rows);
    free_matrix(S, local_rows);
    free_matrix(E, local_rows);
    free_matrix(W, local_rows);
    free_matrix(b, local_rows);
    free_matrix(x, xrows);
    free_matrix(aprox_b, xrows);
    free_matrix(r, xrows);
    free_matrix(p, xrows);
    free_matrix(q, xrows);


    // Finalize the MPI environment.
    MPI_Finalize();
}