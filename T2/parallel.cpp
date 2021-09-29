#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstdlib>
#include <random>
#include <chrono>
#include <time.h>
#include <omp.h>

struct complex_number{
    // Estructura que define un numero complejo
    float re;
    float im;
};

typedef struct complex_number Comp;

float norm(Comp n){
    // Calcula la norma de un numero complejo
    float x = pow(n.re, 2) + pow(n.im, 2);
    return pow(x, 0.5);
}

Comp squared(Comp n){
    // Eleva un numero complejo al cuadrado
    Comp n_squared;
    n_squared.re = pow(n.re, 2) - pow(n.im, 2);
    n_squared.im = 2 * n.re * n.im;
    return n_squared;
} 

Comp add(Comp a, Comp b){
    // Suma dos numeros complejos
    Comp c;
    c.re = a.re + b.re;
    c.im = a.im + b.im;
    return c;
}

bool in_mandelbrot(Comp c, int iters){
    // Calcula el iter-esimo termino de la sucesión de mandelbrot empezando con c
    // y retorna un bool que indica si c pertenece al conjunto o no
    Comp x = {0, 0};
    for (int i = 0; i < iters; i++)
    {
        x = add(squared(x), c); 
        if (norm(x) > 2)
        {
            return false;
        } 
    }
    return true;
}

float mandelbrot_area(int n, float total_area, int in){
    // Calcula el area aproximada del conjunto de mandelbrot en base a la proporcion
    // de numeros aleatorios que quedaron dentro del conjunto
    return total_area * ((float) in / n);
}

void print_comp(Comp n){
    // Imprime un numero complejo
    printf("%f + %fi\n", n.re, n.im);
}

int main() {
    const float a = -2.0;
    const float b = 2.0;
    const float total_area = pow(b - a, 2);

    std::default_random_engine generator;
    std::uniform_real_distribution<float> distribution(a, b);

    const int n = 500000;
    const int iters = 2000;
    Comp* nums = (Comp*) malloc(n * sizeof(Comp));
    auto start = std::chrono::high_resolution_clock::now();
    // Generamos n numeros complejos aleatorios
    for (int i = 0; i < n; i++)
    {
        nums[i].re = distribution(generator);
        nums[i].im = distribution(generator);
    }
    auto end_1 = std::chrono::high_resolution_clock::now();
    auto execution_random = std::chrono::duration_cast<std::chrono::nanoseconds>(end_1 - start);
    printf("\nTime measured generating random numbers: %.3f seconds\n", execution_random.count() * 1e-9);
    float previous_time = execution_random.count();

    // Vemos cuantos de los numeros pertenecen al conjunto
    int num_threads_array[4] = {2, 4, 8, 16};
    for (int j = 0; j < 4; j++)    
    {    
        int num_threads = num_threads_array[j];
        int in = 0;
        int* partial_count = (int*) calloc(num_threads, sizeof(int));

        #pragma omp parallel for num_threads(num_threads)
        for (int i = 0; i < n; i++)
        {
            if (in_mandelbrot(nums[i], iters))
            {
                int id = omp_get_thread_num();
                partial_count[id] += 1;
            }
        }
        for (int i = 0; i < num_threads; i++)
        {
            in += partial_count[i];
        }
        free(partial_count);

        float area = mandelbrot_area(n, total_area, in);
        auto end_thread = std::chrono::high_resolution_clock::now();
        auto execution_area = std::chrono::duration_cast<std::chrono::nanoseconds>(end_thread - start);
        float total_time = execution_area.count();
        printf("El area del conjunto de mandelbrot es %f\n", area);
        printf("Time measured calculating area with %i threads: %.3f seconds\n\n", num_threads, (total_time - previous_time) * 1e-9);
        previous_time = total_time;
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto execution = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    printf("Total time measured: %.3f seconds\n\n", execution.count() * 1e-9);
    free(nums);
    return 0;
}