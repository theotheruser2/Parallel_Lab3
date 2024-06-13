#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <unistd.h>
#include <time.h>

#if defined(_OPENMP)
#include <omp.h>
#define SCHED_TYPE guided
#define CHUNK_SIZE 20
#else
int omp_get_thread_num() { return -1; }
int omp_set_num_threads(int num){return 0;}
int omp_set_nested(int num){return 0;}
double omp_get_wtime(){ return ((double) clock())/CLOCKS_PER_SEC; }
#endif

int i;
void checkProgress(){
    while (i<100){
        sleep(1);
        printf("\nВыполняется эксперимент №%d/100\n", i);
    } 
}
double minArr(double array[], int size) {
    double min = array[0];
    int i;
    #if defined(CHUNK_SIZE)
    #pragma omp parallel for private(i) shared(size, array) reduction(min: min) schedule(SCHED_TYPE, CHUNK_SIZE)
    #elif defined (SCHED_TYPE)
    #pragma omp parallel for private(i) shared(size, array) reduction(min: min) schedule(SCHED_TYPE)
    #endif
    for (i = 0; i < size - 1; i++) {
        if (array[i] < min){
            min = array[i];
        }
    }
    return min;
}

void swap(double *a, double *b) {
    double temp = *a;
    *a = *b;
    *b = temp;
}

void selectionSort(double array[], int size) {
    int step, i, min_idx;
    #if defined(CHUNK_SIZE)
    #pragma omp parallel for private(step) shared(size, array, min_idx) schedule(SCHED_TYPE, CHUNK_SIZE) 
    #elif defined (SCHED_TYPE)
    #pragma omp parallel for private(step) shared(size, array, min_idx) schedule(SCHED_TYPE)
    #endif
    for (step = 0; step < size - 1; step++) {
        min_idx = step;
        for (i = step + 1; i < size; i++) {
            if (array[i] < array[min_idx])
                min_idx = i;
        }
        swap(&array[min_idx], &array[step]);
    }
}

int main(int argc, char *argv[]) {
    int N, j;
    unsigned int seed = 2;
    double delta_ms;
    N = atoi(argv[1]); // N равен первому параметру командной строки
    double start, end; 
    start = omp_get_wtime(); 
    omp_set_num_threads(16);
    omp_set_nested(1);
    #if defined(_OPENMP)
    #pragma omp parallel
    {
    #pragma omp master
        {
            checkProgress();
        }
    #endif
    //#pragma omp single
    if(omp_get_thread_num() != 0){
    for (i = 0; i < 100; i++) { // 100 экспериментов
        srand(seed^2);
        // Этап 1
        double M1[N], M2[N/2], M2Copy[N/2];
        #if defined(CHUNK_SIZE)
        #pragma omp parallel for private(j) shared(seed, N, M1) schedule(SCHED_TYPE, CHUNK_SIZE)
        #elif defined (SCHED_TYPE)
        #pragma omp parallel for private(j) shared(seed, N, M1) schedule(SCHED_TYPE)
        #endif

        for (j = 0; j < N; j++) {
            M1[j] = rand_r(&seed) % 240 + 1;
        }
        #if defined(CHUNK_SIZE)
        #pragma omp parallel for private(j) shared(seed, N, M2) schedule(SCHED_TYPE, CHUNK_SIZE)
        #elif defined (SCHED_TYPE)
        #pragma omp parallel for private(j) shared(seed, N, M2) schedule(SCHED_TYPE)
        #endif
        for (j = 0; j < N/2; j++) {
            M2[j] = rand_r(&seed) % 2160 + 240;
        }
        // Этап 2
        #if defined(CHUNK_SIZE)
        #pragma omp parallel for private(j) shared(N, M1)  schedule(SCHED_TYPE, CHUNK_SIZE)
        #elif defined (SCHED_TYPE)
        #pragma omp parallel for private(j) shared(N, M1)  schedule(SCHED_TYPE)
        #endif
        for (j = 0; j < N; j++) {
            M1[j] = cbrt(M1[j]/ exp(1));
        }
        #if defined(CHUNK_SIZE)
        #pragma omp parallel for private(j) shared(N, M2, M2Copy) schedule(SCHED_TYPE, CHUNK_SIZE)
        #elif defined (SCHED_TYPE)
        #pragma omp parallel for private(j) shared(N, M2, M2Copy) schedule(SCHED_TYPE)
        #endif
        for (j = 0; j < N/2; j++) {
            M2Copy[j] = M2[j];
        }
        M2[0] = pow(log10(M2[0]), exp(1));
        #if defined(CHUNK_SIZE)
        #pragma omp parallel for private(j) shared(N, M2, M2Copy) schedule(SCHED_TYPE, CHUNK_SIZE)
        #elif defined (SCHED_TYPE)
        #pragma omp parallel for private(j) shared(N, M2, M2Copy) schedule(SCHED_TYPE)
        #endif
        for (j = 1; j < N/2; j++) {
            //printf("\nThread Num: %d\n", omp_get_thread_num());
            M2[j] = pow(log10(M2[j] + M2Copy[j-1]), exp(1));
        }
        // Этап 3
        #if defined(CHUNK_SIZE)
        #pragma omp parallel for private(j) shared(N, M2, M1) schedule(SCHED_TYPE, CHUNK_SIZE)
        #elif defined (SCHED_TYPE)
        #pragma omp parallel for private(j) shared(N, M2, M1) schedule(SCHED_TYPE)
        #endif
        for (j = 0; j < N/2; j++) {
            M2[j] = fabs(M1[j] - M2[j]);
        }
        // Этап 4
        double firstHalf[N/4];
        double secondHalf[N/2 - N/4];
        #if defined(CHUNK_SIZE)
        #pragma omp parallel for private(j) shared(N, M2, firstHalf) schedule(SCHED_TYPE, CHUNK_SIZE)
        #elif defined (SCHED_TYPE)
        #pragma omp parallel for private(j) shared(N, M2, firstHalf) schedule(SCHED_TYPE)
        #endif
        for (j = 0; j < N/4; j++) {
            firstHalf[j] = M2[j];
        }
        #if defined(CHUNK_SIZE)
        #pragma omp parallel for private(j) shared(N, M2, secondHalf) schedule(SCHED_TYPE, CHUNK_SIZE)
        #elif defined (SCHED_TYPE)
        #pragma omp parallel for private(j) shared(N, M2, secondHalf) schedule(SCHED_TYPE)
        #endif
        for (j = 0; j < N/2 - + N/4; j++) {
            secondHalf[j] = M2[j + N/4];
        }
        #if defined(_OPENMP)
        #pragma omp parallel sections
        #else
            selectionSort(firstHalf, N/4);
            selectionSort(secondHalf, N/2 - N/4);
        #endif
        {
            #pragma omp section
            {
                selectionSort(firstHalf, N/4);
            }
            #pragma omp section
            {
                selectionSort(secondHalf, N/2 - N/4);
            }

        }

        #if defined(CHUNK_SIZE)
        #pragma omp parallel for private(j) shared(N, M2, firstHalf, secondHalf) schedule(SCHED_TYPE, CHUNK_SIZE)
        #elif defined (SCHED_TYPE)
        #pragma omp parallel for private(j) shared(N, M2, firstHalf, secondHalf) schedule(SCHED_TYPE)
        #endif
        for (j = 0; j < N/2; j++) {
            if (j < N/4){
                M2[j] = firstHalf[j];
            } else{
                M2[j] = secondHalf[j - N/4];
            }
        }
        // Этап 5
        double minM2 = minArr(M2, N/2);
        double X = 0;
        j = 0;
        #if defined(CHUNK_SIZE)
        #pragma omp parallel for private(j) shared(N, M2, minM2) reduction(+: X) schedule(SCHED_TYPE, CHUNK_SIZE)
        #elif defined (SCHED_TYPE)
        #pragma omp parallel for private(j) shared(N, M2, minM2) reduction(+: X) schedule(SCHED_TYPE)
        #endif
        for (j = 0; j < N/2; j++) {
            if (fmod(floor(M2[j]/minM2), 2) == 0){
                X += sin(M2[j]);
            }
        }
        if (i < 100)
        printf("%f, ", X);
        else i = 99;
    }
    }
    #if defined(_OPENMP)
    }
    #endif
    end = omp_get_wtime(); 
    delta_ms = (end - start) * 1000;
    printf("%f\n", delta_ms);
    return 0;
}


