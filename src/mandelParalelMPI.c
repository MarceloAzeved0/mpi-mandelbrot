#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define MAXITER 2000

struct complex
{
    double real;
    double imag;
};

int main(int argc, char *argv[])
{
    int rank;
    int size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    for (int NPOINTS = 500; NPOINTS <= 5000; NPOINTS += 500)
    {
        int i, j, iter, numoutside, value = 0;
        double area, error, ztemp;
        double start, finish;
        struct complex z, c;

        start = MPI_Wtime();

        for (i = rank; i < NPOINTS; i += size)
        {
            for (j = 0; j < NPOINTS; j++)
            {
                c.real = -2.0 + 2.5 * (double)(i) / (double)(NPOINTS) + 1.0e-7;
                c.imag = 1.125 * (double)(j) / (double)(NPOINTS) + 1.0e-7;
                z = c;
                for (iter = 0; iter < MAXITER; iter++)
                {
                    ztemp = (z.real * z.real) - (z.imag * z.imag) + c.real;
                    z.imag = z.real * z.imag * 2 + c.imag;
                    z.real = ztemp;
                    if ((z.real * z.real + z.imag * z.imag) > 4.0e0)
                    {
                        value++;
                        break;
                    }
                }
            }
        }

        MPI_Reduce(&value, &numoutside, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        finish = MPI_Wtime();

        if (rank == 0)
        {
            area = 2.0 * 2.5 * 1.125 * (double)(NPOINTS * NPOINTS - numoutside) / (double)(NPOINTS * NPOINTS);
            error = area / (double)NPOINTS;

            printf("NPOINTS: %d | Area of Mandlebrot set = %12.8f +/- %12.8f\n", NPOINTS, area, error);
            printf("Time = %12.8f seconds\n", finish - start);
        }
    }

    MPI_Finalize();
    return 0;
}
