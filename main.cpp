#include <iostream>
#include<string>
#include<string.h>
#include<mpi.h>
#include<chrono>

using namespace std;


int* getMatrix(int n, int m)
{
    return new int[n*m];
}

void showMatrix(int *matrix, int n, int m)
{
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < m; j++)
        {
            cout<<matrix[i*m +j]<<" ";
        }
        cout<<endl;
    }
}

void initMatrix(int *matrix,int m, int n)
{
    for(int i = 0; i < m*n; i++)
    {
        matrix[i] = std::rand() %20 +1;
    }
}

void initMatrixZero(int *matrix,int m, int n)
{
    for(int i = 0; i < m*n; i++)
    {
        matrix[i] = 0;
    }
}

void multMatrix2x2(int *matrixA, int *matrixB, int *matrixC)
{
    for(int i = 0; i < 2; i++)
    {
        for(int j =0; j < 2; j++)
        {
            for(int k =0; k < 2; k++)
            {
                matrixC[i*2 + j] += matrixA[i*2 + k] * matrixB[k*2 + j];
            }
        }
    }
}

void multMatrix(int *matrixA, int *matrixB, int *matrixC, int m1, int n1,int m2, int n2)
{
    for(int i = 0; i < m1; i++)
    {
        for(int j = 0; j < n2; j++)
        {
            for(int k =0; k < n1; k++)
            {
                matrixC[i*n2 + j] += matrixA[i*n1 + k] * matrixB[k*n2 + j];
            }
        }
    }
}

void copyMatrix(int *matrixA, int* matrixB, int n, int m)
{
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < m; j++)
        {
            matrixA[i*n +j] = matrixB[i*n +j];
        }
    }
}

void addMat(int *matR,int *matC,int offset,int size)
{
    for(int i = 0; i < size; i++)
    {
        matR[offset + i] = matC[i];
    }
}

bool compareMatrix(int *matA,int *matB,int m,int n)
{
    int fl = 0;
    for(int i = 0; i < m*n; i++)
    {
            if(matA[i] != matB[i])
            {
                fl = 1;
            }
    }
    if(fl == 0)
    {
        return true;
    }
    else
    {
        return false;
    }

}

int main(int argc, char* argv[])
{
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int m = 2000;
    int n = 1000;
    int r = m/(size-1);

    if(rank == 0)
    {
        int *matA = getMatrix(m,n);
        initMatrix(matA,m,n);

        int *matB = getMatrix(n,m);
        initMatrix(matB,n,m);

        int *matG = getMatrix(m,m);
        initMatrixZero(matG,m,m);

        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        multMatrix(matA,matB,matG,m,n,n,m);
        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
        cout <<"time one process = "<<time_span.count()<<"seconds"<<endl;
        std::chrono::high_resolution_clock::time_point t3 = std::chrono::high_resolution_clock::now();
        for(int i = 1; i < size; i++)
        {
            MPI_Send(matA + ((i-1)*r*n),r*n,MPI_INT,i,0,MPI_COMM_WORLD);
            MPI_Send(matB,n*m,MPI_INT,i,0,MPI_COMM_WORLD);
        }
        int *matC = getMatrix(r,m);
        int *matR = getMatrix(m,m);
        MPI_Status status;
        for(int i = 1; i < size; i++)
        {
            memset(matC,0,r*n*sizeof(int));
            MPI_Recv(matC,r*m,MPI_INT,i,0,MPI_COMM_WORLD,&status);
            addMat(matR,matC,(i-1)*r*m,r*m);
        }
        std::chrono::high_resolution_clock::time_point t4 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time_span2 = std::chrono::duration_cast<std::chrono::duration<double>>(t4 - t3);
        cout <<"time MPI = "<<time_span2.count()<<"seconds"<<endl;
        if(compareMatrix(matG,matR,m,m))
        {
            cout << "eq" <<endl;
        }
        else
        {
            cout << "no eq "<<endl;
        }

    }
    else
    {
        int *matA = new int[r*n];
        int *matB = new  int[n*m];
        MPI_Status status;
        MPI_Recv(matA,r*n,MPI_INT,0,0,MPI_COMM_WORLD,&status);
        MPI_Recv(matB,n*m,MPI_INT,0,0,MPI_COMM_WORLD,&status);
        int *matC = getMatrix(r,m);
        initMatrixZero(matC,r,m);
        multMatrix(matA,matB,matC,r,n,n,m);
        MPI_Send(matC,r*m,MPI_INT,0,0,MPI_COMM_WORLD);
    }
    MPI_Finalize();
    return 0;
}
