// Copyright 2020 Chirkov Roman
#include <mpi.h>
#include <iostream>
#include <vector>
#include <ctime>
#include <random>
#include "../../../modules/task_2/chirkov_r_vertical_gauss/vertical_gauss.h"

double **generateMatrix(int size)
{
    double **matrix;
    matrix = new double * [size];
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    for (int i = 0; i < size; i++)
    {
        matrix[i] = new double[size];
        for (int j = 0; j < size; j++)
        {
            matrix[i][j] = gen();
        }
    }
    return matrix;
}

//Добавить проверку длины вектора!!!

std::vector<double> sequentialGauss(std::vector<std::vector<double> >  matrix, std::vector<double> vector, int size) {
	std::vector<double>  x(size);
    double max;
    int k, rowMax;
	const double accuracy = 0.00001;  // точность
	k = 0;
	while (k < size) {
		// Поиск строки с максимальным a[i][k] среди оставшихся строк
		max = std::abs(matrix[k][k]);
		rowMax = k;
		for (int i = k + 1; i < size; i++) 
		{
			if (std::abs(matrix[i][k]) > max)
			{
				max = std::abs(matrix[i][k]);
				rowMax = i;
			}
		}
		// Перестановка строк
		if (max < accuracy) {
			//добавить выброс исключения!!!!!!!!!!!!!!!!
			// нет ненулевых диагональных элементов
			std::cout << "Решение получить невозможно из-за нулевого столбца ";
			std::cout << rowMax << " матрицы A" << std::endl;
			std::vector<double>  x(0);
			return x;
		}
		for (int j = 0; j < size; j++) {//перестановка строк матрицы
			double temp = matrix[k][j];
			matrix[k][j] = matrix[rowMax][j];
			matrix[rowMax][j] = temp;
		}
		double temp = vector[k];// перестановка соответствующего элемента вектора
		vector[k] = vector[rowMax];
		vector[rowMax] = temp;

		// Нормализация уравнений
		//переделать в короткий вариант
		for (int i = k; i < size; i++) {
			double temp = matrix[i][k];
			if (std::abs(temp) < accuracy) continue; // для нулевого коэффициента пропустить
			for (int j = 0; j < size; j++) {
				matrix[i][j] = matrix[i][j] / temp;
			}
			vector[i] = vector[i] / temp;
			if (i == k)  continue; // уравнение не вычитать само из себя
			for (int j = 0; j < size; j++)
				matrix[i][j] = matrix[i][j] - matrix[k][j];
			vector[i] = vector[i] - vector[k];
		}
		k++;
	}
  // обратная подстановка
	for (k = size - 1; k >= 0; k--) {
		x[k] = vector[k];
		for (int i = 0; i < k; i++) {
			vector[i] = vector[i] - matrix[i][k] * x[k];
		}
	}
	return x;

}
/*
int parallelFind(std::vector<int> matrix, int rows, int cols)
{
    int size, rank;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Bcast(&rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&cols, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int *sendcounts = new int[size];
    for (int i = 0; i < size; i++)
    {
        sendcounts[i] = rows * cols / size;
    }
    sendcounts[size - 1] += rows * cols % size;

    int *displs = new int[size];
    int sum = 0;
    for (int i = 0; i < size; i++)
    {
        displs[i] = sum;
        sum += sendcounts[i];
    }

    int recvcount = rows * cols / size;
    if (rank == size - 1)
    {
        recvcount += rows * cols % size;
    }

    int *recvbuf = new int[recvcount];

    MPI_Scatterv(rank == 0 ? &matrix[0] : 0, sendcounts, displs, MPI_INT,
                 recvbuf, recvcount, MPI_INT, 0, MPI_COMM_WORLD);

    int partMax = -2147483648;
    for (int i = 0; i < recvcount; i++)
    {
        if (recvbuf[i] > partMax)
        {
            partMax = recvbuf[i];
        }
    }

    delete[] sendcounts;
    delete[] displs;
    delete[] recvbuf;

    int globalMax;
    MPI_Allreduce(&partMax, &globalMax, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    return globalMax;
}*/
