// Copyright 2020 Chirkov Roman
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <ctime>
#include <iostream>
#include <vector>
#include <random>
#include "./vertical_gauss.h"

TEST(max_value, Test_sequential) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        const int size = 3;
		std::vector<double> x(size);
		std::vector<std::vector<double> > matrix = {
			{2, -3, 7},
			{-4, 6, -2},
			{1, 4, -1}
		};
		//double * vector;
		std::vector<double> vector = {17, 2, 6};
		x = sequentialGauss(matrix, vector, size);
		
		for(int i = 0; i < size; i++){
			std::cout << x[i] << " ";
		}
		std::cout << std::endl;
		
		std::vector<double> staticResult = {1, 2, 3};
        ASSERT_EQ(x, staticResult); //переделать
    }
}
/*
TEST(max_value, Test_oneElement) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        const int rows = 1;
        const int cols = 1;
        std::vector<int> matrix = {1};
        int max = parallelFind(matrix, rows, cols);

        ASSERT_EQ(max, 1);
    } else {
        parallelFind({{}}, 0, 0);
    }
}

TEST(max_value, Test_matrix10x10) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        const int rows = 10;
        const int cols = 10;
        std::vector<int> matrix = generateMatrix(rows, cols);
        int parallelMax = parallelFind(matrix, rows, cols);
        int sequentialMax = sequentialFind(matrix, rows, cols);
        ASSERT_EQ(parallelMax, sequentialMax);
    } else {
        parallelFind({{}}, 0, 0);
    }
}

TEST(max_value, Test_matrix100x100) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        const int rows = 100;
        const int cols = 100;
        std::vector<int> matrix = generateMatrix(rows, cols);
        int parallelMax = parallelFind(matrix, rows, cols);
        int sequentialMax = sequentialFind(matrix, rows, cols);
        ASSERT_EQ(parallelMax, sequentialMax);
    } else {
        parallelFind({{}}, 0, 0);
    }
}

TEST(max_value, Test_matrix1000x1000) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        const int rows = 100;
        const int cols = 100;
        std::vector<int> matrix = generateMatrix(rows, cols);
        int parallelMax = parallelFind(matrix, rows, cols);
        int sequentialMax = sequentialFind(matrix, rows, cols);
        ASSERT_EQ(parallelMax, sequentialMax);
    } else {
        parallelFind({{}}, 0, 0);
    }
}
*/
int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);

    ::testing::AddGlobalTestEnvironment(new GTestMPIListener::MPIEnvironment);
    ::testing::TestEventListeners& listeners =
        ::testing::UnitTest::GetInstance()->listeners();

    listeners.Release(listeners.default_result_printer());
    listeners.Release(listeners.default_xml_generator());

    listeners.Append(new GTestMPIListener::MPIMinimalistPrinter);
    return RUN_ALL_TESTS();
}
