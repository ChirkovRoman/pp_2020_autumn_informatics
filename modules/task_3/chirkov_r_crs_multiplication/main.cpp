// Copyright 2020 Chirkov Roman
#include <mpi.h>
#include <gtest/gtest.h>
#include <gtest-mpi-listener.hpp>
#include <iostream>
#include <vector>

#include "./crs_multiplication.h"

void test(const int& acols, const int& arows, const int& bcols) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    crsMatrix A(acols, arows), B(bcols, acols);
    double time;

    if (rank == 0) {
        std::vector<double> vectorA(generate(acols, arows, 1, 50));
        std::vector<double> vectorB(generate(bcols, acols, 1, 50));

        A = crsMatrix(vectorA, acols, arows);
        B = crsMatrix(vectorB, bcols, bcols);

        vectorA.clear();
        vectorB.clear();

        time = MPI_Wtime();
    }
	
	std::vector<double> reality;
    reality = multiply(&A, &B);
	
    if (rank == 0) {
        time = MPI_Wtime() - time;
        std::cout << "Parallel time: " << time << " s" << std::endl;

        time = MPI_Wtime();
        std::vector<double> expectation = A * B;
        time = MPI_Wtime() - time;
        std::cout << "Seq time: " << time << " s" << std::endl;

        ASSERT_EQ(reality, expectation);
    }
}

TEST(multiplication_tests, sequential) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        crsMatrix A(std::vector<double>{1.0, 2.0, 3.0},
            std::vector<int>{2, 1, 2}, std::vector<int>{0, 1, 3, 3}, 3,
            3),
            B(std::vector<double>{5, 3}, std::vector<int>{0, 1},
                std::vector<int>{0, 1, 2, 2}, 2, 3);

        std::vector<double> reality = A * B;
        std::vector<double> expectation{ 0, 0, 0, 6, 0, 0 };

        ASSERT_EQ(reality, expectation);
    }
}

TEST(multiplication_tests, parallel_10x10_10) { test(20, 10, 20); }
TEST(multiplication_tests, parallel_100x50_100) { test(100, 50, 100); }
TEST(multiplication_tests, parallel_300x20_300) { test(300, 20, 300); }
TEST(multiplication_tests, parallel_10x300_10) { test(10, 300, 10); }

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