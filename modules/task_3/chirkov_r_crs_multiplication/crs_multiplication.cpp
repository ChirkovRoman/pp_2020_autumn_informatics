// Copyright 2020 Chirkov Roman
#include <mpi.h>
#include <algorithm>
#include <ctime>
#include <random>
#include <vector>
#include "../../../modules/task_3/chirkov_r_crs_multiplication/crs_multiplication.h"

crsMatrix::crsMatrix(const std::vector<double>& A, const int& _col_size,
    const int& t_rows) : col_size(_col_size), row_size(t_rows) {
    not_empty = 0;
    pointer.push_back(0);
    for (auto row = 0; row < t_rows; row++) {
        int non_zero = std::count_if(A.begin() + row * _col_size,
            A.begin() + row * _col_size + _col_size,
            [](double val) { return val != 0; });
        pointer.push_back(pointer.back() + non_zero);

        for (auto column = 0; column < _col_size; column++) {
            if (A.at(row * _col_size + column) != 0) {
                columns.push_back(column);
            }
        }

        not_empty += non_zero;
    }

    values.resize(not_empty);
    std::copy_if(A.begin(), A.end(), values.begin(),
        [](double val) { return val != 0; });
}

const std::vector<double> crsMatrix::makeVector() const {
    std::vector<double> reality;
    for (auto row = 0; row < row_size; row++) {
        for (auto col = 0; col < col_size; col++) {
            if (pointer.at(row + 1) - pointer.at(row) == 0) {
                reality.push_back(0.0);
            } else {
                bool pushed = false;
                for (auto i = pointer.at(row); i <= pointer.at(row + 1) - 1; i++) {
                    if (columns.at(i) == col) {
                        reality.push_back(values.at(i));
                        pushed = true;
                    }
                }
                if (!pushed) {
                    reality.push_back(0.0);
                }
            }
        }
    }

    return reality;
}

const std::vector<double> crsMatrix::makeColumn(const int& t_col) {
    std::vector<double> reality;
    for (auto row = 0; row < row_size; row++) {
        if (pointer.at(row + 1) - pointer.at(row) == 0) {
            reality.push_back(0.0);
        } else {
            bool pushed = false;
            for (auto i = pointer.at(row); i <= pointer.at(row + 1) - 1; i++) {
                if (columns.at(i) == t_col) {
                    reality.push_back(values.at(i));
                    pushed = true;
                }
            }
            if (!pushed) {
                reality.push_back(0.0);
            }
        }
    }

    return reality;
}

const std::vector<double> operator*(const crsMatrix& A,
    const crsMatrix& B) {
    if (A.col_size != B.row_size) {
        throw SIZE_ERROR;
    }

    std::vector<double> reality;
    for (auto row = 0; row < A.row_size; row++) {
        for (auto col = 0; col < B.col_size; col++) {
            reality.push_back(0.0);
            for (auto i = A.pointer.at(row); i <= A.pointer.at(row + 1) - 1;
                i++) {
                if (B.pointer.at(A.columns.at(i) + 1) -
                    B.pointer.at(A.columns.at(i)) ==
                    0) {
                    continue;
                }
                if (B.columns.at(B.pointer.at(A.columns.at(i))) == col) {
                    reality.back() +=
                        A.values.at(i) * B.values.at(B.pointer.at(A.columns.at(i)));
                }
            }
        }
    }

    return reality;
}

const std::vector<double> operator*(const crsMatrix& A,
    const std::vector<double>& B) {
    if (A.col_size != (signed)B.size()) {
        throw SIZE_ERROR;
    }

    std::vector<double> reality;
    for (auto row = 0; row < A.row_size; row++) {
        reality.push_back(0.0);

        for (auto i = A.pointer.at(row); i <= A.pointer.at(row + 1) - 1; i++) {
            reality.back() += A.values.at(i) * B.at(A.columns.at(i));
        }
    }

    return reality;
}


std::vector<double> multiply(crsMatrix* A, crsMatrix* B) {
    int rank, comm_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    if (A->row_size < comm_size) {
        if (rank == 0) {
            return *A * *B;
        } else {
            return std::vector<double>();
        }
    }

    if (comm_size == 1) {
        return *A * *B;
    }

    MPI_Bcast(&A->not_empty, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank != 0) {
        A->values.resize(A->not_empty);
        A->columns.resize(A->not_empty);
    }
    MPI_Bcast(&A->values[0], A->not_empty, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&A->columns[0], A->not_empty, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&A->pointer[0], A->row_size + 1, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Bcast(&B->not_empty, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank != 0) {
        B->values.resize(B->not_empty);
        B->columns.resize(B->not_empty);
    }
    MPI_Bcast(&B->values[0], B->not_empty, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&B->columns[0], B->not_empty, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&B->pointer[0], B->row_size + 1, MPI_INT, 0, MPI_COMM_WORLD);


    int rows_per_proc = A->row_size / comm_size;
    int rem = A->row_size % comm_size;
    int left_bound = rank * rows_per_proc;
    int right_bound = (rank + 1) * rows_per_proc;
    if (rank == comm_size - 1) {
        right_bound = A->row_size;
    }


    std::vector<double> local_reality;
    for (auto row = left_bound; row < right_bound; row++) {
        for (auto col = 0; col < B->col_size; col++) {
            local_reality.push_back(0.0);
            for (auto i = A->pointer.at(row); i <= A->pointer.at(row + 1) - 1;
                i++) {
                if (B->pointer.at(A->columns.at(i) + 1) -
                    B->pointer.at(A->columns.at(i)) ==
                    0) {
                    continue;
                }
                if (B->columns.at(B->pointer.at(A->columns.at(i))) == col) {
                    local_reality.back() +=
                        A->values.at(i) * B->values.at(B->pointer.at(A->columns.at(i)));
                }
            }
        }
    }


    std::vector<double> global_reality;
    if (rank == 0) {
        global_reality.resize(A->row_size * B->col_size);
    }

    if (rank == 0) {
        MPI_Gather(&local_reality[0], rows_per_proc * B->col_size, MPI_DOUBLE,
            &global_reality[0], rows_per_proc * B->col_size, MPI_DOUBLE, 0,
            MPI_COMM_WORLD);
    } else {
        MPI_Gather(&local_reality[0], rows_per_proc * B->col_size, MPI_DOUBLE,
            MPI_IN_PLACE, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    if (rank == comm_size - 1) {
        MPI_Send(&local_reality[local_reality.size() - (rem * B->col_size) - 1],
            rem * B->col_size + 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    if (rank == 0) {
        MPI_Recv(&global_reality[(A->row_size - rem) * B->col_size - 1],
            rem * B->col_size + 1, MPI_DOUBLE, comm_size - 1, 0, MPI_COMM_WORLD,
            MPI_STATUS_IGNORE);
    }


    return global_reality;
}

std::vector<double> generate(const int& _col_size, const int& _rows_size,
    const int& left_edge,
    const int& right_edge) {
    if (_rows_size <= 0 || _col_size <= 0) {
        throw SIZE_ERROR;
    }

    std::vector<double> result(_rows_size * _col_size);

    std::random_device r;
    std::seed_seq seed{ r() };
    std::mt19937 randomGen(seed);
    std::uniform_int_distribution<int> distr(left_edge, right_edge);

    std::generate(result.begin(), result.end(),
        [&distr, &randomGen, &left_edge]() {
            auto matrix_elem = distr(randomGen);
            return matrix_elem <= left_edge + 4 ? matrix_elem : 0;
        });

    return result;
}