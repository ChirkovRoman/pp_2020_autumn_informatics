// Copyright 2020 Chirkov Roman
#ifndef MODULES_TASK_3_CHIRKOV_R_CRS_MULTIPLICATION_CRS_MULTIPLICATION_
#define MODULES_TASK_3_CHIRKOV_R_CRS_MULTIPLICATION_CRS_MULTIPLICATION_

#include <vector>

#define SIZE_ERROR -2 //

class crsMatrix {
 public:
	int col_size, row_size, not_empty;
    std::vector<double> values;
    std::vector<int> columns;
    std::vector<int> pointer;
    
    crsMatrix() : col_size(0), row_size(0), not_empty(0) {}

    crsMatrix(const std::vector<double>& _values,
        const std::vector<int>& _columns,
        const std::vector<int>& _pointer, const int& _col_size,
        const int& _row_size)
        : values(values),
        columns(_columns),
        pointer(_pointer),
        col_size(_col_size),
        row_size(_row_size),
        not_empty(values.size()) {}

    crsMatrix(const int& _col_size, const int& t_rows)
        : col_size(_col_size), row_size(t_rows) {
        pointer.resize(t_rows + 1);
    }

    crsMatrix(const std::vector<double>& A, const int& _col_size,
        const int& t_rows);

    const std::vector<double> makeVector() const;
    const std::vector<double> makeColumn(const int& t_col);

    friend const std::vector<double> operator*(const crsMatrix& A,
        const crsMatrix& B);
    friend const std::vector<double> operator*(const crsMatrix& A,
        const std::vector<double>& B);
};


std::vector<double> multiply(crsMatrix* A, crsMatrix* B);
std::vector<double> generate(const int& _col_size, const int& _rows_size,
    const int& left_edge,
    const int& right_edge);


#endif  // MODULES_TASK_3_CHIRKOV_R_CRS_MULTIPLICATION_CRS_MULTIPLICATION_
