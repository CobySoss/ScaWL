#include "mm_helper.hpp"
#include "sparse_representation.hpp"
#include <iostream>
#include <fstream>
#include <numeric>
#include <algorithm>
#include <assert.h>
#include <tuple>
#include <vector>
using namespace std;
COO read_matrix_market_to_COO(const char* fname) {
    COO mat;
    std::ifstream f(fname);

    while(f.peek() == '%') {
        f.ignore(2048, '\n');
    }

    f >> mat.nrows >> mat.ncols >> mat.nnz;

    mat.row_id = (unsigned int*) malloc (mat.nnz * sizeof(unsigned int));
    mat.col_id = (unsigned int*) malloc (mat.nnz * sizeof(unsigned int));
    mat.values = (double*) malloc (mat.nnz * sizeof(double));

    for(unsigned int i = 0; i < mat.nnz; i++) {
        unsigned int m;
        unsigned int n;
        double val;
        f >> m >> n >> val;
        
        mat.row_id[i] = --m;
        mat.col_id[i] = --n;
        mat.values[i] = val;
    }

    return mat;
}

CSR read_matrix_market_to_CSR(const char* fname) {
    COO coo_mat  = read_matrix_market_to_COO(fname);
    unsigned int* idx_arr = (unsigned int*) malloc (coo_mat.nnz * sizeof(unsigned int));

    std::iota(idx_arr, idx_arr + coo_mat.nnz, 0 );
    std::sort(idx_arr, idx_arr + coo_mat.nnz, [&coo_mat](unsigned int i, unsigned int j) {
        if(coo_mat.row_id[i] < coo_mat.row_id[j])
            return true;
        else if(coo_mat.row_id[i] > coo_mat.row_id[j])
            return false;
        else if(coo_mat.col_id[i] < coo_mat.col_id[j])
            return true;
        else
            return false;
    });

    CSR csr_mat;
    csr_mat.nnz = coo_mat.nnz;
    csr_mat.nrows = coo_mat.nrows;
    csr_mat.ncols = coo_mat.ncols;

    csr_mat.row_indx = (unsigned int*) malloc ((csr_mat.nrows + 1) * sizeof(unsigned int));
    csr_mat.col_id = (unsigned int*) malloc (csr_mat.nnz * sizeof(unsigned int));
    csr_mat.values = (double*) malloc (csr_mat.nnz * sizeof(double));

    unsigned int prev_row = 0;
    int cnt = 0;
    csr_mat.row_indx[0] = 0;

    for (unsigned int i = 0; i < csr_mat.nnz; ++i) {
        auto cur_idx  = idx_arr[i];
        auto cur_row = coo_mat.row_id[cur_idx];
        assert(prev_row <= cur_row);
        while(prev_row != cur_row) {
            csr_mat.row_indx[prev_row + 1] = csr_mat.row_indx[prev_row] + cnt;
            cnt = 0;
            prev_row++;
        }
        cnt++;

        csr_mat.col_id[i] = coo_mat.col_id[cur_idx];
        csr_mat.values[i] = coo_mat.values[cur_idx];
    }
    while(prev_row < csr_mat.nrows ) {
        csr_mat.row_indx[prev_row + 1] = csr_mat.row_indx[prev_row] + cnt;
        cnt = 0;
        prev_row++;
    }

    free(coo_mat.row_id);
    free(coo_mat.col_id);
    free(coo_mat.values);

    return csr_mat;


}

CSC read_matrix_market_to_CSC(const char* fname) {
    COO coo_mat  = read_matrix_market_to_COO(fname);
    unsigned int* idx_arr = (unsigned int*) malloc (coo_mat.nnz * sizeof(unsigned int));

    std::iota(idx_arr, idx_arr + coo_mat.nnz, 0 );
    std::sort(idx_arr, idx_arr + coo_mat.nnz, [&coo_mat](unsigned int i, unsigned int j) {
        if(coo_mat.col_id[i] < coo_mat.col_id[j])
            return true;
        else if(coo_mat.col_id[i] > coo_mat.col_id[j])
            return false;
        else if(coo_mat.row_id[i] < coo_mat.row_id[j])
            return true;
        else
            return false;
    });

    CSC csc_mat;
    csc_mat.nnz = coo_mat.nnz;
    csc_mat.nrows = coo_mat.nrows;
    csc_mat.ncols = coo_mat.ncols;

    csc_mat.col_indx = (unsigned int*) malloc ((csc_mat.ncols + 1) * sizeof(unsigned int));
    csc_mat.row_id = (unsigned int*) malloc (csc_mat.nnz * sizeof(unsigned int));
    csc_mat.values = (double*) malloc (csc_mat.nnz * sizeof(double));

    unsigned int prev_col = 0;
    int cnt = 0;
    csc_mat.col_indx[0] = 0;

    for (unsigned int i = 0; i < csc_mat.nnz; ++i) {
        auto cur_idx  = idx_arr[i];
        auto cur_col = coo_mat.col_id[cur_idx];
        assert(prev_col <= cur_col);
        while(prev_col != cur_col) {
            csc_mat.col_indx[prev_col + 1] = csc_mat.col_indx[prev_col] + cnt;
            cnt = 0;
            prev_col++;
        }
        cnt++;

        csc_mat.row_id[i] = coo_mat.row_id[cur_idx];
        csc_mat.values[i] = coo_mat.values[cur_idx];
    }
    while(prev_col < csc_mat.ncols ) {
        csc_mat.col_indx[prev_col + 1] = csc_mat.col_indx[prev_col] + cnt;
        cnt = 0;
        prev_col++;
    }

    free(coo_mat.col_id);
    free(coo_mat.row_id);
    free(coo_mat.values);

    return csc_mat;

}

