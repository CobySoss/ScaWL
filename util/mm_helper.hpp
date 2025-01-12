#ifndef MM_HELPER_HPP
#define MM_HELPER_HPP

#include "sparse_representation.hpp"
#include <string>
using namespace std;
COO read_matrix_market_to_COO(const char* fname);
CSR read_matrix_market_to_CSR(const char* fname);
CSC read_matrix_market_to_CSC(const char* fname);
#endif /* !MM_HELPER_HPP */
