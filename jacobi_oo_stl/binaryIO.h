#pragma once

#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

/*
  Matrices are in CRS format: The "dsp" (displacement) vector stores the starting index of each row. So for example
  the third row starts at dsp[2] and ends at (dsp[3] - 1). The "col" vector stores the column indices of each non-zero
  matrix entry. The "ele" vector stores the value af each non-zero matrix entry.
  The redundant "cnt" vector stores the number of elements per row and we have the relations:
  cnt[k] := dsp[k+1]-dsp[k], resp.  dsp[k+1]:=cumsum(cnt[0:k]) with dsp[0]:=0
*/


//! \brief A sparse matrix in CRS format (counter, column index, value) is read from a binary file.
//!
//!        The binary file has to store 4 Byte integers and 8 Byte doubles and contains the following data:
//!        - Number of rows
//!        - Number of non-zero elements/blocks
//!        - Number of non-zero matrix elements (= previous number * dofs per block)
//!        - [#elements per row]
//!        - [column indices]
//!        - [matrix elements]
//!
//! \param[in]   file name of binary file
//! \param[out]  cnt  number of non-zero elements per row
//! \param[out]  col  column indices of non-zero elements,  C-Style
//! \param[out]  ele  non-zero elements of matrix
//!
//!

void read_binMatrix(const std::string& file, std::vector<int> &cnt, std::vector<int> &col, std::vector<double> &ele);

//! \brief A sparse matrix in CRS format (counter, column index, value) is written to a binary file.
//!
//!        The binary file has to store 4 Byte integers and 8 Byte doubles and contains the following data:
//!        - Number of rows
//!        - Number of non-zero elements
//!        - Number of non-zero elements
//!        - [#elements per row]
//!        - [column indices]
//!        - [elements]
//!
//! \param[in]  file name of binary file
//! \param[in]  cnt  number of non-zero elements per row
//! \param[in]  col  column indices of non-zero elements,  C-Style
//! \param[in]  ele  non-zero elements of matrix
//!
//!
void write_binMatrix(const std::string& file, const std::vector<int> &cnt, const std::vector<int> &col, const std::vector<double> &ele);


//! \brief Output of an STL vector.
//!
//! \param[in,out] s    output stream
//! \param[in]     rhs  STL vector
//! \tparam        T    type of vector elements
//!
//! \return modified stream
//!
//!
template <class T>
std::ostream& operator<<(std::ostream &s, std::vector<T> const & rhs)
{
    copy (rhs.cbegin(), rhs.cend(), std::ostream_iterator<T>(s, " "));
    return s;
}


// See item 2 in  http://codeforces.com/blog/entry/15643
//! \brief Macro that prints a variable name together with its contest
//!
//! \param[in]     x  a variable
#define what_is(x) cerr << #x << " is \n" << x << endl;
