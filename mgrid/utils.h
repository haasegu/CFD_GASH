#pragma once

#include <algorithm>
#include <cassert>
#include <iomanip>       // setw
#include <iostream>
#include <numeric>
#include <map>
#include <utility>        // swap()
#include <vector>

///**
 //* Output of a vector of @p v.
 //*
 //* @param[in,out] s    output stream
 //* @param[in]     v    vector
 //* @return             output stream 
 //*/
//template <class T>
//std::ostream& operator<<(std::ostream &s, const std::vector<T> &v)
//{
    //for (auto it: v) { std::cout << "  " << std::setw(5) << it; }
    //return s;
//}


/**
 * Determines the permutation vector for sorting @p v ascending with operator< .
 *
 * @param[in] v    vector to permute
 * @return         index vector for permutation, 
 *                 @p v[idx[0]] denotes the smallest element of original data.
 */
template <typename T>
std::vector<int> sort_indexes(const std::vector<T> &v)
{
    // initialize original index locations
    std::vector<int> idx(v.size());
    iota(begin(idx),end(idx),0);

    // sort indexes based on comparing values in v
    sort(idx.begin(), idx.end(),
           [&v](int i1, int i2) -> bool
            { return v[i1] < v[i2]; }
        );

    return idx;
}


/**
 * Determines the permutation vector for sorting @p v descending with operator> .
 *
 * @param[in] v    vector to permute
 * @return         index vector for permutation, 
 *                 @p v[idx[0]] denotes the smallest element of original data.
 */
template <typename T>
std::vector<int> sort_indexes_desc(const std::vector<T> &v)
{
    // initialize original index locations
    std::vector<int> idx(v.size());
    iota(begin(idx),end(idx),0);

    // sort indexes based on comparing values in v
    sort(idx.begin(), idx.end(),
           [&v](int i1, int i2) -> bool
            { return v[i1] > v[i2]; }
        );

    return idx;
}

/**
 * Generates the inverse permutation vector of @p idx.
 *
 * @param[in] idx  permutation vector
 * @return         inverse permutation, 
 */
template <typename T>
std::vector<T> inverse_indexes(const std::vector<T> &idx)
{
    static_assert(std::is_integral<T>(),"Vector elements have to be integral numbers.");
    return sort_indexes(idx);
}


/**
 * Resort elements of vector @p x according to permutation @p old2new.
 *
 * @param[in]     old2new    permutation vector
 * @param[in,out] x          vector[]
 */
template <class T>
void permute(std::vector<int> const& old2new, std::vector<T> & x)
{
    assert(x.size()==old2new.size());

    auto old(x);
    for (size_t k=0; k<old2new.size(); ++k)
    {
        const size_t &k_old = old2new[k];
        x[k] = old[k_old];
    }
}

/**
 * Resort a vector @p x storing two entries per element 
 * according to permutation @p old2new.
 *
 * @param[in]     old2new    permutation vector
 * @param[in,out] x          vector[][2]
 */
template <class T>
void permute_2(std::vector<int> const& old2new, std::vector<T> & x)
{
    assert(x.size()==2*old2new.size());

    auto old(x);
    for (size_t k=0; k<old2new.size(); ++k)
    {
        const size_t &k_old = old2new[k];
        x[2*k  ] = old[2*k_old  ];
        x[2*k+1] = old[2*k_old+1];
    }
}

/**
 * Changes the entries in @p x accoring to renumbering @p  old2new.
 *
 * @param[in]     old2new    renumbering vector
 * @param[in,out] x          vector
 */
template <class T>
void reNumberEntries(std::vector<T> const& old2new, std::vector<T> & x)
{
    static_assert(std::is_integral<T>(),"Vector elements have to be integral numbers.");
    assert( *max_element(cbegin(x),cend(x)) < static_cast<int>(old2new.size()) );
    auto old(x);
    auto const n2o = inverse_indexes(old2new);
    for (size_t k=0; k<x.size(); ++k)
    {
        T const k_old = old[k];
        x[k] = n2o[k_old];        
    }
}

/**
 * Changes the entries in @p x accoring to renumbering @p  old2new.
 *
 * @param[in]     old2new    renumbering vector
 * @param[in,out] x          vector
 */
template <class T>
void reNumberEntries(std::vector<T> const& old2new, std::map<T,T> & x)
{
    static_assert(std::is_integral<T>(),"Vector elements have to be integral numbers.");
    assert( max_element(cbegin(x),cend(x))->second < static_cast<int>(old2new.size()) );
    auto old(x);
    auto const n2o = inverse_indexes(old2new);   
    for ( auto &[key,val] : x)
    {
        val = n2o[val];
    }
}

/**
 * Sort the two entries per element of vector @p x ascending. 
 *
 * @param[in] x          vector[][2]
 */
template <class T>
void sortAscending_2(std::vector<T> & x)
{
    static_assert(std::is_integral<T>(),"Vector elements have to be integral numbers.");
    //for (auto it=begin(x); it!=end(x); it += 2)    // general solution
    //{
        //sort(it,it+1);
    //}
    
    for (size_t k=0; k<x.size(); k+=2)               // for only 2 entries
    {
        if (x[k]>x[k+1]) { std::swap(x[k],x[k+1]); }
    }
}
