// functions for arrays
#ifndef ARRAY_H
#define ARRAY_H

#include <iostream>

namespace Ticon {

// index computation (to make 1D array look like multi-d-array
// Note: make multi-D arrays column-major (to agree with Julia)

inline int ind(const int i, const int j, const int si, const int sj)
{
  return i + j*si;
}

inline int ind(const int i, const int j, const int k,
               const int si, const int sj, const int sk)
{
  return i + j*si + k*si*sj;
}


/*
 * print an array to a given IO stream.
 * prints to the current position of the stream, and an std::endl; at the end
 */
template <typename T>
void printArray(std::ostream& fout, const T* A, const int si)
{
  fout << "[";
  for (int i=0; i < si; ++i)
  {
    fout << A[i];
    if (i != (si-1))
      fout << ", ";
  }

  fout << "]" << std::endl;
} // function printArray


/*
 * 2D array printing.  Prints each row of the array to a line
 * TODO: add indent argument for each subsequent line
 */
template <typename T>
void printArray(std::ostream& fout, const T* A, const int si, const int sj)
{
  fout << "[";

  for (int i=0; i < si; ++i)
  {
    for (int j=0; j < sj; ++j)
    {
      fout << A[ind(i, j, si, sj)];
      if (j != (sj-1))
        fout << ", ";
    }

    if (i != (si-1))
      fout << std::endl;
  }

  fout << "]" << std::endl;

}

}  // namespace Ticon
#endif
