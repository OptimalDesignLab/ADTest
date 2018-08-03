// functions for arrays
#ifndef ARRAY_H
#define ARRAY_H

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

}  // namespace Ticon
#endif
