// complex step functions

#ifndef COMPLEXIFY_H
#define COMPLEXIFY_H

#include <cmath>
#include <complex>
#include <iostream>
#include "ADTestConfig.h"


namespace Ticon {


template<typename T>
T absvalue(const T& val)
{
  return val;
}

template <typename T>
std::complex<T> absvalue(const std::complex<T>& val)
{
  auto val2 = val;
  if (val.real() < 0)
    val2 *= -1;

  return val2;
}
/*
template <typename T>
bool operator<(const std::complex<T>& val1, const std::complex<T>& val2)
{
  return val1.real() < val2.real();
}
*/

template <typename T>
T max(const T& val1, const T& val2)
{
  return std::max(val1, val2);
}

template <typename T>
std::complex<T> max(const std::complex<T>&val1, const std::complex<T>& val2)
{
  if (val1.real() > val2.real())
    return val1;
  else
    return val2;
}

} // namespace ticon

#endif
