// complex step functions

#ifndef COMPLEXIFY_H
#define COMPLEXIFY_H

#include <cmath>
#include <complex>
#include <iostream>

template<typename T>
T absvalue(const T& val)
{
  std::cout << "calling regular version" << std::endl;
  return val;
}

template <typename T>
std::complex<T> absvalue(const std::complex<T>& val)
{
  std::cout << "calling complex version" << std::endl;
  auto val2 = val;
  if (val.real() < 0)
    val2 *= -1;

  return val2;
}

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


#endif
