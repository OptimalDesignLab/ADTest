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
T absvalue(const std::complex<T>& val)
{
  std::cout << "calling complex version" << std::endl;
  auto val2 = val;
  if (val.real() < 0)
    val2 *= -1;

  return val2;
}


#endif
