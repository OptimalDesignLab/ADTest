
#include <vector>
#include "codi.hpp"

template <typename T>
void foo_codi(T* in_vec, T* out_vec)
{
  auto fac = 1/in_vec[0];
  for (int i=0; i < 4; ++i)
    out_vec[i] = fac/in_vec[i];

  // problem occurs on re-assignment of fac
  // renaming to fac2 makes the problem go away
  fac = 1/in_vec[1];
  for (int i=0; i < 4; ++i)
    out_vec[i] += fac/in_vec[i];

}

int main (int argc, char* argv[])
{

  codi::RealForwardVec<4> in_vec[4];
  codi::RealForwardVec<4> out_vec[4];

  for (int i=0; i < 4; ++i)
    in_vec[i] = i;

  foo_codi(in_vec, out_vec);

  std::vector<double> jac(16);
  for (int i=0; i < 4; ++i)
    for (int j=0; j < 4; ++j)
      jac[i + 4*j] = out_vec[i].getGradient()[j];

  return 0;
}
