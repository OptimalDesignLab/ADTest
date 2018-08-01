#include <cstdlib>
#include <iostream>
#include "ADTestConfig.h"
#include "euler.h"

int main(int argc, char* argv[])
{
  std::cout << "hello world" << std::endl;
  std::cout << "ADTest version " << ADTest_VERSION_MAJOR << "." << ADTest_VERSION_MINOR << std::endl;

  double qL[] = {1.0, 2.0, 3.0, 7.0};
  double qR[] = {1.0, 2.0, 3.0, 7.0};
  double aux_vars[0];
  double nrm[] = {1.0, 1.0};
  double flux[4];
  Ticon::RoeSolver(qL, qR, aux_vars, nrm, flux);

  std::cout << "flux = [";
  for (int i=0; i < 4; ++i)
    std::cout << flux[i] << ", ";
  std::cout << "]" << std::endl;

  std::complex<double> qLc[] = {1.0, 2.0, 3.0, 7.0};
  std::complex<double> qRc[] = {1.0, 2.0, 3.0, 7.0};
  std::complex<double> aux_varsc[0];
  std::complex<double> fluxc[4];

  Ticon::RoeSolver(qLc, qRc, aux_varsc, nrm, fluxc);

  std::cout << "fluxc = [";
  for (int i=0; i < 4; ++i)
    std::cout << fluxc[i] << ", ";
  std::cout << "]" << std::endl;


  return 0;
}
