#include <cstdlib>
#include <iostream>
#include "ADTestConfig.h"
#include "euler.h"
#include "utils.h"
#include "test.h"

using namespace Ticon;
int main(int argc, char* argv[])
{

  std::cout << "hello world" << std::endl;
  std::cout << "ADTest version " << ADTest_VERSION_MAJOR << "." << ADTest_VERSION_MINOR << std::endl;

  double qL[] = {1.0, 2.0, 3.0, 7.0};
  double qR[] = {1.0, 2.0, 3.0, 7.0};
  double aux_vars[0];
  double nrm[] = {1.0, 1.0};
  double flux[4];
  std::cout << "elapsed time = " << TIME(Ticon::RoeSolver(qL, qR, aux_vars, nrm, flux)) << std::endl;

  std::cout << "flux = ";
  Ticon::printArray(std::cout, flux, 4);
//  std::cout << "elapsed time = " << std::chrono::duration<double, std::milli>(etime).count() << std::endl;

//  std::cout << "elapsed time = " << etime << std::endl;
//  Ticon::printEtime(std::cout, etime);

/*
  std::cout << "flux = [";
  for (int i=0; i < 4; ++i)
    std::cout << flux[i] << ", ";
  std::cout << "]" << std::endl;
*/

  std::complex<double> qLc[] = {1.0, 2.0, 3.0, 7.0};
  std::complex<double> qRc[] = {1.0, 2.0, 3.0, 7.0};
  std::complex<double> aux_varsc[0];
  std::complex<double> fluxc[4];

  Ticon::RoeSolver(qLc, qRc, aux_varsc, nrm, fluxc);
  std::cout << "fluxc = ";
  Ticon::printArray(std::cout, fluxc, 4);

  /*
  std::cout << "fluxc = [";
  for (int i=0; i < 4; ++i)
    std::cout << fluxc[i] << ", ";
  std::cout << "]" << std::endl;
  */

  // differentiated version
  double flux_dotL[4*4];
  double flux_dotR[4*4];
  Ticon::RoeSolver_diff(qL, qR, aux_vars, nrm, flux_dotL, flux_dotR);
  std::cout << "flux_dotL =\n";
  Ticon::printArray(std::cout, flux_dotL, 4, 4);

  // big test
  TestAD<double, double, double> testdata(10, 100000);
  testRoeSolver(testdata);
  testRoeSolver_diff(testdata);

  return 0;
}



