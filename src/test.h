#ifndef TEST_H
#define TEST_H

#include <vector>
#include <iostream>

#include "utils.h"
#include "euler.h"

namespace Ticon {


template<typename Tsol, typename Tmsh, typename Tres>
class TestAD
{
  public:
    TestAD() = default;
    TestAD(int _numNodesPerElement, int _numEl);

    // data members
    int numDofPerNode = 0;
    int numNodesPerElement = 0;
    int numEl = 0;

    std::vector<Tsol> qL;
    std::vector<Tsol> qR;
    std::vector<Tmsh> nrm;
    std::vector<Tres> flux;
    std::vector<Tres> flux_jacL;
    std::vector<Tres> flux_jacR;
    Tres aux_vars[0];  // unused


}; // class TestAD


// constructor
template <typename Tsol, typename Tmsh, typename Tres>
TestAD<Tsol, Tmsh, Tres>::TestAD(int _numNodesPerElement, int _numEl)
: numDofPerNode(4), numNodesPerElement(_numNodesPerElement), numEl(_numEl)
{
  int len = numDofPerNode * numNodesPerElement * numEl;
  qL = std::vector<Tsol>(len);
  qR = std::vector<Tsol>(len);
  nrm = std::vector<Tmsh>(2 * numNodesPerElement * numEl);

  // zero initialize the output arrays
  flux = std::vector<Tres>(len, 0);
  flux_jacL = std::vector<Tres>(4*len, 0);
  flux_jacR = std::vector<Tres>(4*len, 0);

  // write values to the arrays
  for (int i=0; i < numEl; ++i)
    for (int j=0; j < numNodesPerElement; ++j)
    {
      // set qL and qR to slightly different states
      int index = ind(0, j, i, 4, numNodesPerElement, numEl);
      qL[index] = 1.0; qL[index+1] = 2.0; qL[index+2] = 3.0; qL[index+3] = 7.0;
      qR[index] = 1.1; qR[index+1] = 2.1; qR[index+2] = 3.1; qR[index+3] = 7.1;

      // set nrm
      index = ind(0, j, i, 2, numNodesPerElement, numEl);
      nrm[index]   = 1.0;
      nrm[index+1] = 1.0;
    }

} // constructor


// regular Roe Solver benchmark
template <typename Tsol, typename Tmsh, typename Tres>
void testRoeSolver(TestAD<Tsol, Tmsh, Tres>& testdata, std::ostream& fout = std::cout);

// hand differentiated Roe Solver benchmark
template <typename Tsol, typename Tmsh, typename Tres>
void testRoeSolver_diff(TestAD<Tsol, Tmsh, Tres>& testdata, std::ostream& fout = std::cout);

// CoDiPack differentiated Roe Solver becnhmark
template <typename Tsol, typename Tmsh, typename Tres>
void testRoeSolver_diff_codi(TestAD<Tsol, Tmsh, Tres>& testdata, std::ostream& fout= std::cout);

// Adept differentiated Roe Solver becnhmark
template <typename Tsol, typename Tmsh, typename Tres>
void testRoeSolver_diff_adept(TestAD<Tsol, Tmsh, Tres>& testdata, std::ostream& fout= std::cout);


//-----------------------------------------------------------------------------
// defintions

template <typename Tsol, typename Tmsh, typename Tres>
void testRoeSolver(TestAD<Tsol, Tmsh, Tres>& testdata, std::ostream& fout /* = std::cout */)
{
  auto aux_vars = testdata.aux_vars;

  // call the flux function a bunch of times
  auto t1 = Clock::now();
  for (int i=0; i < testdata.numEl; ++i)
    for (int j=0; j < testdata.numNodesPerElement; ++j)
    {
      int indexq = 0;
      int indexnrm = 0;

      auto qL = &(testdata.qL[indexq]);
      auto qR = &(testdata.qR[indexq]);
      auto nrm = &(testdata.nrm[indexnrm]);
      auto flux = &(testdata.flux[indexq]);

      RoeSolver(qL, qR, aux_vars, nrm, flux);

      indexq += testdata.numDofPerNode;
      indexnrm += 2;
    }

  auto t2 = Clock::now();

  fout << "Roe solver elapsed time = " << t2 - t1 << std::endl;  
}  // function testRoeSolver




template <typename Tsol, typename Tmsh, typename Tres>
void testRoeSolver_diff(TestAD<Tsol, Tmsh, Tres>& testdata, std::ostream& fout /* = std::cout */)
{
  auto aux_vars = testdata.aux_vars;

  // call the flux function a bunch of times
  auto t1 = Clock::now();
  for (int i=0; i < testdata.numEl; ++i)
    for (int j=0; j < testdata.numNodesPerElement; ++j)
    {
      int indexq = 0;
      int indexnrm = 0;
      int indexfluxjac = 0;


      auto qL = &(testdata.qL[indexq]);
      auto qR = &(testdata.qR[indexq]);
      auto nrm = &(testdata.nrm[indexnrm]);
      auto fluxjacL = &(testdata.flux_jacL[indexfluxjac]);
      auto fluxjacR = &(testdata.flux_jacR[indexfluxjac]);

      RoeSolver_diff(qL, qR, aux_vars, nrm, fluxjacL, fluxjacR);

      indexq += testdata.numDofPerNode;
      indexfluxjac += testdata.numDofPerNode*testdata.numDofPerNode;
      indexnrm += 2;
    }

  auto t2 = Clock::now();

  fout << "Roe solver diff elapsed time = " << t2 - t1 << std::endl;  
  fout << "first jacL = ";
  printArray(std::cout, &(testdata.flux_jacL[0]), 4, 4);
}  // function testRoeSolver


#ifdef ENABLE_CODIPACK

template <typename Tsol, typename Tmsh, typename Tres>
void testRoeSolver_diff_codi(TestAD<Tsol, Tmsh, Tres>& testdata, std::ostream& fout /* = std::cout */)
{
  auto aux_vars = testdata.aux_vars;

  // call the flux function a bunch of times
  auto t1 = Clock::now();
  for (int i=0; i < testdata.numEl; ++i)
    for (int j=0; j < testdata.numNodesPerElement; ++j)
    {
      int indexq = 0;
      int indexnrm = 0;
      int indexfluxjac = 0;


      auto qL = &(testdata.qL[indexq]);
      auto qR = &(testdata.qR[indexq]);
      auto nrm = &(testdata.nrm[indexnrm]);
      auto fluxjacL = &(testdata.flux_jacL[indexfluxjac]);
      auto fluxjacR = &(testdata.flux_jacR[indexfluxjac]);

      RoeSolver_diff_codi(qL, qR, aux_vars, nrm, fluxjacL, fluxjacR);

      indexq += testdata.numDofPerNode;
      indexfluxjac += testdata.numDofPerNode*testdata.numDofPerNode;
      indexnrm += 2;
    }

  auto t2 = Clock::now();

  fout << "Roe solver CoDiPack diff elapsed time = " << t2 - t1 << std::endl;  
}  // function testRoeSolver

#endif // end if ENABLE_CODIPACK

#ifdef ENABLE_ADEPT

#include "adept.h"

template <typename Tsol, typename Tmsh, typename Tres>
void testRoeSolver_diff_adept(TestAD<Tsol, Tmsh, Tres>& testdata, std::ostream& fout /* = std::cout */)
{
//  auto aux_vars = testdata.aux_vars;

  // call the flux function a bunch of times

  adept::Stack stack;
  stack.preallocate_statements(101);
  stack.preallocate_operations(235);
  auto t1 = Clock::now();
  for (int i=0; i < testdata.numEl; ++i)
  {
    for (int j=0; j < testdata.numNodesPerElement; ++j)
    {
      int indexq = 0;
      int indexnrm = 0;
      int indexfluxjac = 0;


      auto qL = &(testdata.qL[indexq]);
      auto qR = &(testdata.qR[indexq]);
      auto nrm = &(testdata.nrm[indexnrm]);
      auto fluxjacL = &(testdata.flux_jacL[indexfluxjac]);
      auto fluxjacR = &(testdata.flux_jacR[indexfluxjac]);

      stack.new_recording();
      RoeSolver_diff_adept(stack, qL, qR, nrm, fluxjacL, fluxjacR);

      indexq += testdata.numDofPerNode;
      indexfluxjac += testdata.numDofPerNode*testdata.numDofPerNode;
      indexnrm += 2;
    }
  }

  auto t2 = Clock::now();

  fout << "Roe solver Adept diff elapsed time = " << t2 - t1 << std::endl;
  fout << "first jacL = ";
  printArray(std::cout, &(testdata.flux_jacL[0]), 4, 4);
  std::cout << "Stack memory usage: " << stack.memory() << " bytes" << std::endl;
  std::cout << "Stack num statements: " << stack.n_statements() << std::endl;
  std::cout << "Stack num operatations: " << stack.n_operations() << std::endl;
  std::cout << "ADEPT_MULITPASS_SIZE = " << ADEPT_MULTIPASS_SIZE << std::endl;
}  // function testRoeSolver
#endif


} // namespace Ticon
#endif
