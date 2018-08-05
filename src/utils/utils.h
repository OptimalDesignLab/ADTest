// main header file for the Utils library

#ifndef UTILS_H
#define UTILS_H

#include <chrono>
#include "array.h"
#include "complexify.h"


namespace Ticon {

/*
 * Use this clock for timing
 */
typedef std::chrono::high_resolution_clock Clock;


/*
 * Returns the elapsed time of  a given statement, returning a std::chrono::duration object
 */
#define TIME(stmt) ({\
  auto t1 = Ticon::Clock::now();\
  stmt;\
  auto t2 = Ticon::Clock::now();\
  t2 - t1;\
})

/*
 * operator<< for a std::chrono::duration object, with decent formatting
 */
template <typename R, typename T>
std::ostream& operator<<(std::ostream& fout, const std::chrono::duration<R, T>& etime)
{
  // figure out the right unit of measure, largest unit = seconds
  if (std::chrono::duration_cast<std::chrono::nanoseconds>(etime).count() < 1000)
    fout << std::chrono::duration<double, std::nano>(etime).count() << " nanoseconds";
  else if (std::chrono::duration_cast<std::chrono::microseconds>(etime).count() < 1000)
    fout << std::chrono::duration<double, std::micro>(etime).count() << " microseconds";
  else if (std::chrono::duration_cast<std::chrono::milliseconds>(etime).count() < 1000)
    fout << std::chrono::duration<double, std::milli>(etime).count() << " milliseconds";
  else //if (std::chrono::duration_cast<std::chrono::seconds>(etime).count() < 1000)
    fout << std::chrono::duration<double, std::ratio<1,1>>(etime).count() << " seconds";

  return fout;
} // function printEtime

} // namespace Ticon

#endif
