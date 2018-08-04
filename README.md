# Code Organization

 * `src`: a single file to run example cases
 * `src/utils`: some utility functions, including `complexify.h` and `ind` for arrays
 * `src/euler`: currently the Roe solver, both primal and forward-mode differentiated version

## Code Facilities

Use the `TIME` macro for timing things.  Note that the macro returns an object
which can be printed to an `std::ostream`.

```c++
using namespace Ticon;  // important!
std::cout << "elapsed time = " << TIME(setCourseAndSpeed(heading, speed)) << std::endl;
```

# C++ Notes

 * Put *everything* in a namespace
 * Do not indent code inside a namespace
 * Currently we use raw pointers rather than `std::vector` or `std::array`.  This helps with Julia compatability, but should be reconsidered in the future (boundschecking would be nice).
 * use `auto` whenever possible
 * Always initialize data members of classes (especially those for which the default initialization is undefined)

# CMake Notes

prefer `target_*` functions to `*`, ie. `target_include_directories` to `include_directories`.
The idea is for each target to tell anyone who tries to link to it what they need to do (what library to link to, where the header files are found etc.)
The top level `CMakeLists.txt` should use `include` to get the `CMakeLists.txt`
of all components of the package.

Try to avoid lateral includes whenever possible.  The main `CMakeLists.txt`
should include every major component, and therefore make names available
to each other.

`CMAKE_CURRENT_LIST_DIR` gives the directory of the current `CMakeLists.txt`
file being run.  Use this to set paths whenver possible.
