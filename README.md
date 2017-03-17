myMCS: a tool for finding the maximum common subgraph between two molecules and the "extra" fragment.

32/64bit --- How to compile in Linux without any IDE:

1. create a "build" directory --> mkdir build
1. goto to build directory --> cd build
2. execture cmake to generate the makefile --> cmake ../src
3. execture the command make to compile both tool and shared library --> make 

32/64bit --- How to compile in Linux for using with Ecplise CDT

1. create a "build" directory --> mkdir build
1. goto to build directory --> cd build
2. execture cmake to generate the makefile --> cmake -G"Eclipse CDT4 - Unix Makefiles" -D CMAKE_BUILD_TYPE=Debug -D_ECLIPSE_VERSION=4.3 ../src
3. open Ecplise CDT and import the project --> File, Import, General, Existing project into workspace and select the root directory of the project. If desired, check "Copy project into workspace"
4. Buid project from Ecplise; both tool and shared library will be generated.

64 bit How to compile in Linux for using with Visual Studio (xy version, yyyy year)

1. create a "build" directory --> mkdir build
1. goto to build directory --> cd build
2. execture cmake to generate the makefile -->  cmake -G "Visual Studio xx yyyy Win64" ../src
(e.g.: for visual studio 2013 x64 --> cmake -G "Visual Studio 12 2013 Win64" ../src)
3. open the solution with Visual Studio
4. Select the prefered build type from the IDE (it is set in Debug mode by default, but you can easily switch to Release)
5. IMPORTANT STEP FOR BUILDING THE LIBRARY: select the "libMCS-x86/libMCS-x64" project--> right click--> properties-->"C/C++"-->Code Generation-->Runtime Library--> select "/MT" (for release) or "/MTd" (for debug). This ensures the compiler to statically link the required libraries.
6. Buid the project; you can build automatically both the two subproject (fmcsWrap for the executable and myMCS for the .dll), or build only the one you need.


BUGS:

TODO:
. add alert in final SDF if the two MCSs are different
. for each atom in dissimilarities in common with the MCS, add tag for specifying aromatic/alifatic/"chain". <Aromatic>


