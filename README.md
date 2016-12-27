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

64 bit How to compile in Linux for using with Visual Studio

TBD.




