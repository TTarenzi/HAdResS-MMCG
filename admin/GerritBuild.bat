<<<<<<< HEAD
REM Requires Windows SDK and CMake to be in the path
REM Run Windows SDK SetEnv script
SetEnv /Release && ^
cmake -G "Visual Studio 10 Win64" -D GMX_DOUBLE=%GMX_DOUBLE% -D GMX_MPI=%GMX_MPI% -D GMX_OPENMP=%GMX_OPENMP% -DGMX_DEFAULT_SUFFIX=off . && ^
msbuild /m:2 /p:Configuration=MinSizeRel All_Build.vcxproj && ^
ctest -D ExperimentalTest -C MinSizeRel -V
=======
mkdir releng
cd releng
git init && git fetch git://git.gromacs.org/releng.git refs/heads/4.6.0 && git checkout -q -f FETCH_HEAD && cd .. && python -u releng/GerritBuild.py
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
