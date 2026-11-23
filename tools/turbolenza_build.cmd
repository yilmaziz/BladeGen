@echo off
rem tools\turbolenza_build.cmd <build-dir-name> [target] [test-filter]
rem Builds one target of TURBOlenza in build\<build-dir-name> (MSVC Release, the bladegen conda
rem env) and, for test_turbolenza, runs it.  Separate build dirs let parallel work not collide.
setlocal
set BDIR=%1
if "%BDIR%"=="" set BDIR=turbolenza
set TARGET=%2
if "%TARGET%"=="" set TARGET=test_turbolenza
set FILTER=%~3
if "%BLADEGEN_ENV%"=="" set BLADEGEN_ENV=%USERPROFILE%\miniconda3\envs\bladegen
set VSWHERE=%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe
for /f "usebackq tokens=*" %%i in (`"%VSWHERE%" -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath`) do set VSINSTALL=%%i
call "%VSINSTALL%\VC\Auxiliary\Build\vcvars64.bat" >nul || exit /b 1
set PATH=%BLADEGEN_ENV%\Library\bin;%BLADEGEN_ENV%;%BLADEGEN_ENV%\Scripts;%PATH%
cd /d "%~dp0.."
if not exist build\%BDIR%\build.ninja (
  cmake -S . -B build\%BDIR% -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=cl -DBLADEGEN_BUILD_TESTS=ON -DBLADEGEN_WARNINGS_AS_ERRORS=ON -DCMAKE_PREFIX_PATH="%BLADEGEN_ENV%\Library" || exit /b 1
)
cmake --build build\%BDIR% --target %TARGET% || exit /b 1
if /i "%TARGET%"=="test_turbolenza" (
  build\%BDIR%\bin\test_turbolenza.exe %FILTER% || exit /b 1
)
exit /b 0
