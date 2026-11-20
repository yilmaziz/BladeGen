@echo off
rem tools\win_dev.cmd <preset> [configure|build|test|all]
rem Windows developer build: finds Visual Studio with vswhere, enters its x64
rem developer environment and uses the conda environment BLADEGEN_ENV (default:
rem %USERPROFILE%\miniconda3\envs\bladegen) for OCCT, Eigen, Catch2, Ninja, CMake.
setlocal
set PRESET=%1
if "%PRESET%"=="" set PRESET=windows-msvc-release
set STEP=%2
if "%STEP%"=="" set STEP=all
if "%BLADEGEN_ENV%"=="" set BLADEGEN_ENV=%USERPROFILE%\miniconda3\envs\bladegen
set VSWHERE=%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe
for /f "usebackq tokens=*" %%i in (`"%VSWHERE%" -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath`) do set VSINSTALL=%%i
if "%VSINSTALL%"=="" (echo Visual Studio with C++ tools not found & exit /b 1)
call "%VSINSTALL%\VC\Auxiliary\Build\vcvars64.bat" >nul || exit /b 1
set PATH=%BLADEGEN_ENV%\Library\bin;%BLADEGEN_ENV%;%BLADEGEN_ENV%\Scripts;%PATH%
cd /d "%~dp0.."
if "%STEP%"=="configure" goto configure
if "%STEP%"=="all" goto configure
if "%STEP%"=="build" goto build
if "%STEP%"=="test" goto test
echo unknown step %STEP% & exit /b 2
:configure
cmake --preset %PRESET% -DCMAKE_PREFIX_PATH="%BLADEGEN_ENV%\Library" -DPython_EXECUTABLE="%BLADEGEN_ENV%\python.exe" || exit /b 1
if "%STEP%"=="configure" exit /b 0
:build
cmake --build --preset %PRESET% || exit /b 1
if "%STEP%"=="build" exit /b 0
:test
ctest --preset %PRESET% --output-on-failure || exit /b 1
exit /b 0
