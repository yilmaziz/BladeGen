@echo off
rem tools\designer.cmd [--rebuild] [design.json]
rem Opens the BladeGen Designer (C++ / Qt, MDI) on Windows.
rem  * uses the conda environment BLADEGEN_ENV (default %USERPROFILE%\miniconda3\envs\bladegen;
rem    create it once with:  conda env create -f environment.yml   (or update: conda env update -f environment.yml)
rem  * builds bladegen_designer (preset windows-msvc-release) when it is missing, or always with --rebuild
setlocal
cd /d "%~dp0.."
if "%BLADEGEN_ENV%"=="" set BLADEGEN_ENV=%USERPROFILE%\miniconda3\envs\bladegen
if not exist "%BLADEGEN_ENV%\Library\lib\cmake\Qt6" (
  echo Qt 6 not found in %BLADEGEN_ENV%
  echo create / update the environment with:  conda env update -f environment.yml
  exit /b 1
)
set REBUILD=0
set DESIGN=
:args
if "%~1"=="" goto run
if "%~1"=="--rebuild" (set REBUILD=1) else (set DESIGN=%~1)
shift
goto args
:run
set EXE=%CD%\build\windows-msvc-release\bin\bladegen_designer.exe
if "%REBUILD%"=="1" goto build
if exist "%EXE%" goto start
:build
echo building the BladeGen Designer (first time: a few minutes) ...
call tools\win_dev.cmd windows-msvc-release configure || exit /b 1
call tools\win_dev.cmd windows-msvc-release build || exit /b 1
:start
set PATH=%BLADEGEN_ENV%\Library\bin;%BLADEGEN_ENV%;%PATH%
start "" "%EXE%" %DESIGN%
