@echo off
rem tools\gui.cmd [--rebuild] [design.json]
rem Opens the BladeGen section designer window (Windows).
rem  * uses the conda environment BLADEGEN_ENV (default %USERPROFILE%\miniconda3\envs\bladegen;
rem    create it once with:  conda env create -f environment.yml)
rem  * builds the Python module (preset 'python') when it is missing or older than this GUI,
rem    or always with --rebuild
setlocal
cd /d "%~dp0.."
if "%BLADEGEN_ENV%"=="" set BLADEGEN_ENV=%USERPROFILE%\miniconda3\envs\bladegen
if not exist "%BLADEGEN_ENV%\python.exe" (
  echo conda environment not found: %BLADEGEN_ENV%
  echo create it with:  conda env create -f environment.yml
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
set PATH=%BLADEGEN_ENV%\Library\bin;%BLADEGEN_ENV%;%BLADEGEN_ENV%\Scripts;%PATH%
set PYTHONPATH=%CD%\build\python\python;%PYTHONPATH%
if "%REBUILD%"=="1" goto build
"%BLADEGEN_ENV%\python.exe" -c "import bladegen, bladegen.gui; bladegen.Naca65Section" >nul 2>&1 && goto gui
:build
echo building the BladeGen Python module (first time: a few minutes) ...
call tools\win_dev.cmd python build || exit /b 1
:gui
"%BLADEGEN_ENV%\python.exe" -m bladegen.gui %DESIGN%
