@echo off
rem tools\jupyter.cmd [--rebuild] [--port N]
rem Opens the BladeGen notebooks in JupyterLab (Windows).
rem  * uses the conda environment BLADEGEN_ENV (default %USERPROFILE%\miniconda3\envs\bladegen;
rem    create it once with:  conda env create -f environment.yml)
rem  * builds the Python module (preset 'python') when it is missing, or always with --rebuild
rem  * starts JupyterLab in python\notebooks with the module on PYTHONPATH
setlocal
cd /d "%~dp0.."
if "%BLADEGEN_ENV%"=="" set BLADEGEN_ENV=%USERPROFILE%\miniconda3\envs\bladegen
if not exist "%BLADEGEN_ENV%\python.exe" (
  echo conda environment not found: %BLADEGEN_ENV%
  echo create it with:  conda env create -f environment.yml
  exit /b 1
)
set PORT=8888
set REBUILD=0
:args
if "%~1"=="" goto run
if "%~1"=="--rebuild" set REBUILD=1
if "%~1"=="--port" (set PORT=%~2& shift)
shift
goto args
:run
set MOD=build\python\python\bladegen
if "%REBUILD%"=="1" goto build
if exist "%MOD%\_core*.pyd" goto lab
:build
echo building the BladeGen Python module (first time: a few minutes) ...
call tools\win_dev.cmd python build || exit /b 1
:lab
set PATH=%BLADEGEN_ENV%\Library\bin;%BLADEGEN_ENV%;%BLADEGEN_ENV%\Scripts;%PATH%
set PYTHONPATH=%CD%\build\python\python;%PYTHONPATH%
echo opening JupyterLab on http://localhost:%PORT%  (Ctrl+C in this window stops it)
"%BLADEGEN_ENV%\python.exe" -m jupyter lab --notebook-dir=python\notebooks --port=%PORT%
