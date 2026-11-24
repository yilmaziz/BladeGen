@echo off
rem tools\turbolenza_lint.cmd [file ...]: the TURBOlenza sources under the GCC/Clang CI warnings (clang++ -fsyntax-only).
setlocal
if "%BLADEGEN_ENV%"=="" set BLADEGEN_ENV=%USERPROFILE%\miniconda3\envs\bladegen
set VSWHERE=%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe
for /f "usebackq tokens=*" %%i in (`"%VSWHERE%" -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath`) do set VSINSTALL=%%i
call "%VSINSTALL%\VC\Auxiliary\Build\vcvars64.bat" >nul || exit /b 1
"%BLADEGEN_ENV%\python.exe" "%~dp0turbolenza_lint.py" %*
