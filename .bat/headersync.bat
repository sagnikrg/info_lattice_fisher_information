@echo off
setlocal




rem ====== USER CONFIG ======
set "local_path=C:\Users\Ritam\OneDrive\Documents\Codespaces\Codespace_Th"
rem Source header repo log path (you said this exists already)
set "SOURCE_LOG=%local_path%\.header\log\headersync.log"
rem =========================




rem Derive repo_name and repo_dir from this script’s path:
set "script_dir=%~dp0"
for %%A in ("%script_dir%\..") do (
    set "repo_dir=%%~fA"
    set "repo_name=%%~nA"
)



echo Detected repo: "%repo_name%"
echo Repo dir: "%repo_dir%"


rem Paths
set "SOURCE=%local_path%\.header\.header\.src"
set "DEST=%repo_dir%\code\.header\.src"
set "REPO=%repo_dir%"
set "MAIN_LOG=%repo_dir%\.bat\headersync.log"
set "TMP_LOG=%repo_dir%\.bat\headersync.tmp"
set "CLEAN_LOG=%repo_dir%\.bat\headersync_clean.tmp"
set "SRC_REPO_DIR=%local_path%\.header"


rem Prompt for a commit message
set /p commitMessage="Enter a commit message: "

rem Ensure destination exists
if not exist "%DEST%" mkdir "%DEST%"

rem === Run robocopy to TEMP log ===
@echo off
robocopy "%SOURCE%" "%DEST%" /MIR /FFT /R:1 /W:1 /COPY:DAT /DCOPY:T /XJ /MT:8 /TEE /LOG:"%TMP_LOG%"



rem === Clean up robocopy output (shorten paths) ===
setlocal enabledelayedexpansion
rem Build both with-trailing-\ and without, so we match all occurrences
set "src_prefix_a=%SOURCE%\"
set "src_prefix_b=%SOURCE%"
set "dst_prefix_a=%DEST%\"
set "dst_prefix_b=%DEST%"

(
  for /f "usebackq delims=" %%L in ("%TMP_LOG%") do (
    set "line=%%L"
    rem Replace DEST first (shows up in *EXTRA File lines)
    set "line=!line:%dst_prefix_a%=%repo_name%\code\.header\.src\!"
    set "line=!line:%dst_prefix_b%=%repo_name%\code\.header\.src!"
    rem Replace SOURCE (shows up in 'Source :' and 'Newer' lines)
    set "line=!line:%src_prefix_a%=.header\.src\!"
    set "line=!line:%src_prefix_b%=.header\.src!"
    echo(!line!
  )
) > "%CLEAN_LOG%"
endlocal





rem === Append run header + cleaned log + divider ===
(
  echo(
  echo ==== [%date% %time%] Header mirror run ====
  echo Repo        : %repo_name%
  echo Source      : %SOURCE%
  echo Destination : %DEST%
  echo Commit Msg  : %commitMessage%
  echo Commit Hash : %GIT_HASH%
)>>"%MAIN_LOG%"

type "%CLEAN_LOG%" >> "%MAIN_LOG%"

(
  echo(
  echo ==== End of run ====
  echo(
  echo(
  echo(
  echo(
  echo(
  echo(
  echo(
)>>"%MAIN_LOG%"





rem === Append run header + cleaned log + divider ===
(
  echo(
  echo ==== [%date% %time%] Header mirror run ====
  echo Repo        : %repo_name%
  echo Source      : %SOURCE%
  echo Destination : %DEST%
  echo Commit Msg  : %commitMessage%
  echo Commit Hash : %GIT_HASH%
)>>"%SOURCE_LOG%"

type "%CLEAN_LOG%" >> "%SOURCE_LOG%"

(
  echo(
  echo ==== End of run ====
  echo(
  echo(
  echo(
  echo(
  echo(
  echo(
  echo(
)>>"%SOURCE_LOG%"





rem Cleanup temp logs
del "%TMP_LOG%" "%CLEAN_LOG%" >nul 2>&1




rem === Git operations ===
call :git_operations "%REPO%" "[header] %commitMessage%"
set "GIT_RESULT=%ERRORLEVEL%"

if "%GIT_RESULT%"=="0" (
    if not defined GIT_HASH set "GIT_HASH=[header] %commitMessage%"
) else (
    set "GIT_HASH=GIT_ERROR_%GIT_RESULT%"
)


rem === Git operations: SOURCE (.header) repo ===
call :git_operations "%SRC_REPO_DIR%" "[header] %commitMessage%"



echo .header mirrored and logged. Commit hash: %GIT_HASH%
goto :eof

echo .header mirrored and logged.
echo   Dest repo hash: %GIT_HASH_DEST%
echo   Src  repo hash: %GIT_HASH_SRC%
goto :eof


:git_operations
rem %~1 = repo_dir, %~2 = commit_tag (e.g., [header]), %~3 = commit message
pushd "%~1"
git add -A
git diff --cached --quiet
if %errorlevel%==0 (
    echo No changes to commit in %~1.
    popd
    exit /b 0
) else (
    git commit -m "%~2-%~3"
    if errorlevel 1 (
        echo Git commit failed in %~1.
        popd
        exit /b 2
    )
    for /f "delims=" %%H in ('git rev-parse --short HEAD') do set "GIT_HASH=%%H"
    git push
    if errorlevel 1 (
        echo Git push failed in %~1.
        popd
        exit /b 3
    )
)
popd
exit /b 0
