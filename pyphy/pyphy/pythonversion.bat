@echo off
cmd /c python -c "import sys,re; print ''.join(re.match(r'(\d+)\.(\d+)', sys.version).groups())"  2> NUL
if errorlevel 1 goto bad
goto end
:bad
echo nopython
:end
