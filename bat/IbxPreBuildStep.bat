@echo
cd ../../bat
Echo Current dir: "%CD%"

echo Generating VS include paths
start python GenerateVSIncludeProps.py
echo Done

cd ..\
Echo Current dir: "%CD%"
rem del /s *.suo
rem del /s *.sdf
rem del /s *vc*.pdb
rem del /s *.ipch
rem del /s *.user
rem del /s *.ncb
REM del /s *.sbr
REM del /s *.log
cd ./proj
Echo Current dir: "%CD%"
rem del /s *.pdb
rem del /s *.exe
rem del /s *.obj
rem del /s *.cu.cache
rem del /s *.cu.obj
rem del /s *.tlog
rem del /s *.pch
rem del /s *.lastbuildstate
rem del /s *.log
rem del /s *.cu.deps

cd ../bat


EXIT /B 0