import os, fnmatch,re, shutil,sys

propsFileDir = "..\\proj\\PropertySheets\\"
propsFileName = "IncludeFolders.props"

propsFileContent = """<?xml version="1.0" encoding="utf-8"?>
<Project ToolsVersion="4.0" xmlns="http://schemas.microsoft.com/developer/msbuild/2003">
  <ImportGroup Label="PropertySheets" />
  <PropertyGroup Label="UserMacros" />
  <PropertyGroup />
  <ItemDefinitionGroup>
    <ClCompile>
      <AdditionalIncludeDirectories>~~~REPLACE ME~~~%(AdditionalIncludeDirectories)</AdditionalIncludeDirectories>
    </ClCompile>
  </ItemDefinitionGroup>
  <ItemGroup />
</Project>
"""

def GenerateIncludeList(directory):
    includeList = ["$(SolutionDir)..\\src\\"]
    for root, dirs, files in os.walk(directory):
        for basename in dirs:
            include = os.path.join("$(SolutionDir)" + root, basename + "\\")
            if "__pycache__" in include: continue
            includeList.append(include)
    return includeList

def GenerateIncludeListStr(directory):
    includeList = GenerateIncludeList(directory)
    result = ""
    for include in includeList:
        print("IncludeDir: " + include)
        result += include + ";"
    return result

def GeneratePropsFile(targetName, includeList):
    target = open(targetName, 'w')
    target.write(propsFileContent.replace("~~~REPLACE ME~~~", includeList))
    target.close()
    print("\nSaved to: " + propsFileDir + propsFileName + "\n")

os.chdir("..\\sln\\")
includeList = GenerateIncludeListStr("..\\src\\")

os.chdir(propsFileDir)
GeneratePropsFile(propsFileName, includeList)

sys.exit(0)
