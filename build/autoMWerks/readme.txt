This directory will hold the metrowerks xml file that can be generated from python by the command.
	python $PHYCAS_ROOT/python/xcode_to_mwerks.py $PHYCAS_ROOT/build/xcode_proj/Phycas.xcode/project.pbxproj > $PHYCAS_ROOT/build/autoMWerks/testPhycas.mcp.xml

(or if you put the $PHYCAS_ROOT/python directory in your path, make the script executable, and navigate to the build directory:
	xcode_to_mwerks.py xcode_proj/Phycas.xcode/project.pbxproj > autoMWerks/testPhycas.mcp.xml
)