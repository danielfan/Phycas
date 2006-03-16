In order to create the auto-generated java files, you need to have an XSLT processor (i.e. Xalan) and python installed.
Here is an example of how to do the translations with xalan.

From the gui\phycasGUI\swixml\ directory, run:

java org.apache.xalan.xslt.Process -IN [command desc xml file] -XSL translate\mainJavaTransform.xsl -OUT PhycasMain.java
java org.apache.xalan.xslt.Process -IN [command desc xml file] -XSL translate\mainXMLTransform.xsl -OUT ..\..\xml\PhycasMain.xml
java org.apache.xalan.xslt.Process -IN [command desc xml file] -XSL translate\commandsXMLTransform.xsl -OUT ..\..\xml\Commands.xml
python translate\commandsJavaTransform.py [command desc xml file]

If you are working on a system with shell scripting, you can cd to  gui\phycasGUI\swixml and invoke the script autotransform.sh
to perform the transformations.

In order for the python script to work, the environment variable PYTHONPATH must be set to the python directory under phycasdev.