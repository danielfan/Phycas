23-March-2008
Miscellaneous notes on creating a Mac installer for Phycas
The final installer will be automated and probably end up just being
a metapackage that does everything, but for now these notes might be
useful for figuring out how some steps along the way can be accomplished.

To hide a folder named "scripts" inside a mounted disk image named "phycas-installer":
--------------------------------------------------------------------------------------
/Developer/Tools/SetFile -a V /Volumes/phycas-installer/scripts

To set a background image for a mounted disk image named "phycas-installer":
----------------------------------------------------------------------------
With phycas-installer mounted and in the foreground, press Command-J. 
Check the "This window only" radio button
Check the "Picture..." radio button
Select the image to serve as the background.
Note: the selected image file should be placed in a hidden folder inside 
the disk image, otherwise it will not be available to users when the dmg 
is distributed.

Applescript script for opening the site-packages directory:
-----------------------------------------------------------
set posixdir to (do shell script "export PATH=/usr/local/bin:$PATH;echo `/Volumes/phycas-installer/Phycas/scripts/opensitepackages.py`")
set macfolder to POSIX file posixdir
tell application "Finder" to open folder macfolder

The do shell script command runs the opensitepackages.py python script inside
the hidden folder scripts. The opensitepackages.py file looks like this:

#!/usr/bin/env python
import sys
print sys.prefix + '/lib/python%s/site-packages' % sys.version[:3]

The opensitepackages.py file needs to be made executable.

The initial PATH modification is needed because the
do shell script command runs sh with a very limited PATH, one which doesn't in 
fact include the path where python is normally installed.

The "POSIX file posixdir" command is needed to convert a POSIX-style path
e.g. /Library/Frameworks/Python.framework/Versions/2.5/lib/python2.5/site-packages
into a Mac path
e.g. Macintosh HD:Library:Frameworks:Python.framework:Versions:2.5:lib:python2.5:site-packages

Changing the icon for something:
--------------------------------
Do a Get Info (Command-I) on the object that has an icon you want to steal.
Click on the icon in the upper-left-hand-corner and use Command-C to copy it
to the clipboard. Now Get Info for the object that you want to have this icon.
Click on the icon in the upper left corner of the Get Info dialog box of this
object and press Command-V to paste the stolen icon there.
