#!/usr/bin/python
import os, sys

print 'Generating ctags...'
print 'PHYCAS_ROOT =',os.environ['PHYCAS_ROOT']

TAGS_TEMP_FILE = '/tmp/tags'
projectDir = os.environ['PHYCAS_ROOT']
if len(projectDir) == 0:
    print '  WARNING: Could not find environmental variable named PHYCAS_ROOT'
    projectDir = '/Users/plewis/Documents/Projects/pdev/trunk'
    print '  Set projectDir to',projectDir
projectName = 'Phycas'
ctagsExecutablePath = "/Users/swofford/Applications/text_editors/BBEdit_8/BBEdit.app/Contents/MacOS/ctags"
baseArgs = '--excmd=number --tag-relative=no --fields=+a+m+n+S -f /tmp/tags -R'
appendArg = '--append'
os.chdir('/')
sourceDir = os.path.join(projectDir, 'phycas', 'src')
tagsFile = os.path.join(projectDir, 'tags')

# create the project's tags in '/tmp'
if os.access(sourceDir, os.F_OK):
    buildTagsCommand = ''''%s' %s '%s' ''' % (ctagsExecutablePath, baseArgs, sourceDir)
    print "buildTagsCommand=",buildTagsCommand
    output = os.popen(buildTagsCommand).read()

# move it where it goes
os.rename(TAGS_TEMP_FILE, tagsFile)