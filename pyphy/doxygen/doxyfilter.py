import sys, re, string

# This Python file is called by Doxygen before any file is processed. This hook
# allows the use of some stylistic features not supported by Doxygen. These
# features are made Doxygen-compatible before Doxygen sees them. The following
# code comment illustrates the features currently handled by convert.py:
#
# /*-----------------------------------------------------------------------------------------------------------------------
# | Note first the first line of the comment (line above this one), which would not be recognized as the start of a
# | Doxygen comment. This line has hyphens extending out to column 120. The middle lines of the comment begin with a 
# | vertical bar (|), which is also not a Doxygen feature. Within the comment, variables to be set off in Courier font
# | can be enclosed in apostrophes like `this' (note leading apostrophe is a leaning apostrophe). Code snippets can be
# | indicated by enclosing the snippet in a block starting with "|>" and ending with "|<"; see below:
# |>
# | #include <iostream>
# |
# | int main()
# | {
# | std::cerr << "This program doesn't do a bloody thing!" << std::endl;
# | return 0;
# | }
# |<
# | These modifications make the comment much more readable and do not seem to slow down the documentation generation
# | much.
# */
#

patStart        = re.compile(r'^/\*---')
patCode         = re.compile(r'^\|>')
patMiddle       = re.compile(r'^\|[^><]')
patDataMember   = re.compile(r"`(.+?)\'")

inf = open(sys.argv[1], 'r')
lines = inf.readlines()
inside_code = False
for line in lines:
    line = patDataMember.sub("<tt>\g<1></tt>", line)

    matchStart = patStart.match(line)
    if matchStart:
        #x = line[3:]
        #string.replace(x, '-', '')
        print '/**',
        continue
    
    matchCodeBegin = patCode.match(line) and not inside_code
    if matchCodeBegin:
        print string.replace(line, '|>', '\code'),
        inside_code = True
        continue
    
    matchCodeEnd = patCode.match(line) and inside_code
    if matchCodeEnd:
        print string.replace(line, '|>', '\endcode'),
        inside_code = False
        continue
    
    matchMiddle = patMiddle.match(line)
    if matchMiddle:
        print string.lstrip(line, '|'),
        continue
    
    print line,
    
inf.close()
