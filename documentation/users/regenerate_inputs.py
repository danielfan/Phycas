from phycas import *
o = dict(locals())
c = []
for k, v in o.iteritems():
    if isinstance(v, Utilities.PhycasCommand.PhycasCommand):
        c.append((k,v))
c.sort()
for cmd in c:
    cmd[1].manual()

o = open("cmdsettings.tex", "w")
o.write("""This section lists all currently available Phycas settings. 
You will find more if you look at the \code{Phycas.\_\_init\_\_} method in the \pathname{Phycas.py} file, but be forewarned that settings found in \pathname{Phycas.py} but not listed here are experimental and not fully tested --- use at your own risk.
Please do not ask for help with undocumented settings: we will document them here when we feel they are ready to be used.
""")

for cmd in c:
    n = cmd[0].lower()
    o.write("""
\subsection{Settings used by \code{%s}}\label{subsec:%ssettings}
\input{%s}
""" % (n,n,n))
