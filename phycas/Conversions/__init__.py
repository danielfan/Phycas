from _Conversions import *

import os
from phycas import release_version
if not release_version:
    print 'importing Conversions from ',os.path.abspath(__file__)
