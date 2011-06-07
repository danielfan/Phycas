from _ConversionsExt import *

import os
from phycas import release_version
if not 'NO_PHYCAS_SPLASH' in os.environ:
    print 'release_version is',release_version
    if not release_version:
        print 'importing Conversions from ',os.path.abspath(__file__)
