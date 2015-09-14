#
# Copyright John Reid 2015
#

"""
Code to choose whether to import debug or release build module.
"""

import sys
import logging
logger = logging.getLogger(__name__)

#
# decide whether to import debug or release stempy C++ extension
#
# only available in python debug build
_python_debug_build = hasattr(sys, "gettotalrefcount")
if _python_debug_build:
    from ._debug_build import *
    from . import _debug_build as hmm_build
else:
    from ._release_build import *
    from . import _release_build as hmm_build

logger.info('Loaded HMM C++-python interface from %s', hmm_build.__name__)
