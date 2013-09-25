#
# Copyright John Reid 2007
#

"""
Test FPU status
"""


import hmm; fpu = hmm.fpu
#import hmm._fpu as fpu

print fpu.__doc__
print fpu.summarise_status()
