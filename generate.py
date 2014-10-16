# Normal shear
#template = """mkdir -p %(path)s && argvtolines.sh %(path)s 1 1 %(num)i 2 | bin/jam2D -screen && argvtolines.sh %(path)s 3 1 %(num)i 3 3 | bin/jam2D -screen"""

# Large shear (25 ccs)
template = """mkdir -p %(path)s && argvtolines.sh %(path)s 1 1 %(num)i 2 | bin/jam2D -screen && argvtolines.sh %(path)s 3 1 %(num)i 25 3 | bin/jam2D -screen"""

# Just make 
#template = """argvtolines.sh %(path)s 3 1 %(num)i 25 3 | bin/jam2D -screen"""

#pressures = ["P%se%i" % (s,i) for s in ["1468", "2154", "4642", "6813"] for i in [-7,-6,-5,-4,-3]]
#nums = ["N%i" % i for i in [16,22,32,64,128,256,512,1024]]

pressures = ["P%02i0e-1" % i for i in [1,3,5,7,9,13,15,17,19]]
#pressures = ["P3162e%i" % x for x in [-7, -6, -5, -4, -3]]
nums = ["N128"]

# NO NO TRAILING /
paths = ["Packings/%(N)s~%(P)s" % {'N': N, 'P': P} for N in nums for P in pressures]

paths = sorted(paths, key = lambda x: (len(x.split("/")[1].split("~")[0]), x))
packrange = range(8000, 8100)

todo = []

for path in paths:
    for num in packrange:
        todo.append(template % locals())

import sys
if '--shuffle' in sys.argv:
    import random
    random.shuffle(todo)

for t in todo:
    print t
