template = """mkdir -p %(path)s && argvtolines.sh %(path)s 1 1 %(num)i 2 | bin/jam2D -screen && argvtolines.sh %(path)s 3 1 %(num)i 3 3 | bin/jam2D -screen"""
#template = """mkdir -p %(path)s && argvtolines.sh %(path)s 1 1 %(num)i 2 | bin/jam2D -screen && argvtolines.sh %(path)s 3 1 %(num)i 25 3 | bin/jam2D -screen"""
#template = """argvtolines.sh %(path)s 3 1 %(num)i 25 3 | bin/jam2D -screen"""

#pressures = ["P%se%i" % (s,i) for s in ["1468", "2154", "4642", "6813"] for i in [-7,-6,-5,-4,-3]]
#nums = ["N%i" % i for i in [16,22,32,64,128,256,512,1024]]

pressures = ["P3162e%i" % x for x in [-7, -6, -5, -4, -3]]
nums = ["N16", "N22", "N32", "N1024"]

# NO NO TRAILING /
paths = ["Packings/%(N)s~%(P)s" % {'N': N, 'P': P} for N in nums for P in pressures]

paths += ["Packings/N22~P1e-3", "Packings/N22~P1e-2"]

paths = sorted(paths, key = lambda x: (len(x.split("/")[1].split("~")[0]), x))

packrange = range(8000, 8100)

for path in paths:
    for num in packrange:
        print template % locals()
