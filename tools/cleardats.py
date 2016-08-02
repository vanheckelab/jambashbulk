import glob, os
datfiles = glob.glob("*.dat")
errfiles = glob.glob("*.err")

all = set(x.split(".")[0] for x in datfiles)
running = set(x.split(".")[0] for x in errfiles)

print "checking .err files for irrelevant errors...",
i = 0
for f in running:
   err = open(f + ".err").readlines()
   if "ParserWarning: Falling back to the 'python' engine because the 'c' engine does not support regex separators; you can avoid this warning by specifying engine='python'" in err[0] and len(err) == 2:
      os.unlink(f+".err")
      i += 1

print i, "found."
errfiles = glob.glob("*.err")
running = set(x.split(".")[0] for x in errfiles)

removable = all - running
print "%i .dat files, of which we can remove %i (%i still running)" % (len(all), len(removable), len(running))

for f in removable:
   os.unlink(f + ".dat")


