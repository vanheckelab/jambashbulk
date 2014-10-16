import glob, os
datfiles = glob.glob("*.dat")
errfiles = glob.glob("*.err")

all = set(x.split(".")[0] for x in datfiles)
running = set(x.split(".")[0] for x in errfiles)

removable = all - running
print "%i .dat files, of which we can remove %i (%i still running)" % (len(all), len(removable), len(running))

for f in removable:
   os.unlink(f + ".dat")
