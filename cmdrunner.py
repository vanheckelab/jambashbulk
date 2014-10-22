# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 15:40:53 2013

@author: Merlijn van Deen
"""
import time, subprocess, psutil

basepath = r"C:\Users\Merlijn van Deen\Documents\GitHub\jambashbulk"
cmds = pickle.load(open(basepath + "\\cmds"))

def cleanup(procs):
    time.sleep(0.01)
    for i, (proc, logfile) in enumerate(procs):
        if proc.poll() is not None:
            del procs[i]
            logfile.write("\n\n")
            logfile.write("="*80 + "\n")
            logfile.write("OUTPUT FINISHED\n")
            logfile.write(time.strftime("%a, %d %b %Y %H:%M:%S\n"))
            logfile.close()
            break

procs = []
for (N,lP,num),cmdlines in sorted(cmds.items(), key=lambda x: x[0]):
    cmdlines = [str(s) for s in cmdlines]
    while(len(procs) > 2):
        cleanup(procs)
            
    logfile = os.path.join(basepath, "logs",  "~".join(str(x) for x in (N,lP,num)) + ".log")
    print "Starting", (N,lP,num), "using commands", cmdlines
    print "Logfile: ", logfile
    
    stdlog = open(logfile, "w")
    stdlog.write(time.strftime("%a, %d %b %Y %H:%M:%S\n"))
    stdlog.write("Shearing packing N=%i, log10P=%.1f, num=%i, using commands: \n" % (N,lP,num))
    for line in cmdlines:
        stdlog.write(line + "\n")
        
    stdlog.write("\n\nOUTPUT FOLLOWS\n" + "="*80 + "\n\n")
    stdlog.flush()
    
    proc = subprocess.Popen(
        [os.path.join(basepath, "bin\\jam2D"), "-screen"],
        stdin=subprocess.PIPE,
        stdout=stdlog,
        stderr=subprocess.STDOUT,
        cwd=basepath,
    )
    print proc.poll()
    psutil.Process(proc.pid).set_nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
    procs.append((proc, stdlog))

    for line in cmdlines:
        proc.stdin.write(line + "\n")
        
print "Waiting for processes to finish..."
while(procs):
    cleanup(procs)