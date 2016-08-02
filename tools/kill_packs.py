#!/usr/bin/env python
import sys

if len(sys.argv) == 2:
  okuntil = sys.argv[1]
else:
  okuntil = None

import subprocess, os

fns = subprocess.Popen("ls -s *.err | grep -e '0 ' | cut -d' ' -f2 | sort", stdout=subprocess.PIPE, shell=True).communicate()[0].strip().split('\n')

for fn in fns:
    basename = fn.split(".")[0]
    hostname = basename.split("_")[1]
    if hostname == "maris063":
        continue
    print basename + " on " + hostname,
    if not basename >okuntil:
        continue
        if raw_input("?") == "n":
            continue
    else:
        print "(automatic)"
    data = subprocess.Popen(["ssh", hostname, "ps aux | grep -v grep | grep jam2D | grep " + basename], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    try:
        pid = data[0].split()[1].strip()
        print "Parent (bash) PID is: ", pid
    except IndexError:
       print data
       print "No PID found; touching", basename + ".err"
       open(basename + ".err", 'a').write("dead?")
       continue
    data = subprocess.Popen(["ssh", hostname, "ps --ppid=%s | grep jam2D" % pid], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    try:
        pid = data[0].split()[0].strip()
        print "Child (jam2D) PID is: ", pid
    except IndexError:
       print data
       print "No PID found; touching", basename + ".err"
       open(basename + ".err", 'a').write("dead?")
       continue
    print subprocess.Popen(["ssh", hostname, "kill %(pid)s && echo %(pid)s killed" % locals()], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0]
