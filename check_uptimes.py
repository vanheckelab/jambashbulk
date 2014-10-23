import subprocess, sys, itertools, getpass

username = getpass.getuser()

usage_error = False
allservers = ["maris%03i" % i for i in range(2, 68)]
if len(sys.argv) == 2:
    if sys.argv[1].lower() == "default":
        servers = ["maris%03i" % i for i in range(25,33) + range(61,68)]
    elif sys.argv[1].lower() == "all":
        servers = allservers
    else:
        usage_error = True

if len(sys.argv) == 1 or usage_error:
    print "Usage: %s [DEFAULT|ALL]" % sys.argv[0]
    print "Assuming ALL"
    print ""
    servers = allservers

if "-v" in sys.argv:
  verbose = True
else:
  verbose = False

justown = False
if "--own" in sys.argv:
  justown = True
  verbose = True

thekilling = False
if "--kill" in sys.argv:
  thekilling = True
  verbose = True
  justown = True

processes = [subprocess.Popen(["ssh", server, "cat /proc/cpuinfo | grep processor | wc -l && uptime && ps -eo pcpu,pid,user,args | awk '{if ($1 > 95) print $0}'"], stdout=subprocess.PIPE, stderr=open("/dev/null", "w")) for server in servers]

total_cpu_available = 0
total_cpu = 0
total_serv_available = 0
for server,process in zip(servers, processes):
    data = process.stdout.read().strip().split("\n")
    if not data or not data[0]:
        print server, "[ offline ]"
        continue
    total_serv_available += 1
    uptime = data[1].strip()
    nprocs = int(data[0].strip())
    load = round(float(data[1].split()[-3][:-1]))
    available = int(nprocs - load)
    if available < 0:
        available = 0
    total_cpu += nprocs
    total_cpu_available += available
    total_serv_available += 1
    print server, "[% 2i/% 2i processors available]\t" % (available, nprocs),

    processes = [d.split(None, 3) for d in data[2:]]
    if not verbose:
       # summarize processes
       key = lambda x: x[2]
       users = []
       for user, procs in itertools.groupby(sorted(processes, key=key), key=key):
           totalcpu = sum(float(p[0]) for p in procs)
           users.append((user, totalcpu))
       print "\t".join("%s (%.1f%%)" % (user, totalcpu) for (user, totalcpu) in users)
    else:
       print "" # newline
       tokill = []
       for process in sorted(processes):
           if process[2] == username or not justown:
               print "\t".join(process)
               tokill.append(process)
       if thekilling and tokill:
           kill = raw_input("Kill %i processes? [y/n]" % len(tokill)).upper()
           if kill == "Y":
               subprocess.Popen(["ssh", server, "kill " + " ".join(str(pid) for pcpu,pid,user,args in tokill)], stderr=open("/dev/null", "w"))
       
           

print total_cpu_available, "/", total_cpu, "processors available on", total_serv_available, "hosts (%i%%)" % (((total_cpu_available) * 100) / total_cpu)
print total_cpu - total_cpu_available, "processors in use (%i%%)" % (((total_cpu - total_cpu_available) * 100) / total_cpu)
