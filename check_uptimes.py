import subprocess, sys

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


processes = [subprocess.Popen(["ssh", server, "cat /proc/cpuinfo | grep processor | wc -l && uptime"], stdout=subprocess.PIPE, stderr=open("/dev/null", "w")) for server in servers]

total_cpu_available = 0
total_serv_available = 0
for server,process in zip(servers, processes):
    data = process.stdout.read().strip().split("\n")
    if not data or not data[0]:
        print server, "[ offline ]"
        continue
    total_serv_available += 1
    uptime = data[1].strip()
    nprocs = int(data[0].strip())
    load = round(float(data[1].split()[-1]))
    available = int(nprocs - load)
    if available < 0:
        available = 0
    total_cpu_available += available
    total_serv_available += 1
    print server, "[% 2i/% 2i processors available]" % (available, nprocs), uptime

print total_cpu_available, "processors available on", total_serv_available, "hosts"
