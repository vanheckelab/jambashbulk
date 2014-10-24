import os, sys, csv
sys.path.append("/home/merlijn/phd-library")
sys.path.append(os.path.join(os.path.split(__file__)[0], "src"))

sys.path.append('fakepackages')
from packing_tools import V_harm
from hdf_tools import pytables_import_shear, getcc
pytables_import_shear.insert_particles = False

from packing_tools.parser.parser import read_packings

import jamBashbulk
print jamBashbulk.__version__


import numpy as np
path,fn = os.path.split(sys.argv[1])
fn = fn.split("log")[1]

packings_reader = read_packings(os.path.join(path, 'particles' + fn))
particle_positions = [pack for pack in packings_reader]

comments, packings, log = pytables_import_shear.load_measurement(path, fn)
data = packings.to_records()

def getAlpha(pack):
    return pack["L2"][0] / pack["L"]

def getDelta(pack):
    return np.sqrt(pack["L2"][1] / pack["L1"][0]) - 1.0;

initial_contact_network = V_harm.get_contacts(particle_positions[0])["connmatrix"]
zeroAlpha = getAlpha(particle_positions[0])
zeroDelta = getDelta(particle_positions[0])
zeroL = particle_positions[0]["L"]
f1dataO = []
f1dataS = []
f2dataO = []

cols = ["step#", "cc_strained", "cc_relaxed", "gamma_strained", "gamma_relaxed", "s_xy_strained-", "s_xy_strained+", "U_strained-", "U_strained+"]
f = open(os.path.join(path, 'cc_relax_data_' + fn), 'wb')
f.write(",".join(cols) + "\n")
w = csv.DictWriter(f, cols)

for (before, after) in getcc.get_multi_ccs(None, data):
    entry = {}
    
    afterid = after["step#"]
    entry["step#"] = afterid
    
    pack = particle_positions[afterid]
    dcc = np.sum(initial_contact_network ^ V_harm.get_contacts(pack)["connmatrix"])
    entry["cc_strained"] = dcc/2
    
    jamBashbulk.import_packing(pack)
    jamBashbulk.relax_packing(True, False, False)
    relaxpack = jamBashbulk.export_packing()
    
    dcc = np.sum(initial_contact_network ^ V_harm.get_contacts(relaxpack)["connmatrix"])
    entry["cc_relaxed"] = dcc/2
    
    entry["gamma_relaxed"] = getAlpha(relaxpack) - zeroAlpha
    
    entry["gamma_strained"] = before["gamma"]
    entry["s_xy_strained-"] = before["s_xy"]
    entry["s_xy_strained+"] = after["s_xy"]
    
    entry["U_strained-"] = before["U"]
    entry["U_strained+"] = after["U"]
    
    w.writerow(entry, )
    
f.close()
