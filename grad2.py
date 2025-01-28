#!/usr/bin/python3

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Copyright 2008-2012 Peter Larsson
#
# Changes made by Ryuichi Suehiro:
# - [Describe change 1]
# - [Describe change 2]

import subprocess
import os
import sys
import math
import re

from optparse import OptionParser

def get_number_of_atoms(where):
    result = subprocess.run(["grep", "NIONS", where], stdout=subprocess.PIPE)
    return int(result.stdout.decode().split()[11])

def get_ediff(where):
    result = subprocess.run(["grep", "EDIFF", where], stdout=subprocess.PIPE)
    return float(result.stdout.decode().split()[2])

def get_selectivedynamics(where):
    result = subprocess.run(['awk', 'NR > 9 && NF > 3 {print $4, $5, $6}', where], stdout=subprocess.PIPE)
    return [[c == 'T' for c in s.split()] for s in result.stdout.decode().split('\n') ]

# Some ANSI colors
OKGREEN = '\033[92m'
FAIL = '\033[91m'
ENDC = '\033[0m'

# *****************************************
#             Main program
# *****************************************
    
parser = OptionParser()
parser.add_option("-v","--verbose",dest="verbose",help="Print extra info",action="store_true")
(options,args) = parser.parse_args()


with open(args[0], 'r') as outcar:
    outcarfile = args[0]
    outcarlines = outcar.readlines()
    
#Find max iterations
nelmax = int(subprocess.run(["grep", "NELM", outcarfile], 
                            stdout=subprocess.PIPE).stdout.decode().split()[2][0:-1])
natoms = get_number_of_atoms(outcarfile)
ediff = math.log10(float(get_ediff(outcarfile)))

# Check if selective dynamics is specified
result = subprocess.run(["grep", "selective dynamics as specified on POSCAR", outcarfile], 
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
selective_dynamics = result.returncode == 0
if selective_dynamics:
    parent_dir = os.path.dirname(outcarfile)
    poscarfile = os.path.join(parent_dir, 'POSCAR')
    if os.path.exists(poscarfile):
        sd = get_selectivedynamics(poscarfile)
    else:
        sys.stderr.write(f'{outcarfile} uses selective dynamics. But POSCAR can not be found.')
        raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), poscarfile)
else:
    sd = [[True, True, True] for i in range(natoms)]


re_energy = re.compile("free  energy")
re_iteration = re.compile("Iteration")
re_timing = re.compile("LOOP:")
re_force = re.compile("TOTAL-FORCE")
re_mag = re.compile("number of electron")
re_volume = re.compile("volume of cell")

lastenergy = 0.0
energy = 0.0
steps = 1
iterations = 0
cputime = 0.0
totaltime = 0.0
dE = 0.0
magmom = 0.0
spinpolarized = False
volume = 0.0
#average = 0.0
#maxforce = 0.0

i = 0
for line in outcarlines:
        
    if re_iteration.search(line):
        iterations = iterations + 1
        
    if re_force.search(line):
        # Calculate forces here...
        forces = []
        magnitudes = []
        for j in range(0,natoms):
            parts = outcarlines[i+j+2].split()
            x = float(parts[3]) if sd[j][0] else 0.0
            y = float(parts[4]) if sd[j][1] else 0.0
            z = float(parts[5]) if sd[j][2] else 0.0
            forces.append([x,y,z])
            magnitudes.append(math.sqrt(x*x + y*y + z*z))
            
        average = sum(magnitudes)/natoms
        maxforce = max(magnitudes)

    if re_mag.search(line):
         parts = line.split()
         if len(parts) > 5 and parts[0].strip() != "NELECT":
             spinpolarized = True
             magmom = float(parts[5])

    if re_timing.search(line):
        cputime = cputime + float(line.split()[6])/60.0
                    
    if re_volume.search(line):
        parts = line.split()
        if len(parts) > 4:
            volume = float(parts[4])

    if re_energy.search(line):
        lastenergy = energy
        energy = float(line.split()[4])
        dE = math.log10(abs(energy-lastenergy+1.0E-12))

        # Construct output string
        try:
            stepstr = str(steps).rjust(4)
            energystr = "Energy: " + ("%3.6f" % (energy)).rjust(12)

            logdestr = "Log|dE|: " + ("%1.3f" % (dE)).rjust(6)                    
            iterstr = "SCF: " + ("%3i" % (iterations))
            avgfstr="Avg|F|: " + ("%2.3f" % (average)).rjust(6)
            maxfstr="Max|F|: " + ("%2.3f" % (maxforce)).rjust(6)
            timestr="Time: " + ("%3.2fm" % (cputime)).rjust(6)
            volstr="Vol.: " + ("%3.1f" % (volume)).rjust(5)
        except NameError:
            print("Cannot understand this OUTCAR file...try to read ahead")
            continue

        #if iterations == nelmax:
        #    sys.stdout.write(FAIL)
            #print("         ^--- SCF cycle reached NELMAX. Check convergence!")

        #if (dE < ediff):
        #    sys.stdout.write(OKGREEN)

        if spinpolarized:
            # sys.stdout.write(("Step %3i  Energy: %+3.6f  Log|dE|: %+1.3f  Avg|F|: %.6f  Max|F|: %.6f  SCF: %3i  Mag: %2.2f  Time: %03.2fm") % (steps,energy,dE,average,maxforce,iterations,magmom,cputime))
            magstr="Mag: " + ("%2.2f" % (magmom)).rjust(6)
            print("%s  %s  %s  %s  %s  %s  %s  %s  %s" % (stepstr,energystr,logdestr,iterstr,avgfstr,maxfstr,volstr,magstr,timestr))
        else:
            print("%s  %s  %s  %s  %s  %s  %s  %s" % (stepstr,energystr,logdestr,iterstr,avgfstr,maxfstr,volstr,timestr))
            # sys.stdout.write(("Step %3i  Energy: %+3.6f  Log|dE|: %+1.3f  Avg|F|: %.6f  Max|F|: %.6f  SCF: %3i  Time: %03.2fm") % (steps,energy,dE,average,maxforce,iterations,cputime))

#        sys.stdout.write(ENDC)

        steps = steps + 1
        iterations = 0
        totaltime = totaltime + cputime
        cputime = 0.0

    i = i + 1
