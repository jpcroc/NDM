#!/usr/bin/env python3
#-*-coding:utf-8-*-

# From LAMMPS to NDM; for the moment, limited to eam/alloy potentials with one species
# TJ, 03/04/2020


import numpy as np
from scipy.interpolate import CubicSpline
import sys


# parameters for NDM
Nstep_rho_ndm = 10000 # for discretization of the embedding energy
Nstep_r_ndm = 10000 # for discretization of the electron density and pair potential
Mat = 26.9815
Z = 13
type = 'Al'
ziegler = [-0.1,-0.2]




# --------------------
# 1 - Read the LAMMPS file
# --------------------
lammps_file = sys.argv[1]
fich = open(lammps_file,'r')
fich.readline()
fich.readline()
fich.readline()
fich.readline()
linesplit = fich.readline().split()
print(linesplit)
Nrho_lammps = int(linesplit[0])
drho_lammps = float(linesplit[1])
Nr_lammps = int(linesplit[2])
dr_lammps = float(linesplit[3])
rcut_lammps = float(linesplit[4])
fich.readline()
nlines = Nrho_lammps//5 # 5 entries per line

# Embedding energy
Xrho_lammps = []
Femb_lammps = []
jj = 0
for ii in range(nlines):
    linesplit = fich.readline().split()
    for elt in linesplit:
        Xrho_lammps.append(jj*drho_lammps)
        jj += 1
        Femb_lammps.append(float(elt))

# Electron density
#fich.readline()
nlines = Nr_lammps//5
Xr_lammps = []
Rho_lammps = []
jj = 0
for ii in range(nlines):
    linesplit = fich.readline().split()
    for elt in linesplit:
        Xr_lammps.append(jj*dr_lammps)
        jj += 1
        Rho_lammps.append(float(elt))

# Pair potential: for eam.alloy files, it is given as r*V (in eV.A)
#fich.readline()
Vpot_lammps = []
jj = 0
for ii in range(nlines):
    linesplit = fich.readline().split()
    for elt in linesplit:
        if jj > 0:
            Vpot_lammps.append(float(elt)/Xr_lammps[jj])
        else:
            Vpot_lammps.append(float(elt))
        jj += 1

fich.close()

# --------------------
# 2 - Fit with splines
# --------------------

fich = open('eamtab.potin','w')
fich.write('1 ! number of types\n')
fich.write('{0:19.16f}    ! cut-off radius in Ang\n'.format(rcut_lammps))
fich.write("{0:14.7} {1:d}  'Al'\n".format(Mat,Z))
fich.write('{0:14.7f} {1:14.7f} ! raccordement Ziegler (en Ang)\n'.format(ziegler[0],ziegler[1]))
fich.write('{0:d}\n'.format(Nstep_r_ndm))
fich.write('1  ! type 1 GLUE/EAM \n')


# Embedding energy

xrho_max = Xrho_lammps[-1]

print('xrho_max = ',xrho_max)

cs_emb = CubicSpline(Xrho_lammps,Femb_lammps)
cs_emb_d = cs_emb.derivative(nu=1) # first derivative

drho_ndm = xrho_max/Nstep_rho_ndm

fich.write('{0:d}   {1:19.16f}\n'.format(Nstep_rho_ndm,drho_ndm))

for ii in range(Nstep_rho_ndm):
    rho_ndm = ii*drho_ndm
    femb = cs_emb(rho_ndm)
    # derivative
    femb_d = cs_emb_d(rho_ndm)
    fich.write('{0:34.15f} {1:34.15f} {2:34.15f}\n'.format(rho_ndm,femb,femb_d))



    
# Electron density

cs_rho = CubicSpline(Xr_lammps,Rho_lammps)
cs_rho_d = cs_rho.derivative(nu=1) # first derivative

dr_ndm = rcut_lammps/Nstep_r_ndm

fich.write('1\n')
fich.write('{0:d}   {1:19.16f}\n'.format(Nstep_r_ndm,dr_ndm))
for ii in range(Nstep_r_ndm):
    r_ndm = ii*dr_ndm
    rho = cs_rho(r_ndm)
    rho_d = cs_rho_d(r_ndm)
    fich.write('{0:34.15f} {1:34.15f} {2:34.15f}\n'.format(r_ndm,rho,rho_d))

# Pair potential

cs_pair = CubicSpline(Xr_lammps,Vpot_lammps)
cs_pair_d = cs_pair.derivative(nu=1) # first derivative

dr_ndm = rcut_lammps/Nstep_r_ndm

fich.write('1\n')
fich.write('{0:d}   {1:19.16f}\n'.format(Nstep_r_ndm,dr_ndm))
for ii in range(Nstep_r_ndm):
    r_ndm = ii*dr_ndm
    vpot = cs_pair(r_ndm)
    vpot_d = cs_pair_d(r_ndm)
    fich.write('{0:34.15f} {1:34.15f} {2:34.15f}\n'.format(r_ndm,vpot,vpot_d))
    

fich.close()
    


