#!/home/suncc/anaconda3/envs/quacc/bin/python
# -*- coding: utf-8 -*-
"""
Gaussian external xtb script, From ADCR Program.

Usage: 
add the following line to the Gaussian input file:
external='Gaussian_xtb.py'

"""
import os
import sys

import numpy as np

Method = "GFN2-xTB"
chemical_symbols = [
    "X",
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
    "Xe",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
    "Lr",
    "Rf",
    "Db",
    "Sg",
    "Bh",
    "Hs",
    "Mt",
    "Ds",
    "Rg",
    "Cn",
    "Nh",
    "Fl",
    "Mc",
    "Lv",
    "Ts",
    "Og",
]

# process command line
SplitTag = 1
while sys.argv[SplitTag] != "R":
    SplitTag += 1
InputFileName, OutputFileName, MsgFileName = sys.argv[SplitTag + 1 : SplitTag + 4]

# build work environment
Label = InputFileName[InputFileName.find("Gau-") + 4 : InputFileName.find(".E")]
WorkDir = f"""{InputFileName[0:InputFileName.find("Gau-") - 1]}/xtb{Label}"""
if os.path.exists(WorkDir):
    Command = f"""rm {WorkDir}"""
    ierr = os.system(Command)
os.mkdir(WorkDir)

# read input information
with open(InputFileName, "r") as InputFile:
    Text = InputFile.readlines()
splitline = (" ".join(Text[0].strip("\n").split())).split()
NAtoms = int(splitline[0])
Derives = int(splitline[1])
Charge = int(splitline[2])
Spin = int(splitline[3])
Symbol = np.zeros(NAtoms, dtype="<U3")
X = np.zeros(NAtoms, dtype="d")
Y = np.zeros(NAtoms, dtype="d")
Z = np.zeros(NAtoms, dtype="d")
EmbedCharge = np.zeros(NAtoms, dtype="d")
Tag = np.zeros(NAtoms, dtype=bool)
ONIOM = False
RealNAtoms = 0
for i in range(NAtoms):
    splitline = (" ".join(Text[i + 1].strip("\n").split())).split()
    Symbol[i] = chemical_symbols[int(splitline[0])]
    X[i] = float(splitline[1])
    Y[i] = float(splitline[2])
    Z[i] = float(splitline[3])
    EmbedCharge[i] = float(splitline[4])
    if len(splitline) > 5 and splitline[5].find("#") != -1:
        Tag[i] = True
        if not ONIOM and EmbedCharge[i] != 0.0:
            ONIOM = True
    else:
        RealNAtoms += 1

# write xtb input file
with open(f"""{WorkDir}/input.coord""", "w") as CoordFile:
    CoordFile.write("""$coord\n""")
    for i in range(NAtoms):
        if not Tag[i]:
            CoordFile.write(
                f""" {X[i]:20.14f} {Y[i]:20.14f} {Z[i]:20.14f}    {Symbol[i]}\n"""
            )
    CoordFile.write("""$end\n""")
if ONIOM:
    with open(f"""{WorkDir}/embedcharge""", "w") as ChargeFile:
        ChargeFile.write(f"""{NAtoms - RealNAtoms}""")
        for i in range(NAtoms):
            if Tag[i]:
                ChargeFile.write(
                    f"""{EmbedCharge[i]:13.8f} {X[i]:20.14f} {Y[i]:20.14f} """
                    f"""{Z[i]:20.14f}    {Symbol[i]}\n"""
                )
    with open(f"""{WorkDir}/embed.input""", "w") as InputFile:
        InputFile.write("""$embedding\n    input=embedcharge\n$end\n""")

Command = f"cd {WorkDir} && xtb input.coord --chrg {Charge} --uhf {Spin - 1} --gfn 2 --acc 0.01 --parallel 1 --alpb Acetone"
if Derives > 1:
    Command = f"""{Command} --hess"""
if Derives > 0:
    Command = f"""{Command} --grad"""
if ONIOM:
    Command = f"""{Command} --input embed.input"""
Command = f"""{Command} > xtbout 2>&1"""
ierr = os.system(Command)
if ierr != 0:
    exit(1)

# process xtb output files
with open(f"""{WorkDir}/xtbout""", "r") as InputFile:
    Text = InputFile.readlines()
Log = open(MsgFileName, "w")
HOMO = LUMO = 1.0e8
for i in range(len(Text) - 1, -1, -1):
    if Text[i].find("TOTAL ENERGY") != -1:
        splitline = (" ".join(Text[i].strip("\n").split())).split()
        Energy = float(splitline[3])
        Log.write(
            f"""SCF Done:  E({Method}) = {Energy:16.9f}     A.U. after 0 cycles\n"""
        )
        Log.write("""           Population analysis using the SCF density.\n""")
    if Text[i].find("(LUMO)") != -1:
        splitline = (" ".join(Text[i].strip("\n").split())).split()
        LUMO = float(splitline[-3])
    if Text[i].find("(HOMO)") != -1:
        splitline = (" ".join(Text[i].strip("\n").split())).split()
        HOMO = float(splitline[-3])
        if LUMO == 1.0e8:
            LUMO = HOMO
        for j in range(Spin):
            if j % 5 == 0:
                Log.write("""\nAlpha  occ. eigenvalues -- """)
            Log.write(f"""{HOMO:10.5f}""")
        Log.write(f"""\nAlpha virt. eigenvalues -- {LUMO:10.5f}\n""")
        if Spin > 1:
            Log.write(f""" Beta  occ. eigenvalues -- {HOMO:10.5f}\n""")
            Log.write(f""" Beta virt. eigenvalues -- {LUMO:10.5f}\n""")
        break
with open(f"""{WorkDir}/charges""", "r") as ChgFile:
    Text = ChgFile.readlines()
Log.write("""Mulliken charges:\n""")
Log.write("""               1\n""")
for i in range(len(Text)):
    Text[i] = Text[i].rstrip("\n").strip()
    Log.write(f"""{i + 1:5d}  {Symbol[i]:<3s} {float(Text[i]):9.6f}\n""")
Log.close()

with open(OutputFileName, "w") as OutputFile:
    OutputFile.write(f"""{Energy:20.12f}{0.0:20.12f}{0.0:20.12f}{0.0:20.12f}\n""")

if Derives > 0:
    with open(f"""{WorkDir}/gradient""", "r") as GradientFile:
        Text = GradientFile.readlines()
    with open(OutputFileName, "a") as OutputFile:
        j = 0
        for i in range(NAtoms):
            if Tag[i]:
                OutputFile.write(f"""{0.0:20.12f}{0.0:20.12f}{0.0:20.12f}\n""")
            else:
                splitline = (
                    " ".join(Text[j + RealNAtoms + 2].strip("\n").split())
                ).split()
                OutputFile.write(
                    f"""{float(splitline[0]):20.12f}"""
                    f"""{float(splitline[1]):20.12f}{float(splitline[2]):20.12f}\n"""
                )
                j += 1

if Derives > 1:
    Hessian = np.zeros(9 * NAtoms**2, dtype="d")
    RealHessian = np.zeros(9 * RealNAtoms**2, dtype="d")
    with open(f"""{WorkDir}/hessian""", "r") as HessianFile:
        Text = HessianFile.readlines()
    cnt = 0
    for line in Text[1:]:
        splitline = (" ".join(line.strip("\n").split())).split()
        for buf in splitline:
            RealHessian[cnt] = float(buf)
            cnt += 1
    k = l = cnt = 0
    for i in range(3 * NAtoms):
        for j in range(3 * NAtoms):
            if not Tag[i // 3] and not Tag[j // 3]:
                Hessian[cnt] = RealHessian[k * 3 * RealNAtoms + l]
                l += 1
                if l == 3 * RealNAtoms:
                    k += 1
                    l = 0
            cnt += 1
    with open(OutputFileName, "a") as OutputFile:
        # polarizability
        OutputFile.write(f"""{0.0:20.12f}{0.0:20.12f}{0.0:20.12f}\n""")
        OutputFile.write(f"""{0.0:20.12f}{0.0:20.12f}{0.0:20.12f}\n""")
        # dipole derivatives
        for i in range(3 * NAtoms):
            OutputFile.write(f"""{0.0:20.12f}{0.0:20.12f}{0.0:20.12f}\n""")
        # hessian
        cnt = 0
        for i in range(3 * NAtoms):
            for j in range(i + 1):
                OutputFile.write(f"""{Hessian[i * 3 * NAtoms + j]:20.12f}""")
                cnt += 1
                if cnt == 3:
                    cnt = 0
                    OutputFile.write("""\n""")
        OutputFile.write("""\n""")

Command = f"""rm -r {WorkDir}"""
ierr = os.system(Command)
