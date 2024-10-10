#!/usr/bin/env python
"""orca_cp_writer.py: Produces an Orca CP correction input file from an xyz"""

import argparse
import sys

def orca_cp_write(args):
    """Generates an Orca compound job file that performs 5 calculations to
        compute the counterpoise corrected interaction energy"""
    readfile = args.xyz_file
    outfile = args.orca_input_file
    split_atom_no = args.split_atom_no
    memory = args.ram
    
    method_keywords = "! " + args.method_keywords + "\n"
    print("Using the method keywords:", method_keywords)
    #Check if the method includes empirical dispersion, may need additional terms adding to this list
    disp_list = ['D3ZERO', 'd3zero', 'D3BJ', 'd3bj', 'D4', 'd4']
    dispersion = any(ele in method_keywords for ele in disp_list)
    if dispersion:
        print("Empirical dispersion correction detected.")
    memory_string = "%maxcore " + memory + "\n"
    
    #Offset the split by the two compulsory lines in an xyz file
    splitter = int(split_atom_no) + 2
    
    xyz = open(readfile, "r")
    lines = xyz.readlines()
    dimer = ''.join(lines[2:])
    monomer1 = '\t\t'.join(lines[2:splitter])
    monomer2 = '\t\t'.join(lines[splitter:])
    
    #Setup geometry for monomer1 with ghost atoms
    mon1_ghost = []
    for count, line in enumerate(lines):
        if count >= splitter:
            newline=line.split(' ')
            newline[3] = ":"
            newline = ' '.join(newline)
            mon1_ghost.append(newline)
        elif count > 1:
            mon1_ghost.append(line)
    mon1_ghost = '\t\t'.join(mon1_ghost)

    #Setup geometry for monomer2 with ghost atoms
    mon2_ghost = []
    for count, line in enumerate(lines):
        if count > 1 and count < splitter:
            newline=line.split(' ')
            newline[3] = ":"
            newline = ' '.join(newline)
            mon2_ghost.append(newline)
        elif count > 1:
            mon2_ghost.append(line)
    mon2_ghost = '\t\t'.join(mon2_ghost)

    inp=open(outfile,'w+')
    
    inp.write(memory_string)
    inp.write("* xyz 0 1\n")
    inp.write(dimer)
    inp.write("*\n\n")
    
    inp.write("%Compound\n")
    if dispersion:
        inp.write("\tvariable DIMER, SP1, SP2, SP3, SP4;\n\tvariable D0, D1, D2, D3, D4;\n")
    else:
        inp.write("\tvariable DIMER, SP1, SP2, SP3, SP4;\n")
    inp.write("\tvariable IE_KCALMOL;\n\tvariable BSSE_AU;\n\tvariable BSSE_KCALMOL;\n\tvariable CP_IE;\n\n")
    
    #Calculation 1 (the dimer)
    inp.write("\t# Calculation 1: dimer calculation\n\tNew_Step\n\t\t")
    inp.write(method_keywords)
    if dispersion:
        inp.write("\tStep_End\n\tRead DIMER = SCF_ENERGY[1] ;\n\tRead D0 = VDW_CORRECTION[1] ;\n\n")
    else:
        inp.write("\tStep_End\n\tRead DIMER = SCF_ENERGY[1] ;\n\n")
    
    #Calculation 2
    inp.write("\t# Calculation 2: fragment 1 @ complex geom with fragment 1 basis\n")
    inp.write("\tNew_Step\n\t\t")
    inp.write(method_keywords)
    inp.write("\t\t* xyz 0 1\n\t\t")
    inp.write(monomer1)
    if dispersion:
        inp.write("\t\t*\n\tStep_End\n\tRead SP1 = SCF_ENERGY[2] ;\n\tRead D1 = VDW_CORRECTION[2] ;\n\n")
    else:
        inp.write("\t\t*\n\tStep_End\n\tRead SP1 = SCF_ENERGY[2] ;\n\n")
    
    #Calculation 3
    inp.write("\t# Calculation 3: fragment 1 @ complex geom with complex basis\n")
    inp.write("\tNew_Step\n\t\t")
    inp.write(method_keywords)
    inp.write("\t\t* xyz 0 1\n\t\t")
    inp.write(mon1_ghost)
    if dispersion:
        inp.write("\t\t*\n\tStep_End\n\tRead SP2 = SCF_ENERGY[3] ;\n\tRead D2 = VDW_CORRECTION[3] ;\n\n")
    else:
        inp.write("\t\t*\n\tStep_End\n\tRead SP2 = SCF_ENERGY[3] ;\n\n")
    
    #Calculation 4
    inp.write("\t# Calculation 4: fragment 2 @ complex geom with fragment 2 basis\n")
    inp.write("\tNew_Step\n\t\t")
    inp.write(method_keywords)
    inp.write("\t\t* xyz 0 1\n\t\t")
    inp.write(monomer2)
    if dispersion:
        inp.write("\t\t*\n\tStep_End\n\tRead SP3 = SCF_ENERGY[4] ;\n\tRead D3 = VDW_CORRECTION[4] ;\n\n")
    else:
        inp.write("\t\t*\n\tStep_End\n\tRead SP3 = SCF_ENERGY[4] ;\n\n")
    
    #Calculation 5
    inp.write("\t# Calculation 5: fragment 2 @ complex geom with complex basis\n")
    inp.write("\tNew_Step\n\t\t")
    inp.write(method_keywords)
    inp.write("\t\t* xyz 0 1\n\t\t")
    inp.write(mon2_ghost)
    if dispersion:
        inp.write("\t\t*\n\tStep_End\n\tRead SP4 = SCF_ENERGY[5] ;\n\tRead D4 = VDW_CORRECTION[5] ;\n\n")
    else:
        inp.write("\t\t*\n\tStep_End\n\tRead SP4 = SCF_ENERGY[5] ;\n\n")
    
    #Ask Orca to work with the variables
    if dispersion:
        inp.write("\tIE_KCALMOL = ((DIMER + D0) - (SP1 + D1 + SP3 + D3)) * 627.5096 ;\n")
        inp.write("\tBSSE_AU = (SP1 - SP2) + (SP3 - SP4) + (D1 - D2) + (D3 - D4) ;\n")
    else:
        inp.write("\tIE_KCALMOL = ((DIMER) - (SP1 + SP3)) * 627.5096 ;\n")
        inp.write("\tBSSE_AU = (SP1 - SP2) + (SP3 - SP4) ;\n")
    inp.write("\tBSSE_KCALMOL = 627.5096 * BSSE_AU ;\n")
    inp.write("\tCP_IE = IE_KCALMOL + BSSE_KCALMOL ;\nEnd\n")

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(
        description="Creates an Orca CP-correction input file from an xyz"
    )
    parser.add_argument("xyz_file", help="The XYZ file to be processed")
    parser.add_argument("orca_input_file", help="The name of the input file to be created")
    parser.add_argument("split_atom_no", help="Number of the atom within the xyz file where the first monomer block ends")
    parser.add_argument(
         "-m",
         "--method_keywords",
         help="string that contains the method and basis sets keywords for Orca, default is 'M062X D3ZERO def2-TZVPD TIGHTSCF defgrid3'",
         default="M062X D3ZERO def2-TZVPD TIGHTSCF defgrid3",
    )
    parser.add_argument(
         "-r",
         "--ram",
         help="Amount of memory (maxcore), in MB, to be requested, default is 4000'",
         default="4000",
    )

    args = parser.parse_args()
    orca_cp_write(args)

