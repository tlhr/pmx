import sys, os, time, math
from copy import deepcopy
from pmx import *
from PositionRestraints import *
from pmx.parser import *
from pmx.options import *
from pmx import library as pmxlib

def main(argv):

    version = "1.0"

    options = [
        Option( "-num", "int", -1, "Position restraints by atom number"),
        Option( "-resnum", "int", -1, "Can also provide residue number"),
        Option( "-name", "string", '', "Position restraints by atom name"),
        Option( "-fcx", "float", 1000, "force constant x"),
        Option( "-fcy", "float", 1000, "force constant y"),
        Option( "-fcz", "float", 1000, "force constant z"),
        ]

    files = [
        FileOption("-i", "r",["itp"],"input.itp", "read itp"),
        FileOption("-ipdb", "r/o",["pdb"],"mol.pdb", "optionally read pdb"),
        FileOption("-o", "w",["itp"],"output.itp", "write itp"),
        FileOption("-residfile", "r/o",["dat"],"resid.dat", "residue id file"),
        FileOption("-atomidfile", "r/o",["dat"],"atomid.dat", "atom id file"),
        ]

    help_text = ( ' Add position restraints for an atom to an .itp',
                  '',
                  ' If -ipdb is set, will calculate COM and search ',
                  ' for the closest heavy atom.',
		)


    cmdl = Commandline( argv, options = options,
                        fileoptions = files,
                        program_desc = help_text,
                        check_for_existing_files = False, version = version)

    residfile = None
    atomidfile = None
    pdbfile = None
    if cmdl.opt['-residfile'].is_set:
        residfile = cmdl['-residfile']
    if cmdl.opt['-atomidfile'].is_set:
        atomidfile = cmdl['-atomidfile']
    if cmdl.opt['-ipdb'].is_set:
        pdbfile = cmdl['-ipdb']

    PositionRestraints(num=cmdl['-num'], resnum=cmdl['-resnum'], name=cmdl['-name'], fcx=cmdl['-fcx'], fcy=cmdl['-fcy'], fcz=cmdl['-fcz'], itpname=cmdl['-i'], outname=cmdl['-o'],residfile=residfile,atomidfile=atomidfile,pdbfile=pdbfile)

main( sys.argv )

