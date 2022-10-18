import sys, os
import copy as cp
import numpy as np
from pmx import *
from pmx import geometry
from pmx import ndx
from pmx import library
from pmx.forcefield2 import *
#from AbsRestraints import *

def find_atom_inY( ind, strX, strY ):
    aX = strX.atoms[ind-1]
    aY = 0
    minDist = 999.99
    for a in strY.atoms:
        d = a-aX
        if d<minDist:
            minDist = d
            aY = a
    return(aY.id)

def gen_new_ii( inpfname, outfname, strX, strY ):
    fp = open(inpfname,'r')
    lines = fp.readlines()
    fp.close()

    fp = open(outfname,'w')
    for l in lines:
        if '[' in l:
            fp.write(l)
        elif l=='\n':
            fp.write(l)
        else:
            l = l.rstrip()
            l = l.lstrip()
            foo = l.split()
#             print(foo,len(foo))

            # go over the array in reverse order to get atomID
            out = []
            for i in range(1,len(foo)+1):
                if i>5:
                    atomIndY = find_atom_inY( int(foo[-i]), strX, strY )
                    out.insert(0,atomIndY)
                else:
                    out.insert(0,foo[-i])

            # go over the elements again to output them
            fp.write('   ')
            for i in range(0,len(foo)):
                fp.write('%s   ' % str(out[i]))
            fp.write('\n')
             
    fp.close()

def main(argv):

    desc=('Transfer relative restraints from stateX to stateY.',
          'The structures in stateX and stateY do not need to match.',
          '',)

# define input/output files

    files= [
        FileOption("-fX", "r",["pdb","gro"],"structX.pdb", "input structure in stateX"),
        FileOption("-fY", "r",["pdb","gro"],"structY.pdb", "input structure in stateY"),
        FileOption("-iiX", "r",["itp"],"iiX.itp", "input restraint itp"),
        FileOption("-iiY", "o",["itp"],"iiY.itp", "output restraint itp"),
#        FileOption("-nX", "r",["ndx"],"ndxX.ndx", "input index file stateX"),
#        FileOption("-nY", "r",["ndx"],"ndxY.ndx", "input index file stateY"),
        ]

# define options

    options=[
       Option( "-fit", "bool", False, "fit"),
        ]

    help_text = ('Read structure in stateX and stateY.',
          'Read itp for restraints in stateX.',
          'Generate itp for stateY.',
          'TODO: fitting not available yet',)

# pass options, files and the command line to pymacs

    cmdl = Commandline( argv, options = options, fileoptions = files, program_desc = help_text, check_for_existing_files = False, version = "0.0" )

    strX = Model(cmdl['-fX'])
    strY = Model(cmdl['-fY'])

    gen_new_ii( cmdl['-iiX'], cmdl['-iiY'], strX, strY  )
    

main( sys.argv )

