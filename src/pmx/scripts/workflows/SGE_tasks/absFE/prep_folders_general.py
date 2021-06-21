#!/usr/bin/env python

import glob
import luigi
import os
import pmx
import shutil as sh
import math
from luigi.parameter import ParameterVisibility
from pmx import ndx
from pmx.model import Model
from pmx.scripts.workflows.utils import check_file_ready
from pmx.scripts.workflows.SGE_tasks.SGETunedJobTask import SGETunedJobTask #tuned for the owl cluster

# ==============================================================================
#                         Derivative Task Classes
# ==============================================================================
class Gather_Inputs_folder(SGETunedJobTask):
    #run on the login node
    run_locally = luigi.BoolParameter(
        visibility=ParameterVisibility.HIDDEN,
        default = True,
        significant=False,
        parsing=luigi.BoolParameter.EXPLICIT_PARSING,
        description="run locally instead of on the cluster")

    #Parameters:
    folder_path = luigi.Parameter(significant=False,
                      visibility=ParameterVisibility.HIDDEN,
                      description='Path to the protein+ligand folder to set up')
    p = luigi.Parameter(description='Protein name')
    l = luigi.Parameter(description='Ligand name')

    study_settings = luigi.DictParameter(significant=False,
                 visibility=ParameterVisibility.HIDDEN,
                 description='Dict of study stettings '
                 'used to propagate settings to dependencies')

    prot_src_override = luigi.Parameter(significant=False,
              visibility=ParameterVisibility.HIDDEN,
              default="",
              description='Path to protein stucture to use. Use only for Apo simulations. If empty will use default Holo structure.')

    #avoid Prameter not a string warnings
    job_name_format = luigi.Parameter(
        visibility=ParameterVisibility.HIDDEN,
        significant=False, default="pmx_{task_family}_p{p}_l{l}",
        description="A string that can be "
        "formatted with class variables to name the job with qsub.")

    #request 1 cores
    n_cpu = luigi.IntParameter(default=1, significant=False,
                               visibility=ParameterVisibility.HIDDEN)

    #variables to be overwriten in sub class' __init__()
    srctop=""
    posre=True

    def find_chains(self, path):
        chains=[]
        with open(path) as itp:
            inmoltype=False
            for line in itp:
                line.rstrip()
                if(not line): continue
                elif(not inmoltype and "[ moleculetype ]" in line):
                    inmoltype=True
                elif(inmoltype):
                    if(line[0]==';'): continue
                    elif(line[0]=='['): inmoltype=False #next block
                    else:
                        s=line.split()
                        chains.append(s[0])
                        inmoltype=False

        if(len(chains)<1):
            raise(Exception("Error: No chains found in %s!"%path))
        return(chains)


    def work(self):

        #make folder
        os.makedirs(self.folder_path, exist_ok=True)
        os.chdir(self.folder_path)

        #initial coordinates
        if(self.p and self.l): #P+L
            sh.copy(self.study_settings['top_path']+"/proteins/"+self.p+"/prot_"+self.l+".pdb",
                    self.folder_path+"/init.pdb")
            #append crystalographic solvent, if any
            if(os.path.isfile(self.study_settings['top_path']+"/proteins/"+self.p+"/water.pdb")):            
                os.system("cat {} >> {}".format(
                    self.study_settings['top_path']+"/proteins/"+self.p+"/water.pdb",
                    self.folder_path+"/init.pdb")
                    )
        elif(self.p and not self.l): #ApoP
            if(self.prot_src_override):
                #sh.copy(self.prot_src_override, self.folder_path+"/init.pdb")
                raise(ValueError("prot_src_override option is no longer supported. Please use a prot_apo.pdb file in the data/proteins/{}/ folder. A matching prot_apo.itp file is now also required.".format(self.p)))
            elif(os.path.isfile(self.study_settings['top_path']+"/proteins/"+self.p+"/prot_apo.pdb")):
                #Replaces prot_src_override functionality and supports missing residues in apo or holo.
                #Residue numbering needs to match between the two.
                sh.copy(self.study_settings['top_path']+"/proteins/"+self.p+"/prot_apo.pdb", self.folder_path+"/init.pdb")
                
                #append crystalographic solvent, if any
                water_pdb_fn=self.study_settings['top_path']+"/proteins/"+self.p+"/water_apo.pdb"
                if(os.path.isfile(water_pdb_fn)):
                    #sanity check: crystalographic water should start with water,
                    #as ligand gets inserted right before first water molecule
                    with open(water_pdb_fn, 'r') as f:
                        firstLine = f.readline()
                        s=firstLine.split()
                        if(s[3]!="HOH" and s[3]!="SOL"):
                            raise(RuntimeError("{} starts wither with a header or a molecule that is not water.\n"
                                "Please remove the header and place any crystalographic ions into prot_apo.pdb "
                                "or at the end of water_apo.pdb instead.".format(water_pdb_fn)))
                                
                    os.system("cat {} >> {}".format(water_pdb_fn, self.folder_path+"/init.pdb"))
                    
                #same check needs to be present for holo, as ligand resid is defined relative to final resid of protein in holo state
                water_pdb_fn=self.study_settings['top_path']+"/proteins/"+self.p+"/water.pdb"
                with open(water_pdb_fn, 'r') as f:
                    firstLine = f.readline()
                    s=firstLine.split()
                    if(s[3]!="HOH" and s[3]!="SOL"):
                        raise(RuntimeError("{} starts wither with a header or a molecule that is not water.\n"
                            "Please remove the header and place any crystalographic ions into prot.pdb "
                            "or at the end of water.pdb instead.".format(water_pdb_fn)))
                
                
            else:  #use holo structure/itp for apo
                sh.copy(self.study_settings['top_path']+"/proteins/"+self.p+"/prot.pdb",
                        self.folder_path+"/init.pdb")
                        
                #append crystalographic solvent, if any
                water_pdb_fn=self.study_settings['top_path']+"/proteins/"+self.p+"/water.pdb"
                if(os.path.isfile(water_pdb_fn)):
                    #sanity check: crystalographic water should start with water,
                    #as ligand gets inserted right before first water molecule
                    with open(water_pdb_fn, 'r') as f:
                        firstLine = f.readline()
                        s=firstLine.split()
                        if(s[3]!="HOH" and s[3]!="SOL"):
                            raise(RuntimeError("{} starts wither with a header or a molecule that is not water.\n"
                                "Please remove the header and place any crystalographic ions into prot.pdb "
                                "or at the end of water.pdb instead.".format(water_pdb_fn)))
                                
                    os.system("cat {} >> {}".format(water_pdb_fn, self.folder_path+"/init.pdb"))
        elif(self.l): #L
            sh.copy(self.study_settings['top_path']+"/ligand/"+self.l+"/ligand.pdb",
                    self.folder_path+"/init.pdb")
        
        
        
        #topology
        sh.copy(self.study_settings['top_path']+"/"+self.srctop, self.folder_path+"/topol.top")
        if(self.l):
            sh.copy(self.study_settings['top_path']+"/ligand/"+self.l+"/lig.itp",self.folder_path+"/lig.itp")
        if(self.p):
            #Both ApoP and P+L both need the apo itp as
            #C->A TI topology will point be importing "prot_apo.itp"
            if(self.prot_src_override):
                raise(ValueError("prot_src_override option is no longer supported. Please use a prot_apo.pdb file in the data/proteins/{}/ folder. A matching prot_apo.itp file is now also required.".format(self.p)))
            elif(os.path.isfile(self.study_settings['top_path']+"/proteins/"+self.p+"/prot_apo.pdb")):
                if(not os.path.isfile(self.study_settings['top_path']+"/proteins/"+self.p+"/prot_apo.itp")):
                    raise(RuntimeError("Detected prot_apo.pdb but no matching prot_apo.itp is available in {} .".format(self.study_settings['top_path']+"/proteins/"+self.p+"/")))
                sh.copy(self.study_settings['top_path']+"/proteins/"+self.p+"/prot_apo.itp",self.folder_path+"/prot_apo.itp")
            else: #use holo structure/itp for apo
                sh.copy(self.study_settings['top_path']+"/proteins/"+self.p+"/prot.itp",self.folder_path+"/prot_apo.itp")
                
                
        if(self.p and self.l): #P+L also needs the holo itp
            sh.copy(self.study_settings['top_path']+"/proteins/"+self.p+"/prot.itp",self.folder_path+"/prot.itp")
            
            ##if present, copy over the ndx to exclude protein atoms eligible for restraint anchors
            #if(os.path.isfile(self.study_settings['top_path']+"/proteins/"+self.p+"/prot_apo_exclude.ndx")):
                #sh.copy(self.study_settings['top_path']+"/proteins/"+self.p+"/prot_apo_exclude.ndx",self.folder_path+"/prot_apo_exclude.ndx")
            

        #generate temporary index file
        os.system("echo 'q\n' | gmx make_ndx -f init.pdb "
                  "-o index.ndx > setup.log 2>&1")
        #clean dublicate entries
        ndx_file = ndx.IndexFile("index.ndx", verbose=False)
        ndx_file.write("index.ndx")
        check_file_ready("index.ndx")

        if(self.posre):
            #generate restraints for equillibration
            #TODO: rewrite this to use the pmx Topology class
            if(self.p):
                #find out how many chains (molecule types) there are in the protein
                source_itp="prot.itp"
                if(not self.l): #apoP
                    source_itp="prot_apo.itp"
                chains=self.find_chains(source_itp)

                #Position restraints are only needed before NPT,
                #so TI itp files don't need them
                #and we can skip generating them for state C->A TI in case of L+P

                if(len(chains)==1):
                    #only one chain; its safe to make a single restraint file
                    #call it prot
                    os.system("echo 'Backbone\n' | gmx genrestr -f init.pdb "
                              "-fc 1000 1000 1000 -o prot_posre.itp "
                              "-n index.ndx >> setup.log 2>&1")
                    check_file_ready("prot_posre.itp")
                    os.system("echo 'Backbone\n' | gmx genrestr -f init.pdb "
                              "-fc 500 500 500 -o prot_posre_soft.itp "
                              "-n index.ndx >> setup.log 2>&1")
                    check_file_ready("prot_posre_soft.itp")

                elif(len(chains)>1):
                    #multiple chains, need separate files that are included from prot.itp

                    #split init.pdb into chains
                    chaindict={}
                    an=1
                    with open("init.pdb","r") as initstruct:
                        for line in initstruct:
                            if("ATOM" in line or "HETATM" in line):
                                s=line.split()
                                if(s[4][0].isalpha()): #forth column of pdb is a letter if atom in chain
                                    key=s[4][0]
                                    if key in chaindict.keys():
                                        chaindict[key].append(an)
                                    else:
                                        chaindict.update({key:[an]})
                                an+=1

                    with open("chains.ndx","w") as chainndx:
                        for ch in chaindict.keys():
                            chainndx.write("[ %s ]\n"%ch)
                            nentries=0;
                            width=int(math.log10(max(chaindict[ch]))+2)
                            perrow=int(79/(width))
                            for ind in chaindict[ch]:
                                if(nentries>0 and nentries%perrow==0):
                                    chainndx.write("\n")
                                chainndx.write("{ind:>{width}}".format(ind=ind, width=width))
                                nentries+=1;
                            chainndx.write("\n")
                    check_file_ready("chains.ndx")


                    for ch in chains:
                        os.system("echo '{key}\n' | gmx editconf -f init.pdb "
                              "-o {ch}.pdb "
                              "-n chains.ndx >> setup.log 2>&1".format(ch=ch, key=ch[-1]))
                        check_file_ready("{ch}.pdb".format(ch=ch))
                        os.system("echo 'Backbone\n' | gmx genrestr -f {ch}.pdb "
                              "-fc 1000 1000 1000 -o {ch}_posre.itp "
                              ">> setup.log 2>&1".format(ch=ch, key=ch[-1]))
                        check_file_ready("{ch}_posre.itp".format(ch=ch))
                        os.system("echo 'Backbone\n' | gmx genrestr -f {ch}.pdb "
                              "-fc 500 500 500 -o {ch}_posre_soft.itp "
                              ">> setup.log 2>&1".format(ch=ch, key=ch[-1]))
                        check_file_ready("{ch}_posre_soft.itp".format(ch=ch))


            if(self.l):
                os.system("echo 'MOL\n' | gmx editconf -f init.pdb "
                          "-o lig.pdb -n index.ndx >> setup.log 2>&1")
                check_file_ready("lig.pdb")
                os.system("echo '2\n' | gmx genrestr -f lig.pdb "
                          "-fc 1000 1000 1000 "
                          "-o lig_posre.itp >> setup.log  2>&1")
                check_file_ready("lig_posre.itp")
                os.system("echo '2\n' | gmx genrestr -f lig.pdb "
                          "-fc 500 500 500 "
                          "-o lig_posre_soft.itp >> setup.log 2>&1")
                check_file_ready("lig_posre_soft.itp")

        #clean overwritten files
        cleanList = glob.glob(self.folder_path+'/#*')
        for filePath in cleanList:
            try:
                os.unlink(filePath)
            except:
                raise OSError("Error while deleting file: "+filePath)

        #Return to basepath
        os.chdir(self.study_settings['base_path'])

    def output(self):
        files=["topol.top", "init.pdb" , "index.ndx"]
        if(self.p):
            files.extend(["prot_apo.itp"])
            if(self.l):
                files.extend(["prot.itp"])
            if(self.posre):
                chains=[]
                if(self.l):
                    chains=self.find_chains(self.study_settings['top_path']+"/proteins/"+self.p+"/prot.itp")
                else:
                    chains=self.find_chains(self.study_settings['top_path']+"/proteins/"+self.p+"/prot_apo.itp")
                if(len(chains)==1):
                    files.extend(["prot_posre.itp", "prot_posre_soft.itp"])
                else:
                    for ch in chains:
                        files.extend(["{ch}_posre.itp".format(ch=ch),
                                      "{ch}_posre_soft.itp".format(ch=ch)])

        if(self.l):
            files.extend(["lig.itp"])
            if(self.posre):
                files.extend(["lig_posre.itp", "lig_posre_soft.itp"])
        
        return [luigi.LocalTarget(os.path.join(self.folder_path, f)) for f in files]


class Prep_folder(SGETunedJobTask):
    #run on the login node
    run_locally = luigi.BoolParameter(
        visibility=ParameterVisibility.HIDDEN,
        default = True,
        significant=False,
        parsing=luigi.BoolParameter.EXPLICIT_PARSING,
        description="run locally instead of on the cluster")

    #Parameters:
    folder_path = luigi.Parameter(significant=False,
                         visibility=ParameterVisibility.HIDDEN,
                         description='Path to the protein+ligand folder to set up')
    p = luigi.Parameter(description='Protein name')
    l = luigi.Parameter(description='Ligand name')

    study_settings = luigi.DictParameter(significant=False,
                 visibility=ParameterVisibility.HIDDEN,
                 description='Dict of study stettings '
                 'used to propagate settings to dependencies')

    #avoid Prameter not a string warnings
    job_name_format = luigi.Parameter(
        visibility=ParameterVisibility.HIDDEN,
        significant=False, default="", description="A string that can be "
        "formatted with class variables to name the job with qsub.")

    #request 1 cores
    n_cpu = luigi.IntParameter(default=1, significant=False,
                               visibility=ParameterVisibility.HIDDEN)
                               
    _found_xray_SOL=False;

    def solvate(self):
        """Solvates the system.

        Subclasses can override this to set a custom number of water molecules

        Returns
        -------
        None.

        """
        os.system("gmx solvate -scale 1.0 -cp box.pdb -o water.pdb "\
                  "-cs spc216.gro -p topol_solvated.top >> prep.log 2>&1")

    def gen_ions(self,top_ions,pdb_ions):
        """Generates ions in the system.
        """
        os.system("echo 'SOL' | gmx genion -s tpr.tpr "
                          "-p %s -conc %f "
                          "-neutral -nname Cl -pname Na "
                          "-o %s >> genion.log 2>&1" %(
                              top_ions, self.study_settings['salt_conc'], pdb_ions) )



    def work(self):
        os.chdir(self.folder_path)

        os.system("gmx editconf -f init.pdb -o box.pdb -bt %s -d %f "\
                  "> prep.log 2>&1"%(self.study_settings['bt'], self.study_settings['d']))
        check_file_ready("box.pdb")
        sh.copy("topol.top","topol_solvated.top")
        
        
        #Check input for sanity:
            #number of protein chains
        m_box = Model("box.pdb", bPDBTER=True)
        n_prot_chains=0 #count number of protein chains. Assume that they are all the ones before the ligand.
        for c in m_box.chains:
            if("-MOL-" in c.get_sequence()):
                break;
            else:
                n_prot_chains+=1
        if(n_prot_chains<1):
            raise(RuntimeError("There should be at least one protein chain before the ligand! Check if the ligand has the same chain id as the protein!"))
        

        #solvate using old VdW radii (found in mutff). This needs pmx's mutff
        #to be the first entry in GMXLIB. It contains the old version of
        #vdwradii.dat which produces better water density.
        if "GMXLIB" in os.environ.keys():
            orig_GMXLIB = os.environ["GMXLIB"]
        else:
            orig_GMXLIB = ""
        os.environ["GMXLIB"] = os.path.join(os.path.dirname(pmx.__file__),
                  "data/mutff/") + os.pathsep + orig_GMXLIB
        self.solvate()
        os.environ["GMXLIB"] = orig_GMXLIB

        check_file_ready("water.pdb")
        
        #find solvent that is not crystalographic and replace only it with ions
        n = ndx.IndexFile("index.ndx")
        if("SOL" in n.names):
            self._found_xray_SOL = True
            xray_SOL_gid=n.names.index("SOL")
            n_gids=len(n.names)
            os.system("echo \"name {xr} xraySOL\n\nr HOH | r SOL & ! {xr}\n\nname {last} SOL\n\nq\n\" | ".format(
                        xr=xray_SOL_gid, last=n_gids) +
                        "gmx make_ndx -f water.pdb -n index.ndx -o gen_ion.ndx >> prep.log 2>&1")
            check_file_ready("gen_ion.ndx")
                        
        os.system("gmx grompp -p topol_solvated.top -c water.pdb -o tpr.tpr "\
                  "-f {mdp} {ndx} -v -maxwarn 3 "\
                  ">> prep.log 2>&1".format(mdp=self.init_mdp, ndx="" if not self._found_xray_SOL else "-n gen_ion.ndx"))
        check_file_ready("tpr.tpr")

        #generate ions for each
        #independent repeat (multiple for confidence estimate)
        for i in range(self.study_settings['n_repeats']):
            #sampling simulations in each repeat
            for m in range(self.study_settings['n_sampling_sims']):
                top_ions="topol_ions%d_%d.top"%(i,m)
                pdb_ions="ions%d_%d.pdb"%(i,m)
                if(os.path.isfile(pdb_ions)): #skip if it already exists
                    continue
                sh.copy("topol_solvated.top", top_ions)
                self.gen_ions(top_ions,pdb_ions)
                check_file_ready(pdb_ions)

        #generate index for alignment
        if(self.p and self.l): #P+L
            os.system("echo \"\nq\n\" | "
                  "gmx make_ndx -f ions0_0.pdb "
                  "-o index_prot_mol.ndx > /dev/null 2>&1")

            #clean duplivates and find correct group indeces
            prot_mol_ndx=ndx.IndexFile("index_prot_mol.ndx", verbose=False)
            prot_mol_ndx.write("index_prot_mol.ndx")
            prot_id = prot_mol_ndx.get_group_id("Protein")
            mol_id = prot_mol_ndx.get_group_id("MOL")
            os.system("echo \"{}|{}\n\nq\n\" | ".format(prot_id, mol_id) +
                  "gmx make_ndx -f ions0_0.pdb -n index_prot_mol.ndx "
                  "-o index_prot_mol.ndx >> setup.log 2>&1")

        cleanList = glob.glob(self.folder_path+'/#*')
        for filePath in cleanList:
            try:
                os.unlink(filePath)
            except:
                print("Error while deleting file: ", filePath)

        #Return to basepath
        os.chdir(self.study_settings['base_path'])

    def output(self):
        files=[]
        #independent repeat (multiple for confidence estimate)
        for i in range(self.study_settings['n_repeats']):
            #sampling simulations in each repeat
            for m in range(self.study_settings['n_sampling_sims']):
                files.append("topol_ions%d_%d.top"%(i,m))
                files.append("ions%d_%d.pdb"%(i,m))

        if(self.p and self.l): #P+L
            files.append("index_prot_mol.ndx")

        return [luigi.LocalTarget(os.path.join(self.folder_path, f)) for f in files]
