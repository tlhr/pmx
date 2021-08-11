#!/usr/bin/env python

import luigi
import MDAnalysis as md
import os
import sys
import glob
from io import StringIO
from luigi.parameter import ParameterVisibility
from pmx import ndx
from pmx.forcefield import Topology
from pmx.ndx import IndexGroup
from pmx.scripts.workflows.SGE_tasks.SGETunedJobTask import SGETunedJobTask #tuned for the owl cluster
from pmx.scripts.workflows.find_anchors_and_write_ii import find_restraints
from pmx.scripts.workflows.find_anchors_and_write_ii_single_traj import find_restraints_align2crystal
from pmx.scripts.workflows.find_avg import find_avg_struct
from pmx.scripts.workflows.SGE_tasks.absFE.LinP.equil_sims import Sim_PL_NPT
from pmx.scripts.workflows.SGE_tasks.absFE.LinP.alignment import Task_PL_align
from pmx.scripts.workflows.SGE_tasks.absFE.LinP.morphes import Task_PL_gen_morphes
from pmx.scripts.workflows.utils import check_file_ready, readii_util, writeii_util
from pmx.scripts.workflows.postHoc_restraining_python3 import main as main_postHock_restr

# ==============================================================================
#                            Helper Functions
# ==============================================================================
def clean_virtual_sites_from_ndx(ndx_fn, mol_sel, filteree, itp):
    my_ndx=ndx.IndexFile(ndx_fn, verbose=False)
    whole_mol_ids = my_ndx[mol_sel].ids
    filteree_ids = my_ndx[filteree].ids
    first_at = whole_mol_ids[0]

    top = Topology(itp, is_itp=True, assign_types=False)
    vsites=[] #indeces of virtual sites in ligand
    for l in [top.virtual_sites2, top.virtual_sites3, top.virtual_sites4]:
        for vs in l:
            vsites.append(vs[0].id - 1 + first_at)

    new_filteree_ids = []
    for i in filteree_ids:
        if not (i in vsites):
            new_filteree_ids.append(i)
    g = IndexGroup(ids=new_filteree_ids, name=filteree+"_&_!vsites")
    my_ndx.add_group(g)
    my_ndx.write(ndx_fn)

   

# ==============================================================================
#                         Derivative Task Classes
# ==============================================================================
class Task_PL_gen_restraints(SGETunedJobTask):

    #Parameters:
    p = luigi.Parameter(description='Protein name')
    l = luigi.Parameter(description='Ligand name')
    i = luigi.IntParameter(description='Repeat number')

    folder_path = luigi.Parameter(significant=False,
                 visibility=ParameterVisibility.HIDDEN,
                 description='Path to the protein+ligand folder to set up')

    study_settings = luigi.DictParameter(significant=False,
                 visibility=ParameterVisibility.HIDDEN,
                 description='Dict of study stettings '
                 'used to propagate settings to dependencies')

    restr_scheme = luigi.Parameter(significant=True,
                 description='Restraint scheme to use. '
                 'Aligned, Fitted or Fixed')

    #request 1 cores
    n_cpu = luigi.IntParameter(visibility=ParameterVisibility.HIDDEN,
                               default=1, significant=False)

    #avoid Prameter not a string warnings
    job_name_format = luigi.Parameter(
        visibility=ParameterVisibility.HIDDEN,
        significant=False, default="pmx_{task_family}_p{p}_l{l}",
        description="A string that can be "
        "formatted with class variables to name the job with qsub.")

    #debug output
    debug = luigi.BoolParameter(
        visibility=ParameterVisibility.HIDDEN,
        significant=False, default=False,
        description="show debug output in a log.")

    min_ang_K = luigi.FloatParameter(visibility=ParameterVisibility.HIDDEN,
                               default=0.0, significant=True,
                               description="Minimal value for an angular force constant (kJ/(mol*rad^2)).")

    restr_source = luigi.Parameter(
        visibility=ParameterVisibility.HIDDEN,
        significant=False, default="aligned",
        description="Which trajectory to base restraints on. Options: aligned, coupled.")
    
    #verbose_completeness=luigi.BoolParameter(
        #visibility=ParameterVisibility.HIDDEN,
        #significant=False, default=False,
        #description="Should we list missing outputfiles?")

    extra_packages=[md]


    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        #set variables
        self.base_path = self.study_settings['base_path']

        if(self.restr_scheme=="Aligned"):
            self.states=["C"]
        elif(self.restr_scheme=="Fitted"):
            self.states=["A", "B"]
        elif(self.restr_scheme=="Fixed"):
            raise(Exception("Fixed restraints not yet implemened."))
            self.states = self.study_settings['states']
        else:
            raise(Exception("Unrecognized restraint scheme '%s'"%self.restr_scheme))

    def work(self):
        #generate morphs for A state
        os.makedirs(self.folder_path, exist_ok=True)
        os.chdir(self.folder_path)


        if(self.debug):
            print("debug: restr_scheme={}".format(self.restr_scheme))
        #generate the restraints
        if(self.restr_scheme=="Aligned"):            
            #make index files for TI
            #A->C
            base_ndx_A=ndx.IndexFile("index_prot_mol.ndx", verbose=False)
            fn_ndx_A="index_prot_mol_noH_A_{i}.ndx".format(i=self.i)
                #add C-alpha_common
            fn_ndx_common_A_src=self.folder_path+"/stateC/repeat{i}/morphes{m}/PL_w_chains.ndx".format(i=self.i, m=0)
            ndx_common_A_src=ndx.IndexFile(fn_ndx_common_A_src, verbose=False)
            #base_ndx_A.add_group(ndx.IndexGroup(ids=ndx_common_A_src["C-alpha_common"].ids, name="C-alpha_common"))
            base_ndx_A.add_group(ndx_common_A_src["C-alpha_common"])
            base_ndx_A.write(fn_ndx_A)
            
                #clean hydrogens and virtual sites
            base_ndx_A=ndx.IndexFile(fn_ndx_A, verbose=False) #reload
            mol_id = base_ndx_A.get_group_id("MOL")
            sys_id = base_ndx_A.get_group_id("System")
            #For some reason naming is detected differently in pdb and gro files.
            #4-symbol names can get frame shifted by one.
            #So need to filter for 'H' in the second slot too.
            os.system("echo \"{mol_id} & ! a H* & ! a ?H*\n{sys_id} & ! {mol_id}\n\nq\n\" | "
                      "gmx make_ndx -f ions{i}_0.pdb -n {inndx} "
                      "-o {outndx} > noH_make_ndx_{i}.log 2>&1".format(
                          i=self.i, outndx=fn_ndx_A, mol_id=mol_id, sys_id=sys_id, inndx=fn_ndx_A
                          ))
            clean_virtual_sites_from_ndx(fn_ndx_A, "MOL", "MOL_&_!H*_&_!?H*", "lig.itp")
            check_file_ready(os.path.join(fn_ndx_A))
            if(self.debug):
                print("debug: made {}".format(fn_ndx_A))
                
            #C->A (decorelation will reuse this)
            fn_frame_src=self.folder_path+"/stateC/repeat{i}/morphes{m}/frame0.gro".format(i=self.i, m=0)
            fn_ndx_C="index_prot_mol_noH_C_{i}.ndx".format(i=self.i)
            os.system("echo \"\nq\n\" | gmx make_ndx -f {} -o {} > /dev/null 2>&1".format(fn_frame_src, fn_ndx_C))

                #add C-alpha_common
            base_ndx_C=ndx.IndexFile(fn_ndx_C, verbose=False)
            fn_ndx_common_C_src=self.folder_path+"/stateC/repeat{i}/morphes{m}/P_w_chains.ndx".format(i=self.i, m=0)
            ndx_common_C_src=ndx.IndexFile(fn_ndx_common_C_src, verbose=False)
            base_ndx_C.add_group(ndx_common_C_src["C-alpha_common"])
                #filter out flexible loops, if exluded ndx is available
            #fn_ndx_exclude=self.folder_path+"/prot_apo_exclude.ndx"
            fn_ndx_exclude=self.study_settings['top_path']+"/proteins/"+self.p+"/prot_apo_exclude.ndx"
            if(os.path.isfile(fn_ndx_exclude)):
                exclude_ndx=ndx.IndexFile(fn_ndx_exclude, verbose=False)
                if(not "exclude" in exclude_ndx.names):
                    raise(Exception("Could not find the [ exclude ] group in {}".format(fn_ndx_exclude)))
                filtered_ids=[j for j in ndx_common_C_src["C-alpha_common"].ids if (not j in exclude_ndx["exclude"].ids)]
                filtered_g = IndexGroup(ids=filtered_ids, name="C-alpha_common_filtered")
                base_ndx_C.add_group(filtered_g)
            base_ndx_C.write(fn_ndx_C)
            
                #Add Protein_MOL group; clean hydrogens and virtual sites
            prot_id = base_ndx_C.get_group_id("Protein")
            mol_id  = base_ndx_C.get_group_id("MOL")
            sys_id  = base_ndx_C.get_group_id("System")
            os.system("echo \"{prot_id}|{mol_id}\n{mol_id} & ! a H* & ! a ?H*\n{sys_id} & ! {mol_id}\n\nq\n\" | ".format(
                            prot_id=prot_id, mol_id=mol_id, sys_id=sys_id) +
                      "gmx make_ndx -f {struct} -n {inndx} -o {outndx} > noH_make_ndx_{i}.log 2>&1".format(
                            struct=fn_frame_src, inndx=fn_ndx_C, outndx=fn_ndx_C, i=self.i))
            clean_virtual_sites_from_ndx(fn_ndx_C, "MOL", "MOL_&_!H*_&_!?H*", "lig.itp")
            check_file_ready(os.path.join(fn_ndx_C))
            if(self.debug):
                print("debug: made {}".format(fn_ndx_C))


            ##if decorelating, make an index file from the aligned structure. It can have a different number of atoms if Apo differs from holo
            #if('decor_decoupled' in self.study_settings and  self.study_settings['decor_decoupled']):
                #os.system("echo \"q\n\" | "
                      #"gmx make_ndx -f stateC/repeat{i}/morphes0/frame0.gro "
                      #"-o decor_{i}.ndx > decor_make_ndx_{i}.log 2>&1".format(i=self.i))

            
            frame_path=self.folder_path+"/state{2}/repeat{3}/{5}{4}/"
            fn_ref_ndx=fn_ndx_C
            
            frame_trjs=""
            for m in range(self.study_settings['n_sampling_sims']):
                frame_trjs+=frame_path.format(self.p,self.l,"C",self.i,m,"morphes")+"/frame*.gro"

            if(not os.path.isfile("ii_C_{i}.itp".format(i=self.i)) or  not os.path.isfile("out_dg_{i}.dat".format(i=self.i))):
                oldstdin = sys.stdin
                oldstdout = sys.stdout
                oldstderr = sys.stderr

                my_ndx=ndx.IndexFile(fn_ref_ndx, verbose=False)
                prot_id = my_ndx.get_group_id("C-alpha_common")
                if("C-alpha_common_filtered" in my_ndx.names): #if we have filtered the indeces to exclude flexible loops, use that group
                    prot_id = my_ndx.get_group_id("C-alpha_common_filtered")
                mol_id = my_ndx.get_group_id("MOL_&_!H*_&_!?H*_&_!vsites")
                sys.stdin = StringIO( "{}\n{}\n".format(prot_id, mol_id) )
                with open("gen_restr{i}.log".format(i=self.i), 'w') as logf:
                    sys.stdout = logf
                    sys.stderr = logf

                    g=glob.glob(frame_trjs)
                    argv = ["postHoc_restraining_python3.py", "-f", *g, "-n", fn_ref_ndx,
                                "-oii", "ii_C_{i}.itp".format(i=self.i),
                                "-odg", "out_dg_{i}.dat".format(i=self.i),
                                "-min_K", "%f"%self.min_ang_K]

                    if(self.debug):
                        print("debug: starting postHoc_restraining")
                    main_postHock_restr(argv)
                    if(self.debug):
                        print("debug: after postHoc_restraining")

                sys.stdin = oldstdin
                sys.stdout = oldstdout
                sys.stderr = oldstderr
            check_file_ready("ii_C_{i}.itp".format(i=self.i))
            
            #create ii_A_#.itp by remaping the indeces from ii_C_#.itp
            if(not os.path.isfile("ii_A_{i}.itp".format(i=self.i))):
                A_ndx=ndx.IndexFile(fn_ndx_A, verbose=False)
                C_ndx=ndx.IndexFile(fn_ndx_C, verbose=False)
                common_A_ids=A_ndx["C-alpha_common"].ids
                common_C_ids=C_ndx["C-alpha_common"].ids
                mol_A_ids=A_ndx["MOL_&_!H*_&_!?H*_&_!vsites"].ids
                mol_C_ids=C_ndx["MOL_&_!H*_&_!?H*_&_!vsites"].ids
                relevant_A_ids=common_A_ids+mol_A_ids
                relevant_C_ids=common_C_ids+mol_C_ids
                C_to_A_dict = {relevant_C_ids[a]: relevant_A_ids[a] for a in range(len(relevant_C_ids))}
                
                lig_C_ids, pro_C_ids, means, ks=readii_util("ii_C_{i}.itp".format(i=self.i))
                lig_A_ids=[C_to_A_dict[a] for a in lig_C_ids]
                pro_A_ids=[C_to_A_dict[a] for a in pro_C_ids]
                
                writeii_util("ii_A_{i}.itp".format(i=self.i), lig_A_ids, pro_A_ids, means, ks)
            check_file_ready("ii_A_{i}.itp".format(i=self.i))
        

        elif(self.restr_scheme=="Fitted"):
            raise(RuntimeError("restr_scheme = Fitted is no longer supported."))
       

        #create a C state topology that holds ligand in place
        #with the restraint from the ii.itp files

        for m in range(self.study_settings['n_sampling_sims']):
            for s in ["A","C"]:
                top_ions="topol_ions%d_%d.top"%(self.i,m)
                if(s=="C"): # Apo sim can have different number of waters (if the whole protocol wasn't rerun from scratch after including Apo crystal structure)
                    top_ions=self.folder_path+"/../apoP/topol_ions%d_%d.top"%(self.i,m)

                check_file_ready(top_ions)
                topAC_ions="topolTI_ions%s%d_%d.top"%(s,self.i,m)
                already_inserted=False
                with open(topAC_ions, 'w') as top:
                    with open(top_ions, 'r') as reftop:
                        for l in reftop:
                            #ligand goes before the first water molecule, even when xray water is present.
                            #xray ions, if present, should come in the same chain as protein inside protein.pdb/protein_apo.pdb
                            #or after the xray water in water.pdb/water_apo.pdb
                            if (s=="C" and "SOL " in l and not already_inserted):
                                top.write("MOL    1\n") #add ligand into the apo topology
                                already_inserted=True
                                top.write(l)
                                
                            elif (s=="C" and ("#include \"prot_apo.itp\"" in l or "#include \"prot.itp\"" in l)):
                                top.write("#include \"lig.itp\"\n") #add ligand into the apo topology
                                top.write(l)
                            else:
                                top.write(l)

                    top.write("\n; Include intermolecular restraints\n")
                    top.write("#include \"ii_{s}_{i}.itp\"\n".format(s=s, i=self.i))

                check_file_ready(topAC_ions)

        #restore base path
        os.chdir(self.base_path)


    def requires(self):
        reqs=[]
        #sampling simulations in each repeat
        for m in range(self.study_settings['n_sampling_sims']):
            #states of equilibrium sims
            if(self.restr_scheme=="Aligned"):
                reqs.append(Task_PL_align(p=self.p, l=self.l, i=self.i, m=m, sTI='C',
                                  study_settings=self.study_settings,
                                  folder_path=self.folder_path,
                                  parallel_env=self.parallel_env,
                                  restr_scheme=self.restr_scheme)
                        )
                if(self.restr_source=="coupled"):
                    reqs.append(Task_PL_gen_morphes(p=self.p, l=self.l,
                              i=self.i, m=m, sTI='A',
                              study_settings=self.study_settings,
                              folder_path=self.folder_path,
                              parallel_env=self.parallel_env,
                              restr_scheme=self.restr_scheme)
                        )
            elif(self.restr_scheme=="Fitted"):
                for s in self.states:
                    reqs.append(
                        Sim_PL_NPT(p=self.p, l=self.l, i=self.i, m=m, s=s,
                                  study_settings=self.study_settings,
                                  folder_path=self.folder_path,
                                  parallel_env=self.parallel_env)
                        )


        return(reqs)

    def output(self):
        targets=[
            luigi.LocalTarget(os.path.join(self.folder_path, "ii_A_{i}.itp".format(i=self.i))),
            luigi.LocalTarget(os.path.join(self.folder_path, "ii_C_{i}.itp".format(i=self.i))),
            luigi.LocalTarget(os.path.join(self.folder_path, "out_dg_{i}.dat".format(i=self.i))),
            luigi.LocalTarget(os.path.join(self.folder_path, "index_prot_mol_noH_A_{i}.ndx".format(i=self.i))),
            luigi.LocalTarget(os.path.join(self.folder_path, "index_prot_mol_noH_C_{i}.ndx".format(i=self.i)))
            ]
        if(self.restr_scheme=="Fitted"):
            targets.append([luigi.LocalTarget(os.path.join(self.folder_path, "averageA_{i}.gro".format(i=self.i)))])
        # for s in self.states:
        #     targets.append([
        #         luigi.LocalTarget(os.path.join(self.folder_path, "dump{s}_{i}.gro".format(s=s,i=self.i))),
        #         luigi.LocalTarget(os.path.join(self.folder_path, "all_eq{s}_{i}_fit.xtc".format(s=s,i=self.i)))
        #         ])

        for m in range(self.study_settings['n_sampling_sims']):
            for s in ["A","C"]:
                targets.append(luigi.LocalTarget(os.path.join(self.folder_path,
                                          "topolTI_ions%s%d_%d.top"%(s,self.i,m))))

        if('decor_decoupled' in self.study_settings and  self.study_settings['decor_decoupled'] and
           (self.study_settings['decor_method']=="md" or self.study_settings['decor_method']=="sampling")):
            targets.append(luigi.LocalTarget(os.path.join(self.folder_path,
                                          "decor_%d.ndx"%(self.i))))
            
        #if(self.verbose_completeness):
            #for t in targets:
                #if(not os.path.isfile(t.path)):
                    #print("Incomplete becasue missing file: ", t.path)

        return targets
