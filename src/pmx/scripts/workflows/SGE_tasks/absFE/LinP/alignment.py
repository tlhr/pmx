#!/usr/bin/env python

import luigi
import numpy as np
import os
from pmx.scripts.workflows.SGE_tasks.SGETunedJobTask import SGETunedJobTask #tuned for the owl cluster
from pmx import ndx
from pmx.model import Model
from luigi.parameter import ParameterVisibility
from pmx.scripts.workflows.utils import check_file_ready
from pmx.scripts.workflows.fit_ligs_multiframes_python3 import fit,rotate_velocities_R, find_last_protein_atom
from pmx.scripts.workflows.utils import read_from_mdp
from pmx.scripts.workflows.SGE_tasks.absFE.LinW.equil_sims import Sim_WL_NPT
from pmx.scripts.workflows.SGE_tasks.absFE.LinP.equil_sims import Sim_PL_NPT
from pmx.scripts.workflows.SGE_tasks.absFE.ApoP.equil_sims import Sim_ApoP_NPT
from pmx.xtc import Trajectory


# ==============================================================================
#                         Derivative Task Classes
# ==============================================================================
class Task_PL_align(SGETunedJobTask):
    #Parameters:
    p = luigi.Parameter(description='Protein name')
    l = luigi.Parameter(description='Ligand name')
    i = luigi.IntParameter(description='Repeat number')
    m = luigi.IntParameter(description='Sampling sim number')
    sTI = luigi.Parameter(description='Coupling state for TI')

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

    stage="morphes"

    reuse_existing_pbc_fixes = luigi.BoolParameter(
        visibility=ParameterVisibility.HIDDEN,
        significant=False, default=False,
        description="If the trj_?.trr exist, reuse them instead of making new ones.")

    #request 1 cores
    n_cpu = luigi.IntParameter(visibility=ParameterVisibility.HIDDEN,
                               default=1, significant=False)

    #avoid Prameter not a string warnings
    job_name_format = luigi.Parameter(
        visibility=ParameterVisibility.HIDDEN,
        significant=False, default="pmx_{task_family}_p{p}_l{l}_{sTI}{i}_{m}",
        description="A string that can be "
        "formatted with class variables to name the job with qsub.")

    extra_packages=[np]

    #debug output
    debug = luigi.BoolParameter(
        visibility=ParameterVisibility.HIDDEN,
        significant=False, default=False,
        description="show debug output in a log.")


    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._setupState()

        #set variables
        self.base_path = self.study_settings['base_path']
        self.sim_path = self.folder_path+"/state%s/repeat%d/%s%d"%(
            self.sTI, self.i, self.stage, self.m)
        self.mdp = self.study_settings['mdp_path'] +\
            "/protein/eq_npt_{0}.mdp".format(
                self.study_settings['TIstates'][self.sTI])

        if(not (self.restr_scheme=="Aligned")):
            raise(Exception("{} should only be used if restraint "
                            "scheme is {}".format(
                                self.__class__.__name__ ,"Aligned")))

    def _setupState(self):
        if(self.sTI != "C"):
            raise(ValueError("Aligning morphes for TI state{}. "
                     "{} should only be done on TI stateC.".format(
                         self.sTI, self.__class__.__name__)))
            exit(1);
        # no need to set self.s as it is never used

    def requires(self):
        #restraints require both state A & B for all repeats and sampling sims
        tasks=[]
        #for i in range(self.study_settings['n_repeats']):
            #for m in range(self.study_settings['n_sampling_sims']):

        #Currently work() is only susing self.i and self.m
        #TODO: all the m's should be concatenated and so align should depend
        #on all of them at once
        i=self.i
        m=self.m
        tasks.append(Sim_WL_NPT(l=self.l, i=i, m=m, s='B',
                  study_settings=self.study_settings,
                  folder_path=self.base_path+"/water/lig_{}".format(self.l),
                  parallel_env=self.parallel_env))
        tasks.append(Sim_PL_NPT(p=self.p, l=self.l, i=i, m=m, s='A',
                  study_settings=self.study_settings,
                  folder_path=self.folder_path,
                  parallel_env=self.parallel_env))
        tasks.append(Sim_ApoP_NPT(p=self.p, i=i, m=m,
                  study_settings=self.study_settings,
                  folder_path=self.base_path+"/prot_{}/apoP".format(self.p),
                  parallel_env=self.parallel_env))

        return(tasks)

    def output(self):
        #find nframes by reading the mdp file
        end_time, dtframe = read_from_mdp(self.mdp)
        nframes=int(end_time/dtframe) - int(self.study_settings['b']/dtframe) +1 #first frame counts

        targets=[]
        for nf in range(nframes):
            targets.append(luigi.LocalTarget(
                os.path.join(self.sim_path, 'frame%d.gro'%nf)) )

        targets.append(luigi.LocalTarget(os.path.join(self.sim_path, 'aligned.trr')) )

        return targets
        
    #check completeness of outputs.
    #frame*.gro files can be incomplete if this task got terminated early (like due to exceeding storage quota).
    def complete(self):
        """
        Check if all dHdl files exist and are of correct length.
        """
        reqs_complete = all(r.complete() for r in luigi.task.flatten(self.requires()))
        if(reqs_complete):
            outputs = luigi.task.flatten(self.output())
            exist = list(map(lambda output: output.exists(), outputs))
            if not all(exist):
                return (False)
            
            finished_gros = []
            for o in outputs:
                if(o.path[-4:]==".gro"): #we only want to check gros (trrs are harder)
                    finished_gros.append(self._check_gro_finished(o.path))
            
            return(all(finished_gros))
        else:
            return(reqs_complete)

    #helper function for .gro completeness check
    def _check_gro_finished(self, fn):
        """Checks if a gro file is complete.

        Parameters
        ----------
        fn: filename

        Returns
        -------
        Boolean: True for finished.
        """
        ret=False
        with open(fn, 'rb') as f:
            #lines = f.read().splitlines()
            #box_size_line = lines[-2] # empty line after this
            
            #Faster version based on https://openwritings.net/pg/python/python-read-last-line-file
            f.seek(-2, os.SEEK_END)
            while f.read(1) != b'\n': # find start of last line
                f.seek(-2, os.SEEK_CUR) 
            while f.read(1) != b'\n': # find start of second to last line
                f.seek(-2, os.SEEK_CUR) 
            box_size_line=f.readline().decode()
            
            if(len(box_size_line.split())==9): # box info line should have 9 columns; atom lines only have 8
                ret=True
        
        return(ret)

    def gen_ndx_w_chains(self, struct_path, ndx_path, nchains):
        #create default
        os.system("echo 'q' | gmx make_ndx -f {gro} -o {out} "
                  ">> align.log 2>&1".format(
                    gro=struct_path, out=ndx_path ) )
        #add chains
        n = ndx.IndexFile(ndx_path)
        m = Model(struct_path, bPDBTER=True)
        for i in range(nchains):
            ch = m.chains[i]
            grp = ndx.make_index_group(ch.atoms, "chain_" + ch.id)
            n.add_group(grp)
        n.write(ndx_path)
    
    
    def gen_ndx_common_Calpha(self, apo_struct_path, holo_struct_path, apo_ndx_path, holo_ndx_path):
        """Adds an index group of common C-alpha atoms to ndx files of both apo and holo states.

        Args:
            apo_struct_path  (str): apo structure file.
            holo_struct_path (str): holo structure file.
            apo_ndx_path     (str): apo index file.
            holo_ndx_path    (str): holo index file.
            
        """
        
        n_a = ndx.IndexFile(apo_ndx_path)
        n_h = ndx.IndexFile(holo_ndx_path)
        m_a = Model(apo_struct_path, bPDBTER=True)
        m_h = Model(holo_struct_path, bPDBTER=True)
        
        ready_a = "C-alpha_common" in n_a.names
        ready_h = "C-alpha_common" in n_h.names
        if(ready_a and ready_h):
            return; #this was alredy done before, so skip
        
        common_a=[]
        common_h=[]
        
        res_arr_a=[r_a for r_a in m_a.residues if r_a.moltype=='protein']
        res_arr_h=[r_h for r_h in m_h.residues if r_h.moltype=='protein']
        
        for r_a in res_arr_a:
            for r_h in res_arr_h:
                if(r_a.orig_id == r_h.orig_id and r_a.resname == r_h.resname): #residue matches
                    for a_a in r_a.atoms:
                        for a_h in r_h.atoms:
                            if(a_a.name == a_h.name and (a_a.id in n_a["C-alpha"].ids) and (a_h.id in n_h["C-alpha"].ids)):
                                common_a.append(a_a.id)
                                common_h.append(a_h.id)
                                
        common_a.sort()
        common_h.sort()
        
        #add the index groups and write to files
        if(not ready_a):
            g_a = ndx.IndexGroup(ids=common_a, name="C-alpha_common")
            n_a.add_group(g_a)
            n_a.write(apo_ndx_path)
        else:
            if common_a != n_a["C-alpha_common"].ids:
                raise(Exception("Generated common protein index group does not match the one that already exists: [C-alpha_common] in %s".format(apo_ndx_path)))
            
        if(not ready_h):
            g_h = ndx.IndexGroup(ids=common_h, name="C-alpha_common")
            n_h.add_group(g_h)
            n_h.write(holo_ndx_path)
        else:
            if common_h != n_h["C-alpha_common"].ids:
                raise(Exception("Generated common protein index group does not match the one that already exists: [C-alpha_common] in %s".format(holo_ndx_path)))


    def work(self):

        mylog=open("align_py.log","w")

        #find number of protein chains
        m_init = Model(self.folder_path+"/init.pdb", bPDBTER=True)
        n_prot_chains=0 #count number of protein chains. Assume that they are all the ones before the ligand
        for c in m_init.chains:
            if("-MOL-" in c.get_sequence()):
                break;
            else:
                n_prot_chains+=1
        if(n_prot_chains<1):
            raise(RuntimeError("There should be at least one protein chain before the ligand! Check if the ligand has the same chain id as the protein!"))

        os.makedirs(self.sim_path, exist_ok=True)
        os.chdir(self.sim_path)
        if(os.path.isfile("align.log")): #clean old partial log if present
            os.unlink("align.log")

        trj_A_src=self.folder_path+"/state{2}/repeat{3}/npt{4}/".format(
            self.p, self.l, 'A', self.i, self.m)  #P+L
        trj_B_src=self.base_path+"/prot_{0}/apoP/repeat{3}/npt{4}/".format(
            self.p, self.l, None, self.i, self.m) #apoP
        trj_C_src=self.base_path+"/water/lig_{1}/state{2}/repeat{3}/npt{4}/".format(
            self.p, self.l, 'B', self.i, self.m)  #vacL

        #make the ndxs for the chains
        self.gen_ndx_w_chains(self.folder_path+"/ions%d_%d.pdb"%(self.i, self.m),
                              "PL_w_chains.ndx", n_prot_chains) #P+L
        check_file_ready("PL_w_chains.ndx")
        self.gen_ndx_w_chains(self.base_path+"/prot_{0}/apoP/ions{1}_{2}.pdb".format(
                                self.p, self.i, self.m),
                              "P_w_chains.ndx", n_prot_chains) #P
        check_file_ready("P_w_chains.ndx")
        
        #Update the PL and P ndx files with a group containing common C-alpha atoms.
        #Used for alignment of structures with missing residues.
        self.gen_ndx_common_Calpha(
            self.base_path+"/prot_{0}/apoP/ions{1}_{2}.pdb".format(self.p, self.i, self.m),
            self.folder_path+"/ions%d_%d.pdb"%(self.i, self.m),
            "P_w_chains.ndx", "PL_w_chains.ndx")
                
        #LW index
        os.system("echo 'q' | gmx make_ndx -f {gro} -o {out} "
                  ">> align.log 2>&1".format(
                  gro=self.base_path+"/water/lig_{1}/ions{3}_{4}.pdb".format(
                    self.p, self.l, 'B', self.i, self.m),
                  out="LW.ndx" ) )
        check_file_ready("LW.ndx")

        #Cut the begining off of trjs and center them
        names = ["A","B","C"]
        sources = [trj_A_src, trj_B_src, trj_C_src]
        ndxs = ["-n PL_w_chains.ndx", "-n P_w_chains.ndx", ""]
        for i in range(3):
            src=sources[i]
            if(self.reuse_existing_pbc_fixes and os.path.exists("trj_{}.trr".format(names[i]))):
                continue;

            #wrap mol centers
            os.system("echo System | gmx trjconv -s {tpr} -f {trj} -o {out} "
                  "-b {b} -ur compact -pbc mol"
                  ">> align.log 2>&1".format(
                      tpr=src+"tpr.tpr", trj=src+"traj.trr",
                      out="trj_{}_temp_cut_pbc.trr".format(names[i]),
                      b=self.study_settings['b']) )

            #center on chain_A and rewrap mol centers
            c = "chain_A"
            if(i==2): #LW
                c = "Other"
            os.system("echo {c} System | gmx trjconv -s {tpr} -f {trj} -o {out} "
                  "{ndx} -ur compact -center -pbc mol"
                  ">> align.log 2>&1".format(
                      tpr=src+"tpr.tpr", ndx=ndxs[i], c=c,
                      trj="trj_{}_temp_cut_pbc.trr".format(names[i]),
                      out="trj_{}.trr".format(names[i]) ) )

            check_file_ready("trj_{}.trr".format(names[i]))

            #clean temp
            os.unlink("trj_{}_temp_cut_pbc.trr".format(names[i]))

        #make the C state
        m_A = Model(self.folder_path+"/ions%d_%d.pdb"%(self.i,self.m),bPDBTER=False)
        m_B = Model(self.base_path+"/prot_{0}/apoP/ions{3}_{4}.pdb".format(
            self.p, self.l, None, self.i, self.m),bPDBTER=False) #apoP
        m_C = Model(self.base_path+"/water/lig_{1}/ions{3}_{4}.pdb".format(
            self.p, self.l, 'B', self.i, self.m),bPDBTER=False) #vacL
        m_A.a2nm()
        m_B.a2nm()
        m_C.a2nm()


        trj_A = Trajectory("trj_A.trr") #P+L
        trj_B = Trajectory("trj_B.trr") #apoP
        trj_C = Trajectory("trj_C.trr") #vacL



        #ndx_file_A = ndx.IndexFile(self.folder_path+"/index_prot_mol.ndx", verbose=False)
        ndx_file_A = ndx.IndexFile("PL_w_chains.ndx", verbose=False)
        ndx_file_B = ndx.IndexFile("P_w_chains.ndx", verbose=False)
        ndx_file_C = ndx.IndexFile("LW.ndx", verbose=False)
        #p_ndx = np.asarray(ndx_file_A["C-alpha"].ids)-1 # as in Vytas' alignment script
        pA_ndx = np.asarray(ndx_file_A["C-alpha_common"].ids)-1 # as in Vytas' alignment script
        pB_ndx = np.asarray(ndx_file_B["C-alpha_common"].ids)-1 # as in Vytas' alignment script
        linA_ndx = np.asarray(ndx_file_A["MOL"].ids)-1
        l_ndx = np.asarray(ndx_file_C["MOL"].ids)-1


        #find chain and resID of the last residue of the protein
        mol_first_atom = m_A.atoms[linA_ndx[0]]
        #chID = mol_first_atom.chain_id
        resID = mol_first_atom.resnr #internal numbering (residue.id, not residue.orig_id)
        A_mol_res_index=-1;
        A_prot_end_res_index=-1;
        for i,r in enumerate(m_A.residues):
            if(r.moltype=='protein'):
                A_prot_end_res_index = m_A.residues.index(r)
            if(r.id==resID):
                A_mol_res_index = m_A.residues.index(r)
                break;
        if(A_mol_res_index<0):
            raise("Could not find residue with resID %d in protein+ligand."%(resID))
        if(A_prot_end_res_index<0):
            raise("Could not find the last protein residue in protein+ligand.")
        
        mol_res_index_shift=A_mol_res_index-A_prot_end_res_index
        
        B_prot_end_res_index=-1;
        for i,r in enumerate(m_B.residues):
            if(r.moltype=='protein'):
                B_prot_end_res_index = m_B.residues.index(r)
        if(B_prot_end_res_index<0):
            raise("Could not find the last protein residue in ApoP.")
        B_mol_res_index=B_prot_end_res_index+mol_res_index_shift #this is where the ligand will go

        

        num_aligned_atoms = len(m_B.atoms) + l_ndx.shape[0]
        #if(len(m_A.atoms) != num_aligned_atoms):
            #print("\nWARNING: number of atoms in Apo + ligand ({}) does not match that of Holo ({})!\n".format(num_aligned_atoms,len(m_A.atoms)))
            #mylog.write("\nWARNING: number of atoms in Apo + ligand ({}) does not match that of Holo ({})!\n".format(num_aligned_atoms,len(m_A.atoms)))
            #mylog.flush()
        trj_out = Trajectory("aligned.trr", mode='Out',
                             atomNum = num_aligned_atoms) #aligned output

        #Frames are not acessible individually, just in sequence.
        #pmx.xtc.Trajectory is based on __iter__, so we need a custom
        #"for" loop to simultaneously go through both trajectories.
        #Based on https://www.programiz.com/python-programming/iterator
        iter_A = iter(trj_A)
        iter_B = iter(trj_B)
        iter_C = iter(trj_C)
        fridx=0
        mylog.write("\n\n\n\t\t###################\nStarting looping through frames:\n")
        mylog.flush()
        while True:
            try:
                frame_A = next(iter_A)
                frame_B = next(iter_B)
                frame_C = next(iter_C)
            except StopIteration:
                break

            #don't want frames while still equilibrating
            if(frame_A.time<self.study_settings['b']):
                continue


            frame_A.update(m_A)
            #m_b needs to be reloaded to have correct # of atoms next iteration
            m_B = Model(self.base_path+"/prot_{0}/apoP/ions{3}_{4}.pdb".format(
                self.p, self.l, None, self.i, self.m),bPDBTER=True) #apoP
            m_B.a2nm()
            frame_B.update(m_B, uv=True)
            frame_C.update(m_C, uv=True)
            mylog.write("frame %d: read models from pdbs\n"%fridx)
            #mylog.flush()

            if(self.debug):
                print("debug: \tframe {} before alignment".format(fridx))

                # # load x and v for output
                # x = np.zeros(len(m_A.atoms)*3)
                # v = np.zeros(len(m_A.atoms)*3)
                # for i, atom in enumerate(m_A.atoms):
                #     x[i*3:(i+1)*3]=atom.x
                #     v[i*3:(i+1)*3]=atom.v
                # print("Lig positions: \n", x[linA_ndx*3])
                # print("Lig velocities: \n", v[linA_ndx*3])


            # m_A.write("m_A1.gro")
            # step1: fit prot from prot+lig onto apo protein
            (v1,v2,R) = fit( m_B, m_A, pB_ndx, pA_ndx )
            # rotate velocities
            # not needed. We aren't saving m_A
            mylog.write("\t\tFit A on B\n")
            #mylog.flush()

            # step2: ligand onto the ligand from prot+lig structure
            (v1,v2,R) = fit( m_A, m_C, linA_ndx, l_ndx )
            mylog.write("\t\tFit C on A\n")
            #mylog.flush()
            # rotate velocities
            rotate_velocities_R( m_C, R )
            mylog.write("\t\tRotated C velocities\n")
            #mylog.flush()

            #insert vac ligand into B
            #do the insertion explicitly without relying on chains
            mol = m_C.residues[0]
            mol.model = m_B
            m_B.residues.insert(B_mol_res_index, mol)
            #don't add to chain, it doesn't matter
            m_B.al_from_resl()
            m_B.renumber_atoms()
            m_B.al_from_resl()


            mylog.write("\t\tInserted ligand from C into B\n")
            #mylog.flush()

            if(self.debug):
                print("debug: \tframe {} after alignment".format(fridx))

                # # load x and v for output
                # x = np.zeros(len(m_B.atoms)*3)
                # v = np.zeros(len(m_B.atoms)*3)
                # for i, atom in enumerate(m_B.atoms):
                #     x[i*3:(i+1)*3]=atom.x
                #     v[i*3:(i+1)*3]=atom.v
                # print("Lig positions: \n", x[linA_ndx*3])
                # print("Lig velocities: \n", v[linA_ndx*3])
                # print("Lig velocities in atomsB:")
                # for i, atom in enumerate(m_B.atoms):
                #     if(i in linA_ndx):
                #         print(atom.v)


            # output
            frame_fn="frame%d.gro"%fridx
            if(not os.path.isfile(frame_fn) or not self._check_gro_finished(frame_fn)):
                m_B.write(frame_fn)
                mylog.write("\t\tWrote B as %s\n"%frame_fn)
            mylog.flush()

            x = np.zeros(len(m_B.atoms)*3)
            v = np.zeros(len(m_B.atoms)*3)
            for i, atom in enumerate(m_B.atoms):
                x[i*3:(i+1)*3]=atom.x
                v[i*3:(i+1)*3]=atom.v

            mylog.write("\t\t\tSet atomx.x &.v; ready for writing trj frame to aligned.trr\n")
            trj_out.write_xtc_frame(step=frame_B.step, time=frame_B.time,
                                    lam=1.0, box=frame_B.box, x=x, v=v,
                                    units=m_B.unity, bTrr=True )

            mylog.write("\t\tWrote B to aligned trajectory\n")
            mylog.flush()

            # if(self.debug):
            #     raise(Exception("Debugging"))

            fridx+=1

        trj_out.close()
        os.system("rm -f \\#*")
        mylog.close()

        #restore base path
        os.chdir(self.base_path)
