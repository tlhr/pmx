#!/usr/bin/env python

import luigi
import os
from luigi.parameter import ParameterVisibility
from pmx.scripts.workflows.SGE_tasks.Sim import SGE_Sim
from pmx.scripts.workflows.SGE_tasks.absFE.LinP.prep_folders import Prep_PL_folder

# ==============================================================================
#                         Derivative Task Classes
# ==============================================================================
class Sim_PL_EM(SGE_Sim):

    #Parameters:
    p = luigi.Parameter(description='Protein name')
    l = luigi.Parameter(description='Ligand name')
    i = luigi.IntParameter(description='Repeat number')
    m = luigi.IntParameter(description='Sampling sim number')
    s = luigi.Parameter(description='Coupling state')

    folder_path = luigi.Parameter(significant=False,
                 visibility=ParameterVisibility.HIDDEN,
                 description='Path to the protein+ligand folder to set up')

    study_settings = luigi.DictParameter(significant=False,
                 visibility=ParameterVisibility.HIDDEN,
                 description='Dict of study stettings '
                 'used to propagate settings to dependencies')

    stage="em"
    #request 2 cores
    n_cpu = luigi.IntParameter(visibility=ParameterVisibility.HIDDEN,
                               default=2, significant=False)

    job_name_format = luigi.Parameter(
        visibility=ParameterVisibility.HIDDEN,
        significant=False, default="pmx_{task_family}_p{p}_l{l}_{s}{i}_{m}",
        description="A string that can be "
        "formatted with class variables to name the job with qsub.")

    restr_to_EM = luigi.BoolParameter(
        visibility=ParameterVisibility.HIDDEN,
        significant=False,
        default=False,
        description="restrain to EM or to crystal structure.")
        
    use_dbl_precision = luigi.BoolParameter(
        visibility=ParameterVisibility.HIDDEN,
        significant=False,
        default=False,
        description="Use double precision mdrun instead of single precision.")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        #set required file names
        self.base_path = self.study_settings['base_path']
        self.sim_path = self.folder_path+"/state%s/repeat%d/%s%d"%(
            self.s, self.i, self.stage, self.m)
        self.mdp = self.study_settings['mdp_path'] +\
            "/protein/em_posre_{0}.mdp".format(self.study_settings['states'][self.s])
        self.top = self.folder_path+"/topol_ions{3}_{4}.top".format(
            self.p, self.l, self.s, self.i, self.m)
        self.struct = self.folder_path+"/ions{3}_{4}.pdb".format(
            self.p, self.l, self.s, self.i, self.m)
        self.posre = self.folder_path+"/ions{3}_{4}.pdb".format(
            self.p, self.l, self.s, self.i, self.m)
        if(self.posre_ref_override):
            self.posre = self.posre_ref_override
        
        self.mdrun = self.study_settings['mdrun']
        if(self.use_dbl_precision):
            self.mdrun = self.study_settings['mdrun_double']
        self.mdrun_opts = self.study_settings['mdrun_opts']

    #work(): same as SGE_Sim

    def requires(self):
        return( Prep_PL_folder(p=self.p, l=self.l,
                               study_settings=self.study_settings,
                               folder_path=self.folder_path) )
                                #no need to pass parallel_env as
                                #Prep_PL_folder runs on the login node
    def output(self):
        return luigi.LocalTarget(os.path.join(self.sim_path, 'confout.gro'))


class Sim_PL_NVT_posre(Sim_PL_EM):
    stage="nvt_posre"
    #request 4 cores
    n_cpu = luigi.IntParameter(visibility=ParameterVisibility.HIDDEN,
                               default=4, significant=False)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        if(self.restr_to_EM):
            self.posre = self.folder_path+"/state{s}/repeat{i}/em{m}/confout.gro".format(
                s=self.s, i=self.i, m=self.m)
        if(self.posre_ref_override):
            self.posre = self.posre_ref_override

        #override relevant file names
        self.mdp = self.study_settings['mdp_path'] +\
            "/protein/eq_nvt_posre_{0}.mdp".format(self.study_settings['states'][self.s])
        self.struct = self.folder_path+"/state{2}/repeat{3}/em{4}/confout.gro".format(
            self.p, self.l, self.s, self.i, self.m)

    def requires(self):
        return( Sim_PL_EM(p=self.p, l=self.l, i=self.i, m=self.m, s=self.s,
                          study_settings=self.study_settings,
                          folder_path=self.folder_path,
                          parallel_env=self.parallel_env) )

class Sim_PL_NVT_posre_soft(Sim_PL_EM):
    stage="nvt_posre_soft"
    #request 2 cores
    n_cpu = luigi.IntParameter(visibility=ParameterVisibility.HIDDEN,
                               default=2, significant=False)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        if(self.restr_to_EM):
            self.posre = self.folder_path+"/state{s}/repeat{i}/em{m}/confout.gro".format(
                s=self.s, i=self.i, m=self.m)
        if(self.posre_ref_override):
            self.posre = self.posre_ref_override

        #override relevant file names
        self.mdp = self.study_settings['mdp_path'] +\
            "/protein/eq_nvt_posre_soft_{0}.mdp".format(self.study_settings['states'][self.s])
        self.struct = self.folder_path+"/state{2}/repeat{3}/nvt_posre{4}/confout.gro".format(
            self.p, self.l, self.s, self.i, self.m)

    def requires(self):
        return( Sim_PL_NVT_posre(p=self.p, l=self.l, i=self.i, m=self.m, s=self.s,
                          study_settings=self.study_settings,
                          folder_path=self.folder_path,
                          parallel_env=self.parallel_env) )

class Sim_PL_NPT(Sim_PL_EM):
    stage="npt"
    #request 4 cores
    n_cpu = luigi.IntParameter(visibility=ParameterVisibility.HIDDEN,
                               default=4, significant=False)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        #override relevant file names
        self.mdp = self.study_settings['mdp_path'] +\
            "/protein/eq_npt_{0}.mdp".format(self.study_settings['states'][self.s])
        self.struct = self.folder_path+"/state{2}/repeat{3}/nvt_posre_soft{4}/confout.gro".format(
            self.p, self.l, self.s, self.i, self.m)
        self.posre = None
        if(self.posre_ref_override):
            self.posre = self.posre_ref_override

    def requires(self):
        return( Sim_PL_NVT_posre_soft(p=self.p, l=self.l, i=self.i, m=self.m, s=self.s,
                          study_settings=self.study_settings,
                          folder_path=self.folder_path,
                          parallel_env=self.parallel_env) )
