#!/usr/bin/env python

import luigi
from luigi.parameter import ParameterVisibility
from pmx.scripts.workflows.SGE_tasks.absFE.LinW.prep_folders import Prep_WL_folder
from pmx.scripts.workflows.SGE_tasks.absFE.LinP.equil_sims import Sim_PL_EM, Sim_PL_NVT_posre, Sim_PL_NPT

# ==============================================================================
#                         Derivative Task Classes
# ==============================================================================
class Sim_WL_EM(Sim_PL_EM):

    #Parameters:
    p = None #disables base class' p

    folder_path = luigi.Parameter(significant=False,
                 visibility=ParameterVisibility.HIDDEN,
                 description='Path to the protein+ligand folder to set up')

    #request 1 cores
    n_cpu = luigi.IntParameter(visibility=ParameterVisibility.HIDDEN,
                               default=1, significant=False)

    job_name_format = luigi.Parameter(
        visibility=ParameterVisibility.HIDDEN,
        significant=False, default="pmx_{task_family}_l{l}_{s}{i}_{m}",
        description="A string that can be "
        "formatted with class variables to name the job with qsub.")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        #override relevant file names
        self.mdp = self.study_settings['mdp_path'] +\
            "/water/em_{0}.mdp".format(self.study_settings['states'][self.s])
        self.posre = None
        if(self.posre_ref_override):
            self.posre = self.posre_ref_override

    def requires(self):
        return( Prep_WL_folder(l=self.l,
                               study_settings=self.study_settings,
                               folder_path=self.folder_path) )
                                #no need to pass parallel_env as
                                #Prep_WL_folder runs on the login node


class Sim_WL_NVT(Sim_PL_NVT_posre):
    stage="nvt"

    #Parameters:
    p = None #disables base class' p

    #request 1 cores
    n_cpu = luigi.IntParameter(visibility=ParameterVisibility.HIDDEN,
                               default=1, significant=False)

    job_name_format = luigi.Parameter(
        visibility=ParameterVisibility.HIDDEN,
        significant=False, default="pmx_{task_family}_l{l}_{s}{i}_{m}",
        description="A string that can be "
        "formatted with class variables to name the job with qsub.")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        #override relevant file names
        self.mdp = self.study_settings['mdp_path'] +\
            "/water/eq_nvt_{0}.mdp".format(self.study_settings['states'][self.s])
        self.posre = None
        if(self.posre_ref_override):
            self.posre = self.posre_ref_override

    def requires(self):
        return( Sim_WL_EM(l=self.l, i=self.i, m=self.m, s=self.s,
                          study_settings=self.study_settings,
                          folder_path=self.folder_path,
                          parallel_env=self.parallel_env) )


class Sim_WL_NPT(Sim_PL_NPT):
    stage="npt"

    #Parameters:
    p = None #disables base class' p

    #request 1 core
    n_cpu = luigi.IntParameter(visibility=ParameterVisibility.HIDDEN,
                               default=1, significant=False)

    job_name_format = luigi.Parameter(
        visibility=ParameterVisibility.HIDDEN,
        significant=False, default="pmx_{task_family}_l{l}_{s}{i}_{m}",
        description="A string that can be "
        "formatted with class variables to name the job with qsub.")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        #override relevant file names
        self.mdp = self.study_settings['mdp_path'] +\
            "/water/eq_npt_{0}.mdp".format(self.study_settings['states'][self.s])
        self.struct = self.folder_path+"/state{2}/repeat{3}/nvt{4}/confout.gro".format(
            self.p, self.l, self.s, self.i, self.m)
        self.posre = None
        if(self.posre_ref_override):
            self.posre = self.posre_ref_override

    def requires(self):
        return( Sim_WL_NVT(l=self.l, i=self.i, m=self.m, s=self.s,
                          study_settings=self.study_settings,
                          folder_path=self.folder_path,
                          parallel_env=self.parallel_env) )
