#!/usr/bin/env python

import glob
import luigi
import os
import sys
import pmx.scripts.analyze_dhdl
from luigi.parameter import ParameterVisibility
from pmx.scripts.workflows.SGE_tasks.SGETunedJobTask import SGETunedJobTask #tuned for the owl cluster
from pmx.scripts.workflows.SGE_tasks.absFE.LinP.TI import Task_PL_TI_simArray


# ==============================================================================
#                         Derivative Task Classes
# ==============================================================================
class Task_PL_analysis_aligned(SGETunedJobTask):

    #Parameters:
    p = luigi.Parameter(description='Protein name')
    l = luigi.Parameter(description='Ligand name')
    i = luigi.IntParameter(description='Repeat number')

    folder_path = luigi.Parameter(significant=False,
        visibility=ParameterVisibility.HIDDEN,
        description='Path to the protein+ligand folder to set up')

    study_settings = luigi.DictParameter(significant=True,
        visibility=ParameterVisibility.HIDDEN,
        description='Dict of study stettings '
        'used to propagate settings to dependencies')

    #request 1 core
    n_cpu = luigi.IntParameter(visibility=ParameterVisibility.HIDDEN,
                               default=1, significant=False)
                               
    job_name_format = luigi.Parameter(
        visibility=ParameterVisibility.HIDDEN,
        significant=False, default="pmx_{task_family}_p{p}_l{l}_{i}",
        description="A string that can be "
        "formatted with class variables to name the job with qsub.")
    
    n_bins = luigi.IntParameter(visibility=ParameterVisibility.HIDDEN,
                               default=10, significant=False,
                               description="Number of histogram bins in plot.")
                               
    n_bootstrap = luigi.IntParameter(default=0, significant=True,
                               description="Number of times estimators are bootstrapped.")
                               
    T = luigi.FloatParameter(default=298.0, significant=True,
        description="Simulation temperature.")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        #set variables
        self.base_path = self.study_settings['base_path']
        self.ana_folder=self.folder_path+"/analysis/repeat%d"%self.i
        self.TIpath=self.folder_path+'/state{s}/repeat{i}/morphes{m}/dHdl*.xvg'


    def work(self):
        dHdlA=[]
        dHdlB=[]

        os.makedirs(self.ana_folder, exist_ok=True)
        os.chdir(self.ana_folder)

        TIstate_list=[*self.study_settings['TIstates']]
        for m in range(self.study_settings['n_sampling_sims']):
            dHdlA.extend(glob.glob(self.TIpath.format(
                            s=TIstate_list[0], i=self.i, m=m)))
            dHdlB.extend(glob.glob(self.TIpath.format(
                            s=TIstate_list[1], i=self.i, m=m)))

        #set script params and call analyze_dhdl
        orig_argv = sys.argv
        orig_stdout=sys.stdout
        orig_stderr=sys.stderr
        sys.argv = [['analyze_dhdl.py'],
                    ['-fA'], dHdlA,
                    ['-fB'], dHdlB,
                    ['--nbins', str(int(self.n_bins))],
                    ['-b', str(int(self.n_bootstrap))],
                    ['-t', str(float(self.T))]]
        sys.argv = [item for sublist in sys.argv for item in sublist] #flatten argv

        with open("analysis.log", "w") as f:
            sys.stdout = f
            sys.stderr = f
            pmx.scripts.analyze_dhdl.entry_point()

        sys.argv = orig_argv #reset argv
        sys.stdout=orig_stdout #reset output
        sys.stderr=orig_stderr

        os.chdir(self.base_path)#reset cwd
        pass

    def output(self):
        files=['wplot.png', 'analysis.log', 'integA.dat', 'integB.dat', 'results.txt']
        return([luigi.LocalTarget(os.path.join(self.ana_folder, f)) for f in files])
        
    def complete(self):
        """
        Check if Bootstraping was done (if required).
        """
        c=super().complete()
        if(c and self.n_bootstrap>0):
            found_boots=0
            with open(os.path.join(self.ana_folder, 'analysis.log'), "r") as anaf:
                content = anaf.readlines()
                for l in content:
                    if("Bootstrap (Std Err): iteration" in l):
                        s=l.split('/')
                        found_boots=int(s[-1])
                        break;
            if(found_boots!=self.n_bootstrap):
                c=False
            
        return(c)

    def requires(self):
        tasks=[]
        for sTI in self.study_settings['TIstates']:
            for m in range(self.study_settings['n_sampling_sims']):
                tasks.append(Task_PL_TI_simArray(
                      p=self.p, l=self.l, i=self.i, m=m, sTI=sTI,
                      study_settings=self.study_settings,
                      folder_path=self.folder_path,
                      parallel_env=self.parallel_env,
                      restr_scheme="Aligned") )

        return(tasks)


class Task_PL_analysis_aligned2crystal(Task_PL_analysis_aligned):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        #overwrite variables
        self.ana_folder=self.folder_path+"/analysis/aligned2crystal_repeat%d"%self.i
        self.TIpath=self.folder_path+'/state{s}/repeat{i}/aligned2crystal_morphes{m}/dHdl*.xvg'


    def requires(self):
        tasks=[]
        for sTI in self.study_settings['TIstates']:
            for m in range(self.study_settings['n_sampling_sims']):
                tasks.append(Task_PL_TI_simArray(
                      p=self.p, l=self.l, i=self.i, m=m, sTI=sTI,
                      study_settings=self.study_settings,
                      folder_path=self.folder_path,
                      parallel_env=self.parallel_env,
                      restr_scheme="Aligned_crystal") )

        return(tasks)
