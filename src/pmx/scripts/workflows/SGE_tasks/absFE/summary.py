#!/usr/bin/env python

import copy
import luigi
import os
import numpy as np
#import matplotlib as plt
from luigi.parameter import ParameterVisibility
from pmx.scripts.workflows.SGE_tasks.SGETunedJobTask import SGETunedJobTask #tuned for the owl cluster
from pmx.scripts.workflows.SGE_tasks.absFE.LinP.analysis import Task_PL_analysis_aligned,Task_PL_analysis_aligned2crystal
from pmx.scripts.workflows.SGE_tasks.absFE.LinW.analysis import Task_WL_analysis
from pmx.scripts.workflows.SGE_tasks.absFE.LinP.decorrelate_energetically import Task_PL_energetic_decorrelation


# ==============================================================================
#                         Derivative Task Classes
# ==============================================================================
class Task_summary_aligned(SGETunedJobTask):
    #run on the login node
    run_locally = luigi.BoolParameter(
        default = True,
        significant=False,
        parsing=luigi.BoolParameter.EXPLICIT_PARSING,
        description="run locally instead of on the cluster")

    #Parameters:
    hosts = luigi.ListParameter(description='list of protein names to evaluate')
    ligands = luigi.ListParameter(description='list of ligand names to evaluate')
    n_repeats = luigi.IntParameter(description='Number of repeats',
                                   default=3)
    n_sampling_sims = luigi.IntParameter(description='Number of sampling simulations',
                                         default=1)

    show_incomplete = luigi.BoolParameter(default=False, significant=True,
                               visibility=ParameterVisibility.HIDDEN)

    show_sep_repeats = luigi.BoolParameter(default=False, significant=True,
                               visibility=ParameterVisibility.HIDDEN)

    only_LinW = luigi.BoolParameter(default=False, significant=True,
                               visibility=ParameterVisibility.HIDDEN)

    decor_decoupled = luigi.Parameter(default="energetic", significant=True)

    #TODO: add default
    study_settings = luigi.DictParameter(significant=False,
        visibility=ParameterVisibility.HIDDEN,
        description='Dict of study stettings '
        'used to propagate settings to dependencies')

    #request 1 core
    n_cpu = luigi.IntParameter(default=1, significant=False,
                               visibility=ParameterVisibility.HIDDEN)

    #avoid Prameter not a string warnings
    job_name_format = luigi.Parameter(
        visibility=ParameterVisibility.HIDDEN,
        significant=False, default="pmx_{task_family}",
        description="A string that can be "
        "formatted with class variables to name the job with qsub.")


    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)



        self.WL_settings=copy.deepcopy(self.study_settings.get_wrapped())
        self.WL_settings['TIstates']=self.WL_settings['states']
        self.WL_settings['n_repeats']=self.n_repeats
        self.WL_settings['n_sampling_sims']=self.n_sampling_sims
        self.PL_settings=copy.deepcopy(self.study_settings.get_wrapped())
        self.PL_settings['n_repeats']=self.n_repeats
        self.PL_settings['n_sampling_sims']=self.n_sampling_sims

        if(self.decor_decoupled.casefold()=="false" or self.decor_decoupled.casefold()=="no"):
            self.PL_settings['decor_decoupled']=False
            self.PL_settings['decor_method']="none"
        elif(self.decor_decoupled.casefold()=="md"):
            self.PL_settings['decor_decoupled']=True
            self.PL_settings['decor_method']="md"
        elif(self.decor_decoupled.casefold()=="gaussian" or self.decor_decoupled.casefold()=="sampling"):
            self.PL_settings['decor_decoupled']=True
            self.PL_settings['decor_method']="sampling"
        elif(self.decor_decoupled.casefold()=="energetic" or self.decor_decoupled.casefold()=="MCI"):
            self.PL_settings['decor_decoupled']=True
            self.PL_settings['decor_method']="MCI"
        else:
            raise(Exception("Invalid value for decor_decoupled detected (%s)"%self.decor_decoupled))


        self.base_path = self.study_settings['base_path']
        self.outname="summary_aligned.txt"
        if(self.show_incomplete):
            self.outname="summary_aligned_incomplete.txt"
        if(self.PL_settings['decor_method']=="MCI"):
            self.outname=self.outname[:-4]+"_energ_decor"+self.outname[-4:]
        self.restrname="out_dg_{i}.dat"

        self.anafolderfmt_P="/analysis/repeat{i}"
        if(self.PL_settings['decor_method']=="MCI"):
            self.anafolderfmt_P="/decor_analysis/repeat{i}"
        self.anafolderfmt_W="/analysis/repeat{i}"

    def work(self):

        def read_results(inP=False):
            rs=np.ndarray(self.n_repeats)
            nfound=0
            for i in range(self.n_repeats):
                if(inP):
                    ana_folder=folder_path+self.anafolderfmt_P.format(i=i)
                else:
                    ana_folder=folder_path+self.anafolderfmt_W.format(i=i)

                #discard missing pieces
                if(self.show_incomplete):
                    if(not os.path.isfile(ana_folder+"/results.txt") or (inP and not os.path.isfile(folder_path+"/"+self.restrname.format(i=i)))):
                        continue;
                with open(ana_folder+"/results.txt", 'r') as f:
                    for line in f:
                        if "BAR: dG" in line:
                            s = line.split()
                            rs[nfound]=float(s[3])
                            break

                if(inP): #analytical correction for this repeat
                    with open(folder_path+"/"+self.restrname.format(i=i), 'r') as f:
                        for line in f:
                            if("Restraint contribution to free energy" in line and
                               "kJ/mol" in line):
                                s = line.split()
                                rs[nfound]+=float(s[-2])
                                break
                nfound+=1

            if(nfound>0):
                dGpart = np.mean(rs[:nfound])
                std = np.std(rs[:nfound])
                return([dGpart,std,nfound], rs)
            else:
                return([np.nan,np.nan,nfound], rs)


        #dG in water
        inws={}
        p="water"
        w_reps={}
        for l in self.ligands:
            key=l
            folder_path = self.base_path+'/'+p+'/lig_'+l
            tup, rs=read_results(inP=False)
            inws.update({key:tup})
            w_reps.update({key:rs})

        #dG in protein
        inps={}
        p_reps={}
        for p in self.hosts:
            for l in self.ligands:
                key=p+' '+l
                if(not self.only_LinW):
                    folder_path = self.base_path+'/prot_'+p+'/lig_'+l
                    tup, rs=read_results(inP=True)
                    inps.update({key:tup})
                    p_reps.update({key:rs})
                else:
                    inps.update({key:[np.nan,np.nan,0]})

        #read analytical corrections for restraints
        anacorrs={}
        corrs_reps={}
        for p in self.hosts:
            for l in self.ligands:
                key=p+' '+l
                if(self.only_LinW):
                    anacorrs.update({key:[np.nan, np.nan, 0]})
                    continue;
                folder_path = self.base_path+'/prot_'+p+'/lig_'+l
                cors=np.zeros(self.n_repeats)
                nfound=0
                for i in range(self.n_repeats):
                    if(self.show_incomplete): #skip missing
                        if(not os.path.isfile(folder_path+"/"+self.restrname.format(i=i))):
                            continue;
                    with open(folder_path+"/"+self.restrname.format(i=i), 'r') as f:
                        for line in f:
                            if("Restraint contribution to free energy " in line and
                               "kJ/mol" in line):
                                s=line.split()
                                cors[nfound]=float(s[-2])
                                nfound+=1

                if(nfound>0):
                    anacorrs.update({key:[np.mean(cors[:nfound]), np.std(cors[:nfound]), nfound]})
                else:
                    anacorrs.update({key:[np.nan, np.nan, nfound]})

                corrs_reps.update({key:cors})


        #print summary table
        if(self.show_sep_repeats):
            rep_file=open("sep_repeats_"+self.outname, 'w')

        with open(self.outname, 'w') as sf:

            print("{:^20s} \t{:^20s}   {:^20s}   {:^20s}   {:^20s}".format(
                            "host guest","ddG (kJ/mol)","dG in prot+restr" ,
                            "dG in water", "ana. corr."))
            sf.write("{:^20s} \t{:^20s}   {:^20s}   {:^20s}   {:^20s}\n".format(
                            "host guest","ddG (kJ/mol)","dG in prot+restr" ,
                            "dG in water", "ana. corr."))
            for p in self.hosts:
                for l in self.ligands:
                    key=p+' '+l
                    if(self.show_sep_repeats):
                        rep_file.write("{:<20s}:\t".format(key))

                    if(np.isfinite(inws[l][0]) and np.isfinite(inps[key][0]) and np.isfinite(anacorrs[key][0]) or
                       (np.isfinite(inws[l][0]) and self.only_LinW) ):
                        ddG = inws[l][0] - inps[key][0] #water - (protein + restr corr.)
                        sigma = np.sqrt(inps[key][1]**2 + inws[l][1]**2) #standard dev.
                        print("{:<20s}:\t{:>8.2f} +- {:<8.2f}   {:>8.2f} +- {:<8.2f}   {:>8.2f} +- {:<8.2f}   {:>8.2f} +- {:<8.2f}".format(
                            key, ddG, sigma,
                            inps[key][0], inps[key][1],
                            inws[l][0], inws[l][1],
                            anacorrs[key][0], anacorrs[key][1]) )
                        sf.write("{:<20s}:\t{:>8.2f} +- {:<8.2f}   {:>8.2f} +- {:<8.2f}   {:>8.2f} +- {:<8.2f}   {:>8.2f} +- {:<8.2f}\n".format(
                            key, ddG, sigma,
                            inps[key][0], inps[key][1],
                            inws[l][0], inws[l][1],
                            anacorrs[key][0], anacorrs[key][1]) )

                        if(self.show_sep_repeats):
                            for r in range(self.n_repeats):
                                dG = w_reps[l][r] - p_reps[key][r]
                                rep_file.write("{:>8.2f}:\t".format(dG))
                            rep_file.write("\n")

        if(self.show_sep_repeats):
            rep_file.close()

    def output(self):
        files=[self.outname]
        if(self.show_sep_repeats):
            files.append("sep_repeats_"+self.outname)
        return([luigi.LocalTarget(os.path.join(self.base_path, f)) for f in files])

    def requires(self):
        tasks=[]
        if(not self.show_incomplete):
            #Ligand in Water
            p="water"
            for l in self.ligands:
                folder_path = self.base_path+'/'+p+'/lig_'+l
                for sTI in self.WL_settings['states']: #uses equil states for TI
                    for i in range(self.WL_settings['n_repeats']):
                        tasks.append(Task_WL_analysis(
                            l = l, i = i,
                            study_settings = self.WL_settings,
                            folder_path = folder_path,
                            parallel_env=self.parallel_env))

            #Ligand in Protein
            if(not self.only_LinW):
                for p in self.hosts:
                    for l in self.ligands:
                        folder_path = self.base_path+'/prot_'+p+'/lig_'+l
                        for sTI in self.PL_settings['TIstates']:
                            for i in range(self.PL_settings['n_repeats']):
                                if(self.PL_settings['decor_method']=="MCI"):
                                    tasks.append(Task_PL_energetic_decorrelation(
                                        p = p, l = l, i = i,
                                        study_settings = self.PL_settings,
                                        folder_path = folder_path,
                                        parallel_env=self.parallel_env))
                                else:
                                    tasks.append(Task_PL_analysis_aligned(
                                        p = p, l = l, i = i,
                                        study_settings = self.PL_settings,
                                        folder_path = folder_path,
                                        parallel_env=self.parallel_env))

        return(tasks)

class Task_summary_aligned2crystal(Task_summary_aligned):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        #overwrite values
        self.outname="summary_aligned2crystal.txt"
        self.restrname="out_dg_aligned2crystal_{i}.dat"
        self.anafolderfmt_P="/analysis/aligned2crystal_repeat{i}"

    def requires(self):
        tasks=[]

        if(not self.show_incomplete):
            #Ligand in Water
            p="water"
            for l in self.ligands:
                folder_path = self.base_path+'/'+p+'/lig_'+l
                for sTI in self.WL_settings['states']: #uses equil states for TI
                    for i in range(self.WL_settings['n_repeats']):
                        tasks.append(Task_WL_analysis(
                            l = l, i = i,
                            study_settings = self.WL_settings,
                            folder_path = folder_path,
                            parallel_env=self.parallel_env))

            #Ligand in Protein
            for p in self.hosts:
                for l in self.ligands:
                    folder_path = self.base_path+'/prot_'+p+'/lig_'+l
                    for sTI in self.PL_settings['TIstates']:
                        for i in range(self.PL_settings['n_repeats']):
                            tasks.append(Task_PL_analysis_aligned2crystal(
                                p = p, l = l, i = i,
                                study_settings = self.PL_settings,
                                folder_path = folder_path,
                                parallel_env=self.parallel_env))

        return(tasks)
