import errno
import logging
import luigi
import os
import shutil
import subprocess
import tempfile
import pmx.scripts.workflows.SGE_tasks.SGETunedRunner as sge_runner
from luigi.contrib.sge import logger
from pmx.scripts.workflows.SGE_tasks.SGETunedJobTask import SGETunedJobTask #tuned for the owl cluster



def _parse_qsub_job_array_id(qsub_out):
    """Parse job id from qsub output string.

    Assume format:

        "Your job-array <job_id>.<min>-<max>:<step> ("<job_name>") has been submitted"

    """
    return int(qsub_out.split()[2].split('.')[0])

def _build_qsub_array_command(cmd, job_name, pe, n_cpu, work_dir, nsubtasks, runtime=None, extra_options=""):
    """Submit array job shell command to SGE queue via `qsub`"""
    h_rt=""
    if(runtime):
        h_rt="-l h_rt="+runtime
    qsub_template = """echo {cmd} | qsub -t 1:{nsubtasks} -V -r y {h_rt} -pe {pe} {n_cpu} -N {job_name} -wd {work_dir} {extra_options}"""
    return qsub_template.format(
        cmd=cmd, job_name=job_name, work_dir=work_dir,
        pe=pe, n_cpu=n_cpu, h_rt=h_rt, nsubtasks=nsubtasks, extra_options=extra_options)



class SGETunedArrayJobTask(SGETunedJobTask):

    unfinished=[] #set this in __init__()

    def _find_unfinished(self):
        """Finds the list of unfinished jobs in the array
        which need to be (re)run.

        Overload this in subclasses.
        Will be called in self._init_local() before pickling
        an instance of this class for execution on the worker nodes
        so the nodes can map SGE_TASK_ID to the correct job.

        Parameters
        ----------
        None.

        Returns
        -------
        List of unfinished jobs in the array.
        """
        return([])

    def _init_local(self):
        #find previously unfinished tasks before pickling this instance
        self.unfinished=self._find_unfinished()
        #pickles the instance
        super()._init_local()


    def run(self):
        #check if we want this job to fail (for debug purposes)
        if(self.debug_exit):
            raise(Exception("Exiting for debug reasons. We don't want this job to actually execute."))
    
        if self.run_locally:
            self.unfinished=self._find_unfinished()
            #loop over each instance and run them separately
            for t in range(len(self.unfinished)):
                with tempfile.TemporaryDirectory() as dirpath:
                    os.environ['SGE_TASK_ID']="%d"%(t+1)
                    os.environ['TMPDIR']=dirpath
                    self.work()
                    os.environ['SGE_TASK_ID']=""
                    os.environ['TMPDIR']=""
        else:
            self._init_local()
            self._run_job()

    def _run_job(self):

        # Build a qsub argument that will run sge_runner.py on the directory we've specified
        runner_path = sge_runner.__file__
        if runner_path.endswith("pyc"):
            runner_path = runner_path[:-3] + "py"
        job_str = 'python {0} "{1}" "{2}"'.format(
            runner_path, self.tmp_dir, os.getcwd())  # enclose tmp_dir in quotes to protect from special escape chars
        if self.no_tarball:
            job_str += ' --no-tarball'

            # #force loading of dependencies by sourcing a custom profile
            # if(os.path.isfile(self.source_conda)):
            #     job_str = '"source {}; '.format(self.source_conda) + job_str+'"'
            # else:
            #     mylogger = logging.getLogger(self.__class__.__name__)
            #     mylogger.error("Tarballing of dependencies is disabled and "
            #                   "{} does not exist. "
            #                   "Will not be able to load all workflow "
            #                   "dependencies without it. Please create it and "
            #                   "within activate a conda environment containing "
            #                   "at least python>3.6, "
            #                   "pmx, luigi, MDanalysis, matplotlib, and numpy."
            #                   "".format(self.source_conda))
            #     raise Exception("Could not source " + self.source_conda)

        #tell runner that this is an array job
        job_str += ' --arrayjob'

        self.errfile=""; #no errorfile. mdrun dumps too much into stderr.
        #Let SGE assign a separate one for each job in array

        # Build qsub submit command
        submit_cmd = _build_qsub_array_command(job_str, self.job_name,
                                               self.parallel_env, self.n_cpu,
                                               self.sim_path,
                                               len(self.unfinished), self.runtime,
                                               self.qsub_extra_options)
        logger.debug('qsub command: \n' + submit_cmd)

        # Submit the job and grab job ID
        output = subprocess.check_output(submit_cmd, shell=True).decode('utf-8')
        logger.debug("Submitted job to qsub with response:\n" + output)
        output=output.strip().split("\n")[-1] # remove AWS SGE warning about LD_LIBRARY_PATH
        self.job_id = _parse_qsub_job_array_id(output)
        #logger.debug("Submitted job to qsub with response:\n" + output)

        self._track_job()

        # Now delete the temporaries, if they're there.
        if (self.tmp_dir and os.path.exists(self.tmp_dir) and not self.dont_remove_tmp_dir):
            logger.info('Removing temporary directory %s' % self.tmp_dir)
            subprocess.call(["rm", "-rf", self.tmp_dir])
