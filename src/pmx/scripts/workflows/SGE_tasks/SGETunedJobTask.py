import logging
import luigi
import os
import pmx
import random
import subprocess
import sys
import time
import pmx.scripts.workflows.SGE_tasks.SGETunedRunner as sge_runner
from luigi.parameter import ParameterVisibility
from luigi.contrib.hadoop import create_packages_archive
from luigi.contrib.sge import SGEJobTask, logger, _parse_qstat_state, _parse_qsub_job_id
try:
    import cPickle as pickle
except ImportError:
    import pickle

def _parse_qstat_state_opt_header(qstat_out, job_id, header=True):
    """Parse "state" column from `qstat` output for given job_id

    Returns state for the *first* job matching job_id. Returns 'u' if
    `qstat` output is empty or job_id is not found.

    If header=False, can also be used with qstat | grep jobid

    """
    if qstat_out.strip() == '':
        return 'u'
    lines = qstat_out.split('\n')
    if(header):
        # skip past header
        while not lines.pop(0).startswith('---'):
            pass
    for line in lines:
        if line:
            job, prior, name, user, state = line.strip().split()[0:5]
            if int(job) == int(job_id):
                return state
    return 'u'

def extended_build_qsub_command(cmd, job_name, outfile, errfile, pe, n_cpu, runtime=None, extra_options="", priority=0):
    """Submit shell command to SGE queue via `qsub`"""
    h_rt=""
    if(runtime):
        h_rt="-l h_rt="+runtime
    qsub_template = """echo {cmd} | qsub -o ":{outfile}" -e ":{errfile}" -V -r y {h_rt} -pe {pe} {n_cpu} -N {job_name} {extra_options}"""
    #handle priority
    priority=int(priority)
    if(priority>0):
        print("WARNING: qsub will not allow a positive priority from a non-admin account. Setting priority to 0.")
        priority=0
    if(priority<0):
        qsub_template += f" -p {priority}"
    #handle extra options
    qsub_template += f" {extra_options}"
    
    return qsub_template.format(
        cmd=cmd, job_name=job_name, outfile=outfile, errfile=errfile,
        pe=pe, n_cpu=n_cpu, h_rt=h_rt, extra_options=extra_options)


class SGETunedJobTask(SGEJobTask):
    run_locally = luigi.BoolParameter(
        default=False,
        significant=False,
        description="run locally instead of on the cluster")

    #temp files
    #Don't archive luigi or pmx. Jobs load them by copying env vars from login node
    no_tarball = luigi.BoolParameter(
        visibility=ParameterVisibility.HIDDEN,
        significant=False,
        default=True,
        description="don't tarball (and extract) the luigi project files")

    shared_tmp_dir = luigi.Parameter(
        visibility=ParameterVisibility.HIDDEN,
        default=os.path.join(os.getenv("HOME"), 'temp'),
        significant=False)

    dont_remove_tmp_dir = luigi.BoolParameter(
        visibility=ParameterVisibility.HIDDEN,
        significant=False,
        default=True,
        description="don't delete the temporary directory used (for debugging)")

    #change default parallel environment
    parallel_env = luigi.Parameter(
        default='*fast',
        visibility=ParameterVisibility.HIDDEN,
        significant=False)
    #poll time in seconds.
    #Large value is better to avoid overloading the login node.
    #Needs to be less than time between queue updates.
    poll_time = luigi.IntParameter(
        significant=False, default=600,
        visibility=ParameterVisibility.HIDDEN,
        description="specify the wait time to poll qstat for the job status")

    source_conda = luigi.Parameter(significant=False,
        default=os.path.join(os.getenv("HOME"),".luigi_profile"),
        visibility=ParameterVisibility.HIDDEN,
        description="File to source to load an appropriate conda environment "
        "in case tarballing of required python packages is disabled "
        "for performance.")

    #avoid Prameter not a string warnings
    job_name_format = luigi.Parameter(
        significant=False,
        visibility=ParameterVisibility.HIDDEN,
        default="pmx_{task_family}", description="A string that can be "
        "formatted with class variables to name the job with qsub.")
    job_name = luigi.Parameter(
        significant=False, default="",
        visibility=ParameterVisibility.HIDDEN,
        description="Explicit job name given via qsub.")
        
    qsub_extra_options = luigi.Parameter(
        significant=False, default="",
        visibility=ParameterVisibility.HIDDEN,
        description="Extra options passed to qsub.")

    #set runtime
    runtime = luigi.Parameter(
        significant=False, default="24:00:00",
        visibility=ParameterVisibility.HIDDEN,
        description="Hard realtime limit SGE will let jobs from this task run for. "
        "Format: hh:mm:ss")

    end_wait_buffer = luigi.IntParameter(
        significant=False, default=60,
        visibility=ParameterVisibility.HIDDEN,
        description="Number of seconds to wait for the (remote) file system "
        "to sync output files at end of task before declaring task as failed.")

    debug_exit = luigi.BoolParameter(
        visibility=ParameterVisibility.HIDDEN,
        significant=False,
        default=False,
        description="Don't execute. Just raise an exception before submitting to queue.")
    
    priority = luigi.IntParameter(
        significant=False, default=0,
        visibility=ParameterVisibility.HIDDEN,
        description="Job priority to pass to qsub. Can be 0 or negative. Used to lower priority of certain tasks in the SGE queue.")

    extra_packages=[] #extra packages to be tarballed. Overloaded by subclasses.

    #override scheduler based retry policy
    _disable_window_seconds = luigi.IntParameter(
        significant=False, default=3600*24*7, # 7 days
        visibility=ParameterVisibility.HIDDEN,
        description="Number of seconds in which to count failures.")
    _retry_count = luigi.IntParameter(
        significant=False, default=0,
        visibility=ParameterVisibility.HIDDEN,
        description="Number of failures to allow in the last _disable_window_seconds.")

    @property
    def retry_count(self):
        """
        Override this positive integer to have different ``retry_count`` at task level
        Check :ref:`scheduler-config`
        """
        return self._retry_count

    @property
    def disable_window_seconds(self):
        """
        Override this positive integer to have different ``disable_window_seconds`` at task level.
        Check :ref:`scheduler-config`
        """
        return self._disable_window_seconds

    def run(self):
        """
        Sometimes writing to the network drives isn't finished before the next
        task starts checking for unfulfilled dependencies. In such cases we
        should wait for the file system to sync before continuing.
        """

        #check if we want this job to fail (for debug purposes)
        if(self.debug_exit):
            raise(Exception("Exiting for debug reasons. We don't want this job to actually execute."))

        super().run()

        #make sure output files actually exist and are acessible for next tasks
        waitbuffer = self.end_wait_buffer;
        while(waitbuffer>0 and not self.complete()):
            time.sleep(1)
            waitbuffer-=1

        if(not self.complete()):
            errmsg = "Task "+self.__class__.__name__
            pars = [(t[0], getattr(self, t[0])) for t in self.get_params() if t[1].significant]
            errmsg+= str(dict((x, y) for x, y in pars))
            errmsg+= " Failed due to being incomplete."
            raise(Exception(errmsg))


    def _init_local(self):

        # Set up temp folder in shared directory (trim to max filename length)
        base_tmp_dir = self.shared_tmp_dir
        random_id = '%016x' % random.getrandbits(64)
        folder_name = self.task_id + '-' + random_id
        self.tmp_dir = os.path.join(base_tmp_dir, folder_name)
        max_filename_length = os.fstatvfs(0).f_namemax
        self.tmp_dir = self.tmp_dir[:max_filename_length]
        logger.info("Tmp dir: %s", self.tmp_dir)
        os.makedirs(self.tmp_dir)

        # Dump the code to be run into a pickle file
        logging.debug("Dumping pickled class")
        self._dump(self.tmp_dir)

        if not self.no_tarball:
            # Make sure that all the class's dependencies are tarred and available
            # This is not necessary if luigi is importable from the cluster node
            logging.debug("Tarballing dependencies")
            # Grab luigi, the whole of pmx, and the module containing the code to be run
            packages = [luigi] + [pmx] + self.extra_packages +\
                [__import__(self.__module__, None, None, 'dummy')]
            create_packages_archive(packages, os.path.join(self.tmp_dir, "packages.tar"))

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

        # Build qsub submit command
        self.outfile = os.path.join(self.tmp_dir, 'job.out')
        self.errfile = os.path.join(self.tmp_dir, 'job.err')
        submit_cmd = extended_build_qsub_command(job_str, self.job_name,
                             self.outfile, self.errfile, self.parallel_env,
                             self.n_cpu, self.runtime, self.qsub_extra_options, self.priority)
        logger.debug('qsub command: \n' + submit_cmd)

        # Submit the job and grab job ID
        output = subprocess.check_output(submit_cmd, shell=True).decode('utf-8')
        logger.debug("Submitted job to qsub with response:\n" + output)
        output=output.strip().split("\n")[-1] # remove AWS SGE warning about LD_LIBRARY_PATH
        self.job_id = _parse_qsub_job_id(output)

        self._track_job()

        # Now delete the temporaries, if they're there.
        if (self.tmp_dir and os.path.exists(self.tmp_dir) and not self.dont_remove_tmp_dir):
            logger.info('Removing temporary directory %s' % self.tmp_dir)
            subprocess.call(["rm", "-rf", self.tmp_dir])

    def _track_job(self):
        while True:
            # Sleep for a little bit
            time.sleep(self.poll_time)

            # See what the job's up to
            # ASSUMPTION
            #qstat_out = subprocess.check_output(['qstat']).decode('utf-8')
            #sge_status = _parse_qstat_state(qstat_out, self.job_id)
            qstat_out = subprocess.run("qstat | grep {}".format(self.job_id), shell=True, check=False, stdout=subprocess.PIPE).stdout.decode('utf-8')
            sge_status = _parse_qstat_state_opt_header(qstat_out, self.job_id, header=False)
            if sge_status == 'r' or sge_status == 'Rr' :
                logger.info('Job is running...')
            elif sge_status == 'qw':
                logger.info('Job is pending...')
            elif sge_status == 't' or sge_status == 'Rt':
                logger.info('Job is transferring...')
            #elif 'E' in sge_status:
                #logger.error('Job has FAILED:\n' + '\n'.join(self._fetch_task_failures()))
                #break
            elif 'd' in sge_status:
                logger.error('Job has been scheduled for deletion...')
            elif sge_status == 'u':
                # Then the job could either be failed or done.
                errors = self._fetch_task_failures()
                if not errors:
                    logger.info('Job is done')
                else:
                    logger.error('Job has FAILED:\n' + '\n'.join(errors))
                break
            else:
                logger.info('Job status is UNKNOWN!')
                logger.info('Status is : %s' % sge_status)
                #raise Exception("job status isn't one of ['r', 'qw', 'E*', 't', 'u', 'Rr', 'Rt', 'd*']: %s" % sge_status)
                raise Exception("job status isn't one of ['r', 'qw', 't', 'u', 'Rr', 'Rt', 'd*']: %s" % sge_status)

    def _dump(self, out_dir=''):
        """Dump instance to file."""
        with self.no_unpicklable_properties():
            self.job_file = os.path.join(out_dir, 'job-instance.pickle')
            if self.__module__ == '__main__':
                d = pickle.dumps(self)
                module_name = os.path.basename(sys.argv[0]).rsplit('.', 1)[0]
                d = d.replace(b'c__main__', b"c" + module_name.encode('utf-8'))
                open(self.job_file, "wb").write(d)
            else:
                pickle.dump(self, open(self.job_file, "wb"))



#     """A local version of SGETunedJobTask, for easier debugging.

#     This version skips the ``qsub`` steps and simply runs ``work()``
#     on the local node, so you don't need to be on an SGE cluster to
#     use your Task in a test workflow.

#     Copied from luigi.contrib.sge.SGELocalJobTask to allow for consistent
#     retry policy management in the current file.
#     """
#     def run(self):
#         self.work()
