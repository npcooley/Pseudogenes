universe = vanilla
executable = OSGWrapperA_v2.sh
log = LogFilesA/out.$(v01).log
output = LogFilesA/out.$(v01).out
error = LogFilesA/out.$(v01).err

# OSG has different requirement syntax than CHTC
# Kernel version comes from major.minor.patch in the format:
# major * 10000 + minor * 1000 + patch
# heron comes from r-base which uses debian:latest, which
# as of 20190722 is "buster" from 4.19.105
# as per advice, 31000 will keep jobs away from RHEL 6

requirements = Arch == "X86_64" && \
  HAS_SINGULARITY == True && \
  OSG_HOST_KERNEL_VERSION >= 31000 && \
  SINGULARITY_CAN_USE_SIF == False

request_cpus = 6
# ifThenElse(Expr, True, False)
# request 12GB initially, double every time the job restarts
request_memory = ifThenElse(MemoryUsage =!= undefined, (20000 * NumJobStarts * 2), 20000)
request_disk = ifThenElse(MemoryUsage =!= undefined, (40000 * NumJobStarts * 2), 40000)

+SingularityImage = "/cvmfs/singularity.opensciencegrid.org/npcooley/synextend:1.18.0"

# IF using a docker container, only send the RScript and any associated data files
transfer_input_files = JobScriptA_v2.R

# If just has been running for more than x hours - 3600 seconds to an hour! - put it on hold
periodic_hold = (JobStatus == 2) && ((CurrentTime - EnteredCurrentStatus) > 54000)

# register as completed
on_exit_hold = (ExitBySignal == True) || (ExitCode != 0) && \
               (NumJobStarts < 5)

# release jobs held for going over memory
periodic_release = (NumJobStarts < 5) && \
                    # ((HoldReasonCode == 26) && (HoldReasonSubCode == 0)) && \
                    ((CurrentTime - EnteredCurrentStatus > 180))
                    
periodic_remove = NumJobStarts >= 5

# release a job that has been held, after it has been held for 10 minutes, up to a maximum of 5 retries
# periodic_release = (NumJobStarts < 10) && \
#                    (((HoldReasonCode == 12) && (HoldReasonSubCode == 0 )) || (HoldReasonCode == 3)) && \
#                    ((CurrentTime - EnteredCurrentStatus) > 600)

# don't overwhelm the submit node
# max_idle = 1000
max_materialize = 10000

arguments = $(v01) $(v02) $(v03) $(v04) $(v05) $(v06) $(v07) $(v08) $(v09) $(v10)

# transfer_output_remaps = "Result*.RData = OutDir/Result$(Process).RData"
# This doesn't seem to work?
# Always test a small number of jobs first!
queue v01, v02, v03, v04, v05, v06, v07, v08, v09, v10 from TestMapA_v2.txt

