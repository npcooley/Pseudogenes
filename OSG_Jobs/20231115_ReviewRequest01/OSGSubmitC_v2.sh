universe = vanilla
executable = OSGWrapperC.sh
log = LogFilesC/out.$(v01).log
# PGAP err and output files are too large to collect at scale once this
# job is debugged
# output = LogFilesB/out.$(v01).out
# error = LogFilesB/out.$(v01).err
output = LogFilesC/out.$(v01).out
error = LogFilesC/out.$(v01).err

# OSG has different requirement syntax than CHTC
# Kernel version comes from major.minor.patch in the format:
# major * 10000 + minor * 1000 + patch
# heron comes from r-base which uses debian:latest, which
# as of 20190722 is "buster" from 4.19.105
# as per advice, 31000 will keep jobs away from RHEL 6

requirements = Arch == "X86_64" && \
  HAS_SINGULARITY == True && \
  OSG_HOST_KERNEL_VERSION >= 31000 && \
  SINGULARITY_CAN_USE_SIF == False && \
  has_sse4_2 == True
request_cpus = 6
request_memory = ifThenElse(MemoryUsage =!= undefined, (24000 * NumJobStarts * 2), 24000)
request_disk = 64GB

+SingularityImage = "/cvmfs/singularity.opensciencegrid.org/ncbi/pgap:2022-04-14.build6021"

# IF using a docker container, only send the RScript and any associated data files
transfer_input_files = YAMLFiles_v2/controller.yaml, \
                        YAMLFiles_v2/$(v06), \
                        $(v02), \
			                  osdf:///ospool/ap20/data/npcooley/Training/input-2022-04-14.build6021.tgz

# release jobs held for going over memory
periodic_release = (NumJobStarts < 4) && \
                    ((HoldReasonCode == 26) && (HoldReasonSubCode == 0)) && \
                    ((CurrentTime - EnteredCurrentStatus > 120))
                    
periodic_remove = NumJobStarts >= 4

# register as completed
# on_exit_hold = (ExitBySignal == True) || (ExitCode != 0) && \
#                (NumJobStarts < 5)

# release a job that has been held, after it has been held for 10 minutes, up to a maximum of 5 retries
# periodic_release = (NumJobStarts < 10) && \
#                    (((HoldReasonCode == 12) && (HoldReasonSubCode == 0 )) || (HoldReasonCode == 3)) && \
#                    ((CurrentTime - EnteredCurrentStatus) > 600)

# don't overwhelm the submit node
# max_idle = 1000
max_materialize = 4000

arguments = $(v01) $(v02) $(v03) $(v04) $(v05) $(v06) $(v07)
# transfer_output_remaps = "Result*.RData = OutDir/Result$(Process).RData"
# This doesn't seem to work?
# Always test a small number of jobs first!
queue v01, v02, v03, v04, v05, v06, v07 from TestMapC_v2.txt

