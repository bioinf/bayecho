import subprocess
import threading
import itertools
import os
import struct
import random
import mmap
import time
import datetime
import tempfile
import sys
import shutil
import gzip
from optparse import OptionParser

from scipy.stats.distributions import norm, poisson
from scipy.special import gammaln
from scipy import floor, ceil
from numpy import array, arange
import numpy
import math

##########################

HashingCmd = "hashing"
NeighborJoinCmd = "NeighborJoin"
NeighborJoinParamCmd = "NeighborJoinParam"
ParallelNeighborJoinCmd = "ParallelNeighborJoin"
ParallelNeighborJoinParamCmd = "ParallelNeighborJoinParam"
VotingCmd = "Voting"

##########################

CompBase = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
maxThread = 1 # Default max number of threads
maxEMIter = 10 # Default max iterations for EM algorithm
KeepAllFiles = False
MsgLogger = None
Verbose = True

##########################

# Default parameters
TmpBASEDIR = "tmp"
OutputDIR = "output"
LogDIR = "log"

OutputFName = "corrected"

RevCompFName = "revdataMMAP.txt"
ConfMatFName = "confMatEst.txt"

ModelSelectionSetSize = 100000
ModelSelectionNHashFile = 1

min_parame = None
max_parame = None
min_paramh = None
max_paramh = None

# Threshold for EM in parameter estimation
ConfMatCovergenceThreshold = 1e-10

##########################

class LogMgr:
    def __init__(self, log_file_name):
        self.logfile = open(log_file_name, "w")
        self.Lock = threading.Lock()

    def __delete__(self):
        self.logfile.close()

    def log(self, msg, console_output=False):
        with self.Lock:
            output_str = "%s: %s"%(time.strftime("%Y_%m_%d-%H_%M_%S", time.localtime()), msg)
            self.logfile.write("%s\n"%output_str)
            self.logfile.flush()
            if Verbose or console_output:
                print output_str + "\n"

##########################

class CmdExecuter:
    def __init__(self, nprocess):
        self.nCurThread = threading.BoundedSemaphore(nprocess)
        self.Lock = threading.Lock()
        self.LockStatus = threading.Lock()
        self.threads = []
        self.status = []

    def __exeCmd(self, cmd, t_id):
        finished = False

        self.nCurThread.acquire()

        cmd_string = " ".join(cmd)
        MsgLogger.log(cmd_string)

        try:
            max_attempts = 5 # additional attempts
            attempts = 0
            while not finished and attempts <= max_attempts:
                try:
                    proc = subprocess.Popen(cmd, bufsize=-1, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    stdout, stderr = proc.communicate()
                    if stdout != "":
                        print stdout
                    returncode = proc.returncode
                    if stderr is not None and len(stderr)>0:
                        MsgLogger.log("ERROR: " + stderr, True)

                    if returncode == 0:
                        finished = True
                    else:
                        print "Return code in executer ", returncode
                        break
                except Exception as e:
                    attempts += 1
                    if attempts <= max_attempts:
                        MsgLogger.log("ERROR:\n" + cmd_string + "\nfailed to execute.\nTrying again. Attempt #: " + str(attempts), True)
                    else:
                        MsgLogger.log("ERROR:\n" + cmd_string + "\nfailed to execute.", True)

            if finished == True: # success
                with self.LockStatus:
                    self.status[t_id] = 0
            else: # failure
                with self.LockStatus:
                    MsgLogger.log("ERROR:\n" + cmd_string + "\nfailed to execute. Aborting.")
                    self.status[t_id] = 1
        finally:
            self.nCurThread.release()

    def exeCmd(self, cmd):
        t = threading.Thread(target=self.__exeCmd, args=(cmd,len(self.threads)))
        with self.Lock:
            self.threads.append(t)
        with self.LockStatus:
            self.status.append(-1)
        t.start()

    def wait(self):
        with self.Lock:
            for t in self.threads:
                t.join()
            self.threads = []

        with self.LockStatus:
            passback_status = list(self.status)
            self.status = []

        # Check status of processes
        for st in passback_status:
            if st != 0:
                MsgLogger.log("ERROR: A command failed to execute successfully. Aborting.")
                sys.exit(1)
        return passback_status

##########################

def poiTotalVar(hist):
    n = sum(hist)
    mode = numpy.argmax(hist[2:]) + 2

    est_mean = mode+1
    est_tv = 1.0
    for m in range(mode-2, mode+3):
        emp_pmf = hist/n
        poi_pmf = poisson.pmf(arange(len(hist)), m)
        residual_mass = 1.0 - sum(poi_pmf)

        total_var = 0.5*sum( abs(p1-p2) for p1, p2 in itertools.izip(emp_pmf, poi_pmf)) + 0.5*residual_mass

        if total_var < est_tv:
            est_tv = total_var
            est_mean = m

    return (est_tv, est_mean)

##########################

def makeHashingCmd(readFName, read_st, read_ed, nhash, ihash_st, ihash_ed, K, fileNamePrefix, fileNameSuffix):
    cmd = [HashingCmd,
           "-read", "%s"%readFName,
           "-st",  "%d"%read_st,
           "-ed",  "%d"%read_ed,
           "-nhash",  "%d"%nhash,
           "-ihash_st", "%d"%ihash_st,
           "-ihash_ed", "%d"%ihash_ed,
           "-K", "%d"%K,
           "-fpre", "%s"%fileNamePrefix,
           "-fsuf", "%s"%fileNameSuffix]
    return cmd

def makeNeighborJoinCmd(readFName, inputFNames, read_st, read_ed, read_st2, read_ed2, blocksize, K, h, e, maxCov, fileNamePrefix, fileNameSuffix, ihash_st, ihash_ed):
    cmd = [NeighborJoinCmd,
           "-read", "%s"%readFName,
           "-st", "%d"%read_st,
           "-ed", "%d"%read_ed,
           "-st2", "%d"%read_st2,
           "-ed2", "%d"%read_ed2,
           "-K", "%d"%K,
           "-h", "%d"%h,
           "-e", "%f"%e,
           "-fpre", "%s"%fileNamePrefix,
           "-fsuf", "%s"%fileNameSuffix,
           "-ihash_st", "%d"%ihash_st,
           "-ihash_ed", "%d"%ihash_ed,
           "-max_cov", "%d"%maxCov
           ]
    for fname in inputFNames:
        cmd += ["-input", "%s"%fname]
    return cmd

def makeNeighborJoinParamCmd(readFName, inputFNames, read_st, read_ed, read_st2, read_ed2, blocksize, K, h, e, maxCov, fileNamePrefix, fileNameSuffix, ihash_st, ihash_ed):
    cmd = [NeighborJoinParamCmd,
           "-read", "%s"%readFName,
           "-st", "%d"%read_st,
           "-ed", "%d"%read_ed,
           "-st2", "%d"%read_st2,
           "-ed2", "%d"%read_ed2,
           "-K", "%d"%K,
           "-h", "%d"%h,
           "-e", "%f"%e,
           "-fpre", "%s"%fileNamePrefix,
           "-fsuf", "%s"%fileNameSuffix,
           "-ihash_st", "%d"%ihash_st,
           "-ihash_ed", "%d"%ihash_ed,
           ]
    for fname in inputFNames:
        cmd += ["-input", "%s"%fname]
    return cmd

def makeVotingCmd(readFName, inputFNames, read_st, read_ed, K, h, e, fileNamePrefix, fileNameSuffix, confMatFName=None, maxCov=None, minCov=None, estCov=None, h_rate=None, saveParam=False):
    cmd = [VotingCmd,
           "-read", "%s"%readFName,
           "-st", "%d"%read_st,
           "-ed", "%d"%read_ed,
           "-K", "%d"%K,
           "-h", "%d"%h,
           "-e", "%f"%e,
           "-fpre", "%s"%fileNamePrefix,
           "-fsuf", "%s"%fileNameSuffix]
    for fname in inputFNames:
        cmd += ["-input", "%s"%fname]

    if confMatFName is not None:
        cmd += ["-confMat", "%s"%confMatFName]
    if maxCov is not None:
        cmd += ["-max_cov", "%d"%maxCov]
    if minCov is not None:
        cmd += ["-min_cov", "%d"%minCov]
    if estCov is not None:
        cmd += ["-cov", "%d"%estCov]
    if h_rate is not None:
        cmd += ["-h_rate", "%f"%h_rate]
    if saveParam:
        cmd += ["-save_stats", "1"]
    return cmd

# Hash the kmers of the reads and create hash files. Hash files are then merged and used as input to NeighborJoin.
def Hashing(cmdexecuter, target_st, target_ed, nData, paramK):
  # Check for file
  tmpf = open(RevCompFName) # Will throw error if file does not open
  tmpf.close()

  # Building hash
  ihash_st = 0
  ihash_ed = 1
  cmdexecuter.exeCmd(makeHashingCmd(RevCompFName, target_st, target_ed, 1, ihash_st, ihash_ed, paramK, TmpDIR, "%d"%target_st))
  cmdexecuter.wait() # This will prevent parallization in hashing

  # Rename final hash and index files from tmp names to permanent names
  cmdexecuter.exeCmd(["mv", os.path.join(TmpDIR, "%d.hash"%(target_st)), os.path.join(TmpDIR, "all.hash")])
  cmdexecuter.exeCmd(["mv", os.path.join(TmpDIR, "%d.index"%(target_st)), os.path.join(TmpDIR, "all.index")])
  cmdexecuter.wait()

def Neighboring(cmdexecuter, target_st, target_ed, nData, paramK, paramh, parame, maxCov, nKmers, param):
  read_st = target_st
  read_ed = target_ed 
  read_st2 = 0
  read_ed2 = nData

  if read_st<=read_st2:
    st, ed, st2, ed2 = read_st, read_ed, read_st2, read_ed2
  else:
    st, ed, st2, ed2 = read_st2, read_ed2, read_st, read_ed

  # Create empty neighbor file
  NeighborFiles = [[]]
  ihash_st = 0
  ihash_ed = nKmers

  inputFNames = [ os.path.join(TmpDIR, "all.hash"), os.path.join(TmpDIR, "all.index") ]
  # NeighborJoin command
  if param == False:
    cmdexecuter.exeCmd(makeNeighborJoinCmd(RevCompFName, inputFNames, st, ed, st2, ed2, nData, paramK, paramh, parame, maxCov, TmpDIR, "%d_%d_%d"%(0, st, st2), ihash_st, ihash_ed))
  else:
    cmdexecuter.exeCmd(makeNeighborJoinParamCmd(RevCompFName, inputFNames, st, ed, st2, ed2, nData, paramK, paramh, parame, maxCov, TmpDIR, "%d_%d_%d"%(0, st, st2), ihash_st, ihash_ed))
  cmdexecuter.wait()
  NeighborFiles[0].append(os.path.join(TmpDIR, "neighbors_%d_%d_%d.list"%(0, st, st2)))
  
  cmdexecuter.exeCmd(["mv", NeighborFiles[-1][0], os.path.join(TmpDIR, "neighbors_%d_%d.list"%(st, st2))])
  cmdexecuter.wait()

def Voting(cmdexecuter, target_st, target_ed, nData, paramh, parame, confMatFName = None, maxCov = None, minCov = None, estCov = None, h_rate = None, saveParam=False):
    # Vote!
  inputFNames = [ os.path.join(TmpDIR, "neighbors_%d_%d.list"%(0,0))]
  cmdexecuter.exeCmd(makeVotingCmd(RevCompFName, inputFNames, target_st, target_ed, paramK, paramh, parame, TmpDIR, "%d_%d_%f"%(target_st, paramh, parame), confMatFName, maxCov, minCov, estCov, h_rate, saveParam))

def loadConfMat(fname):
    ConfMat = []
    with open(fname) as fin:
        nMat = int(fin.readline().strip())
        for i in range(nMat):
            ConfMat.append([[ float(token) for token in fin.readline().split()] for b in range(4)])
            fin.readline()

    return array(ConfMat)

if __name__ == '__main__':
    # Set starting time
    ST_TIME = time.time()

    ##########################
    # Program paths
    ProgramDir = sys.path[0]
    HashingCmd = os.path.join(ProgramDir, HashingCmd)
    NeighborJoinCmd = os.path.join(ProgramDir, NeighborJoinCmd)
    NeighborJoinParamCmd = os.path.join(ProgramDir, NeighborJoinParamCmd)
    VotingCmd = os.path.join(ProgramDir, VotingCmd)

    ##########################
    # Parse Options
    parser = OptionParser(usage = "usage: %prog [options] read_file_name")
    parser.add_option("-o", "--output", action="store", dest="output_filename", type="string", help="Output file name", default=None)
    parser.add_option("-l", "--log", action="store", dest="log_filename", type="string", help="Log file name", default=None)
    parser.add_option("--DD", "--tmp_dir", action="store", dest="tmp_directory", type="string", help="Temporary data directory", default=None)
    parser.add_option("-u", "--ncpu", action="store", dest="ncpu", type="int", help="Number of processes used in training", default=maxThread)

    parser.add_option("-k", "--kmer", action="store", dest="k", type="int", help="k-mer size used for hashing", default=None)
    parser.add_option("-e", "--min_error_tolerance", action="store", dest="e", type="float", help="Minimum error tolerance for parameter searching", default=min_parame)
    parser.add_option("-E", "--max_error_tolerance", action="store", dest="E", type="float", help="Maximum error tolerance for parameter searching", default=max_parame)
    parser.add_option("--hh", "--min_minimum_overlap", action="store", dest="h", type="int", help="Minimum minimum overlap length for parameter searching", default=min_paramh)
    parser.add_option("--hH", "--max_minimum_overlap", action="store", dest="H", type="int", help="Maximum minimum overlap length for parameter searching", default=max_paramh)
    parser.add_option("--h_rate", "--heterozygous_rate", action="store", dest="h_rate", type="float", help="Rate for heterozygous site", default=None)
    parser.add_option("--model_selection_size", action="store", dest="msize", type="int", help="Model selection data set size", default=ModelSelectionSetSize)
    parser.add_option("--keep_all_files", action="store_true", dest="keep_all_files", help="Keep all temporary files. By default, temporary files are deleted automatically.", default=False)
    parser.add_option("--hash", action="store", dest="hash_file", type="string", help="hash file", default=None)
    parser.add_option("--index", action="store", dest="index_file", type="string", help="index_file", default=None)

    parser.add_option("--pn", action="store_true", dest="pn", help="Do parallel neighboring (default: False)", default=False)

    (options, args) = parser.parse_args()
    if len(args)!=1:
        parser.print_help()
        sys.exit(0)
    OrigReadFName = args[0]
    hash_file = options.hash_file
    index_file = options.index_file

    if options.pn:
        NeighborJoinCmd = os.path.join(ProgramDir, ParallelNeighborJoinCmd)
        NeighborJoinParamCmd = os.path.join(ProgramDir, ParallelNeighborJoinParamCmd)
 

    # Check for eistence of input file
    if not os.path.isfile(OrigReadFName):
        print "ERROR: Input file " + OrigReadFName + " does not exist."
        sys.exit(1)

    # Create log directory
    if options.log_filename is None:
        # default
        if not os.path.isdir(LogDIR):
            os.mkdir(LogDIR)
    else:
        # create intermediate dirs for user-specified log file if necessary
        head, tail = os.path.split(options.log_filename)
        assert(tail != "" and tail is not None)
        if head.strip() != "" and not os.path.isdir(head):
            os.makedirs(head)

    if options.log_filename is None:
        # default
        MsgLogger = LogMgr(os.path.join(LogDIR, "run_"+time.strftime("%Y_%m_%d-%H_%M_%S", time.localtime())+".log"))
    else:
        MsgLogger = LogMgr(options.log_filename)

    MsgLogger.log("Run with parameters: " + str(sys.argv))

    # Create tmp directory
    try:
        if options.tmp_directory is None:
            # Set default to "./tmp/(random_name)"
            if not os.path.exists(TmpBASEDIR):
                os.mkdir(TmpBASEDIR)
            elif not os.path.isdir(TmpBASEDIR):
                MsgLogger.log("File with same name as default tmp directory. Please rename file or choose different tmp directory.")
                raise OSError("File name conflict: " + TmpBASEDIR)
            options.tmp_directory = tempfile.mkdtemp(dir=TmpBASEDIR)
        else:
            # User-specified tmp directory
            assert(options.tmp_directory != "" and options.tmp_directory is not None)
            if not os.path.isdir(options.tmp_directory):
                os.mkdir(options.tmp_directory)
            # Create inner tmp directory within the one specified by user
            # The inner tmp directory is given a randomized name that is
            # guaranteed to be unique
            options.tmp_directory = tempfile.mkdtemp(dir=options.tmp_directory)
        assert(options.tmp_directory is not None)
    except Exception as e:
        MsgLogger.log("Error in creating tmp directory.")
        MsgLogger.log(str(e))
        raise

    TmpDIR = os.path.normpath(options.tmp_directory)+'/'
    MsgLogger.log("Temporary Directory: %s"%TmpDIR, True)

    # Create output directory
    try:
        if options.output_filename is None:
            # default
            if os.path.isfile(OutputDIR):
                MsgLogger.log("File with same name as default output directory exists. Please rename file or choose different output directory.")
                raise OSError("File name conflict: " + OutputDIR)
            if not os.path.isdir(OutputDIR):
                os.mkdir(OutputDIR)
            options.output_filename = os.path.join(OutputDIR, OutputFName+"_"+time.strftime("%Y_%m_%d-%H_%M_%S", time.localtime())+".txt")
            MsgLogger.log("Default output file is " + options.output_filename + ".")
        else:
            # create intermediate dirs for user-specified output file if necessary
            head, tail = os.path.split(options.output_filename)
            assert(tail != "" and tail is not None)
            if head.strip() != "" and not os.path.isdir(head):
                os.makedirs(head)
        assert(options.output_filename is not None)
    except Exception as e:
        MsgLogger.log("Error in creating output directory.")
        MsgLogger.log(str(e))
        raise

    # Ensure that output file name does not conflict with a directory
    if(os.path.isdir(options.output_filename)):
        MsgLogger.log(options.output_filename + " is a directory. Please choose a different output file name")
        sys.exit(1)
    if os.path.isdir(options.output_filename + ".qual"):
        MsgLogger.log(options.output_filename + ".qual is a directory. Please choose a different output file name")
        sys.exit(1)
    if os.path.isdir(options.output_filename + ".seq"):
        MsgLogger.log(options.output_filename + ".seq is a directory. Please choose a different output file name")
        sys.exit(1)

    # Ensure that output file name does not conflict with a directory
    if(os.path.isdir(options.output_filename)):
        MsgLogger.log(options.output_filename + " is a directory. Please choose a different output file name")
        sys.exit(1)
    if os.path.isdir(options.output_filename + ".qual"):
        MsgLogger.log(options.output_filename + ".qual is a directory. Please choose a different output file name")
        sys.exit(1)
    if os.path.isdir(options.output_filename + ".seq"):
        MsgLogger.log(options.output_filename + ".seq is a directory. Please choose a different output file name")
        sys.exit(1)

    # It's important to create directories early on so if there is a file
    # system conflict, the script aborts before doing too much work
    #######################################################################

    # Reverse complement file
    RevCompFName = os.path.join(TmpDIR, "revdataMMAP.txt")
    # Confusion Matrix file
    ConfMatFName = os.path.join(TmpDIR, "confMatEst.txt")
    # Output file name
    OutputFName = options.output_filename

    maxThread = options.ncpu

    KeepAllFiles = options.keep_all_files

    ##########################
    if options.k is not None:
        paramK = options.k
    else:
        # If K is not specified, then compute one later
        paramK = None

    ##########################
    # Initialize
    MsgLogger.log("Starting...", True)

    # Set the number of threads cmdexecuter should use
    cmdexecuter = CmdExecuter(maxThread)

    ###############################
    # Reverse Complement and MMAP
    ###############################

    #########################################
    # Generate reverse complements of reads
    with tempfile.NamedTemporaryFile(dir=TmpDIR,delete=False) as TmpFile:
        file_is_fastq = False
        (file_head, file_ext) = os.path.splitext(OrigReadFName)
        if(file_ext == ".gz"):
            file_open_function = gzip.open
            # Check second extension
            (file_head2, file_ext2) = os.path.splitext(file_head)
            if(file_ext2 == ".fastq"):
                file_is_fastq = True
            else:
                file_is_fastq = False
        else:
            file_open_function = open
            if(file_ext == ".fastq"):
                file_is_fastq = True
            else:
                file_is_fastq = False

        if file_is_fastq: # fastq file
            try:
                fin = file_open_function(OrigReadFName, "r")
                try:
                    while True:
                        line0 = fin.next()
                        line1 = fin.next()
                        line2 = fin.next()
                        line3 = fin.next()

                        # Quick sanity checks to ensure FASTQ format
                        assert(line0[0] == "@")
                        assert(line2[0] == "+")

                        line1 = line1.strip().upper()
                        TmpFile.write("%s\n"%line1)
                        TmpFile.write("%s\n"%"".join([ CompBase[b] for b in reversed(line1)]))
                except StopIteration as e:
                    pass
            finally:
               if fin is not None:
                   fin.close()
        else: # regular text file
            try:
                fin = file_open_function(OrigReadFName, "r")
                for line in fin:
                    line = line.strip().upper()
                    TmpFile.write("%s\n"%line)
                    TmpFile.write("%s\n"%"".join([ CompBase[b] for b in reversed(line)]))
            finally:
                if fin is not None:
                    fin.close()

    ######################
    # Create MMAP index
    nData = 0
    ReadLen = []
    ReadSt = [0]
    with open(TmpFile.name, "r") as fin:
        for line in fin:
            nData += 1
            ReadLen.append(len(line.strip()))
            ReadSt.append(ReadSt[-1]+ReadLen[-1]+1+1) # 1 byte for string array NULL end, 1 byte for indicator of original read
    ReadSt = ReadSt[:-1]

    ModelSelectionSetSize = min(options.msize, nData)

    random.seed(1)
    SeqOrdering = range(nData)
    random.shuffle(SeqOrdering)

    with open(TmpFile.name, "r") as fin:
        with open(RevCompFName, "wb") as fout:
            fout.write(struct.pack("Q", nData))
            for idx in SeqOrdering:
                fout.write(struct.pack("Q", ReadSt[idx]))

            for line, id in itertools.izip(fin, itertools.count()):
                line = line.strip()
                if id%2==0:
                    fout.write("Y%s\x00"%line)
                else:
                    fout.write("N%s\x00"%line)
    os.unlink(TmpFile.name)

    MsgLogger.log("Data preprocessing complete.", True)

    maxReadLen = max(ReadLen)

    # Set minparam_h, maxparam_h, minparam_e, maxparam_e
    if options.h is None:
        min_paramh = int(maxReadLen * 0.5)
    else:
        min_paramh = int(options.h)

    if options.H is None:
        max_paramh = int(maxReadLen * 0.66)
    else:
        max_paramh = int(options.H)

    if options.e is None:
        min_parame = 0.05
    else:
        min_parame = options.e

    if options.E is None:
        max_parame = 0.15
    else:
        max_parame = options.E

    assert(min_parame <= max_parame)
    assert(min_paramh <= max_paramh)

    # Initialize K by using number of data
    if options.k is None:
        paramK = max(int(maxReadLen/5.0), int(math.ceil(math.log(nData*max(ReadLen))/math.log(4))))
    assert(paramK > 0)

    MsgLogger.log("Use parameter k = %d"%paramK, True)
    MsgLogger.log("Overlap search range: [%d, %d]"%(min_paramh, max_paramh), True)
    MsgLogger.log("Error tolerance search range: [%f, %f]"%(min_parame, max_parame), True)

    ##########################
    # HASHING STAGE
    ##########################

    # Hashing
    # Use only one thread in this stage to maximize available memory
    # This reduces the number of "read blocks" that are required
    cmdexecuter = CmdExecuter(1)

    # Hashing arguments: (cmdexecuter, startread, endread, stepsize, num reads, K, num hash_blocks)
    if (hash_file is None) or (index_file is None):
    	Hashing(cmdexecuter, 0, nData, nData, paramK)
    else:
      os.system("cp " + hash_file + " "+ index_file + " " + TmpDIR + "/")

    # Read in total number of kmers
    nKmers = -1
    with open(os.path.join(TmpDIR, "all.index"), "rb") as tmpfile:
        nKmers = int(struct.unpack("@I", tmpfile.read(4))[0])

    if nKmers < 0:
        print "ERROR: Number of kmers in " + os.path.join(TmpDIR, "all.index") + " is not valid."
        sys.exit(1)

    ##########################
    # END OF HASHING STAGE
    ##########################

    #################################
    # PARAMETER SELECTION STAGE
    #################################

    # Use multiple threads for parameter selection
    cmdexecuter = CmdExecuter(maxThread)

    #################################
    # Cluster reads for parameter selection
    # Construct Neighbor Sets
    Neighboring(cmdexecuter, 0, ModelSelectionSetSize, nData,  paramK, min_paramh, max_parame, 0, nKmers, True)
    cmdexecuter.wait()

    # Parameter Search
    # Voting has multiple purposes
    # Here it is used to create histograms for a grid search of parameters
    cmdexecuter = CmdExecuter(maxThread)
    for cur_paramh in arange(min_paramh, max_paramh+1, 2):
        for cur_parame in arange(min_parame, max_parame+0.01, 0.05):
            Voting(cmdexecuter, 0, ModelSelectionSetSize, nData, cur_paramh, cur_parame,  h_rate=options.h_rate, saveParam = True)
    cmdexecuter.wait()

    # Histogram Selection
    TotalVariations = []
    for cur_paramh in arange(min_paramh, max_paramh+1, 2):
        for cur_parame in arange(min_parame, max_parame+0.01, 0.05):
            with open(os.path.join(TmpDIR, "histogram_%d_%d_%f.txt"%(0, cur_paramh, cur_parame)), "r") as fin:
                line = fin.readline()
                hist = array([ int(token) for token in line.split()])
                TotalVariations.append(((cur_paramh, cur_parame), poiTotalVar(1.0*hist/sum(hist))))
                MsgLogger.log("Total Variation: %s"%str(TotalVariations[-1]))

    TotalVariations.sort(key=lambda x:x[1][0])
    best_paramh = TotalVariations[0][0][0]
    best_parame = TotalVariations[0][0][1]
    estCov = TotalVariations[0][1][1]
    maxCov = max(3, int(math.ceil(estCov + pow(estCov, 0.5)*5)))
    #minCov = max(3, int(math.ceil(estCov - pow(estCov, 0.5)*5)))
    minCov = 0
    MsgLogger.log("Estimated Coverage: %d"%TotalVariations[0][1][1], True)
    MsgLogger.log("param: h=%s e=%s cov=%s max_cov=%s min_cov=%s"%(best_paramh, best_parame, estCov, maxCov, minCov), True)
    MsgLogger.log("h, e, mu selection complete.", True)

    ######################################
    # Use EM to estimate Confusion Matrix
    oldConfMat = []
    newConfMat = []
    confMatFName = os.path.join(TmpDIR, "confmat_%d_%d_%f.txt"%(0, best_paramh, best_parame))
    predSeqFName = os.path.join(TmpDIR, "output_%d_%d_%f.txt"%(0, best_paramh, best_parame))

    oldConfMat = loadConfMat(confMatFName)
    for em_iter in range(maxEMIter):
        # Voting has many purposes
        # Here it is used to compute the EM estimator for the Confusion Matrix
        Voting(cmdexecuter, 0, ModelSelectionSetSize, nData, best_paramh, best_parame, confMatFName,
               maxCov, minCov, estCov, options.h_rate, saveParam = True)
        cmdexecuter.wait()

        # Convergence Check
        newConfMat = loadConfMat(confMatFName)
        MsgLogger.log("EM Iteration %s: relative error = %s"%(em_iter, numpy.sum(numpy.abs(oldConfMat - newConfMat))/numpy.sum(numpy.abs(oldConfMat))), True)
        if numpy.sum(numpy.abs(oldConfMat - newConfMat))/numpy.sum(numpy.abs(oldConfMat)) <= ConfMatCovergenceThreshold:
            break
        else:
            oldConfMat = newConfMat

    # Backup estimated Confusion Matrix
    cmdexecuter.exeCmd(["cp", confMatFName, ConfMatFName])
    cmdexecuter.wait()
    confMatFName = ConfMatFName
    MsgLogger.log("EM complete.", True)

    ############################################
    # END OF PARAMETER SELECTION STAGE
    ############################################

    ############################################
    # ERROR CORRECTION STAGE
    ############################################

    # Use best parameters to correct all reads!
    BatchSize = nData # Correct all data at once

    BatchSt = 0
    while BatchSt<nData:
        BatchEd = min(nData, BatchSt + BatchSize)
        ######################################
        # Construct adjacency lists of reads
        # Use one thread to maximize available memory
        cmdexecuter = CmdExecuter(1)
        Neighboring(cmdexecuter, BatchSt, BatchEd, nData, paramK, best_paramh, best_parame, maxCov,  nKmers, False)
        cmdexecuter.wait()

        ######################################
        # Voting
        # Voting has multiple purposes
        # Here it is used to correct reads
        Voting(cmdexecuter, BatchSt, BatchEd, nData,  best_paramh, best_parame,  confMatFName, maxCov, minCov, estCov, options.h_rate)
        cmdexecuter.wait()

        BatchSt += BatchSize

        ######################################
        # Delete Neighbor Files
        if not KeepAllFiles:
            MsgLogger.log("Deleting Neighbor Files...")
            cmdexecuter.exeCmd(["find", TmpDIR,  "-name", "neighbor*", "-delete"])
            cmdexecuter.wait()

    MsgLogger.log("Error correction complete.", True)

    ############################################
    # END OF ERROR CORRECTION STAGE
    ############################################

    ######################################
    # Collect and reorder output
    with tempfile.NamedTemporaryFile(dir=TmpDIR,delete=False) as TmpFile:
        with tempfile.NamedTemporaryFile(dir=TmpDIR,delete=False) as QualTmpFile:
            # Recompute ReadSt
            ReadSt = [0]
            for readlen in ReadLen:
                ReadSt.append(ReadSt[-1] + readlen + 1)
            ReadSt = ReadSt[:-1]

            TmpFile.seek(sum(ReadLen)+len(ReadLen)-1)
            TmpFile.write("\n")
            TmpFile.flush()
            QualTmpFile.seek(sum(ReadLen)+len(ReadLen)-1)
            QualTmpFile.write("\n")
            QualTmpFile.flush()

            readmap = mmap.mmap(TmpFile.fileno(), sum(ReadLen)+len(ReadLen))
            qualmap = mmap.mmap(QualTmpFile.fileno(), sum(ReadLen)+len(ReadLen))
            cur_read_order = 0
            for read_st in range(0, nData, nData):
                with open(os.path.join(TmpDIR, "output_%d_%d_%f.txt"%(read_st, best_paramh, best_parame)), "r") as fin:
                    with open(os.path.join(TmpDIR, "quality_%d_%d_%f.txt"%(read_st, best_paramh, best_parame)), "r") as qualfin:
                        for line, qualline in itertools.izip(fin,qualfin):
                            readmap[ReadSt[SeqOrdering[cur_read_order]]:(ReadSt[SeqOrdering[cur_read_order]]+ReadLen[SeqOrdering[cur_read_order]]+1)] = line.strip() + '\n'
                            qualmap[ReadSt[SeqOrdering[cur_read_order]]:(ReadSt[SeqOrdering[cur_read_order]]+ReadLen[SeqOrdering[cur_read_order]]+1)] = qualline.strip() + '\n'
                            cur_read_order += 1
            readmap.flush()
            readmap.close()
            qualmap.flush()
            qualmap.close()

    # Eliminate reverse complement reads
    with open(OutputFName+".seq", "w+b") as fout:
        with open(TmpFile.name, "r") as fin:
            for (line, counter) in itertools.izip(fin, itertools.count()):
                if counter%2==0:
                    fout.write(line)
    os.unlink(TmpFile.name)

    with open(OutputFName+".qual", "w+b") as fout:
        with open(QualTmpFile.name, "r") as fin:
            for (line, counter) in itertools.izip(fin, itertools.count()):
                if counter%2==0:
                    fout.write(line)
    os.unlink(QualTmpFile.name)

    with open(OutputFName, "w+b") as fout:
        with open(OutputFName+".seq", "r") as fseq:
            with open(OutputFName+".qual", "r") as fqual:
                for (seq, qual, id) in itertools.izip(fseq, fqual, itertools.count()):
                    fout.write("@%d\n%s\n+\n%s\n"%(id, seq.strip(), qual.strip()))

    if not KeepAllFiles:
        MsgLogger.log("Removing Temporary Directory...")
        shutil.rmtree(TmpDIR, True)

    MsgLogger.log("Output is in file: " + OutputFName)
    ED_TIME = time.time()
    MsgLogger.log("Total running time (h:m:s): " + str(datetime.timedelta(seconds=(ED_TIME - ST_TIME))))
    MsgLogger.log("End...", True)
