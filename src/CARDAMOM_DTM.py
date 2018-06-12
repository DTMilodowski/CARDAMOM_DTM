# This defines the class of object CARDAMOM.
# The object stores the parameters for a given CARDAMOM project setup and provides methods
# to compile code, write input, execute the MCMC, read output, etc.
# ----------------------------------------------------------------------------------------
# The original version was written by J-F Exbrayat, and I have modified this slightly to
# work with Luke Smallman's versions of DALEC, which will be used during both the BALI and
# Forests2020 projects.
# ----------------------------------------------------------------------------------------
import os, cPickle, time, struct, sys         
import datetime as dt       
import numpy as np
from netCDF4 import Dataset
import shutil, socket

from matplotlib import pyplot as plt
from matplotlib import rcParams
# Set up some basiic parameters for the plots
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['arial']
rcParams['font.size'] = 8
rcParams['legend.numpoints'] = 1
axis_size = rcParams['font.size']+2
colour = ['#46E900','#1A2BCE','#E0007F']

# Now load rerun stuff
#sys.path.append('./rerun/')
#import f2py_dalec_gsi_dfol_cwd_fr as f2py

#from DALEC_f2py import readParsDALEC, readDevelopmentParams, write2netCDF, GR


class CARDAMOM(object):
    
    def __init__(self, **kwargs):
        """
        Defines the CARDAMOM object 
        """
        #set the paths for the current project
        self.set_paths()

        # check if project details have been defined
        if 'project_name' in kwargs:
            self.project_name = kwargs['project_name']
        else:
            self.project_name = raw_input('Enter project name: ')
      
        # load if project exists
        if self.project_name in os.listdir(self.paths["projects"]): 
            if self.project_name+".pData" in os.listdir('/'.join([self.paths["projects"],self.project_name])):
                print "Project \"%s\" found" % self.project_name          
                self.load_project()
            else:
                print "Project directory found but file %s.pData could not be found" % self.project_name
                self.new_project() 
        else:
            os.mkdir('/'.join([self.paths["projects"],self.project_name]))
            self.new_project()

    #-----------------------------------------------------------------------------------------
    # This method creates the variables needed by CARDAMOM but not yet set to definite values:
    # - nPoints.............. the number of points to run CARDAMOM on
    # - latitude............. a sequence of length nPoints with the latitude 
    # - longitude............ a sequence of length nPoints with the longitude
    # - drivers.............. an array of dimensions [nPoints, nTsteps, nDrivers]
    # - observations......... an array of dimensions [nPoints, nTsteps, nObs]
    # NOTE: nPoints, nTsteps and nObs will be defined automatically from the shape of input
    # arrays
    def new_project(self):
               
        # This is initialized here but will contain a dictionary will all the project details
        self.details = {}
        projtype=raw_input("Enter project number of type: \n1\tDALEC_CDEA\n2\tDALEC_GSI_BUCKET\n3\tAT-DALEC\n4\tAT-DALEC_CROP\n5\tDALEC_CDEA_FR\n6\tDALEC_GSI_FR\n7\tDALEC_GSI_FR_DBio\n8\tDALEC_GSI_MFOL_FR\n9\tDALEC_GSI_FR_LABILE\n10\tDALECN_GSI_FR\n11\tDALEC_GSI_DFOL_FR\n12\tDALEC_GSI_DFOL_FROOT_FR\n13\tDALEC_GSI_DFOL_LABILE_FR\n14\tDALECN_GSI_DFOL_LABILE_FR\n15\tDALECN_GSI_DFOL_LABILE_FROOT_FR\n16\tDALEC_GSI_DFOL_CWD_FR\n")
        if projtype != "":       
            self.project_type=projtype
            if self.project_type == "0":
                self.project_type = "ACM"
            elif self.project_type == "1":
                self.project_type = "DALEC_CDEA"
            elif self.project_type == "2":
                self.project_type = "DALEC_GSI_BUCKET"
            elif self.project_type == "3":
                self.project_type = "AT-DALEC"
            elif self.project_type == "4":
                self.project_type = "AT-DALEC_CROP"
            elif self.project_type == "5":
                self.project_type = "DALEC_CDEA_FR"
            elif self.project_type == "6":
                self.project_type = "DALEC_GSI_FR"
            elif self.project_type == "7":
                self.project_type = "DALEC_GSI_FR_DBio"
            elif self.project_type == "8":
                self.project_type = "DALEC_GSI_MFOL_FR"
            elif self.project_type == "9":
                self.project_type = "DALEC_GSI_FR_LABILE"
            elif self.project_type == "10":
                self.project_type = "DALECN_GSI_FR"
            elif self.project_type == "11":
                self.project_type = "DALEC_GSI_DFOL_FR"
            elif self.project_type == "12":
                self.project_type = "DALEC_GSI_DFOL_FROOT_FR"
            elif self.project_type == "13":
                self.project_type = "DALEC_GSI_DFOL_LABILE_FR"
            elif self.project_type == "14":
                self.project_type = "DALECN_GSI_DFOL_LABILE_FR"
            elif self.project_type == "15":
                self.project_type = "DALECN_GSI_DFOL_LABILE_FROOT_FR"
            elif self.project_type == "16":
                self.project_type = "DALEC_GSI_DFOL_CWD_FR" # this is Luke's model that includes CWD

        else:
            print "project type not found, defaulting to CARDAMOM_CDEA"
            self.project_type="CARDAMOM_CDEA"

        self.save_project()

    #-------------------------------------------------------------------------------------
    # Saves the object attributes in a file named <project_name>.pData
    def save_project(self):

        f = open('/'.join([self.paths["projects"],self.project_name,self.project_name+".pData"]),"wb")
        cPickle.dump(self.__dict__,f,2)
        f.close()

    #-------------------------------------------------------------------------------------
    # Copies the project into a <newproject>
    def copy_project(self,newproject):

        keep_name = self.project_name
        self.project_name = newproject
        if newproject not in os.listdir(self.paths["projects"]):
            os.mkdir(self.paths["projects"]+"/"+newproject)
        f = open('/'.join([self.paths["projects"],newproject,newproject+".pData"]),"wb")
        cPickle.dump(self.__dict__,f,2)
        f.close()
        self.project_name = keep_name

    #-------------------------------------------------------------------------------------
    # Loads the object attributes stored in a cPickle file named <project_name>.pData
    def load_project(self):

        f = open('/'.join([self.paths["projects"],self.project_name,self.project_name+".pData"]),"rb")
        tmp_dict = cPickle.load(f)
        f.close()

   #     if "paths" in tmp_dict:
   #         tmp_dict.pop("paths",None)

        self.__dict__.update(tmp_dict)

    #-------------------------------------------------------------------------------------
    # Loads paths to the different directories required by CARDAMOM
    def set_paths(self):

        hostname = socket.gethostname()
        if "paths" in dir(self):
            print "Paths found in project"
        else:
            print "Paths need to be defined"
            if "default_paths_%s.pData" % (hostname) in os.listdir(os.getcwd()):
                usedefault = raw_input("Use default paths <y/n>? ")
                if usedefault == 'y':
                    f = open("default_paths_%s.pData" % (hostname),"r")
                    self.paths = cPickle.load(f)
                    f.close()
                else:
                    self.define_paths()
            else:
                print "No default paths found... defining paths now"
                self.define_paths()

    #-------------------------------------------------------------------------------------
    # Define the paths used by the project and asked whether to use them as default ones
    def define_paths(self):

        hostname = socket.gethostname()
        self.paths = {}
        
        keepcurrent = raw_input("Keep current directory (%s) as root one <y/n>? " % os.getcwd())
        
        if keepcurrent == 'y':
            self.paths["CARDAMOM"] = os.getcwd()+"/"
        else:
            self.paths["CARDAMOM"] = raw_input("Enter CARDAMOM root directory: ")

        self.paths["projects"] = raw_input("Where will this project be saved? (provide full path or leave blank to use subdirectory \"projects\"): ")
        if self.paths["projects"] == "":
            self.paths["projects"] = self.paths["CARDAMOM"]+"/projects/"

        self.paths["library"] = raw_input("Where are source codes kept? (provide full path or leave blank to use subdirectory /CARDAMOM/trunk/LIBRARY/CARDAMOM_F): ")
        if self.paths["library"] == "":
            self.paths["library"] = self.paths["CARDAMOM"]+"/CARDAMOM/trunk/LIBRARY/CARDAMOM_F/"

        usecluster = raw_input("Will you run this project on a cluster <y/n>? ")
        if usecluster == "y":
            self.paths["cluster_username"] = raw_input("Enter your username (leave blank for default): ")
            if self.paths["cluster_username"] == "":
                self.paths["cluster_username"] = "dmilodow"

            self.paths["cluster_address"] = raw_input("Enter address of cluster (leave blank for eddie): ")
            if self.paths["cluster_address"] == "":
                self.paths["cluster_address"] = "eddie3.ecdf.ed.ac.uk"

            self.paths["cluster_directory"] = raw_input("Enter cluster working directory (full path or leave blank for default): ")
            if self.paths["cluster_directory"] == "":
                self.paths["cluster_directory"] = "/exports/csce/eddie/geos/groups/gcel/dmilodow/CARDAMOM/"
            
        savedefault = raw_input("Save current paths as default ones for this machine <y/n>? ")
        
        if savedefault == 'y':
            f = open("default_paths_%s.pData" % hostname,"w")
            cPickle.dump(self.paths,f)
            f.close()


    #-------------------------------------------------------------------------------------
    # This method is a wrapper to the function used to setup the project
    def setup(self,latitude,longitude,drivers,observations,parprior,parpriorunc,otherprior,otherpriorunc,**kwargs):

        # pixel data     
        self.details["latitude"] = latitude
        self.details["longitude"] = longitude
        self.details["drivers"] = drivers
        self.details["observations"] = observations

        # MCMC specific values
        self.details["parprior"] = parprior
        self.details["parpriorunc"] = parpriorunc
        self.details["otherprior"] = otherprior
        self.details["otherpriorunc"] = otherpriorunc

        # technical options
        EDCs = raw_input("Use EDCs <y/n>? ")
        if EDCs == 'y':
            self.details["EDCs"] = True
            
        warning = 0
        #get the number of time steps
        if drivers.ndim == 2:
            self.details["no_pts"] = 1
            self.details["met_fields"] = drivers.shape[1]
            self.details["obs_fields"] = observations.shape[1]
            if drivers.shape[0] != observations.shape[0]:
                print " /!\ Warning: dimensions of drivers and observations arrays do not agree on the number of time steps /!\ "           
                warning += 1
            else:
                self.details["tsteps"] = drivers.shape[0]
                 
        else:
            self.details["no_pts"] = drivers.shape[0]
            self.details["met_fields"] = drivers.shape[2]
            self.details["obs_fields"] = observations.shape[2]
            if drivers.shape[1] != observations.shape[1]:
                print " /!\ Warning: dimensions of drivers and observations arrays do not agree on the number of time steps /!\ "
                warning += 1
            else:
                self.details["tsteps"] = drivers.shape[1]

            if drivers.shape[0] != latitude.shape[0]:
                print " /!\ Warning: dimensions of drivers and latitude arrays do not agree on the number of sites /!\ "  
                warning += 1
    
        #assign the number of met fields
        if warning == 0:
            print "Saving project before writing binary files...   ",
            self.save_project()
            print "OK"

            #write input data
            if "runid" in kwargs:
                self.runid = kwargs["runid"]
            else:
                self.runid = -1 # placeholder  
            self.createInput()
            #copy the source code
            self.backup_source()

            """            
            #create a directory for the executable
            if "exec" not in os.listdir(self.paths["projects"]+self.project_name):
                os.mkdir(self.paths["projects"]+self.project_name+"/exec")
            """

            compile_code = raw_input("Compile local code <y/n>? ")
            if compile_code == "y":
                print "Compiling local code"
                self.backup_source()
                self.compile_local_code()
                

            cluster = raw_input("Send project files on the cluster <y/n>? ")
            
            if cluster == "y":
                self.details["cluster"] = True
                #copy data on cluster
                self.send_to_cluster()

            cluster = raw_input("Compile code on the cluster <y/n>? ")
            if cluster == "y":
                self.compile_cluster()

        else:
            print "Too many warning raised... please check input data according to previous messages" 

    #-------------------------------------------------------------------------------------
    # This method writes the input files in the local directory
    def createInput(self):

        # Defines where the input files will be written        
        path2project = '/'.join([self.paths["projects"],self.project_name])
        path2data = path2project+"/data/"
        path2cardamom_output = path2project+"/cardamom_output/"

        if "data" not in os.listdir(path2project):
            print "Directory \"%s\" not found... creating" % path2data
            os.mkdir(path2data)

        # check if no run already specified
        if self.runid < 0:
            runlist = os.listdir("%s" % (path2data))
            if len(runlist) == 0:
                self.runid = 1
            else:
                runlist.sort()
                self.runid = int(runlist[-1].split("_")[-1])+1
            print "WARNING - no runid specified in setup. Default set to %03i" % self.runid
            
        if ("%03i" % self.runid) not in os.listdir(path2data):
            print "Directory \"%s\%03i\" not found... creating" % (path2data,self.runid)
            os.mkdir("%s/%03i" % (path2data,self.runid))
        if "cardamom_output" not in os.listdir(path2project):
            print "Directory \"%s\" not found... creating" % path2cardamom_output
            os.mkdir(path2cardamom_output)
        if ("%03i" % self.runid) not in os.listdir(path2cardamom_output):
            print "Directory \"%s\%03i\" not found... creating" % (path2cardamom_output,self.runid)
            os.mkdir("%s/%03i" % (path2cardamom_output,self.runid))


        print "Now creating input data in \"%s\" for project \"%s\"" % (path2data,self.project_name)

        if self.project_type == "DALEC_CDEA":
            modelid = 1
        elif self.project_type == "DALEC_GSI_BUCKET":
            modelid = 2
        elif self.project_type == "AT-DALEC":
            modelid = 3 
        elif self.project_type == "AT-DALEC_CROP":
            modelid = 4 
        elif self.project_type == "DALEC_CDEA_FR":
            modelid = 5 
        elif self.project_type == "DALEC_GSI_FR":
            modelid = 6
        elif self.project_type == "DALEC_GSI_FR_DBio":
            modelid = 7
        elif self.project_type == "DALEC_GSI_MFOL_FR":
            modelid = 8
        elif self.project_type == "DALEC_GSI_FR_LABILE":
            modelid = 9
        elif self.project_type == "DALECN_GSI_FR":
            modelid = 10
        elif self.project_type == "DALEC_GSI_DFOL_FR":
            modelid = 11
        elif self.project_type == "DALEC_GSI_DFOL_FROOT_FR":
            modelid = 12
        elif self.project_type == "DALEC_GSI_DFOL_LABILE_FR":
            modelid = 13
        elif self.project_type == "DALECN_GSI_DFOL_LABILE_FR":
            modelid = 14
        elif self.project_type == "DALECN_GSI_DFOL_LABILE_FROOT_FR":
            modelid = 15
        elif self.project_type == "DALEC_GSI_DFOL_CWD_FR":
            modelid = 16
        else:
            print "PROBLEM"
        for ii in xrange(self.details["no_pts"]):
            #create an empty array to store data to be written
            towrite=np.zeros(500+self.details["tsteps"]*(self.details["met_fields"]+self.details["obs_fields"]),dtype="d")-9999.
                    
            #provide fixed values
            towrite[0]=modelid                       # pixel number           
            towrite[2]=self.details["tsteps"]        # no of time steps
            towrite[3]=self.details["met_fields"]    # no of met fields
            towrite[4]=self.details["obs_fields"]    # no of obs fields
            towrite[5]=self.details["EDCs"]          # use EDCs? 1: Yes / 0: No

            #provide priors and met data

            #assign dummies to make code easier to read
            parprior = self.details["parprior"]
            parpriorunc = self.details["parpriorunc"]
            otherprior = self.details["otherprior"]
            otherpriorunc = self.details["otherpriorunc"]

            if self.details["no_pts"] == 1:           
                towrite[1]=self.details["latitude"]  # pixel latitude

                towrite[100:100+len(parprior)]=parprior
                towrite[200:200+len(parpriorunc)]=parpriorunc
                towrite[300:300+len(otherprior)]=otherprior
                towrite[400:400+len(otherpriorunc)]=otherpriorunc

                metobs=np.hstack([self.details["drivers"],self.details["observations"]])

            else:
                towrite[1]=self.details["latitude"][ii]  # pixel latitude

                towrite[100:100+len(parprior[ii])]=parprior[ii]
                towrite[200:200+len(parpriorunc[ii])]=parpriorunc[ii]
                towrite[300:300+len(otherprior[ii])]=otherprior[ii]
                towrite[400:400+len(otherpriorunc[ii])]=otherpriorunc[ii]

                metobs=np.hstack([self.details["drivers"][ii],self.details["observations"][ii]])

            towrite[500:]=metobs.flatten()

            #create binary data
            f=file(path2data+"%03i/%s_%05i.bin" % (self.runid,self.project_name,ii+1),'wb')
            f.write(struct.pack(len(towrite)*'d',*towrite))
            f.close()

        print "Written a total of %i/%i input file" % (ii+1,self.details["no_pts"])

    #-------------------------------------------------------------------------------------
    #  This method copies the source code in a project sub-directory
    def backup_source(self):
        if "src" not in os.listdir(self.paths["projects"]+self.project_name):
            os.mkdir(self.paths["projects"]+self.project_name+"/src")
        if ("%03i" % self.runid) not in os.listdir(self.paths["projects"]+self.project_name+'/src/'):
            os.mkdir(self.paths["projects"]+self.project_name+"/src/%03i" % self.runid)
        os.system("cp -r %s/* %s/%s/src/%03i/"  % (self.paths["library"],self.paths["projects"],self.project_name,self.runid))
        recompile_local = raw_input("Recompile local code <y/n>? ")
        if recompile_local=="y":
            self.compile_local_code()

    #-------------------------------------------------------------------------------------
    # This method backs up the source code in a project sub-directory and sends it to the cluster
    def update_source_cluster(self):
        self.backup_source()
        dest = self.paths["cluster_username"]+"@"+self.paths["cluster_address"]+":"+self.paths["cluster_directory"]+"/"+self.project_name
        os.system("scp -r %s/%s/src %s" % (self.paths["projects"],self.project_name,dest))

        recompile_cluster = raw_input("Recompile on cluster <y/n>? ")
        if recompile_cluster == "y":
            self.compile_cluster()

    #-------------------------------------------------------------------------------------
    # This method compiles the code and saves a backup
    def compile_local_code(self):

        path2model=self.paths["library"]+'model/'+self.project_type+'/src/'+self.project_type
        path2likelihood=self.paths["library"]+'model/'+self.project_type+'/likelihood/'
        path2misc=self.paths["library"]+'misc/'
        path2method=self.paths["library"]+'method/'
        path2general = self.paths["library"]+'general/'
        path2exe = self.paths["library"]+"executable/"

        #compile directly in good directory
        print "\nifort -O2  -xhost -ipo -no-ftz  %smath_functions.f90 \n%soksofar.f90 \n%s.f90 \n%s_CROP.f90 \n%scardamom_structures.f90 \n%sMHMCMC/MCMC_FUN/MHMCMC_STRUCTURES.f90 \n%s_PARS.f90 \n%scardamom_io.f90 \n%sMHMCMC/MCMC_FUN/MHMCMC.f90 \n%sMODEL_LIKELIHOOD.f90 \n%scardamom_main.f90 \n-o %s/cardamom.exe\n" % (path2misc,path2misc,path2model,path2model,path2general,path2method,path2model,path2general,path2method,path2likelihood,path2general,path2exe)

        os.system("ifort -O2  -xhost -ipo -no-ftz  %smath_functions.f90 %soksofar.f90 %s.f90 %s_CROP.f90 %scardamom_structures.f90 %sMHMCMC/MCMC_FUN/MHMCMC_STRUCTURES.f90 %s_PARS.f90 %scardamom_io.f90 %sMHMCMC/MCMC_FUN/MHMCMC.f90 %sMODEL_LIKELIHOOD.f90 %scardamom_main.f90 -o %s/cardamom.exe" % (path2misc,path2misc,path2model,path2model,path2general,path2method,path2model,path2general,path2method,path2likelihood,path2general,path2exe))

        print "\ngfortran -O2 -shared \n%s.f90 %s_CROP.f90 \n%s_R_interface.f90 \n-o %sdalec.so \n-fPIC\n" % (path2model,path2model,path2model,path2exe)
        
        os.system("gfortran -O2 -shared %s.f90 %s_CROP.f90 %s_R_interface.f90 -o %sdalec.so -fPIC" % (path2model,path2model,path2model,path2exe))
        
        os.system("mv *.mod %s" % path2exe)

    #-------------------------------------------------------------------------------------
    # This method compiles the code on the cluster
    def compile_cluster(self):
        
        path2cluster = self.paths["cluster_address"]
        path2source = self.paths["cluster_directory"]+self.project_name+"/src/"
        path2include = "%s/models/%s/likelihood/MODEL_LIKELIHOOD.c" % (path2source,self.project_type)
        path2exe = self.paths["cluster_directory"]+self.project_name+"/exec/"

        #os.system("ssh %s 'gcc -O3 %s/general/cardamom_main.c --include %s -o %s.exe -lm'" % (path2cluster,path2source,path2include,path2exe+self.project_name))
        os.system("ssh %s 'ifort -O2  %s/misc/math_functions.f90 %s/misc/oksofar.f90 %s/model/DALEC_GSI_DFOL_CWD_FR/src/DALEC_GSI_DFOL_CWD_FR.f90 %s/model/DALEC_GSI_DFOL_CWD_FR/src/DALEC_GSI_DFOL_CWD_FR_CROP.f90 %s/general/cardamom_structures.f90 %s/method/MHMCMC/MCMC_FUN/MHMCMC_STRUCTURES.f90 %s/model/DALEC_GSI_DFOL_CWD_FR/src/DALEC_GSI_DFOL_CWD_FR_PARS.f90 %s/general/cardamom_io.f90 %s/method/MHMCMC/MCMC_FUN/MHMCMC.f90 %s/model/DALEC_GSI_DFOL_CWD_FR/likelihood/MODEL_LIKELIHOOD.f90 %s/general/cardamom_main.f90 -o %s/cardamom.exe'" % (path2cluster,path2source,path2source,path2source,path2source,path2source,path2source,path2source,path2source,path2source,path2source,path2source,path2exe))

    #-------------------------------------------------------------------------------------
    # Redo the setup with attributes of the object
    def resetup(self):

        latitude      =  self.details["latitude"]
        longitude     =  self.details["longitude"]
        drivers       =  self.details["drivers"]
        observations  =  self.details["observations"]

        # MCMC specific values
        parprior      =  self.details["parprior"]
        parpriorunc   =  self.details["parpriorunc"]
        otherprior    =  self.details["otherprior"]
        otherpriorunc =  self.details["otherpriorunc"]

        print "Starting setup for project "+self.project_name

        self.setup(latitude,longitude,drivers,observations,parprior,parpriorunc,otherprior,otherpriorunc)

    #-------------------------------------------------------------------------------------
    # Download results from the cluster
    def download_cluster(self, **kwargs):
        
        print "Preparing to download data from \"%s\"" % self.paths["cluster_address"]

        cluster_details = "%s@%s" % (self.paths["cluster_username"],self.paths["cluster_address"])

        src = "%s/%s/output/*" % (self.paths["cluster_directory"],self.project_name)
        
        #create the output folder 
        if "output" not in os.listdir("%s/%s" % (self.paths["projects"],self.project_name)):
            os.mkdir("%s/%s/output" % (self.paths["projects"],self.project_name))

        runlist = os.listdir("%s/%s/output/" % (self.paths["projects"],self.project_name))
        if "runid" in kwargs:
            runid = kwargs["runid"]
        else:
            if len(runlist) == 0:
                runid = 1
            else:
                runlist.sort()
                runid = int(runlist[-1].split("_")[-1])+1


        dst = "%s/%s/output/run_%03i/" % (self.paths["projects"],self.project_name,runid)
        if "run_%03i" % runid not in runlist:
            os.mkdir(dst)
            print "scp -r %s:%s %s" % (cluster_details,src,dst)
            os.system("scp -r %s:%s %s" % (cluster_details,src,dst))
        else:
            if len(os.listdir(dst)) != 0.:
                notempty = raw_input("Destination folder \"%s\" not empty... continue <y/n>?" % dst)
                if notempty == "y":
                    print "Downloading data in \"%s\"" % dst
                    os.system("scp -r %s:%s %s" % (cluster_details,src,dst))
            else:
                os.system("scp -r %s:%s %s" % (cluster_details,src,dst))


    #-------------------------------------------------------------------------------------
    # run CARDAMOM locally for single point
    def run_CARDAMOM_single(self, **kwargs):
    
        if "cardamom_output"  not in os.listdir("%s/%s/" % (self.paths["projects"],self.project_name)):
            os.mkdir("%s/%s/cardamom_output/" % (self.paths["projects"],self.project_name))

        runlist = os.listdir("%s/%s/cardamom_output/" % (self.paths["projects"],self.project_name))
        if "runid" in kwargs:
            runid = kwargs["runid"]
        else:
            if len(runlist) == 0:
                runid = 1
            else:
                runlist.sort()
                runid = int(runlist[-1].split("_")[-1])+1
        if "run_%03i" % runid not in runlist:
            path_to_output = self.paths["projects"]+self.project_name+"/cardamom_output/%03i/" % (runid)
            print path_to_output
            os.system("mkdir %s" % (path_to_output))

        if "executable" in kwargs:
            executable = kwargs["executable"]
        else:
            executable = self.paths["library"]+"executable/cardamom.exe"

        if 'path_to_data' in kwargs:
            path_to_data = kwargs['path_to_data']
        else:
            path_to_data = self.paths["projects"]+self.project_name+"/data/%03i/" % (runid) 

        if 'accepted_params' in kwargs:
            accepted_params= kwargs['accepted_params']
        else:
            accepted_params = 10000000

        if 'printing_freq' in kwargs:
            printing_freq = kwargs['printing_freq']
        else:
            printing_freq = 0

        if 'sample_freq' in kwargs:
            sample_freq = kwargs['sample_freq']
        else:
            sample_freq = 10000

        if 'n_chains' in kwargs:
            n_chains = kwargs['n_chains']
        else:
            n_chains = 3

        
        if "target_point" in kwargs:
            target_point = kwargs["target_point"]
        else:
            target_point = 0

        #----------------------------
        # Run CARDAMOM for each point
        print 'Running CARDAMOM locally'
        print '\t- path to exe: ', self.paths["library"]+"executable/" 
        print '\t- path to data: ', path_to_data
        print '\t- path to output: ', path_to_output
        print '\t- number of accepted parameters: ', str(accepted_params)
        print '\t- printing frequency: ', str(printing_freq)
        print '\t- sample frequency: ', str(sample_freq)

        data_bin=path_to_data+"%s_%05i.bin" % (self.project_name,target_point)
        for cc in range(0,n_chains):
            output_prefix = path_to_output+"%s_%05i_%i_" % (self.project_name,target_point,cc+1)
            print "%s %s %s %s %s %s" % (executable,data_bin,output_prefix,str(accepted_params),str(printing_freq),str(sample_freq))
            os.system("%s %s %s %s %s %s" % (executable,data_bin,output_prefix,str(accepted_params),str(printing_freq),str(sample_freq)))

    #-------------------------------------------------------------------------------------
    # run CARDAMOM locally for all points
    def run_CARDAMOM_local(self, **kwargs):
    
        if "cardamom_output"  not in os.listdir("%s/%s/" % (self.paths["projects"],self.project_name)):
            os.mkdir("%s/%s/cardamom_output/" % (self.paths["projects"],self.project_name))

        runlist = os.listdir("%s/%s/cardamom_output/" % (self.paths["projects"],self.project_name))
        if "runid" in kwargs:
            runid = kwargs["runid"]
        else:
            if len(runlist) == 0:
                runid = 1
            else:
                runlist.sort()
                runid = int(runlist[-1].split("_")[-1])+1
        if "run_%03i" % runid not in runlist:
            path_to_output = self.paths["projects"]+self.project_name+"/cardamom_output/%03i/" % (runid)
            os.system("mkdir %s" % (path_to_output))

        if "executable" in kwargs:
            executable = kwargs["executable"]
        else:
            executable = self.paths["library"]+"executable/cardamom.exe"

        if 'path_to_data' in kwargs:
            path_to_data = kwargs['path_to_data']
        else:
            path_to_data = self.paths["projects"]+self.project_name+"/data/%03i/" % (runid)  

        if 'accepted_params' in kwargs:
            accepted_params= kwargs['accepted_params']
        else:
            accepted_params = 10000000

        if 'printing_freq' in kwargs:
            printing_freq = kwargs['printing_freq']
        else:
            printing_freq = 0

        if 'sample_freq' in kwargs:
            sample_freq = kwargs['sample_freq']
        else:
            sample_freq = 10000

        if 'n_chains' in kwargs:
            n_chains = kwargs['n_chains']
        else:
            n_chains = 3

        #----------------------------
        # Run CARDAMOM for each point
        print 'Running CARDAMOM locally'
        print '\t- path to exe: ', self.paths["library"]+"executable/cardamom.exe"
        print '\t- path to data: ', path_to_data
        print '\t- path to output: ', path_to_output
        print '\t- number of chains: ', str(n_chains)
        print '\t- number of accepted parameters: ', str(accepted_params)
        print '\t- printing frequency: ', str(printing_freq)
        print '\t- sample frequency: ', str(sample_freq)

        for ii in xrange(self.details["no_pts"]):
            data_bin=path_to_data+"%s_%05i.bin" % (self.project_name,ii+1)
            for cc in range(0,n_chains):
                output_prefix = path_to_output+"%s_%05i_%i_" % (self.project_name,ii+1,cc+1)
                os.system("%s %s %s %s %s %s &" % (executable,data_bin,output_prefix,str(accepted_params),str(printing_freq),str(sample_freq)))


    #-------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------
    # Plotting scripts
    #-------------------------------------------------------------------------------------
    # -first, plot up the drivers used to run CARDAMOM.
    #  There are 14 fields
    #  0 -> run day
    #  1 -> min temperature in oC
    #  2 -> max temperature
    #  3 -> ssrd in Mj.m-2.day-1
    #  4 -> atm CO2 in ppm
    #  5 -> day of year
    #  6 -> lagged pptn
    #  7 -> fire burned fraction
    #  8 -> deforestation fraction
    #  9 -> 21d average min temp in K
    # 10 -> 21d average max temp in K
    # 11 -> 21d average vpd in Pa
    # 12 -> forest management
    # 13 -> mean temperature in oC    
    # 
    # The 'site' parameter is an index to specify the site to be plotted up
    def plot_drivers(self,site):
        tstep = self.details["drivers"][site,:,0]

        plt.figure(1, facecolor='White',figsize=[10,10])
        
        # day of year
        ax1a = plt.subplot2grid((3,3),(0,0))
        ax1a.annotate('a - run day', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
        ax1a.set_ylabel('day',fontsize=axis_size)
        ax1a.set_xlabel('tstep',fontsize=axis_size)
        ax1a.plot(tstep,self.details["drivers"][site,:,5],'-',color=colour[0])

        # min, max & average temperature
        ax1b = plt.subplot2grid((3,3),(0,1),sharex = ax1a)
        ax1b.annotate('b - temperature', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
        ax1b.set_ylabel('temperature / $^o$C',fontsize=axis_size)
        ax1b.set_xlabel('tstep',fontsize=axis_size)
        ax1b.plot(tstep,self.details["drivers"][site,:,2],'-',color=colour[2])
        ax1b.plot(tstep,self.details["drivers"][site,:,1],'-',color=colour[1])
        ax1b.plot(tstep,self.details["drivers"][site,:,13],'-',color=colour[0])
    
        # ssrd
        ax1c = plt.subplot2grid((3,3),(0,2),sharex = ax1a)
        ax1c.annotate('c - ssrd', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
        ax1c.set_ylabel('ssrd / MJ.m$^{-2}$.d$^{-1}$',fontsize=axis_size)
        ax1c.set_xlabel('tstep',fontsize=axis_size)
        ax1c.plot(tstep,self.details["drivers"][site,:,3],'-',color=colour[0])

        # atmCO2
        ax2a = plt.subplot2grid((3,3),(1,0),sharex = ax1a)
        ax2a.annotate('d - CO$_2$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
        ax2a.set_ylabel('[CO$_2$] / ppm',fontsize=axis_size)
        ax2a.set_xlabel('tstep',fontsize=axis_size)
        ax2a.plot(tstep,self.details["drivers"][site,:,4],'-',color=colour[0])

        # pptn
        ax2b = plt.subplot2grid((3,3),(1,1),sharex = ax1a)
        ax2b.annotate('e - pptn', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
        ax2b.set_ylabel('pptn / mm',fontsize=axis_size)
        ax2b.set_xlabel('tstep',fontsize=axis_size)
        ax2b.plot(tstep,self.details["drivers"][site,:,6],'-',color=colour[1])

        # disturbance
        ax2c = plt.subplot2grid((3,3),(1,2),sharex = ax1a)
        ax2c.annotate('f - disturbance', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
        ax2c.set_ylabel('fraction loss',fontsize=axis_size)
        ax2c.set_xlabel('tstep',fontsize=axis_size)
        ax2c.plot(tstep,self.details["drivers"][site,:,7],'-',color=colour[2])
        ax2c.plot(tstep,self.details["drivers"][site,:,8],'-',color=colour[1])

        # 21d average T
        ax3a = plt.subplot2grid((3,3),(2,0),sharex = ax1a)
        ax3a.annotate('g - 21d average temp', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
        ax3a.set_ylabel('temperature / $^o$C',fontsize=axis_size)
        ax3a.set_xlabel('tstep',fontsize=axis_size)
        ax3a.plot(tstep,self.details["drivers"][site,:,9],'-',color=colour[1])
        #ax3a.plot(tstep,self.details["drivers"][site,:,10],'-',color=colour[2])

        # 21d average vpd
        ax3b = plt.subplot2grid((3,3),(2,1),sharex = ax1a)
        ax3b.annotate('h - 21d average VPD', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
        ax3b.set_ylabel('VPD / Pa',fontsize=axis_size)
        ax3b.set_xlabel('tstep',fontsize=axis_size)
        ax3b.plot(tstep,self.details["drivers"][site,:,11],'-',color=colour[0])

        # management
        ax3c = plt.subplot2grid((3,3),(2,2),sharex = ax1a)
        ax3c.annotate('i - management', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
        ax3c.set_ylabel('management',fontsize=axis_size)
        ax3c.set_xlabel('tstep',fontsize=axis_size)
        ax3c.plot(tstep,self.details["drivers"][site,:,12],'-',color=colour[0])
        
        plt.tight_layout()
        plt.show()

