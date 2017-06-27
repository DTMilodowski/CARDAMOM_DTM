# This defines the class of object CARDAMOM.
# The object stores the parameters for a given CARDAMOM project setup and provides methods
# to compile code, write input, execute the MCMC, read output, etc.
# ----------------------------------------------------------------------------------------
# The original version was written by J-F Exbrayat, and I have modified this slightly to
# work with Luke Smallman's versions of DALEC, which will be used during both the BALI and
# Forests2020 projects.
# ----------------------------------------------------------------------------------------
import os, cPickle, time, struct         
import datetime as dt       
import numpy as np
from netCDF4 import Dataset
import shutil, socket

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
        projtype=raw_input("Enter project number of type: \n1\tDALEC_GSI\n2\tDALEC_CDEA_LU_FIRES\n3\tDALEC_GSI_newalloc\n4\tDALEC_CDEA_LU_FIRES_LIU\n5\tDALEC_GSI_LIU\n6\tDALEC_CDEA_LU_FIRES_ET\n7\tDALEC_GSI_DFOL_CWD_FR_MHMCMC\n")
        if projtype != "":       
            self.project_type=projtype
            if self.project_type == "1":
                self.project_type = "DALEC_GSI"
            elif self.project_type == "2":
                self.project_type = "DALEC_CDEA_LU_FIRES"
            elif self.project_type == "3":
                self.project_type = "DALEC_GSI_newalloc"
            elif self.project_type == "4":
                self.project_type = "DALEC_CDEA_LU_FIRES_LIU"
            elif self.project_type == "5":
                self.project_type = "DALEC_GSI_LIU"
            elif self.project_type == "6":
                self.project_type = "DALEC_CDEA_LU_FIRES_ET"
            elif self.project_type == "7":
                self.project_type == "DALEC_GSI_DFOL_CWD_FR_MHMCMC" # this is Luke's model that includes CWD and a bucket model(?)

        else:
            self.project_type="CARDAMOM_LU"

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

        self.paths["library"] = raw_input("Where are source codes kept? (provide full path or leave blank to use subdirectory \"library\"): ")
        if self.paths["library"] == "":
            self.paths["library"] = self.paths["CARDAMOM"]+"/library/"

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
    def setup(self,latitude,longitude,drivers,observations,parprior,parpriorunc,otherprior,otherpriorunc):

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
            self.createInput()
            #copy the source code
            self.backup_source()

            #create a directory for the executable
            if "exec" not in os.listdir(self.paths["projects"]+self.project_name):
                os.mkdir(self.paths["projects"]+self.project_name+"/exec")

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

        if "data" not in os.listdir(path2project):
            print "Directory \"%s\" not found... creating" % path2data
            os.mkdir(path2data)

        print "Now creating input data in \"%s\" for project \"%s\"" % (path2data,self.project_name)

        if self.project_type == "DALEC_GSI":
            modelid = 1
        elif self.project_type == "DALEC_CDEA_LU_FIRES":
            modelid = 2
        elif self.project_type == "DALEC_GSI_newalloc":
            modelid = 3 
        elif self.project_type == "DALEC_CDEA_LU_FIRES_LIU":
            modelid = 4 
        elif self.project_type == "DALEC_GSI_LIU":
            modelid = 5 
        elif self.project_type == "DALEC_CDEA_LU_FIRES_ET":
            modelid = 6
        elif self.project_type == "DALEC_GSI_DFOL_CWD_FR_MHMCMC":
            modelid = 7

        for ii in xrange(self.details["no_pts"]):
            #create an empty array to store data to be written
            towrite=np.zeros(300+self.details["tsteps"]*(self.details["met_fields"]+self.details["obs_fields"]),dtype="d")-9999.
                    
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
                towrite[150:150+len(parpriorunc)]=parpriorunc
                towrite[200:200+len(otherprior)]=otherprior
                towrite[250:250+len(otherpriorunc)]=otherpriorunc

                metobs=np.hstack([self.details["drivers"],self.details["observations"]])

            else:
                towrite[1]=self.details["latitude"][ii]  # pixel latitude

                towrite[100:100+len(parprior[ii])]=parprior[ii]
                towrite[150:150+len(parpriorunc[ii])]=parpriorunc[ii]
                towrite[200:200+len(otherprior[ii])]=otherprior[ii]
                towrite[250:250+len(otherpriorunc[ii])]=otherpriorunc[ii]

                metobs=np.hstack([self.details["drivers"][ii],self.details["observations"][ii]])

            towrite[300:]=metobs.flatten()

            #create binary data
            f=file(path2data+"%s_%05i.bin" % (self.project_name,ii+1),'wb')
            f.write(struct.pack(len(towrite)*'d',*towrite))
            f.close()

        print "Written a total of %i/%i input file" % (ii+1,self.details["no_pts"])

    #-------------------------------------------------------------------------------------
    #  This method copies the source code in a project sub-directory
    def backup_source(self):
        if "src" not in os.listdir(self.paths["projects"]+self.project_name):
            os.mkdir(self.paths["projects"]+self.project_name+"/src")
        os.system("cp -r %s/* %s/%s/src" % (self.paths["library"],self.paths["projects"],self.project_name))
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
    # This method compiles the code and saves a backup #### CHANGE TO COMPILE LUKE'S FORTRAN CODE!!!
    def compile_local_code(self):

        path2source=self.paths["projects"]+self.project_name+"/src/"
        path2include="%s/models/%s/likelihood/MODEL_LIKELIHOOD.c" % (path2source,self.project_type)
        path2exe = self.paths["projects"]+self.project_name+"/exec/"

        if "exec" not in os.listdir(self.paths["projects"]+self.project_name):
            os.mkdir(self.paths["projects"]+self.project_name+"/exec")
        #compile directly in good directory
        os.system("gcc -O3 %s/general/cardamom_main.c --include %s -o %s.exe -lm" % (path2source,path2include,path2exe+self.project_name))
                 
    #-------------------------------------------------------------------------------------
    # This method compiles the code on the cluster #### CHANGE TO COMPILE LUKE'S FORTRAN CODE!!!
    def compile_cluster(self):
        
        path2cluster = self.paths["cluster_address"]
        path2source = self.paths["cluster_directory"]+self.project_name+"/src/"
        path2include = "%s/models/%s/likelihood/MODEL_LIKELIHOOD.c" % (path2source,self.project_type)
        path2exe = self.paths["cluster_directory"]+self.project_name+"/exec/"

        os.system("ssh %s 'gcc -O3 %s/general/cardamom_main.c --include %s -o %s.exe -lm'" % (path2cluster,path2source,path2include,path2exe+self.project_name))

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
