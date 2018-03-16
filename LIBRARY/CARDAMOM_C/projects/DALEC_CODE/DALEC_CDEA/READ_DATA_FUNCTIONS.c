
/*READING ALL FILES - FUNCTIONS GROUPED HERE!!!*/








int READ_LAIOBS(char site, double *MODISLAI, int nodays, double *info)
{
char laifilename[]="LAIOBS_";
strcat(laifilename,site);


/*THIS IS A BINARY FILE!!!*/
FILE *fid;
fid = fopen(laifilename, "rb");
printf("MODIS file pointer = %i\n",fid);

/*reading data here*/
int n, nolaipts=0;
fread(MODISLAI,sizeof(double),nodays,fid)

for (n=0;n<(int)nodays;n++){if (MODISLAI[n]>0){nolaipts=nolaipts+1;}}
info[4]=(double)nolaipts;

/*close file*/
fclose(fid);
/*close function*/
return 0;
}












int READ_METDATA(char site, double *MET, int nodays)
{
char metfilename[]="METDATA_";
strcat(metfilename,site);

FILE *fid;
fid = fopen(infofilename, "r");
printf("MET file pointer = %i\n",fid);

/*reading data here*/
int n;
float metval;
for (n=0;n<(int)nodays*6;n++){fscanf(fid,"%f",&metval);MET[n]=(double)metval;}
/*close file*/
fclose(fid);
/*close function*/
return 0;
}



int READ_RUNINFO(char site, double *info)
/*this file is ASCII - in order for user to make easy changes*/

{
char infofilename[]="RUNINFO_";
strcat(infofilename,site);
printf("Filename = %s\n",infofilename);
/*reading all run info, for use later*/
FILE *fid;
fid = fopen(infofilename, "r");
printf("INFO file pointer = %i\n",fid);
/*getting pointer for file*/
float hv1,hv2,hv3,hv4;
fscanf(fid,"%f",&hv1);info[0]=(double)hv1;
fscanf(fid,"%f",&hv2);info[1]=(int)hv2;
fscanf(fid,"%f",&hv3);info[2]=(double)hv3;
fscanf(fid,"%f",&hv4);info[3]=(int)hv4;

fclose(d);
return 0;
}




