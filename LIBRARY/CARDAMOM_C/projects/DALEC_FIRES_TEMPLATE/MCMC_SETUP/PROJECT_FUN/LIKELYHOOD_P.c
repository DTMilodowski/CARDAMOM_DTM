double LIKELYHOOD_P(DATA DATA,double *PARS)
{
/*remember - LOG likelihood*/
double p=0,p_lma,pn;
int n;

/*looping through all priors for P*/
for (n=0;n<50;n++){if (DATA.parpriors[n]>-9999){p=p-0.5*pow(log(PARS[n]/DATA.parpriors[n])/log(DATA.parpriorunc[n]),2);}
/*if (n==22){
printf("OBS.SOM = %f\n",DATA.parpriors[n]);
printf("OBS.SOMunc = %f\n",DATA.parpriorunc[n]);
printf("MOD.SOM = %f\n",PARS[n]);
printf("P = %f\n",-0.5*pow(log(PARS[n]/DATA.parpriors[n])/log(DATA.parpriorunc[n]),2));
}*/
}




/*for any other priors, explicitly define functions based on values in DATA.otherpriors*/


return p;}










