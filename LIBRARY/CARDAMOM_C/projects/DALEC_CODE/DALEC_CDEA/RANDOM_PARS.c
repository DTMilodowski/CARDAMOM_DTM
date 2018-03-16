/*function for sampling random parameters*/






double SAMPLEINIV(const struct  parameterinfo *PI, const int n,double minv,double maxv)
{

if (minv<PI->inimin[n]){minv=PI->inimin[n];}
if (maxv<PI->inimin[n]){maxv=PI->inimin[n];}
if (maxv>PI->inimax[n]){maxv=PI->inimax[n];}
if (minv>PI->inimax[n]){minv=PI->inimax[n];}
/*absurd*/
if (maxv<minv){
printf("iniv = %d - maxv = %f minv = %f\n",n,maxv,minv);
}

double r=(double)random()/RAND_MAX;
double I;
if (PI->inilog[n]==0){I=(maxv-minv)*r+minv;}
else{I=exp((log(maxv)-log(minv))*r+log(minv));}
return I;
}



/*function for sampling random initial values*/

double SAMPLEPARS(const struct parameterinfo *PI, const int n,double minv,double maxv)
{


if (minv<PI->parmin[n]){minv=PI->parmin[n];}
if (maxv>PI->parmax[n]){maxv=PI->parmax[n];}
if (maxv<minv){printf("maxminmismatch --- np = %d",n);}

double P;
double r=(double)random()/RAND_MAX;
if (PI->parlog[n]==0){P=(maxv-minv)*r+minv;}
else{P=exp((log(maxv)-log(minv))*r+log(minv));}
return P;
}




/* must feed it RAND and INIV pointers, as well as parameter ranges*/


int DALEC_CDEA_RANDOM_PARS_RC(const struct parameterinfo *PI, double *P, double *I,double const maxtemp, double const maxrad,double const lat)
{




int ni,np,n;
double r;
double gppmax=100;
double gppmin=0.1;
double laimax=5;

/*SAMPLING PARAMETERS W. Reality Checks*/


/*LMA*/
np=16;
P[np]=SAMPLEPARS(PI,np,PI->parmin[np],PI->parmax[np]);


/*Ceff -> as function of LMA*/
/*sampling centred at LMA/25*/
np=10;
P[np]=SAMPLEPARS(PI,np,PI->parmin[np],PI->parmax[np]);




/*Autotrophic*/
np=1;
P[np]=SAMPLEPARS(PI,np,PI->parmin[np],PI->parmax[np]);


/*Bday*/ 
np=11;
P[np]=SAMPLEPARS(PI,np,PI->parmin[np],PI->parmax[np]);

/*labrelease period*/
np=13;
P[np]=SAMPLEPARS(PI,np,PI->parmin[np],PI->parmax[np]);


/*leaffall period*/
np=15;
P[np]=SAMPLEPARS(PI,np,PI->parmin[np],PI->parmax[np]);


/*Trate*/
np=9;
P[np]=SAMPLEPARS(PI,np,PI->parmin[np],PI->parmax[np]);

/*PARAMETERS AS FUNCTIONS OF OTHERS*/


/*ffoliar*/
np=2;
P[np]=SAMPLEPARS(PI,np,PI->parmin[np],PI->parmax[np]);


/*flab -> as a function of ffol*/
np=12;
P[np]=SAMPLEPARS(PI,np,PI->parmin[np],(5.0 - ffol/(1-fauto-ffol))/6.0);


/*Fday -> AS FUNCTION OF Bday*/
np=14;
P[np]=(int)SAMPLEPARS(PI,np,P[11]+121,P[11]+365-121) % (int)365;





/*Clab -> as a function of Flab*/
ni=0;
I[ni]=SAMPLEINIV(PI,ni,minclab,maxclab);


/*leaf lifespan -> as a function of ffol flab gppmax and Cfol*/
np=4;
P[np]=SAMPLEPARS(PI,np,PI->parmin[np],PI->parmax[np]);

/*C_foliar -> AS FUNCTION OF LMA and LAImax=15*/
ni=1;
I[ni]=SAMPLEINIV(PI,ni,PI->inimin[ni],folcmax);


/*re-deriving LAImax*/
laimax=(I[0]+I[1])*10/P[16];
if (laimax>15){laimax=15;}



/*froot -> as a function of ffol and flab*/
np=3;
P[np]=SAMPLEPARS(PI,np,0.2*(ffol+flab)/(1-fauto-ffol-flab),5*(ffol+flab)/(1-fauto-ffol-flab));


/*troot -> as a function of froot Ciroot*/
np=6;
P[np]=SAMPLEPARS(PI,np,PI->parmin[np],PI->parmax[np]);

/*Croot -> AS FUNCTION OF Cfoliar*/
ni=2;
I[ni]=SAMPLEINIV(PI,ni,PI->inimin[ni],maxcroo);
/*printf("I[5] = %f\n",I[5]);
printf("I[6] = %f\n",I[6]);*/





/*TURNOVERS (as functions of others)*/



/*twood* -> as a function of tfol*/
np=5;
P[np]=SAMPLEPARS(PI,np,PI->parmin[np],torfol);





/*Cwood -> AS FUNCTION OF torwood*/
ni=3;
I[ni]=SAMPLEINIV(PI,ni,0.1*fwood*gppmin/P[5],fwood*gppmax/P[5]);

/*tsom* ->As FUNCTION of twood*/
np=8;
if (mintor>P[6]){mintor=P[6];}
P[np]=SAMPLEPARS(PI,np,PI->parmin[np],mintor/trmed);

/*Mrate*/
np=0;
P[np]=SAMPLEPARS(PI,np,P[8],PI->parmax[np]);

/*tlitter*/
np=7;
P[np]=SAMPLEPARS(PI,np,P[8],PI->parmax[np]);

/*Clitter*/
/*max litter = max production / min decomposition rate*/
/*also, since values for other pools are chosen, now min and max Clit are derivable*/
ni=4;
double minclit=(torfol*(I[0]+I[1])*0.1+P[6]*I[2]*0.1)/(trmax*(P[7]+P[0]));
double maxclit=(torfol*P[16]*laimax+P[6]*I[2]*10)/(trmin*(P[7]+P[0]));
I[ni]=SAMPLEINIV(PI,ni,minclit,maxclit);

/*Csom*/
ni=5;
double mincsom=(I[4]*trmin*P[0]*0.1 + I[3]*0.1*P[5])/(trmax*P[8]);
double maxcsom=(P[0]*1000*trmax + 50000*P[5])/(P[8]*trmin);
I[ni]=SAMPLEINIV(PI,ni,mincsom,maxcsom);





/*
for (n=0;n<PI->numpars;n++){SAMPLEPARS(PI,P,I,n);}


for (n=0;n<PI->numiniv;n++){SAMPLEINIV(PI,P,I,n);}
*/
return 0;

}
