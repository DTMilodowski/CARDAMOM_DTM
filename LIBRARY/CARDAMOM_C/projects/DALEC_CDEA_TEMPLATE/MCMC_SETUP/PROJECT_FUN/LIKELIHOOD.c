double LIKELIHOOD(DATA D)
{
int n,dn;
double P=0;
double tot_exp;
double CPG,PP;


/*GPP LOG likelyhood*/
tot_exp=0;
if (D.ngpp>0){for (n=0;n<D.ngpp;n++){dn=D.gpppts[n];tot_exp+=pow((D.M_GPP[dn]-D.GPP[dn])/2,2);}
P=P-0.5*tot_exp;}

/*LAI likelyhood*/
tot_exp=0;
/*if (D.nlai>0){for (n=0;n<D.nlai;n++){dn=D.laipts[n];tot_exp+=pow((D.M_LAI[dn]-D.LAI[dn])/4,2);}*/
if (D.nlai>0){for (n=0;n<D.nlai;n++){dn=D.laipts[n];tot_exp+=pow(log(D.M_LAI[dn]/D.LAI[dn])/log(2),2);}
P=P-0.5*tot_exp;}

/*NEE likelyhood*/
tot_exp=0;
if (D.nnee>0){for (n=0;n<D.nnee;n++){dn=D.neepts[n];tot_exp+=pow((D.M_NEE[dn]-D.NEE[dn])/2,2);}
P=P-0.5*tot_exp;}

/*Cwood likelyhood*/
/*continue from HERE!!!!*/
/*fraction of woody pool (with uncertainty of 2)*/
tot_exp=0;
if (D.nwoo>0){
for (n=1;n<D.nwoo;n++){dn=D.woopts[n];tot_exp+=pow(log((D.M_POOLS[D.nopools*dn+3]/D.M_POOLS[D.nopools*D.woopts[0]+3])/D.WOO[dn])/((D.WOO[dn]-D.WOO[D.woopts[0]])*log(2)),2);}
/*
printf("CWM = %f\n",D.M_POOLS[D.nopools*3+3]);
printf("CWO = %f\n",D.WOO[dn]);
printf("P = %f\n", pow(log(D.M_POOLS[D.nopools*dn+3]/D.WOO[dn])/log(2),2));
printf("totexp = %f\n",tot_exp);*/
P=P-0.5*tot_exp;}




/*growth rate constraint - useful for large-scale assimilation*/
/*switch off for individual disturbed pools or where not applicable*/
if (D.otherpriors[0]>-9999){
for (n=0;n<D.nopools;n++){
CPG=D.M_POOLS[D.nodays*D.nopools+n]/D.M_POOLS[n];
PP=-0.5*pow(log(CPG/D.otherpriors[0])/log(D.otherpriorunc[0]),2);
P=P+PP;}}



/*LOG likelyhood*/
if (isnan(P)){P=log(0);}
return P;

}

