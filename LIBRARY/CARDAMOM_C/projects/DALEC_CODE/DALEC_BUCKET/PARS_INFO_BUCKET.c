
/*PARAMETER_INFO (typedef struct) must have at least 3 fields
 * npars,
 * parmax
 * parmin*/



int PARS_INFO_BUCKET(PARAMETER_INFO *PI){



PI->npars=26;

/*contains 6 fields with min max log for par and par*/


/*Decomposition rate*/
PI->parmin[0]=0.00001;
PI->parmax[0]=0.01;



/*Fraction of GPP respired*/
PI->parmin[1]=0.3;
PI->parmax[1]=0.7;


/*Fraction of (1-fgpp) to foliage*/
PI->parmin[2]=0.01;
PI->parmax[2]=0.5;


/*Fraction of (1-fgpp) to roots*/
PI->parmin[3]=0.01;
PI->parmax[3]=1;


/*Leaf Lifespan*/
/*Wright et al. 2004*/
PI->parmin[4]=1.001;
PI->parmax[4]=8;


/*TOR wood*/
PI->parmin[5]=0.00001;
PI->parmax[5]=0.005;


/*TOR roots*/
PI->parmin[6]=0.0001;
PI->parmax[6]=0.01;

/*TOR litter*/
PI->parmin[7]=0.0003;
PI->parmax[7]=0.01;

/*TOR SOM*/
PI->parmin[8]=0.00001;
PI->parmax[8]=0.003;


/*Temp factor* = Q10 = 1.2-1.6*/
PI->parmin[9]=0.018;
PI->parmax[9]=0.08;


/*Canopy Efficiency*/
PI->parmin[10]=1;
PI->parmax[10]=80;


/*Bday*/
PI->parmin[11]=365.25;
PI->parmax[11]=365.25*4;


/*Fraction to Clab*/
PI->parmin[12]=0.01;
PI->parmax[12]=0.5;


/*Clab Release period*/
PI->parmin[13]=10;
PI->parmax[13]=100;


/*Fday*/
PI->parmin[14]=365.25;
PI->parmax[14]=365.25*4;


/*Leaf fall period*/
PI->parmin[15]=20;
PI->parmax[15]=150;

/*LMA*/
/*Kattge et al. 2011*/
PI->parmin[16]=10;
PI->parmax[16]=400;




/*INITIAL VALUES DECLARED HERE*/



/*C labile*/
PI->parmin[17]=10.0;
PI->parmax[17]=1000.0;

/*C foliar*/
PI->parmin[18]=10.0;
PI->parmax[18]=1000.0;

/*C roots*/
PI->parmin[19]=10.0;
PI->parmax[19]=1000.0;

/*C_agb*/
PI->parmin[20]=10.0;
PI->parmax[20]=50000.0;

/*C litter*/
PI->parmin[21]=10.0;
PI->parmax[21]=1000.0;

/*C_som*/
PI->parmin[22]=10.0;
PI->parmax[22]=200000.0;

/*Hydrology parameters go here*/

/*water limitation factor*/
PI->parmin[23]=0.01;
PI->parmax[23]=10;

/*inherent water use efficiency*/
PI->parmin[24]=1e-1;
PI->parmax[24]=1e-6;


/*B_water - mm*/
PI->parmin[25]=1;
PI->parmax[25]=10000;



}
