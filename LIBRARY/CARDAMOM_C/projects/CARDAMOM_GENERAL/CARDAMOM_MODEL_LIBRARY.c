/*This function attributes model specific variables based
 * on the ID number stored in DATA->ID*/

int CARDAMOM_MODEL_LIBRARY(DATA *DATA){

/*ID = 1 - DALEC_CDEA*/
if (DATA->ID==1){
/*DALEC_CDEA - 6 pools*/
DATA->nopools=6;
DATA->nopars=23;
}


/*ID = 2 - DALEC_BUCKET*/
if (DATA->ID==2){
/*DALEC_CDEA - 6 pools*/
DATA->nopools=7;
DATA->nopars=26;
}

return 0;

}
