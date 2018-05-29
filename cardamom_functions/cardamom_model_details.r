
###
## Function which contains basic information about the models
###

cardamom_model_details <-function(modelname,specific_pft,ctessel_pft) {
    if (modelname == "ACM") {
	# information contains is
	# The model name
	# Number of met parameters
	# Number of model parameters to be optimised
	nopools=array(2,dim=c(length(ctessel_pft)))
	nopars=array(18,dim=c(length(ctessel_pft)))
	nofluxes=array(3,dim=c(length(ctessel_pft)))
	cardamom_model_details=list(name="ACM",nopools=nopools,nofluxes=nofluxes,nomet=19+4,nopars=nopars)
    } else if (modelname == "DALEC_CDEA") {
	# information contains is
	# The model name
	# Number of met parameters
	# Number of model parameters to be optimised
	nopools=array(6,dim=c(length(ctessel_pft)))
	nopars=array(23,dim=c(length(ctessel_pft)))
	nofluxes=array(16,dim=c(length(ctessel_pft)))
	cardamom_model_details=list(name="DALEC_CDEA",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
   } else if (modelname == "DALEC_GSI_BUCKET") {
	# information contains is 
	# The model name
        # Number of met parameters
        # Number of model parameters to be optimised
        nopools=array(8,dim=c(length(ctessel_pft)))
        nopars=array(40,dim=c(length(ctessel_pft)))
        nofluxes=array(21,dim=c(length(ctessel_pft)))
	if (specific_pft == "pft_specific") {nopars[which(ctessel_pft == 1)]=35+2 ; nofluxes[which(ctessel_pft == 1)]=21 ; nopools[which(ctessel_pft == 1)]=9}
        cardamom_model_details=list(name="DALEC_GSI_BUCKET",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
   } else if (modelname == "DALECN_GSI_BUCKET") {
	# information contains is 
	# The model name
        # Number of met parameters
        # Number of model parameters to be optimised
        nopools=array(8,dim=c(length(ctessel_pft)))
        nopars=array(48,dim=c(length(ctessel_pft)))
        nofluxes=array(25,dim=c(length(ctessel_pft)))
	if (specific_pft == "pft_specific") {nopars[which(ctessel_pft == 1)]=35+2 ; nofluxes[which(ctessel_pft == 1)]=21 ; nopools[which(ctessel_pft == 1)]=9}
        cardamom_model_details=list(name="DALECN_GSI_BUCKET",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
    } else if (modelname == "AT_DALEC") {
	# information contains is
	# The model name
	# Number of met parameters
	# Number of model parameters to be optimised
	# set default parameter number for AT_DALEC
	nopools=array(6,dim=c(length(ctessel_pft)))
	nopars=array(22,dim=c(length(ctessel_pft)))
	nofluxes=array(16,dim=c(length(ctessel_pft)))
	# select for special case for crops
	if (specific_pft == "pft_specific") {nopars[which(ctessel_pft == 1)]=34 ; nofluxes[which(ctessel_pft == 1)]=16 ; nopools[which(ctessel_pft == 1)]=8}
	# output neeeded values
	cardamom_model_details=list(name="AT_DALEC",nopools=nopools,nofluxes=nofluxes,nomet=14,nopars=nopars)
    } else if (modelname == "DALEC_CDEA_FR") {
	# information contains is
	# The model name
	# Number of met parameters
	# Number of model parameters to be optimised
	nopools=array(6,dim=c(length(ctessel_pft)))
	nopars=array(23,dim=c(length(ctessel_pft)))
	nofluxes=array(18,dim=c(length(ctessel_pft)))
	cardamom_model_details=list(name="DALEC_CDEA_FR",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
   } else if (modelname == "DALECN_GSI_FR") {
	# information contains is
	# The model name
	# Number of met parameters
	# Number of model parameters to be optimised
	nopools=array(10,dim=c(length(ctessel_pft)))
	nopars=array(49,dim=c(length(ctessel_pft)))
	nofluxes=array(21,dim=c(length(ctessel_pft)))
	cardamom_model_details=list(name="DALECN_GSI_FR",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
   } else if (modelname == "DALEC_GSI_FR") {
	# information contains is
	# The model name
	# Number of met parameters
	# Number of model parameters to be optimised
	nopools=array(6,dim=c(length(ctessel_pft)))
	nopars=array(33,dim=c(length(ctessel_pft)))
	nofluxes=array(18,dim=c(length(ctessel_pft)))
	cardamom_model_details=list(name="DALEC_GSI_FR",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
   } else if (modelname == "DALEC_GSI_DFOL_FR") {
	# information contains is 
	# The model name
        # Number of met parameters
        # Number of model parameters to be optimised
        nopools=array(6,dim=c(length(ctessel_pft)))
        nopars=array(36,dim=c(length(ctessel_pft)))
        nofluxes=array(18,dim=c(length(ctessel_pft)))
	if (specific_pft == "pft_specific") {nopars[which(ctessel_pft == 1)]=35 ; nofluxes[which(ctessel_pft == 1)]=16 ; nopools[which(ctessel_pft == 1)]=8}
        cardamom_model_details=list(name="DALEC_GSI_DFOL_FR",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
   } else if (modelname == "DALEC_GSI_DFOL_CWD_FR") {
        # information contains is 
        # The model name
        # Number of met parameters
        # Number of model parameters to be optimised
        nopools=array(7,dim=c(length(ctessel_pft)))
        nopars=array(38,dim=c(length(ctessel_pft)))
        nofluxes=array(24,dim=c(length(ctessel_pft)))
	if (specific_pft == "pft_specific") {nopars[which(ctessel_pft == 1)]=35 ; nofluxes[which(ctessel_pft == 1)]=17 ; nopools[which(ctessel_pft == 1)]=8}
        cardamom_model_details=list(name="DALEC_GSI_DFOL_CWD_FR",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
   } else if (modelname == "DALEC_GSI_DFOL_LABILE_FR") {
	# information contains is 
	# The model name
        # Number of met parameters
        # Number of model parameters to be optimised
        nopools=array(9,dim=c(length(ctessel_pft)))
        nopars=array(44,dim=c(length(ctessel_pft)))
        nofluxes=array(19,dim=c(length(ctessel_pft)))
        cardamom_model_details=list(name="DALEC_GSI_DFOL_LABILE_FR",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
   } else if (modelname == "DALECN_GSI_DFOL_LABILE_FR") {
	# information contains is 
	# The model name
        # Number of met parameters
        # Number of model parameters to be optimised
        nopools=array(12,dim=c(length(ctessel_pft)))
        nopars=array(57,dim=c(length(ctessel_pft)))
        nofluxes=array(21,dim=c(length(ctessel_pft)))
        cardamom_model_details=list(name="DALECN_GSI_DFOL_LABILE_FR",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
   } else if (modelname == "DALECN_GSI_DFOL_LABILE_FROOT_FR") {
        # information contains is 
        # The model name
        # Number of met parameters
        # Number of model parameters to be optimised
        nopools=array(12,dim=c(length(ctessel_pft)))
        nopars=array(61,dim=c(length(ctessel_pft)))
        nofluxes=array(21,dim=c(length(ctessel_pft)))
        cardamom_model_details=list(name="DALECN_GSI_DFOL_LABILE_FROOT_FR",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
   } else if (modelname == "DALEC_GSI_DFOL_FROOT_FR") {
	# information contains is 
	# The model name
        # Number of met parameters
        # Number of model parameters to be optimised
        nopools=array(9,dim=c(length(ctessel_pft)))
        nopars=array(46,dim=c(length(ctessel_pft)))
        nofluxes=array(19,dim=c(length(ctessel_pft)))
        cardamom_model_details=list(name="DALEC_GSI_DFOL_FROOT_FR",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
   } else if (modelname == "DALEC_GSI_FR_LABILE") {
	# information contains is
	# The model name
	# Number of met parameters
	# Number of model parameters to be optimised
	nopools=array(7,dim=c(length(ctessel_pft)))
	nopars=array(37,dim=c(length(ctessel_pft)))
	nofluxes=array(20,dim=c(length(ctessel_pft)))
	cardamom_model_details=list(name="DALEC_GSI_FR_LABILE",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
   } else if (modelname == "DALEC_GSI_MFOL_FR") {
	# information contains is
	# The model name
	# Number of met parameters
	# Number of model parameters to be optimised
	nopools=array(7,dim=c(length(ctessel_pft)))
	nopars=array(36,dim=c(length(ctessel_pft)))
	nofluxes=array(18,dim=c(length(ctessel_pft)))
	cardamom_model_details=list(name="DALEC_GSI_MFOL_FR",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
   } else if (modelname == "DALEC_GSI_DBio_FR") {
        # information contains is
        # The model name
        # Number of met parameters
        # Number of model parameters to be optimised
        nopools=array(10,dim=c(length(ctessel_pft)))
        nopars=array(53,dim=c(length(ctessel_pft)))
        nofluxes=array(28,dim=c(length(ctessel_pft)))
        cardamom_model_details=list(name="DALEC_GSI_DBio_FR",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
    }
}
