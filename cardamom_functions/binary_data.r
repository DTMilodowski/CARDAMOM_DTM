###
## Function to create binary input files
###

# /*TEMPLATE FOR ALL DALEC MCMC DATA files*/
# /*Static Elements: 1-100 - use as many as needed*/
# /*Parameter Priors: 101-150*/
# /*Parameter prior uncertainty: 151-200*/
#/*Other priors & uncertainties: 201-300*/
# /*TEMPORAL DRIVERS & DATA: 301-end*/
 
binary_data<-function(met,OBS,file,EDC,latlon_in,ctessel_pft,modelname,parameter_type,nopars) {
    print(paste("writing out binary...",Sys.time(),sep=""))

    # set model ID
    if (modelname == "ACM") {
	modelid=0
    } else if (modelname == "DALEC_CDEA") {
	modelid=1
    } else if (modelname == "DALEC_CDEA_FR") {
	modelid=5
    } else if (modelname == "DALECN_GSI_FR") {
	modelid=10
    } else if (modelname == "DALEC_GSI_FR") {
	modelid=6
    } else if (modelname == "DALEC_GSI_FR_LABILE") {
	modelid=9
    } else if (modelname == "DALEC_GSI_MFOL_FR") {
	modelid=8
    } else if (modelname == "DALEC_GSI_DFOL_FR") {
        modelid=11
    } else if (modelname == "DALEC_GSI_DFOL_FROOT_FR") {
        modelid=12
    } else if (modelname == "DALEC_GSI_DFOL_LABILE_FR") {
        modelid=13
    } else if (modelname == "DALECN_GSI_DFOL_LABILE_FR") {
        modelid=14
    } else if (modelname == "DALECN_GSI_DFOL_LABILE_FROOT_FR") {
        modelid=15
    } else if (modelname == "DALEC_GSI_DFOL_CWD_FR") {
        modelid=16
    } else if (modelname == "DALEC_GSI_DBio_FR") {
	modelid=7
    } else if (modelname == "DALEC_GSI_BUCKET"){
	modelid=2
    } else if (modelname == "DALECN_GSI_BUCKET"){
	modelid=17
    } else if (modelname == "AT_DALEC" & parameter_type == "pft_specific" & ctessel_pft == 1){
	# i.e. crop model
	modelid=4
    } else if (modelname == "AT_DALEC"){
	# i.e. default AT_DALEC
	modelid=3
    }

    # some drivers may be passed as single values assuming this will apply across the whole time series
    # if this is so we need to now adjust that
    if (length(OBS$forest_management) != length(met$run_day)) {
	# we will assume that this is a constant value and will now repeast it
	OBS$forest_management=array(OBS$forest_management, dim=c(length(met$run_day)))
    }


###
## REALLY BAD HACK
###

#if (grepl("Duke",file)) {
#    met$mint=met$mint-1.75
#    met$maxt=met$maxt-2.00
#    met$swrad=met$swrad-0.8698
#} else if (grepl("Harwood",file)){
#    met$mint=met$mint-1.73
#    met$maxt=met$maxt-1.73
#}

#par(mfrow=c(2,3))
#plot(met$mint)
#plot(met$maxt)
#plot(met$swrad)
#plot(met$precip)

    # extract information from list to array
    if (modelname == "ACM") {
	MET=array(NA,dim=c(length(met$run_day),(length(met)+2)))
    } else {
	MET=array(NA,dim=c(length(met$run_day),(length(met)+3)))
    }

    MET[,1] = met$run_day
    MET[,2] = met$mint  ; if (min(met$mint) < -200) {stop('mint error in binary_data')}
    MET[,3] = met$maxt  ; if (min(met$maxt) < -200) {stop('maxt error in binary_data')}
    MET[,4] = met$swrad ; if (min(met$swrad) < 0) {stop('RAD error in binary_data')}
    MET[,5] = met$co2#+200
    MET[,6] = met$doy
    MET[,7] = pmax(0,met$precip)
    MET[,8] = OBS$deforestation
    MET[,9] = OBS$burnt_area
    MET[,10] = met$avgTmin
    MET[,11] = met$photoperiod
    MET[,12] = met$vpd_lagged
    MET[,13] = OBS$forest_management
    if (modelname == "ACM") {
	MET[,14] = met$avgN
	MET[,15] = met$lai
	MET[,16] = met$lat
	MET[,17] = met$wind_spd
	MET[,18] = met$vpd
	MET[,19] = met$Rtot
	MET[,20] = met$top_sand
	MET[,21] = met$bot_sand
	MET[,22] = met$top_clay
	MET[,23] = met$bot_clay
    } else {
	MET[,14] = met$avgTemp
	MET[,15] = met$wind_spd
	MET[,16] = met$vpd
    }

    # TEMPLATE FOR ALL DALEC MCMC DATA files
    # Static Elements: 1-100 - use as many as needed
    # Parameter Priors: 101-150
    # Parameter prior uncertainty: 151-200
    # Other priors & uncertainties: 201-300
    # TEMPORAL DRIVERS & DATA: 301-end
    # construct obs/drivers
    OBSMAT=array(-9999.0,dim=c(length(met$run_day),32))
    #OBSMAT=MET[,1:9]*0-9999.0 # all rows, first 3 columns, which are 1= day of simulation run, 2= mint, 3= maxt, 4 =RAD, 5= CO2, 6=doy
			    # line makes the correct array size but with -9999 in place of all
    OBSMAT[,1]=OBS$GPP      # loads into column 1 GPP, 2 LAI, 3 NEE, 4 Woodinc, 5 Reco obs where they exist, leaving the rest -9999
    OBSMAT[,2]=OBS$LAI
    OBSMAT[,3]=OBS$NEE
    OBSMAT[,4]=OBS$woodinc
    OBSMAT[,5]=OBS$Reco
    OBSMAT[,6]=OBS$Cfol_stock
    OBSMAT[,7]=OBS$Cwood_stock
    OBSMAT[,8]=OBS$Croots_stock
    OBSMAT[,9]=OBS$Clit_stock
    OBSMAT[,10]=OBS$Csom_stock
    OBSMAT[,11]=OBS$Cagb_stock
    # now for the corresponding uncertainty values
    OBSMAT[,12]=OBS$GPP_unc 
    OBSMAT[,13]=OBS$LAI_unc
    OBSMAT[,14]=OBS$NEE_unc
    OBSMAT[,15]=OBS$woodinc_unc
    OBSMAT[,16]=OBS$Reco_unc
    OBSMAT[,17]=OBS$Cfol_stock_unc
    OBSMAT[,18]=OBS$Cwood_stock_unc
    OBSMAT[,19]=OBS$Croots_stock_unc
    OBSMAT[,20]=OBS$Clit_stock_unc
    OBSMAT[,21]=OBS$Csom_stock_unc
    OBSMAT[,22]=OBS$Cagb_stock_unc
    OBSMAT[,23]=OBS$Cstem_stock
    OBSMAT[,24]=OBS$Cstem_stock_unc
    OBSMAT[,25]=OBS$Cbranch_stock
    OBSMAT[,26]=OBS$Cbranch_stock_unc
    OBSMAT[,27]=OBS$Ccoarseroot_stock
    OBSMAT[,28]=OBS$Ccoarseroot_stock_unc
    OBSMAT[,29]=OBS$Cfolmax_stock
    OBSMAT[,30]=OBS$Cfolmax_stock_unc
    OBSMAT[,31]=OBS$Evap
    OBSMAT[,32]=OBS$Evap_unc
    DATA_TEMP=t(cbind(MET,OBSMAT))

    #STATIC DATA
    # Model ID = static_data[1]; DALEC_CDEA, DALEC_BUCKET etc
    # LAT=static_data[2]; Latitude of site(Degrees)
    # nodays=static_data[3]; Number of days (or time steps) in simulation
    # nomet=static_data[4]; Number of met variables
    # noobs=static_data[5]; Number of observation streams
    # EDC=static_data[6]; EDCs on (1) or off (0)
    # pft=static_data[7]; CTESSEL plant functional type, only used by ACM_TESSEL

    # if force_random_search == 1 then CARDAMOM ignores parameter priors even if present in the file during the EDC initialisation
    force_random_search = -9999
    # pass static information
    static_data=rep(-9999.0,length.out=100)
    tmp=c(modelid,latlon_in[1],dim(MET)[1],dim(MET)[2],dim(OBSMAT)[2],EDC,ctessel_pft,OBS$yield_class,OBS$age,nopars,force_random_search,
	  OBS$top_sand[1],OBS$bot_sand[1],OBS$top_clay[1],OBS$bot_clay[1])
    static_data[1:length(tmp)]=tmp

    #ONLY USED FOR LOG NORMALLY PSERIBUTED PARAMETER PRIORS
    PARPRIORS=rep(-9999.0,length.out=100)
    PARPRIORUNC=rep(-9999.0,length.out=100)
    #For all other multiparameter user-defined priors
    OTHERPRIORS=rep(-9999.0,length.out=100)
    OTHERPRIORUNC=rep(-9999.0,length.out=100)

    #introducing some commonly used priors (assumes DALEC_CDEA default parameters)
    #PARPRIORS(2)=0.5;PARPRIORUNC(2)=1.2;%P_AUTO
    #PARPRIORS(10)=0.03;PARPRIORUNC(10)=1.25;%Temp_rate (Mahecha 2010)
    #PARPRIORS(17)=70;PARPRIORUNC(17)=2;%LMA - Kattge 2011
    if (modelname == "DALEC_CDEA") {
        PARPRIORS[11]=20.52048   ; PARPRIORUNC[11]=1.617705 # Ceff
	PARPRIORS[19]=OBS$Cfol_initial 		; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[19]=2.0} # Cfoliar prior
	PARPRIORS[20]=OBS$Croots_initial 	; if (OBS$Croots_initial != -9999) {PARPRIORUNC[20]=2.0} # Croots prior
	PARPRIORS[21]=OBS$Cwood_initial 	; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=2.0} # Cwood prior
	PARPRIORS[22]=OBS$Clit_initial 		; if (OBS$Clit_initial != -9999) {PARPRIORUNC[22]=2.0} # Clitter prior
	PARPRIORS[23]=OBS$SOM			; if (OBS$SOM != -9999) {PARPRIORUNC[23]=2.0} # Csom prior
    } else if (modelname == "AT_DALEC" & parameter_type == "pft_specific" & ctessel_pft == 1) {
#        PARPRIORS[11]=20.52048   ; PARPRIORUNC[11]=1.617705 # Ceff
	PARPRIORS[19]=OBS$Cfol_initial 		; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[19]=2.0} # Cfoliar prior
	PARPRIORS[20]=OBS$Croots_initial 	; if (OBS$Croots_initial != -9999) {PARPRIORUNC[20]=2.0} # Croots prior
	PARPRIORS[21]=OBS$Cwood_initial 	; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=2.0} # Cwood prior
	PARPRIORS[22]=OBS$Clit_initial 		; if (OBS$Clit_initial != -9999) {PARPRIORUNC[22]=2.0} # Clitter prior
	PARPRIORS[23]=OBS$SOM 			; if (OBS$SOM != -9999) {PARPRIORUNC[23]=2.0} # Csom prior
    } else if (modelname == "AT_DALEC") {
#        PARPRIORS[11]=20.52048   ; PARPRIORUNC[11]=1.617705 # Ceff
	PARPRIORS[18]=OBS$Cfol_initial 		; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[18]=2.0} # Cfoliar prior
	PARPRIORS[19]=OBS$Croots_initial 	; if (OBS$Croots_initial != -9999) {PARPRIORUNC[19]=2.0} # Croots prior
	PARPRIORS[20]=OBS$Cwood_initial 	; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[20]=2.0} # Cwood prior
	PARPRIORS[21]=OBS$Clit_initial 		; if (OBS$Clit_initial != -9999) {PARPRIORUNC[21]=2.0} # Clitter prior
	PARPRIORS[22]=OBS$SOM 			; if (OBS$SOM != -9999) {PARPRIORUNC[22]=2.0} # Csom prior
    } else if (modelname == "DALEC_CDEA_FR") {
#        PARPRIORS[11]=20.52048   ; PARPRIORUNC[11]=1.617705 # Ceff
	PARPRIORS[17]=140.0   ; PARPRIORUNC[17]=1.5 # LMA gC.m-2 prior (Duke Forest; Akers et al 2013)
	PARPRIORS[19]=OBS$Cfol_initial 		; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[19]=2.0} # Cfoliar prior
	PARPRIORS[20]=OBS$Croots_initial 	; if (OBS$Croots_initial != -9999) {PARPRIORUNC[20]=2.0} # Croots prior
	PARPRIORS[21]=OBS$Cwood_initial 	; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=2.0} # Cwood prior
	PARPRIORS[22]=OBS$Clit_initial 		; if (OBS$Clit_initial != -9999) {PARPRIORUNC[22]=2.0} # Clitter prior
	PARPRIORS[23]=OBS$SOM 			; if (OBS$SOM != -9999) {PARPRIORUNC[23]=2.0} # Csom prior
    } else if (modelname == "DALEC_GSI_DBio_FR") {
        PARPRIORS[11]=0.2764618		; PARPRIORUNC[11]=0.2014871 # log10 avg foliar N (gN.m-2)
        PARPRIORS[17]=-9999   ; PARPRIORUNC[17]=-9999 # LMA 
        PARPRIORS[19]=OBS$Cfol_initial          ; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[19]=2.0} # Cfoliar prior
        PARPRIORS[20]=OBS$Croots_initial        ; if (OBS$Croots_initial != -9999) {PARPRIORUNC[20]=2.0} # Croots prior
        PARPRIORS[21]=OBS$Cwood_initial         ; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=OBS$Cwood_initial_unc} # Cwood prior
#       PARPRIORS[21]=OBS$Cwood_initial         ; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=2.0} # Cwood prior
        PARPRIORS[22]=OBS$Clit_initial          ; if (OBS$Clit_initial != -9999) {PARPRIORUNC[22]=2.0} # Clitter prior
        PARPRIORS[23]=OBS$SOM                   ; if (OBS$SOM != -9999) {PARPRIORUNC[23]=2.0} # Csom priors
    } else if (modelname == "DALECN_GSI_FR") {
        PARPRIORS[11]=0.2432501   ; PARPRIORUNC[11]=0.2014871 # log10 avg foliar N (gN.m-2)
        PARPRIORS[49]=0.001868948 ; PARPRIORUNC[49]=0.0005951156 # NUE**(1/-2.999299929993) (gC/gN)
#        PARPRIORS[17]=140.0   ; PARPRIORUNC[17]=11.7*3 # (not log) LMA gC.m-2 prior (Duke Forest; Akers et al 2013)
	PARPRIORS[17]=-9999   ; PARPRIORUNC[17]=-9999 # (not log) LMA gC.m-2 prior (Duke Forest; Akers et al 2013)
	PARPRIORS[19]=OBS$Cfol_initial 		; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[19]=2.0} # Cfoliar prior
	PARPRIORS[20]=OBS$Croots_initial 	; if (OBS$Croots_initial != -9999) {PARPRIORUNC[20]=2.0} # Croots prior
	PARPRIORS[21]=OBS$Cwood_initial 	; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=2.0} # Cwood prior
	PARPRIORS[22]=OBS$Clit_initial 		; if (OBS$Clit_initial != -9999) {PARPRIORUNC[22]=2.0} # Clitter prior
	PARPRIORS[23]=OBS$SOM 			; if (OBS$SOM != -9999) {PARPRIORUNC[23]=2.0} # Csom prior
    } else if (modelname == "DALEC_GSI_FR") {
#        PARPRIORS[11]=20.52048   ; PARPRIORUNC[11]=1.617705 # Ceff
	PARPRIORS[17]=-9999   ; PARPRIORUNC[17]=-9999 # (not log) LMA gC.m-2 prior (Duke Forest; Akers et al 2013)
	PARPRIORS[19]=OBS$Cfol_initial 		; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[19]=2.0} # Cfoliar prior
	PARPRIORS[20]=OBS$Croots_initial 	; if (OBS$Croots_initial != -9999) {PARPRIORUNC[20]=2.0} # Croots prior
	PARPRIORS[21]=OBS$Cwood_initial 	; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=2.0} # Cwood prior
	PARPRIORS[22]=OBS$Clit_initial 		; if (OBS$Clit_initial != -9999) {PARPRIORUNC[22]=2.0} # Clitter prior
	PARPRIORS[23]=OBS$SOM 			; if (OBS$SOM != -9999) {PARPRIORUNC[23]=2.0} # Csom prior
    } else if (modelname == "DALEC_GSI_DFOL_FR") {
        PARPRIORS[11]=20.52048   ; PARPRIORUNC[11]=1.2 #1.617705 # Ceff 
	PARPRIORS[17]=-9999   ; PARPRIORUNC[17]=-9999 # LMA 
	PARPRIORS[19]=OBS$Cfol_initial 		; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[19]=2.0} # Cfoliar prior
	PARPRIORS[20]=OBS$Croots_initial 	; if (OBS$Croots_initial != -9999) {PARPRIORUNC[20]=2.0} # Croots prior
	PARPRIORS[21]=OBS$Cwood_initial 	; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=OBS$Cwood_initial_unc} # Cwood prior
#	PARPRIORS[21]=OBS$Cwood_initial 	; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=2.0} # Cwood prior
	PARPRIORS[22]=OBS$Clit_initial 		; if (OBS$Clit_initial != -9999) {PARPRIORUNC[22]=2.0} # Clitter prior
	PARPRIORS[23]=OBS$SOM 			; if (OBS$SOM != -9999) {PARPRIORUNC[23]=2.0} # Csom prior
    } else if (modelname == "DALEC_GSI_DFOL_CWD_FR") {
        #PARPRIORS[11]=20.52048   ; PARPRIORUNC[11]=1.6 #1.617705 # Ceff 
        PARPRIORS[11]=0.2764618		; PARPRIORUNC[11]=0.2014871*0.5 # log10 avg foliar N (gN.m-2)
	if (parameter_type == "pft_specific" & ctessel_pft == 1) {
	      PARPRIORS[12]=OBS$plant   ; PARPRIORUNC[12]=1.1
	      PARPRIORS[15]=OBS$harvest ; PARPRIORUNC[15]=1.1
	}
        PARPRIORS[17]=-9999   ; PARPRIORUNC[17]=-9999 # LMA 
        PARPRIORS[19]=OBS$Cfol_initial          ; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[19]=2.0} # Cfoliar prior
        PARPRIORS[20]=OBS$Croots_initial        ; if (OBS$Croots_initial != -9999) {PARPRIORUNC[20]=2.0} # Croots prior
	if (parameter_type == "pft_specific" & ctessel_pft == 1) {
      
	} else {
	    PARPRIORS[21]=OBS$Cwood_initial         ; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=OBS$Cwood_initial_unc} # Cwood prior
    #       PARPRIORS[21]=OBS$Cwood_initial         ; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=2.0} # Cwood prior
	}
        PARPRIORS[22]=OBS$Clit_initial          ; if (OBS$Clit_initial != -9999) {PARPRIORUNC[22]=2.0} # Clitter prior
        PARPRIORS[23]=OBS$SOM                   ; if (OBS$SOM != -9999) {PARPRIORUNC[23]=2.0} # Csom prior
    } else if (modelname == "DALEC_GSI_DFOL_LABILE_FR") {
        PARPRIORS[11]=20.52048   ; PARPRIORUNC[11]=1.6 #1.617705 # Ceff 
	PARPRIORS[17]=-9999   ; PARPRIORUNC[17]=-9999 # LMA 
	PARPRIORS[19]=OBS$Cfol_initial 		; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[19]=2.0} # Cfoliar prior
	PARPRIORS[20]=OBS$Croots_initial 	; if (OBS$Croots_initial != -9999) {PARPRIORUNC[20]=2.0} # Croots prior
	PARPRIORS[21]=OBS$Cwood_initial 	; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=2.0} # Cwood prior
	PARPRIORS[22]=OBS$Clit_initial 		; if (OBS$Clit_initial != -9999) {PARPRIORUNC[22]=2.0} # Clitter prior
	PARPRIORS[23]=OBS$SOM 			; if (OBS$SOM != -9999) {PARPRIORUNC[23]=2.0} # Csom priors
    } else if (modelname == "DALECN_GSI_DFOL_LABILE_FR") {
        PARPRIORS[11]=0.2432501      ; PARPRIORUNC[11]=0.2014871 # log10 avg foliar N (gN.m-2)
#        PARPRIORS[52]=0.0479264      ; PARPRIORUNC[52]=0.01904211 # NUE**(1/-1.38513851385139) (gC/gN)
        PARPRIORS[52]=0.001868948    ; PARPRIORUNC[52]=0.0005951156 # NUE**(1/-2.999299929993) (gC/gN)
	PARPRIORS[17]=-9999          ; PARPRIORUNC[17]=-9999 # LMA 
	PARPRIORS[19]=OBS$Cfol_initial 		; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[19]=2.0} # Cfoliar prior
	PARPRIORS[20]=OBS$Croots_initial 	; if (OBS$Croots_initial != -9999) {PARPRIORUNC[20]=2.0} # Croots prior
	PARPRIORS[21]=OBS$Cwood_initial 	; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=2.0} # Cwood prior
	PARPRIORS[22]=OBS$Clit_initial 		; if (OBS$Clit_initial != -9999) {PARPRIORUNC[22]=2.0} # Clitter prior
	PARPRIORS[23]=OBS$SOM 			; if (OBS$SOM != -9999) {PARPRIORUNC[23]=2.0} # Csom priors
    } else if (modelname == "DALECN_GSI_DFOL_LABILE_FROOT_FR") {
        PARPRIORS[11]=0.2432501      ; PARPRIORUNC[11]=0.2014871 # log10 avg foliar N (gN.m-2)
#        PARPRIORS[52]=0.0479264      ; PARPRIORUNC[52]=0.01904211 # NUE**(1/-1.38513851385139) (gC/gN)
        PARPRIORS[52]=0.001868948    ; PARPRIORUNC[52]=0.0005951156 # NUE**(1/-2.999299929993) (gC/gN)
        PARPRIORS[17]=-9999          ; PARPRIORUNC[17]=-9999 # LMA 
        PARPRIORS[19]=OBS$Cfol_initial          ; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[19]=2.0} # Cfoliar prior
        PARPRIORS[20]=OBS$Croots_initial        ; if (OBS$Croots_initial != -9999) {PARPRIORUNC[20]=2.0} # Croots prior
        PARPRIORS[21]=OBS$Cwood_initial         ; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=2.0} # Cwood prior
        PARPRIORS[22]=OBS$Clit_initial          ; if (OBS$Clit_initial != -9999) {PARPRIORUNC[22]=2.0} # Clitter prior
        PARPRIORS[23]=OBS$SOM                   ; if (OBS$SOM != -9999) {PARPRIORUNC[23]=2.0} # Csom priors
    } else if (modelname == "DALEC_GSI_DFOL_FROOT_FR") {
        PARPRIORS[11]=20.52048   ; PARPRIORUNC[11]=1.1 # 1.617705 # Ceff 
	PARPRIORS[17]=-9999   ; PARPRIORUNC[17]=-9999 # LMA 
	PARPRIORS[19]=OBS$Cfol_initial 		; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[19]=2.0} # Cfoliar prior
	PARPRIORS[20]=OBS$Croots_initial 	; if (OBS$Croots_initial != -9999) {PARPRIORUNC[20]=2.0} # Croots prior
	PARPRIORS[21]=OBS$Cwood_initial 	; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=2.0} # Cwood prior
	PARPRIORS[22]=OBS$Clit_initial 		; if (OBS$Clit_initial != -9999) {PARPRIORUNC[22]=2.0} # Clitter prior
	PARPRIORS[23]=OBS$SOM 			; if (OBS$SOM != -9999) {PARPRIORUNC[23]=2.0} # Csom prior
    } else if (modelname == "DALEC_GSI_MFOL_FR") {
#        PARPRIORS[11]=20.52048   ; PARPRIORUNC[11]=1.617705 # Ceff
	PARPRIORS[17]=-9999   ; PARPRIORUNC[17]=-9999 # LMA 
#	PARPRIORS[19]=OBS$Cfol_initial 		; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[19]=2.0} # Cfoliar prior
	PARPRIORS[20]=OBS$Croots_initial 	; if (OBS$Croots_initial != -9999) {PARPRIORUNC[20]=2.0} # Croots prior
	PARPRIORS[21]=OBS$Cwood_initial 	; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=2.0} # Cwood prior
	PARPRIORS[22]=OBS$Clit_initial 		; if (OBS$Clit_initial != -9999) {PARPRIORUNC[22]=2.0} # Clitter prior
	PARPRIORS[23]=OBS$SOM 			; if (OBS$SOM != -9999) {PARPRIORUNC[23]=2.0} # Csom prior
    } else if (modelname == "DALEC_GSI_BUCKET") {
        PARPRIORS[11]=0.2764618		; PARPRIORUNC[11]=0.2014871 # log10 avg foliar N (gN.m-2)
	if (parameter_type == "pft_specific" & ctessel_pft == 1) {
	    PARPRIORS[12]=OBS$plant         ; PARPRIORUNC[12]=1.1
	    PARPRIORS[15]=OBS$harvest       ; PARPRIORUNC[15]=1.1      
	    #PARPRIORS[36]=50	            ; PARPRIORUNC[36]=1.4 # Croot to half max depth (gbio.m-2)
	    PARPRIORS[37]=1	            ; PARPRIORUNC[37]=1.4 # Max rooting depth (m)
	} else {
	    #PARPRIORS[39]=150               ; PARPRIORUNC[39]=1.4 # Croot to half max depth (gbio.m-2)
	    PARPRIORS[40]=2	      	    ; PARPRIORUNC[40]=1.4 # Max rooting depth (m)
	    PARPRIORS[21]=OBS$Cwood_initial ; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=OBS$Cwood_initial_unc} # Cwood prior
	}
        PARPRIORS[22]=OBS$Clit_initial          ; if (OBS$Clit_initial != -9999) {PARPRIORUNC[22]=2.0} # Clitter prior
        PARPRIORS[23]=OBS$SOM                   ; if (OBS$SOM != -9999) {PARPRIORUNC[23]=2.0} # Csom prior
    } else if (modelname == "DALECN_GSI_BUCKET") {
        PARPRIORS[11]=0.2764618		; PARPRIORUNC[11]=0.2014871*0.5 # log10 avg foliar N (gN.m-2)
	if (parameter_type == "pft_specific" & ctessel_pft == 1) {
	    PARPRIORS[12]=OBS$plant         ; PARPRIORUNC[12]=1.1
	    PARPRIORS[15]=OBS$harvest       ; PARPRIORUNC[15]=1.1      
	    #PARPRIORS[36]=50	            ; PARPRIORUNC[36]=1.4 # Croot to half max depth (gbio.m-2)
	    #PARPRIORS[37]=1	            ; PARPRIORUNC[37]=1.4 # Max rooting depth (m)
	} else {
	    PARPRIORS[2]=51.70631	    ; PARPRIORUNC[2]=0.53905 # C:N root (gC/gN) Kattge et al., (2011)
	    PARPRIORS[27]=416.6667	    ; PARPRIORUNC[27]=6.277688 # C:N wood (gC/gN) Kattge et al., (2011)
	    PARPRIORS[21]=OBS$Cwood_initial ; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=OBS$Cwood_initial_unc} # Cwood prior
	    #PARPRIORS[39]=150               ; PARPRIORUNC[39]=1.4 # Croot to half max depth (gbio.m-2)
	    #PARPRIORS[40]=2	      	    ; PARPRIORUNC[40]=1.4 # Max rooting depth (m)
	    PARPRIORS[41]=1.639		    ; PARPRIORUNC[41]=0.125 # Rm_leaf N**exponent (gC/gN) Reich et al., (2008)
	    PARPRIORS[43]=1.352		    ; PARPRIORUNC[43]=0.150 # Rm_root N**exponent (gC/gN) Reich et al., (2008)
	    PARPRIORS[45]=1.344		    ; PARPRIORUNC[45]=0.150 # Rm_wood N**exponent (gC/gN) Reich et al., (2008)
	}
        PARPRIORS[22]=OBS$Clit_initial          ; if (OBS$Clit_initial != -9999) {PARPRIORUNC[22]=2.0} # Clitter prior
        PARPRIORS[23]=OBS$SOM                   ; if (OBS$SOM != -9999) {PARPRIORUNC[23]=2.0} # Csom prior
    } else if (modelname == "ACM") {
	
    # For ACM_GPP_ET
    # p(1) = nitrogen use efficiency at optimum temperature (oC)
           #,unlimited by CO2, light and photoperiod (34gC/gN)
#    PARPRIORS[1] = 35.65 ; PARPRIORUNC[1] = 26.19
    # p(2) = maximum temperature at which photosynthesis occurs (oC) (59.04677oC)
#    PARPRIORS[2] = 59.04677 ; PARPRIORUNC[2] = 5.0
    # p(3) = optimum temperature for photosynthesis (oC) (30oC)
#    PARPRIORS[3] = 30.0 ; PARPRIORUNC[3] = 5.0
    # p(4) = kurtosis for temperature response of photosynthesis (0.185912)
#    PARPRIORS[4] = 0.185912 ; PARPRIORUNC[4] = 0.05
    # p(5) = maximum canopy quantum yield (gC/MJ)
    # p(6) = constant on daylength impact
    # p(7) = coefficient on daylength impact
    # p(8) = maximum absorbed radiation (fraction)
    # p(9) = LAI at which radiation absorption is at half saturation (m2/m2)
    # p(10) = Maximum (most negative) leaf-soil WP difference (MPa) (i.e. minLWP)
#    PARPRIORS[10] = -2.0 ; PARPRIORUNC[10] = 0.1

    # p(11) = temperature (oC) at which gc not limited
    # p(12) = sensitivity of gc to temperature
    # p(13) = absorbed shortwave radiation (W.m-2) at which gc is 50 % limited by light
    # p(14) = sensitivity of gc to VPD (kPa)
    # p(15) = soil sw radiation absorption (fraction)

	# for GPP temperature optimum = 30oC (Default) 
        # CO2 compensation point and half saturation
#	PARPRIORS[3] = 45.29614 ; PARPRIORUNC[3] = 20.54073
#	PARPRIORS[4] = 319.4397 ; PARPRIORUNC[4] = 60.01601
	# for GPP temperature optimum = 13oC (Arctic)
        # CO2 compensation point and half saturation
#	PARPRIORS[3] = 4.78419 ; PARPRIORUNC[3] = 4.806224
#	PARPRIORS[4] = 347.4542 ; PARPRIORUNC[4] = 107.3912
    }

    # combine the static data
    DATA_STAT=c(PARPRIORS,PARPRIORUNC,OTHERPRIORS,OTHERPRIORUNC)

    # open the binary file
    zz <- file(file, "wb")
    # write with 8 bit precision (i.e. double)
    writeBin(as.double(static_data), zz)
    writeBin(as.double(DATA_STAT), zz)
    writeBin(as.double(as.vector(DATA_TEMP)), zz)
    close(zz)

}
## Use byte compile
binary_data<-cmpfun(binary_data)
