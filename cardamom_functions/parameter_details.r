###
## Function which contains details of parameter names and units
###

parameter_details<- function (modelname,parameter_type,ctessel_pft) {
    if (modelname == "DALEC_CDEA") {
	output=list(npars=23)
	output$parameter_symbols=c("p1","p2","p3","p4","p5","p6","p7","p8","p9","p10","p11","p12","p13","p14","p15","p16","p17","i1","i2","i3","i4","i5","i6")
	output$parameter_names=c("D_{rate}","F_{Rau}","F_{Fol}","F_{Roots}","L","TOR_{Wood}","TOR_{Roots}","MR_{litter}"
				,"MR_{SOM}","Temp_{par}","C_{eff}","B_{day}","F_{Lab}","R_{L}","F_{day}","R_{F}","LMA"
				,"C_{LABILE}","C_{FOLIAR}","C_{ROOT}","C_{WOOD}","C_{LITTER}","C_{SOM}")
	output$parameter_units=c("day^{-1}","fraction","fraction","fraction","years","day^{-1}","day^{-1}","day^{-1}","day^{-1}"
				,"(unitless)","(unitless)","doy","fraction","days","doy","days","gC m^{-2}","gC m^{-2}","gC m^{-2}"
				,"gC m^{-2}","gC m^{-2}","gC m^{-2}","gC m^{-2}")
    } else if (modelname == "AT_DALEC" & parameter_type == "pft_specific" & ctessel_pft == 1) {
	output=list(npars=33)
	output$parameter_symbols=c("p1","p2","p3","p4","p5","p6","p7","p8","p9","p10","p11","p12","p13","p14","p15","p16","p17","i1","i2","i3","i4","i5","i6","p26","p27","p28","p29","p30","p31","p32","p33")
	output$parameter_names=c("D_rate","F_Rau","DR_pre","DR_post","TOR_FOL","TOR_Wood","RDRSHMAX","VDh"
				,"MR_SOM","MR_LIT","TOR_auto","Sow_Day","RC_lab_trans","PHUem","Harvest_day","Plough_day"
				,"LMA","C_LABILE","C_FOLIAR","C_ROOT","C_WOOD","C_LITTER","C_SOM","C_auto","C_storage_organ"
				,"tmin","tmax","topt","tmin_v","tmax_v","topt_v","PHCR","PHSC")
	output$parameter_units=c("hour^{-1}","fraction","coef","coef","hour^{-1}","hour^{-1}","hour^{-1}","days","hour^{-1}"
				,"hour^{-1}","doy","hour^{-1}","oC/days","doy","doy","gC m^{-2}","gC m^{-2}","gC m^{-2}"
				,"gC m^{-2}","gC m^{-2}","gC m^{-2}","gC m^{-2}","gC m^{-2}","gC m^{-2}"
				,"oC","oC","oC","oC","oC","oC","hours","unitless")
    } else if (modelname == "AT_DALEC") {
	output=list(npars=22)
	output$parameter_symbols=c("p1","p2","p3","p4","p5","p6","p7","p8","p9","p10","p11","p12","p13","p14","p15","p16","i1","i2","i3","i4","i5","i6")
	output$parameter_names=c("D_{rate}","F_{Rau}","F_{Fol}","F_{Roots}","L","TOR_{Wood}","TOR_{Roots}","MR_{litter}"
				,"MR_{SOM}","Temp_{par}","B_{day}","F_{Lab}","R_{L}","F_{day}","R_{F}","LMA"
				,"C_{LABILE}","C_{FOLIAR}","C_{ROOT}","C_{WOOD}","C_{LITTER}","C_{SOM}")
	output$parameter_units=c("day^{-1}","fraction","fraction","fraction","years","day^{-1}","day^{-1}","day^{-1}","day^{-1}"
				,"(unitless)","doy","fraction","days","doy","days","gC m^{-2}","gC m^{-2}","gC m^{-2}"
				,"gC m^{-2}","gC m^{-2}","gC m^{-2}","gC m^{-2}")
    } else {
	print("model name not found in 'parameter_details'")
    }
    return(output)
}

