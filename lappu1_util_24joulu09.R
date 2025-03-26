## pakettien lataus ------------------------------------------------------
libs <- unique(c("data.table","ggplot2","reshape2","fusedlasso","glue","pROC","plyr","dplyr","MASS","KScorrect","glmnet","MASS","DT","gtExtras","gt","gtsummary","stringr","tidyr", "purrr", "Hmisc","CalibrationCurves", "patchwork", "knitr","kableExtra", "ggsurvfit","tibble","rms", "ggtext", "magick"))
lapply(libs, require, character.only=TRUE)




## apufunktioita ------------------------------------------------------

#laskee logitin betoista todennäköisyydet.
beta_to_probs <- function(x, beta){
  prob.logit <- cbind(1,x) %*% beta
  prob <- exp(prob.logit) / (1 + exp(prob.logit))
  return(prob)
}


expit <- function(x){
  return(exp(x) / (1+ exp(x)))
}

logit <- function(x){
  return(log(x/(1-x)))
}

#vuodet kuten 02 ja 07 muuntuu muotoon 2 ja 7, jos ne annetaan numeerisina. tämä funktio laittaa nollat takaisin mukaan
d2 <- function(x){
  ifelse(nchar(x)==1, paste0("0",x), x)
}


## parametrien listoja ----
get_tab1_vars <- function(){
  clinical_vars <- c('(Intercept)','IKA','SUKUP2','ALKI1','TUPAKOINTI', 'TUPAKOINTI1', 'HDL','KOL','TRIG','SYSm','DIASm','VYOTARO','BMI', 'GLC', 'BP_TREAT1', 'LIPID_TREAT1', 'DIAB_FAMILYHIST1')
  return(clinical_vars)
}


get_fresh_clinical_vars_may24 <- function(){
  clinical_vars <- c('(Intercept)','IKA','SUKUP2', 'ALKI1', 'TUPAKOINTI', 'TUPAKOINTI1', 'HDL','KOL','TRIG','SYSm','DIASm','VYOTARO','BMI', 'GLC','BP_TREAT','LIPID_TREAT','DIAB_FAMILYHIST','BP_TREAT1','LIPID_TREAT1','DIAB_FAMILYHIST1')
  return(clinical_vars)
}


## raakadatan kokoaminen ------------------------------------------------------

get_finriskdata <- function(vuosi){
  if(vuosi %in% c("02", "07", "12", "92", "97", "00", "11")) {
    fp <- file.path("F:/profi5/data/isofinrisk", grep(vuosi, list.files("F:/profi5/data/isofinrisk"),value=TRUE))
  }
  if(vuosi == 2) fp <- "F:/profi5/data/isofinrisk/FINRISK_02.txt"
  if(vuosi == 7) fp <- "F:/profi5/data/isofinrisk/FINRISK_07.txt"
  if(vuosi == 12) fp <- "F:/profi5/data/isofinrisk/FINRISK_12.txt"
  if(vuosi == 92) fp <- "F:/profi5/data/isofinrisk/FINRISK_92.txt"
  if(vuosi == 97) fp <- "F:/profi5/data/isofinrisk/FINRISK_97.txt"
  if(vuosi == 11) fp <- "F:/profi5/data/isofinrisk/Health_2011.txt"
  if(vuosi == 0) fp <- "F:/profi5/data/isofinrisk/Health_2000.txt"
  if(vuosi %in% c("E", "e")) fp <- "F:/profi5/data/isofinrisk/other/FR_T2000_T2011_KUU_first_events_forOYmetaboanalysis_2021_05_19.txt"
  
  if(!vuosi %in% c("02", "07", "12", "92", "97", "00", "11","E","e", 02, 07, 12, 92, 97, 00, 11)){
    cat('Valitse joku vuosista c("02", "07", "12", "92", "97", "00", "11", tai "E")')
  }
  
  if(vuosi %in% c("02", "07", "12", "92", "97", "00", "11", "E", "e", 02, 07, 12, 92, 97, 00, 11)){
    dat <- data.frame(fread(fp))
    if("SP2" %in% names(dat)) setnames(dat, "SP2", "SUKUP")
    return(dat)
  }
}

#Hakee ja asettaa yksittäisen finrisk-datan/eventit
load_finrisk <- function(vuosi, returndf=F){
  if(returndf == F) assign(paste0("fr",vuosi), get_finriskdata(vuosi), envir=.GlobalEnv)
  if(returndf == T) return(get_finriskdata(vuosi))
}

#Asettaa kaikki finrisk-datat muistiin nimillä frXX, sekä events
load_all_finrisk <- function(){
  fr92 <- get_finriskdata(92)
  fr97 <- get_finriskdata(97)
  fr02 <- get_finriskdata(02)
  fr07 <- get_finriskdata(07)
  fr12 <- get_finriskdata(12)
  fr00 <- get_finriskdata(00)
  fr11 <- get_finriskdata(11)
  events <- get_finriskdata("e")
  assign("fr92", fr92,envir=.GlobalEnv)
  assign("fr97", fr97,envir=.GlobalEnv)
  assign("fr02", fr02,envir=.GlobalEnv)
  assign("fr07", fr07,envir=.GlobalEnv)
  assign("fr12", fr12,envir=.GlobalEnv)
  assign("fr00", fr00,envir=.GlobalEnv)
  assign("fr11", fr11,envir=.GlobalEnv)
  assign("events", events, envir=.GlobalEnv)
}

#Hakee eventdatan ja asettaa sen muistiin nimellä events
load_eventdata <- function(){
  assign("events", get_finriskdata("e"), envir=.GlobalEnv)
}


#Luo eventdataan sarakkeen followup-rajatulle vasteelle. Eli esim I9_CHD kymmenen vuoden sisään baselinestä. Näin eri keräykset vertailukelpoisia.
create_followup_event_variable <- function(events=events, response_var, fu_length=10, fu_min=0, new_var_name){
  new_var_name <- glue("{response_var}_FU{fu_length}")
  agevar <- glue("{response_var}_AGE")
  
  
  # eventage_over_0 <- events[[agevar]] - events$BL_AGE > 0 #onko event tapahtunut jo ennen baselinea?
  eventage_over_0 <- events[[agevar]] - events$BL_AGE > fu_min #onko event tapahtunut jo ennen baselinea?
  events <- events[eventage_over_0,] #Jos on, niin poistetaan henkilö datasta
  
  eventage_under_fu_max <- (events[[agevar]] - events$BL_AGE) < fu_length #Onko event tapahtunut followupin sisällä? 
  # eventage_over_fu_min <- (events[[agevar]] - events$BL_AGE) > fu_min #Onko event tapahtunut fu_min jälkeen?
  
  # events$response_fuOK <- ifelse(events[[response_var]]==1 & eventage_under_fu_max & eventage_over_fu_min, 1, 0) #jos on, hyväksytään event. muuten ei.
  events$response_fuOK <- ifelse(events[[response_var]]==1 & eventage_under_fu_max, 1, 0) #jos on, hyväksytään event. muuten ei.
  setnames(events, "response_fuOK", new_var_name)
  return(events)
}



## datan kokoamisen apufunktioita ------------------------------------------------------
make_dummy_binary_vars <- function(x){
  x %>% model.matrix(~.-1, .) %>% .[,colnames(.) != "SUKUP1"]
}

load_overlap_names <- function() {
  fread("F:/profi5/data/util/overlap_names_97_12.txt", header = F)$V1 #GLOL poistetty tästä joukosta 1.12.2023
}

get_newvars <- function(overlap_names) {
  newvars <- fread("F:/profi5/data/isofinrisk/FR97_FR12_poiminta0524/FR97_FR12_BP_LIPID_TREAT_DIAB_MI_FAMHIST_hba1c_UID_0524.txt") %>% select(c(UID,any_of(overlap_names)))
}

prep_covariates <- function(fr_set_list) {
  overlap_names <- load_overlap_names()
  
  x_raw <- data.frame()
  for(i in 1:length(fr_set_list)){
    frset <- fr_set_list[[i]]
    cleanset <- frset[,which(names(frset) %in% intersect(names(frset),overlap_names))]
    x_raw <- rbind(x_raw, cleanset)
  }
  
  #uudet muuttujat
  newvars <- get_newvars(overlap_names)
  x_raw <- left_join(x_raw, newvars, by="UID") #%>% rename(ACTV=Q57) %>%rename(ACTT=Q56) #aktiivisuusmuuttujia ei käytetä, mutte ne uudelleennimettäisiin tässä.
  return(x_raw)
}


#siistitään eventabia, otetaan vain rivit joilla studyID on numeerinen. muut oletettavasti roskaa.
filter_events <- function(events_tab, response_var, extra_vars = NULL) {
  vars <- c("studyID", "cohort", response_var, extra_vars)
  events_tab %>%
    filter(!grepl("[A-Za-z]", studyID)) %>%
    select(all_of(vars)) %>%
    mutate(studyID = as.numeric(studyID))
}


clean_dataset <- function(dataset, dropID, extra_drops = NULL) {
  drops <- c("UID", "cohort", "SAMPLE_COLLECTION", "VUOSI", "ALUE", "IKA10", "K66", "K34", grep("Q|KY", names(dataset), value = TRUE), extra_drops)
  dataset <- dataset %>%
    select(-all_of(drops)) %>%
    suppressWarnings(mutate(across(everything(), as.numeric)))
  if(dropID) dataset <- dataset %>% select(-studyID)
  dataset
}


#apufunktio. korjaa FR97:n väärät yksiköt
fix_units <- function(X, yr) {
  if(yr == 97){
    wrong_unit_vars <- c('XXL.VLDL.P','XL.VLDL.P','L.VLDL.P','M.VLDL.P','S.VLDL.P','XS.VLDL.P','IDL.P','L.LDL.P',
                         'M.LDL.P','S.LDL.P','XL.HDL.P','L.HDL.P','M.HDL.P','S.HDL.P')
    X[,wrong_unit_vars] <- X[,wrong_unit_vars] * 10^-6
  }
  return(X)
}

#ottaa verenpainemittausten keskiarvot
take_bp_mean <- function(X) {
  X[,"SYS1"] <- rowMeans(X[,colnames(X) %in% c("SYS1", "SYS2")])
  X[,"DIAS1"] <- rowMeans(X[,colnames(X) %in% c("DIAS1", "DIAS2")])
  X <- X[,!(colnames(X) %in% c("SYS2", "DIAS2"))]
  colnames(X)[colnames(X) == c("SYS1")] <- "SYSm"
  colnames(X)[colnames(X) == c("DIAS1")] <- "DIASm"
  return(X)
}

#korjaa tupakointimuuttujan kaksitasoiseksi
fix_tupakointi <- function(X) {
  TUPAKOINTI <- ifelse(X[,"TUPI3"]==4, 1, 0) 
  #DD: 1=ei tup. 2=lop. yli 1/2-v s. 3=lop. alle 1/2-v s. 4=tupakoi
  #tässä 1=tupakoi, 0=ei tupakoi
  X[,"TUPI3"] <- TUPAKOINTI
  colnames(X)[colnames(X) == "TUPI3"] <- "TUPAKOINTI"
  X <- X[,grep("TUPI", colnames(X), invert=T)]
  return(X)
}


get_variable_order <- function(response_var_name) {
  x_order <- fread("F:/profi5/data/util/elo24_diab_data_muuttuja_jarjestys.txt",header=F)$V1
  var_order <- c(x_order, response_var_name, "studyID")
}


get_scalevars <- function(){
  return(c("SUKUP2","TUPAKOINTI1","Y","BP_TREAT1" , "LIPID_TREAT1", "DIAB_FAMILYHIST1", "LOW_ACTIVITY1"))
}


censor_glucose_atlim <- function(dat, execute=FALSE, fu_min){
  ### Vanha tapa:: glclim = ifelse(cohort=="FR97", 6.633116, 5.6). cohort tuli jostain muistista.
  # paastoglukoosin raja on 5.6, käytetään nyt aluksi sitä, saadaan tosi varovainen raja, katsotaan kuinka sen kanssa käy.
  # > varovainen_paasto_limit <- bind_rows(
  #   +   prop.table(table(fr02$GLC<5.6)) %>% as.data.frame,
  #   +   prop.table(table(fr07$GLC<5.6)) %>% as.data.frame,
  #   +   prop.table(table(fr12$GLC<5.6)) %>% as.data.frame
  #   + ) %>% 
  #   +   filter(Var1 ==TRUE) %>% 
  #   +   pull(Freq) %>%
  #   +   min
  # > quantile(fr97$GLC, varovainen_paasto_limit,na.rm=T)
  # 97.61674% 
  # 6.633116 
  ### ---
  
  ### Uusi tapa: 
  # >   datasetup <- setup_xy_DIAB(followup_length = 1, yr=2, fu_min=0, censor_glucose=F)$fullset %>% 
  #   +     make_dummy_binary_vars(.) %>% 
  #   +     data.frame
  # > datasetup$GLC%>%sd
  # [1] 0.5097442
  # > datasetup$GLC%>%mean
  # [1] 4.049101
  # > (5.6 - 4.049101) / 0.5097442
  # [1] 3.042504
  
    # >   datasetup <- setup_xy_DIAB(followup_length = fu_len, yr=2, fu_min=1, censor_glucose=F)$fullset %>% 
  #   +     make_dummy_binary_vars(.) %>% 
  #   +     data.frame
  # > datasetup$GLC %>% sd
  # [1] 0.4757874
  # > datasetup$GLC %>% mean
  # [1] 4.038504
  # > (5.6 - 4.038504) / 0.4757874
  # [1] 3.28192
  
  #vanha (6.633116, fr97, vs uusi tapa, plus 5.6 raja vertailun vuoksi. FR02 datalla ablinet on samassa kohtaa jakaumaa, muilla seteillä (erityisesti nyt siis fr97) hiukan sivussa. 
  # par(mfrow=c(2,1))
  # hist(dat$GLC, breaks=100)
  # mean(dat$GLC,na.rm=T)
  # abline(v=c(5.6))
  # hist(dat$scaled_glc, breaks=100)
  # mean(dat$scaled_glc,na.rm=T)
  # abline(v=glc_lim)
  
  if(execute){
    
    if(fu_min==0) {glc_lim <- 3.042504; glc_sd <- 0.5097442} # ei minimiä follow-upille
    if(fu_min==1) {glc_lim <- 3.28192; glc_sd <- 0.4757874} # kyllä minimi follow-upille

    dat <- dat %>% 
      mutate(scaled_glc = (GLC - mean(dat$GLC, na.rm=T)) / glc_sd) %>%
      filter(scaled_glc < glc_lim) %>%
      select(-scaled_glc)
  }
  return(dat)
}

## käyttödatan kokoaminen ------------------------------------------------------


make_default_dataset <- function(fr_set_list, events_tab, response_var, censor_glucose=T, event_ages=F, dropID=T, fu_min){
  
  #siistitään selittäjiä
  x_raw <- prep_covariates(fr_set_list)
  
  if(!event_ages) { #tavallinen tilanne
    y_raw <- filter_events(events_tab, response_var) #vastemuuttuja
    prepped_dat <- inner_join(x_raw, y_raw, by = c("studyID", "cohort")) %>%
      censor_glucose_atlim(execute=censor_glucose, fu_min) %>% 
      clean_dataset(., dropID=dropID) #valmis dataset
    
  } else { #erikoistilanne, kun halutaan time-to-event -tietoa: event_ages == TRUE
    response_agevar <- glue("{response_var}_AGE") %>% gsub("_FU[0-9]+", "", .)
    response_yearvar <- glue("{response_var}_YEAR") %>% gsub("_FU[0-9]+", "", .)
    y_raw <- filter_events(events_tab, response_var, extra_vars = c(response_agevar, response_yearvar, "BL_AGE", "BL_YEAR", "FU_END_AGE"))
    prepped_dat <- inner_join(x_raw, y_raw, by = c("studyID", "cohort")) %>%
      censor_glucose_atlim(execute=censor_glucose,fu_min) %>%
      clean_dataset(., dropID=FALSE, extra_drops = "cohort")
    
    default_pick <- make_default_dataset(fr_set_list, events_tab, response_var, event_ages = FALSE,  censor_glucose=censor_glucose, dropID = FALSE, fu_min) #ajetaan funktio uusiksi, mutta nyt event_ages == FALSE. näin saadaan tietoon oikeat studyID:t
    prepped_dat <- prepped_dat[prepped_dat$studyID %in% default_pick$studyID, ] #filtteröidään 
    
  }
  
  
  return(prepped_dat)
}

#kokoaa "lopullisen" datasetin, eli sellaisen jossa oikea määrä muuttujia oikeassa muodossa. skaalaus tehdään muualla.
setup_xy_DIAB <- function(yr, followup_length, fu_min=0, censor_glucose=T, dropID=T){
  vaste <- "E4_DIABETES"
  fu_suffix <- glue("_FU{followup_length}")
  response_var_name <-  paste0(vaste, fu_suffix)
  events <- get_finriskdata("e")
  frSet <- get_finriskdata(yr)
  
  #kokoaa event-dataan uuden muuttujan vaste-followup -parille. fu_min katsotaan täällä.
  eventsFUx <- create_followup_event_variable(events, 
                                              response_var=vaste, 
                                              fu_length=followup_length, 
                                              fu_min=fu_min) 
  #kokoaa selittäjät ja vasteen datasetiksi
  fulldat <- make_default_dataset(fr_set_list=list(frSet), 
                                  events_tab = eventsFUx, 
                                  response_var = response_var_name, 
                                  censor_glucose=censor_glucose,
                                  event_ages = FALSE, 
                                  dropID=dropID,
                                  fu_min=fu_min) 
  
  #järjestetään datasetin sarakkeita
  var_order <- get_variable_order(response_var_name)
  fulldat <- fulldat %>%
    select(any_of(intersect(var_order, names(.)))) %>%
    apply(.,2,as.numeric) %>% #HOX! Tämä on toteutettava nimenomaan näin. mutate(), tai supressWarnings() molemmat rikkoo toimminnon. Idea on että ei-numeeriset vedetään NA:ksi ja myöhemmin otetaan complete.cases(). 
    data.frame 
  names(fulldat) <- gsub(response_var_name, "Y", names(fulldat))
  
  #siistitään ja korjataan datasetin muuttujia
  fulldat <- fulldat %>%
    fix_units(yr) %>%
    take_bp_mean %>%
    fix_tupakointi %>%
    mutate(across(any_of(c("SUKUP", "TUPAKOINTI", "BP_TREAT", "LIPID_TREAT", "DIAB_FAMILYHIST", "LOW_ACTIVITY")), as.factor)) %>%
    filter(complete.cases(.)) #complete caset vasta nyt lopuksi, kun muuttujat on rajattu.
  
  
  #dropID==FALSE -> ajetaan yo. skripti ensin "normaalisti", ja sitten uudelleen niin että event_ages==FALSE. 
  if(!dropID){
    included_ids <- fulldat$studyID
    fulldat <- fulldat %>% select(-studyID)
    eventsFU50 <- create_followup_event_variable(events, response_var="E4_DIABETES", fu_length=50, fu_min=fu_min) # niin pitkä followup että saadaan kaikkien time-to-event -tiedot
    fulldat.t2e <- make_default_dataset(fr_set_list=list(frSet), 
                                        events_tab = eventsFU50, 
                                        response_var = "E4_DIABETES_FU50", 
                                        censor_glucose = censor_glucose,
                                        event_ages = TRUE, 
                                        dropID=dropID,
                                        fu_min) %>%
      filter(studyID %in% included_ids) #filtteröidään "normaaliin" joukkoon rivejä
    fulldat <- fulldat.t2e
  }
  
  return(list(fullset=fulldat))
}


scale_to_SD_units <- function(df, vars_to_scale, train_means=NULL, train_stds=NULL, metab_to_mean0=NULL){
  
  vars_not_to_scale <- setdiff(1:ncol(df), vars_to_scale) #muuttujaindeksit, joita ei skaalata
  df_to_scale <- df[,vars_to_scale] #skaalattava osa datasta
  
  #opetusdatalle: funktio määrittää centerin ja stdn. palauttaa ne listana.
  if(is.null(train_means) | is.null(train_stds)){
    df_scaled <- scale(df_to_scale)
    df_out <- df 
    df_out[,vars_to_scale] <- df_scaled #input-datan skaalattavien sarakkeiden paikoille laitetaan skaalattu data. muut sarakkeen pysyy alkup. kunnossa
    df_out <- data.frame(df_out)
    return(list(data=df_out, center=attr(df_scaled, "scaled:center"), std=attr(df_scaled, "scaled:scale")))
  }
  
  #testidatoille: annetaan opetusdatan center ja std. palautetaan ne yhdenmukaisuuden nimissä.
  if(!(is.null(train_means) | is.null(train_stds))){
    all_vars_to_scale <- names(df_to_scale)
    
    
    #skaalataanko metabolomiikka erikseen? jos kyllä, niin center=T näille ->
    #skaalaus mean=0. kliiniset menee sitten erikseen opetusdatan mukaan. samoin
    #keskihajonnat menee opetusdatan perusteella.
    if(!is.null(metab_to_mean0)){
      metabol_to_scale <- metab_to_mean0
      metabol_inds <- which(all_vars_to_scale %in% metabol_to_scale)
      clin_inds <- which(!all_vars_to_scale %in% metabol_to_scale)
      
      df_scaled_metabol <- scale(df_to_scale[,metabol_inds], 
                                 center=T, 
                                 scale=train_stds[metabol_inds])
      
      df_scaled_clin <- scale(df_to_scale[,clin_inds], 
                              center=train_means[clin_inds], 
                              scale=train_stds[clin_inds])
      
      df_scaled <- df_to_scale
      df_scaled[,metabol_inds] <- df_scaled_metabol
      df_scaled[,-metabol_inds] <- df_scaled_clin
      
      df_out <- df 
      df_out[,vars_to_scale] <- df_scaled 
      df_out <- data.frame(df_out)
    }else{
      df_scaled <- scale(df_to_scale, center=train_means, scale=train_stds)
      df_out <- df 
      #input-datan skaalattavien sarakkeiden paikoille laitetaan skaalattu data. muut sarakkeen pysyy alkup. kunnossa
      df_out[,vars_to_scale] <- df_scaled 
      df_out <- data.frame(df_out)
    }


    
    return(list(data=df_out, center=train_means, std=train_stds))
    
  }
}




build_followupX_datasets <- function(fu_len, fu_min=0, censor_dataset_glucose=T, doNewScale=T, clinical_vars){
  
  #opetusdata ennen skaalausta
  datasetup <- setup_xy_DIAB(followup_length = fu_len, yr=2, fu_min=fu_min, censor_glucose=censor_dataset_glucose)$fullset %>% 
    make_dummy_binary_vars(.) %>% 
    data.frame
  vars_to_scale <- which(names(datasetup) %in% setdiff(names(datasetup), get_scalevars()))
  
  #opetusdatan skaalaus
  datasetup_scale <- scale_to_SD_units(datasetup, vars_to_scale)
  train_center <- datasetup_scale$center
  train_std <- datasetup_scale$std
  traindat <- datasetup_scale$data
  
  #skaalattu opetusdata
  xtrain <- traindat[,-ncol(traindat)]
  ytrain <- traindat[,ncol(traindat)]
  
  #testisetit 97,07,12
  testset_years <- c(97,07,12)
  
  #normaalitapaus, jossa testisetit skaalataan opetusdatan mukaan
  if(doNewScale){
    cat("\nSkaalataan dataa...\n")
    testset_list <- lapply(testset_years, function(yr){
      unscaled_testset <- setup_xy_DIAB(followup_length = fu_len, yr=yr, fu_min=fu_min, censor_glucose=censor_dataset_glucose)$fullset %>% 
        make_dummy_binary_vars(.) %>% 
        data.frame
      
      #mitkä muuttujanimet on a) skaalattavia ja b) ei-kliinisiä => metabolomiikka
      metab_to_mean0_vars <- setdiff(names(unscaled_testset[,vars_to_scale]), setdiff(clinical_vars,"GLC"))

      scaled_testset <- unscaled_testset %>%
        scale_to_SD_units(., vars_to_scale, train_center, train_std, metab_to_mean0=metab_to_mean0_vars) %>%
        pluck(.,"data")
      
      list(dat=scaled_testset, metab_vars=metab_to_mean0_vars)
      
    })
    
    testset_dats <- map(testset_list, "dat")
    names(testset_dats) <- paste0("test",d2(testset_years))
    
    
    return(
      list(
        traindat = traindat,
        xtrain = xtrain,
        ytrain = ytrain,
        test97 = testset_dats$test97,
        test07 = testset_dats$test07,
        test12 = testset_dats$test12,
        vars_to_scale = vars_to_scale,
        scale_train_center = train_center,
        scale_train_std = train_std,
        metab_to_mean0_vars = testset_list[[1]]$metab_vars
      )
    )
  }
  
  #erikoistapaus, datat kootaan ilman skaalausta resultpagerin kliinisten muuttujien taulua varten. 
  if(!doNewScale){
    
    testset_list <- lapply(c(97,07,12), function(yr){
      setup_xy_DIAB(followup_length = fu_len, yr=yr, fu_min=fu_min, censor_glucose=censor_dataset_glucose)$fullset %>% 
        make_dummy_binary_vars(.) %>% 
        data.frame
    })
    names(testset_list) <- paste0("test",d2(testset_years))
    
    unscaled_datasetup <- datasetup
    unscaled_xtrain <- unscaled_datasetup[,-ncol(unscaled_datasetup)]
    unscaled_ytrain <- unscaled_datasetup[,ncol(unscaled_datasetup)]
    return(
      list(
        traindat = unscaled_datasetup,
        xtrain = unscaled_xtrain,
        ytrain = unscaled_ytrain,
        test97 = testset_list$test97,
        test07 = testset_list$test07,
        test12 = testset_list$test12,
        metab_to_mean0_vars = NULL
        )
    )
  }
}


## mallinteko apufunktioita ------------------------------------------------------

get_bonferroni_correction_for_stepwise_mods <- function(fu_min, censor_dataset_glucose){
  datasetup <- setup_xy_DIAB(followup_length = 5, yr=2, fu_min=fu_min, censor_glucose = censor_dataset_glucose)
  traindat <- datasetup$fullset %>%  
    make_dummy_binary_vars(.) %>% 
    data.frame

  pc <- prcomp(apply(traindat[,-ncol(traindat)],2,as.numeric),center = TRUE,scale. = TRUE)
  pval_divide_by_X <- as.numeric(which(cumsum(summary(pc)$importance[3,]>0.99)==1)) #false+false+...+false+true -> cumsum=0,0,...,0,1,2,3,... -> missä cumsum==1, siellä ensimmäinen true, eli cumprop>0.99 ekaa kertaa.
  return(pval_divide_by_X) 
}



get_penaltyfac_ja_kiinteet <- function(xtrain,clinical_vars, model_vars){
  pfac <- rep(1,ncol(xtrain))
  if(model_vars=="fix_none"){
    fixed_vars  <- NULL
  }  
  if(model_vars=="fix_clinical"){
    pfac[which(colnames(xtrain) %in% clinical_vars)]<-0
    fixed_vars  <-  clinical_vars
  } 
  return(list(penalty_factor=pfac, fixed_clinicals=fixed_vars))
}



fix_factor_names <- function(x){
  x<-gsub("SUKUP2", "SUKUP",x)
  x<-gsub("TUPAKOINTI1", "TUPAKOINTI",x)
  x<-gsub("DIAB_BL1", "DIAB_BL",x)
  return(x)
}

refactor_plain_names <- function(x){
  x<-gsub("SUKUP", "SUKUP2", x)
  x<-gsub("TUPAKOINTI", "TUPAKOINTI1",x)
  x<-gsub("DIAB_BL", "DIAB_BL1",x)
  return(x)
}


get_clinical_only_set <- function(datalist,clinical_vars, clin_only_list_names=c("train02","test97","test07","test12")){
  
  dataset_names <- c("traindat", "test97", "test07", "test12")
  
  clol <- lapply(datalist[dataset_names], function(df){select(df, all_of(intersect(c(clinical_vars, "Y"), colnames(df))))})
  names(clol) <- clin_only_list_names
  return(clol)
}


## MALLIT -------------------

build_logreg_model <- function(traindat){
  glm(Y~., family="binomial",  data=data.frame(traindat))
}


# Elanet mallin teko
build_elanet_model <- function(xtrain, ytrain, alphaval, pfac, ...){
  set.seed(2024)
  cv.glmnet(as.matrix(xtrain), as.matrix(ytrain), family="binomial", alpha = alphaval, penalty.factor=pfac, trace.it=F, ...)
}




# Create a custom stepwise logistic regression function based on p-values
stepwise_logistic_regression <- function(significance_level = 0.05, fixed_vars=NULL, traindat=traindat,do_forward=T) {
  set.seed(2024)
  full_model <- glm(Y~., family="binomial",  data=data.frame(traindat))
  latest_model <- full_model
  
  #backward stepwise
  cat("\nPudotetaan muuttujat ")
  while (TRUE) {
    summary_model_coefs <- summary(latest_model)$coefficients
    summary_model_coefs <- summary_model_coefs[rownames(summary_model_coefs) != "(Intercept)", ,drop=F]
    
    if(!is.null(fixed_vars)){summary_model_coefs <- summary_model_coefs[!rownames(summary_model_coefs) %in% fixed_vars,, drop=F]} 
    
    z_values <- summary_model_coefs[, "z value", drop=F]
    p_values <- summary_model_coefs[, "Pr(>|z|)", drop=F]
    estimates <- summary_model_coefs[,"Estimate", drop=F]
    bigBtag <- ifelse(nrow(estimates)==0, FALSE, max(abs(estimates)) > 10^10) #onko malli "valmis", eli onko mukana isoja betoja, tai onko viimeinen kierros fixed-mallia?
    
    # suurin p-arvo == pienin abs(z-arvo)
    # poistettava muuttuja voidaan siis yhtä hyvin katsoa z-arvon perusteella kuin p-arvonkin.
    variable_to_remove <- rownames(z_values)[which.min(abs(z_values))] 
    
    #p-arvon perusteella katsotaan poistetaanko MITÄÄN, z-arvon perusteella katsotaan MITÄ poistetaan. tulos on sama, jos p-arvot ei mene nollaksi, mikä pitäisi jäädä kiinni bigBtagiin.
    if(max(p_values) < significance_level & !bigBtag){ 
      break #poistettavan muuttujan p-arvo pienempi kuin bonferroni-ehto, ja malli on valmis -> break.
    } else if(is.null(variable_to_remove)){
      break #viimeinen kierros kliinisillä, mitään ei saa tiputtaa mutta edellinen ei mennyt vielä läpi seulasta. -> break
    }
    else {
      cat(glue("{variable_to_remove}{ifelse(bigBtag, '* ', ' ')}"))
      latest_model <- update(latest_model, as.formula(paste(". ~ . - ", variable_to_remove)))
    }
  }
  
  forward_model <- latest_model
  
  if(do_forward){
    
    candidate_vars <- setdiff(names(full_model$coefficients), names(forward_model$coefficients))
    cat(glue("\n\nKoitetaan lisata muuttujia. Listan pituus {length(candidate_vars)}\n\n"))
    
    
    # Perform forward stepwise logistic regression
    while(length(candidate_vars) > 0) {
      best_p_value <- Inf
      best_variable <- NULL
      
      for(variable_to_add in candidate_vars) {
        temp_model <- update(forward_model, as.formula(paste(". ~ . + ", variable_to_add)))
        temp_summary <- summary(temp_model)
        temp_p_value <- temp_summary$coefficients[variable_to_add, "Pr(>|z|)"]
        if (temp_p_value < best_p_value) {
          best_p_value <- temp_p_value
          best_variable <- variable_to_add
          best_model <- temp_model
        }
      }
      
      if (best_p_value < significance_level) {
        cat(best_variable," ")
        forward_model <- best_model
        candidate_vars <- setdiff(candidate_vars, best_variable)
      } else {
        cat("\nEi muuta lisattavaa\n")
        break
      }
    }
    
    cat("Valmis malli:", names(forward_model$coefficients))
    
    final_model <- forward_model
  }
  
  return(final_model)
}



## tulosten kokoamisfunktioita ------------------------------------------------------

#Yleistetty tulostenkeräysfunktio. Pitäisi toimia myös C6646! (ei toimi, vaatii koko model-objektin. Tehdään tästä hiukan muokattu versio joka käyttää beta_to_probs():ia.)
gather_model_results <- function(model, data_list, lambda.type = "lambda.min") {
  
  if(inherits(model, "cv.glmnet")) { #glmnet (elanet)
    problist <- lapply(data_list,function(x) do.call(predict, list(object=model, newx=as.matrix(x[,-ncol(x)]), type="response", s=model[[lambda.type]])))
    beta <- predict(model, type = "coefficients", s = model$lambda.min)[,1]
    betatab <- data.frame(beta) %>% 
      mutate(var =rownames(.)) %>%       
      arrange(match(var, names(data_list[[1]]))) %>%       
      arrange(desc(var == "(Intercept)")) %>%
      select(beta)
  }else if(inherits(model,"glm")) { #glm (logreg, stepwise)
    problist <- lapply(data_list, function(x) do.call(predict, list(object = model, newdata = x, type = "response")))
    beta <- summary(model)$coefficients[,1]
    betatab <- data.frame(beta) %>% 
      mutate(var=rownames(.)) %>% 
      bind_rows(data.frame(var=setdiff(names(data_list[[1]]),names(beta)))) %>% 
      mutate(beta=replace_na(beta,0)) %>% 
      arrange(match(var, names(data_list[[1]]))) %>%
      arrange(desc(var == "(Intercept)")) %>%
      filter(var!="Y")
    rownames(betatab) <- betatab$var
    betatab <- betatab %>% select(beta)
  }else{
    stop("Model not supported")
  }
  

  names(problist) <- toupper(names(data_list))
  roclist <- lapply(seq_along(problist), function(i) roc(response = data_list[[i]]$Y, predictor = as.numeric(problist[[i]]), quiet = TRUE))
  names(roclist) <- toupper(names(data_list))
  auclist <- unlist(lapply(roclist, function(roclist_item) as.numeric(roclist_item$auc)))
  names(auclist) <- toupper(names(data_list))
  
  return(list(probs = problist, rocs = roclist, aucs = auclist, betas = betatab))
}



get_datasizetab <- function(datasets){
  
  dataset_names <- c("traindat", "test97", "test07", "test12")
  dataset_dims <- do.call(rbind, lapply(datasets[dataset_names], dim))
  dataset_classes <- do.call(rbind, lapply(datasets[dataset_names], function(df) table(df$Y)))
  datasizetab <- cbind(dataset_dims, dataset_classes)
  
  colnames(datasizetab) <- c("rows","columns","Y0","Y1")
  rownames(datasizetab) <- c("fr02","fr97","fr07","fr12")
  
  return(data.frame(datasizetab))
}


