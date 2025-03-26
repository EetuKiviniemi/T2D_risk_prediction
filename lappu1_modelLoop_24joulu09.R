run_modelLoop <- function(output_dir,fu_lengths,model_vars, fu_min, censor_dataset_glucose) {
  load_eventdata()
  
  indeksit <- expand.grid(model_vars, fu_lengths)
  indeksit$tunniste <- 1:nrow(indeksit)
  colnames(indeksit)<-c("model_vars","followup", "tunniste")
  
  clinical_vars <- get_fresh_clinical_vars_may24() #kliiniset muuttujat listattuina
  pval_divide_by_X <- get_bonferroni_correction_for_stepwise_mods(fu_min, censor_dataset_glucose) #bonferronikorjauskerroin :: sama 38 kaikilla ekskluusioilla. otetaan kuitenkin nyt tähän mukaan.
  
  
  modlist <- list()
  reslist <- list()
  last_fu_len <- 0

  for(i in 1:nrow(indeksit)){
    cat("\n\n kierroksen ", i, " alustukset\n")
    #kierroksen parametrit
    print(indeksit[i,])
    model_vars <- as.character(indeksit[i,1])
    fu_len <- as.numeric(indeksit[i,2])
    tunniste <- as.numeric(indeksit[i,3])   
    

    #Haetaan datasetit uudelleen vain kierroksilla, joilla followupin pituus muuttuu. eli sama data kaksi kierrosta peräkkäin.
    if(fu_len != last_fu_len){
      last_fu_len <- fu_len
      datasetname <- glue("modeldat_24joulu09_folloupw{fu_len}.RDS")
      
      if(recompile_model_data){
        datasets <- build_followupX_datasets(fu_len, fu_min, censor_dataset_glucose, doNewScale=T, clinical_vars)
        saveRDS(datasets, file.path("F:/profi5/foxgroup_rscript/paperi1/pipeline_24joulu09/prepped_data", datasetname))
      }else{
        datasets <- readRDS(file.path("F:/profi5/foxgroup_rscript/paperi1/pipeline_24joulu09/prepped_data", datasetname))
      }
      
      traindat <- datasets$traindat
      xtrain <- datasets$xtrain
      ytrain <- datasets$ytrain
      
      vars_to_scale <- datasets$vars_to_scale
      train_center <- datasets$scale_train_center
      train_std <- datasets$scale_train_std
      metab_to_mean0_vars <- datasets$metab_to_mean0_vars
      
      #datan koot talteen
      datasizetab <- get_datasizetab(datasets)
      cat("\n")
      print(datasizetab)
      cat("\n")
      
    }
    
    
    #fixaukset
    pfac_fixlist <- get_penaltyfac_ja_kiinteet(datasets$xtrain, clinical_vars, model_vars)
    pfac <- pfac_fixlist$penalty_factor
    fixed_vars <- pfac_fixlist$fixed_clinicals
    
    #testidatat
    clinonly_dats <- get_clinical_only_set(datasets,clinical_vars)
    fulltest_dats <- list(train02=datasets$traindat, 
                          test97=datasets$test97, 
                          test07=datasets$test07, 
                          test12=datasets$test12)
    
    ## KLIINISET SELITTÄJÄT ONLY -Malli
    cat("\nAjetaan logistinen regressio")
    logregmod <- build_logreg_model(clinonly_dats$train02)
    logregres <- gather_model_results(model = logregmod, data_list=clinonly_dats)
    
    # STEPWISE Malli
    cat("\n\nAjetaan stepwise")
    stepmod <- stepwise_logistic_regression(significance_level = 0.05/pval_divide_by_X, fixed_vars = fixed_vars, traindat=traindat, do_forward=T)
    stepres <- gather_model_results(model = stepmod, data_list = fulltest_dats)
    
    ## RIDGE  HOX! cv.glmnet: Note also that the results of cv.glmnet are random, since the folds are selected at random. siksi build_elanet_model sisältää set.seed(2024)
    cat("\n\nAjetaan ridge regressio")
    elamod_ridge <- build_elanet_model(xtrain, ytrain, alphaval=0, pfac)
    elanet_ridge_res <- gather_model_results(model = elamod_ridge, fulltest_dats)
    
    ## ELA05
    cat("\n\nAjetaan elastic net 0.5")
    elamod_05 <- build_elanet_model(xtrain, ytrain, alphaval=0.5, pfac)
    elanet_05_res <- gather_model_results(model = elamod_05, fulltest_dats)
    
    ## LASSO
    cat("\n\nAjetaan LASSO regressio")
    elamod_lasso <- build_elanet_model(xtrain, ytrain, alphaval=1, pfac)
    elanet_lasso_res <- gather_model_results(model = elamod_lasso, fulltest_dats)
    
    
    # Tallennetaan listolle
    listobj_name <- paste0("round", tunniste, "_", model_vars, "_", fu_len)
    reslist[[listobj_name]] <- list(ridge=elanet_ridge_res, 
                                    ela05=elanet_05_res, 
                                    lasso=elanet_lasso_res,  
                                    stepwise=stepres, 
                                    logreg=logregres, 
                                    datasizetab=datasizetab, 
                                    scaleparams=list(vars_to_scale=vars_to_scale, 
                                                     center=train_center, 
                                                     scale=train_std, 
                                                     metab_to_mean0_vars=metab_to_mean0_vars))
  }
  
  reslist_fp <- glue("{output_dir}/allRounds_resList.RDS")
  saveRDS(list(results=reslist, indeksit=indeksit), reslist_fp)
}
