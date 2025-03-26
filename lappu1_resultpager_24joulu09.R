#funktioon on hardkoodattu tiettyjä seikkoja, kuten että followupit on aina 5,8,10, ja datasetit on fr02 - fr12 + c6646. 
aja_resultpager <- function(output_dir, muka_kapseli_fp, sisaan_kapseliin_data_name, 
                            timeskip, fu_values, censor_glucose, build_c66_data, anal_c66_data) {
  
  start <- Sys.time()
  
  cat("\nLuetaan mallien tulostiedosto...")
  #kiinteät parametrit jotka on samat joka kierroksella. =~ Kliiniset muuttujat.
  clinical_vars <- get_fresh_clinical_vars_may24()
  reslist_fr <- readRDS(file.path(output_dir, "allRounds_resList.RDS"))
  indeksit <- reslist_fr[[2]] %>% data.frame
  
  ## prep----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  if(anal_c66_data){
    prep_c66_kapseli_metadata(reslist_fr, clinical_vars, 
                              censor_glucose = censor_glucose, 
                              kapseliin_fp=file.path(muka_kapseli_fp, sisaan_kapseliin_data_name))
  }
  
  cat("\nLuetaan C6646 datan Kapselista tuotuja tuloksia...")
  
  c66_results <- list.files("F:/profi5/foxgroup_rscript/muka_kapseli/Tulokset_THL_1230_14.05.00_2022_a89_21022025",pattern=".csv", full.names=T) %>%
    map(fread) %>%
    setNames(c("auc_ci", "calstat.tab", "clinvartab", "datameans", "datasizetab", "rawdata_n", "roctest"))
  
  ## reslist----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  reslist<- reslist_fr
  all_tags <- kokoa_alltags() 
  # dummydata <- kokoa_resultpager_ylist(timeskip=timeskip, fulens=c(5,8), censor_glucose=censor_glucose, clinical_vars)
  # saveRDS(dummydata, file.path(output_dir, "resultpager_dummydata.RDS"))
  dummydata <- readRDS(file.path(output_dir, "resultpager_dummydata.RDS"))
  
  ## tab1 ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  cat("\nKootaan taulua 1...")
  fulltbl1 <- kokoa_valmis_table1_taulu(dummydata, reslist, c66_results, clinical_vars)
  saveRDS(fulltbl1, glue("{output_dir}/clinvartab_all.RDS"))
  
  ## auctab ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  cat("\nKootaan AUC-lukuja...")
  auctablist <- kokoa_valmis_AUC_taulu(reslist, c66_results)
  prepped_auctabs <- beautify_auctabs(auctablist)
  saveRDS(prepped_auctabs$fu5, glue("{output_dir}/auctab_fu5.RDS"))
  saveRDS(prepped_auctabs$all, glue("{output_dir}/auctab_fu810.RDS"))
  
  ## kalibraatiot ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  cat("\nKootaan kalibraatioita...")
  calibrations <- kokoa_valmis_calibration_taulu(ylist=dummydata$ylist, c66_results, all_tags, reslist)
  calstats_allmods <- kokoa_calstats_tagged(calstat.tab=calibrations$caltab, all_tags)
  cleaned_caltab <- beautify_caltab(calstats_allmods)
  saveRDS(cleaned_caltab$fu5, glue("{output_dir}/caltab_fu5.RDS"))
  saveRDS(cleaned_caltab$fuall, glue("{output_dir}/caltab_fuall.RDS"))
  
  ## roc.test taulut ------------------------------------------------------
  cat("\nKootaan roctest-tauluja...")
  problist <- get_full_problist(reslist, dummydata)
  roctest_tabs <- kokoa_roctest_taulu(problist=problist, c66_results, dummydata)
  prepped_roctesttabs <- beautify_roctesttab(roctest_tabs)
  saveRDS(prepped_roctesttabs$fu5, glue("{output_dir}/roctesttab_fu5.RDS"))
  saveRDS(prepped_roctesttabs$fu810, glue("{output_dir}/roctesttab_fu810.RDS"))
  
  ## betat ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  cat("\nKootaan beta-tauluja ja kuvia...")
  betalist <- kokoa_betalist(reslist)
  betataulut <- kokoa_betavals_by_fu(betalist, clinical_vars)
  saveRDS(betataulut, glue("{output_dir}/fulen_betataulut.RDS"))
  
  betakuvat <- kokoa_valmiit_betaval_kuvat(betalist)
  saveRDS(betakuvat, glue("{output_dir}/betaplotit.RDS"))
  
  ## datameans  ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  datameansd <- build_meanSD_tab(dummydata, c66meansd=c66_results$datameans)
  saveRDS(datameansd, glue("{output_dir}/data_means_sds.RDS"))
  
  ## calib tallennus ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  cat("\nTallennetaan kuvia...")
  prepped_calplots <- prep_calibplots(fu_values, calibrations)
  combine_4by5_calplot(prepped_calplots, 
                       c66_plot_fp = "F:/profi5/foxgroup_rscript/muka_kapseli/Tulokset_THL_1230_14.05.00_2022_a89_21022025/c66plot.pdf",
                       output_dir = output_dir,
                       target_filename = "doublecal_plots_4by5.pdf")
  

  ## Valmista! ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  end <- Sys.time()
  cat("\nValmista tuli, ajassa ", round(end-start,2), "min")
}
