#Lopullinen (toivottavasti) skriptipaketti. rajataan tämä vain tilanteisiin,
#joita lapussa oikeasti käytetään, eli fu=c(5,8), exclusions=main/normal.
#Nopetutaa ajoa ja turha työ jää tekemättä.

source("F:/profi5/foxgroup_rscript/paperi1/pipeline_24joulu09/lappu1_util_24joulu09.R")
source("F:/profi5/foxgroup_rscript/paperi1/pipeline_24joulu09/lappu1_modelLoop_24joulu09.R")
source("F:/profi5/foxgroup_rscript/paperi1/pipeline_24joulu09/resultpager_util_24joulu09.R")
source("F:/profi5/foxgroup_rscript/paperi1/pipeline_24joulu09/lappu1_resultpager_24joulu09.R")

muka_kapseli_fp <- "F:/profi5/foxgroup_rscript/muka_kapseli/elo09_eteenpain"
source(file.path(muka_kapseli_fp,"c66_kapselifunktiot.R"))
source(file.path(muka_kapseli_fp,"c66_datafunktiot.R"))
source(file.path(muka_kapseli_fp,"c66_analyze.R"))
sisaan_kapseliin_data_name <- "sisaan_kapseliin.RDS"
ulos_kapselista_fp <- "F:/profi5/foxgroup_rscript/paperi1/pipeline_24joulu09/from_kapseli"


run_analysis <- function(model_output_dir, censor_glucose, fu_min, recompile_model_data = F,
                         doModel=T, doResults=T, doKnit=T, 
                         build_c66_data=T, anal_c66_data=T, knitpdf=T) {
  
  #malli
  if(doModel){
    run_modelLoop(output_dir=model_output_dir,
                  fu_lengths=c(5,8),
                  model_vars=c("fix_none","fix_clinical"),
                  fu_min=fu_min, 
                  censor_dataset_glucose=censor_glucose,
                  recompile_model_data = recompile_model_data)
  }
  
  #resultpager
  if(doResults){
    aja_resultpager(output_dir=model_output_dir, 
                    muka_kapseli_fp,
                    sisaan_kapseliin_data_name, 
                    ulos_kapselista_fp, 
                    timeskip=as.list(rep(fu_min,2)),
                    fu_values=c(5,8), 
                    censor_glucose=censor_glucose,
                    build_c66_data=build_c66_data, 
                    anal_c66_data=anal_c66_data)
  }
  
  #results
  if(doKnit){
    rmarkdown::render(input = file.path(dirname(model_output_dir), "results_24joulu09.Rmd"),
                      params = list(pipeline_fp = model_output_dir, title = "results_25maalis11"), 
                      output_format = "pdf_document",
                      output_file = "results_25maalis11.pdf")
  }
}


# GLUKOOSISENSUROINTI KYLLÄ, EI AIKASENSUROINTIA (PÄÄMALLI)
run_analysis(model_output_dir="F:/profi5/foxgroup_rscript/paperi1/pipeline_24joulu09/model_results", 
             censor_glucose=T, fu_min=0, recompile_model_data = F,
             F,T,T, 
             build_c66_data=F, anal_c66_data=F, knitpdf=T)



