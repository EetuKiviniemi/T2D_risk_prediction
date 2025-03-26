# Yleisiä apufunktioita ---------
# kokoaa vektorin mahdollisista modeltageistä.
kokoa_alltags <- function(aikapisteet=c("FR02","FR97","FR07","FR12","C66"), 
                          mods=c("stepwise","lasso","ridge","ela05","logreg"), 
                          fu_lengths=c(5,8), 
                          kiinteet=c("fixnone","fixclinical")){

  combined_params <- expand.grid(kiinteet = kiinteet, followup = fu_lengths, ap = aikapisteet, mt = mods,stringsAsFactors = F)
  combined_params <- combined_params[!(combined_params$ap == "C66" & combined_params$followup == 10),] #c66 ei ole tehty 10v followup malleja
  combined_params <- combined_params[!(combined_params$mt == "logreg" & combined_params$kiinteet == "fixnone"),] 
  combined_params$kiinteet[combined_params$mt == "logreg"] <- "onlyclinicals"
  #logregmallit tehtiin sekä none että clinfix kierroksillla. Eroa näiden välillä ei ole, otetaan vain toiset niistä.
  combined_params <- unique(combined_params)
  
  combined_params <- combined_params %>% 
    mutate(tag= paste(ap, paste0("FU",followup), kiinteet, mt, sep="_"))
  
  return(combined_params$tag)
}



get_model_problist <- function(baselist, z, inds=1:6, indeksit){
  lapply(inds, function(x) baselist[[x]][[1]][[z]]) %>%
    do.call(cbind,.) %>%
    data.frame %>%
    setNames(glue("{indeksit$model_vars[inds]}_FU{trimws(indeksit$followup[inds])}"))
}

get_full_problist <- function(reslist, dummydata){
  mods <- c("ridge", "ela05", "lasso", "stepwise", "logreg")
  frsets <- c("fr02","fr97","fr07","fr12")
  dataset_inds <- 1:4
  ylist <-dummydata$ylist
  indeksit <- reslist[[2]]%>% data.frame
  #probsit listana
  problist <- lapply(dataset_inds, function(x) lapply(mods, function(y){
    mod <- map(reslist[[1]], y)
    get_model_problist(mod, x, inds=1:length(reslist[[1]]), indeksit=indeksit)
  }))
  names(problist) <- frsets
  problist <- lapply(problist, function(x) setNames(x,mods))
  return(problist)
}



# Taulu 1 kokoaminen ---------
# Kokoaa ja palauttaa listan y-arvoista sekä yhden datasets-listan
kokoa_resultpager_ylist <- function(timeskip=0,fulens=c(5,8,10),censor_glucose=T, clinical_vars) {
  cat("\nKootaan listaa Y-arvoista eri datoilla ja followupeilla \n")
  y02tab <- data.frame()
  y97tab <- data.frame()
  y07tab <- data.frame()
  y12tab <- data.frame()
  
  #looppi haluaa lista-muotoisen timeskipin
  if(typeof(timeskip) != "list") timeskip <- as.list(rep(timeskip, length(fulens)))
  
  #kootaan ylist. loopissa rbindillä koska cbind ei onnistu tyhjään tauluun. käännetään ja siistitään.
  for(i in seq_along(fulens)){
    cat("\t kootaan followupia", fulens[i], "\n")
    datasets <- build_followupX_datasets(fu_len=fulens[i],  
                                         fu_min=timeskip[[i]], 
                                         censor_dataset_glucose=censor_glucose,
                                         doNewScale=F,
                                         clinical_vars=clinical_vars)
    data02 <- datasets$traindat
    data97 <- datasets$test97
    data07 <- datasets$test07
    data12 <- datasets$test12
    
    #kootaan ylist resultpageria varten
    y02tab <- rbind(y02tab, data02$Y)
    y97tab <- rbind(y97tab, data97$Y)
    y07tab <- rbind(y07tab, data07$Y)
    y12tab <- rbind(y12tab, data12$Y)
    
  }
  
  ylist <- list(fr02_5810=t(y02tab)%>%matrix(ncol=length(fulens))%>%data.frame()%>%setNames(rep("Y",, length(fulens))),
                fr97_5810=t(y97tab)%>%matrix(ncol=length(fulens))%>%data.frame()%>%setNames(rep("Y",, length(fulens))),
                fr07_5810=t(y07tab)%>%matrix(ncol=length(fulens))%>%data.frame()%>%setNames(rep("Y",, length(fulens))),
                fr12_5810=t(y12tab)%>%matrix(ncol=length(fulens))%>%data.frame()%>%setNames(rep("Y",, length(fulens)))) 
  
  
  return(list(datasets=datasets, ylist=ylist))
}


kokoa_valmis_clinvar_taulu <- function(dummydata, c66_results) {
  clinical_vars <- get_tab1_vars() 
  #kootaan kliinisten taulu viimeisen kierroksen ajoista. Followupilla kun ei ole väliä. Otetaan mukaan myös C6646 data.
  table1dat <- rbind(cbind(dummydata$datasets$traindat %>% select(any_of(clinical_vars)), aikapiste="FR02"),
                     cbind(dummydata$datasets$test97 %>% select(any_of(clinical_vars)), aikapiste="FR97"),
                     cbind(dummydata$datasets$test07 %>% select(any_of(clinical_vars)), aikapiste="FR07"),
                     cbind(dummydata$datasets$test12 %>% select(any_of(clinical_vars)), aikapiste="FR12"))
  table1dat$aikapiste <- factor(table1dat$aikapiste, levels=c("FR02","FR97","FR07","FR12"))
  
  clinvartbl.fr <- table1dat %>% 
    tbl_summary(by="aikapiste",
                statistic=list(
                  all_continuous() ~ "{mean} ({sd})")) 
  # ALKI1 ~ "{median} ({p25}, {p75})")) #6.3.2025 // alkoholista mean (sd).
  
  #kliinisten taulu KAPSELISSA (input clinical_vars)
  if(!is.null(c66_results)){
    clinvartbl.c66 <- c66_results$clinvartab
    datamean_alk1 <- c66_results$datameans %>%
      filter(var == "ALKI1") %>%
      mutate(
        alkimean = str_split_i(c66, " ", 1) %>%as.numeric %>% round(.,0),
        alkisd = str_split_i(c66, " ", 2) %>%gsub("\\(|\\)", "", .) %>%as.numeric,
        alki_meansd = paste0(alkimean, " (", alkisd, ")")
      )
    clinvartbl.c66$`**C6646**, N = 4,785`[clinvartbl.c66$`**Characteristic**` == "ALKI1"] <- datamean_alk1$alki_meansd
    clinvartbl <- bind_cols(clinvartbl.fr%>%as.data.frame, clinvartbl.c66%>%select(-1))
  }else clinvartbl <- clinvartbl.fr %>% as.data.frame
  
  return(clinvartbl)
}


#datasizetab finriskeille
kokoa_valmis_datasize_taulu <- function(reslist, c66_results) {
  
  res_inds <- which((1:length(reslist$results)) %% 2 == 1)
  reslist_fulens <- str_split(names(reslist$results), "_",simplify=T)[,4] %>% unique
  possible_names <- c("FU5_cases","FU8_cases","FU10_cases")
  used_names <- possible_names[grep(paste0(reslist_fulens,collapse="|"), possible_names)]
  
  dst <-  map(reslist$results[res_inds], pluck, "datasizetab")
  
  fakeYlist <- lapply(1:length(res_inds), function(y){
    apply(dst[[y]], 1, function(x){
      c(rep(0, x[3]), rep(1, x[4]))
    }) %>% 
      Map(cbind, names(.), . ) %>% 
      do.call(rbind, .) %>% 
      data.frame
  }) %>% 
    do.call(cbind,.) %>%
    data.frame() %>%
    select(2,4,1) %>%
    setnames(., c(used_names, "aikapiste"))
  
  fakeYtab <- fakeYlist %>%
    data.frame %>%
    mutate(aikapiste=toupper(aikapiste),
           aikapiste=factor(aikapiste, levels=c("FR02","FR97","FR07","FR12")),
           across(any_of(used_names), as.numeric))
  
  datasizetbl.fr <- fakeYtab %>%
    tbl_summary(by=aikapiste) %>%
    remove_row_type(type="level", level_value="0")
  
  # datasizetab
  if(!is.null(c66_results)){
    datasizetab.c66 <- c66_results$datasizetab %>% slice(1:2)
    datasizetbl <-   bind_cols(datasizetbl.fr%>%as.data.frame, datasizetab.c66%>%select(-1))
  }else datasizetbl <- datasizetbl.fr %>% as.data.frame
  
  return(datasizetbl)
}


kokoa_valmis_table1_taulu <- function(dummydata, reslist, c66_results, clinical_vars) {
  
  clinvartbl <- kokoa_valmis_clinvar_taulu(dummydata, c66_results)
  datasizetbl <- kokoa_valmis_datasize_taulu(reslist, c66_results)
  
  #ännät ennen mitään filttereitä. (tallenetaan rds, nopeuttaa ajoa)
  # raw_N_counts <- data.frame(rawN=sapply(c(02,97,07,12), FUN= function(x) nrow(get_finriskdata(x))), fr=c(02,97,07,12))
  # saveRDS(raw_N_counts, "F:/profi5/data/util/raw_N_count_cache.RDS")
  raw_N_counts <- readRDS("F:/profi5/data/util/raw_N_count_cache.RDS")
  raw_N_counts$rawN <- as.character(raw_N_counts$rawN)
  
  #yhdistetty datasize ja clinvar tab
  fulltbl1 <- bind_rows(list(clinvartbl,datasizetbl))
  if(!is.null(c66_results)){
    fulltbl1 <- rbind(fulltbl1, c("rawN", raw_N_counts$rawN, as.numeric(c66_results$rawdata_n)))
  }else{
    fulltbl1 <- rbind(fulltbl1, c("rawN", raw_N_counts$rawN))
  }
  return(fulltbl1)
}



# C66 datan käsittely ----
# Yhdistää reslistit FR datoille ja C66 datalle kapselista
prep_c66_kapseli_metadata <- function(reslist_fr, clinical_vars, censor_glucose=T, fu_min=0, kapseliin_fp) {
  kapseliin_viemiset <- list()
  
  #indeksit
  kapseliin_viemiset[["indeksit"]] <- reslist_fr[[2]] %>% as.data.frame
  
  #muuttujanimet
  suppressWarnings(
  datasetup <- setup_xy_DIAB(followup_length = 5, 
                              yr=2,
                              censor_glucose = censor_glucose, 
                              fu_min = fu_min)$fullset %>% 
    make_dummy_binary_vars(.) %>% 
    data.frame
  )
  nameFR <- names(datasetup)
  kapseliin_viemiset[["nameFR"]] <- nameFR
  
  #kliiniset muuttujat
  kapseliin_viemiset[["clinical_vars"]] <- clinical_vars
  kapseliin_viemiset[["clinical_vars_old"]] <- get_tab1_vars()
  
  #skaalaus parametrit
  scale_params <- reslist_fr[["results"]][["round1_fix_none_5"]][["scaleparams"]] #skaalausparametrit on aina samat, eivät riiput datasta eikä followupista tai muustakaan
  kapseliin_viemiset[["scale_params"]]<-scale_params
  
  
  #kierrosten ja mallien nimet
  roundnames <- names(reslist_fr[[1]])
  kapseliin_viemiset[["roundnames"]] <- roundnames
  modnames <- names(reslist_fr[[1]][[1]][1:5])
  kapseliin_viemiset[["modnames"]] <- modnames
  
  #betojen lista
  betalist.raw <- lapply(roundnames, function(x) lapply(modnames, function(y){
    bdf <- as.vector(reslist_fr[["results"]][[x]][[y]][["betas"]][["beta"]])
    names(bdf)<-rownames(reslist_fr[["results"]][[x]][[y]][["betas"]])
    bdf
  }))
  
  names(betalist.raw) <- roundnames
  for(i in 1:length(roundnames)){ names(betalist.raw[[roundnames[i]]]) <- modnames}
  kapseliin_viemiset[["betalist"]] <- betalist.raw
  
  
  #tallennetaan rds-muodossa
  saveRDS(kapseliin_viemiset, kapseliin_fp)
}

merge_nested_lists <- function(main_list, c66list) {
  names_to_merge <- intersect(names(c66list), 
                              names(main_list))
  
  for (first_level_name in names_to_merge) {
    second_level_names <- intersect(names(c66list[[first_level_name]]), 
                                    names(main_list[[first_level_name]]))
    
    for (second_level_name in second_level_names) {
      third_level_names <- intersect(names(c66list[[first_level_name]][[second_level_name]]),
                                     names(main_list[[first_level_name]][[second_level_name]]))
      
      for (third_level_name in third_level_names) {
        main_list[[first_level_name]][[second_level_name]][[third_level_name]] <-
          c(main_list[[first_level_name]][[second_level_name]][[third_level_name]], 
            c66list[[first_level_name]][[second_level_name]][[third_level_name]])
      }
    }
  }
  
  return(main_list)
}



# AUC taulut -------------
# Hakee followupin kierrosnimestä
get_roundname_fulen <- function(strings) {
  numbers <- sapply(strings, function(string) {
    split_string <- str_split(string, "_")[[1]]
    last_item <- tail(split_string, 1)
    number <- str_extract(last_item, "\\d+")
    as.numeric(number)
  })
  return(numbers)
}

#hakee käytettyjen nonzero-betojen lukumäärät
get_betacounts <- function(modelname, reslist, indeksit){
  lapply(1:length(reslist[[1]]), function(x) sum(as.numeric((map(reslist[[1]], modelname)[[x]][[4]][,1]))!=0)) %>%
    unlist %>% 
    cbind(tunniste=indeksit$tunniste, beta=.) %>%
    apply(.,2,as.numeric) %>%
    cbind(.,model=modelname) %>%
    data.frame
}


#kokoaa malli-kohtaiset AUC-luvut + betacountin
kokoa_model_aucit <- function(model, reslist, indeksit, c66_results) {
  #rajaa indeksit niin että mukaan tulee vain ne followupit joista meillä on aucceja (kts. c6646)
  indeksit <- indeksit[as.numeric(indeksit$followup) %in% get_roundname_fulen(names(reslist[[1]])),] 
  if(model=="logreg") indeksit <- indeksit[indeksit$model_vars == "fix_clinical",]
  range <- as.numeric(indeksit$tunniste)
  roclist <- lapply(range, function(x) map(reslist[[1]], model)[[x]][[2]]) %>%
    lapply(., function(x) x[-5])

  auc_cilist <- lapply(roclist, function(sublist) {
    lapply(sublist, function(roc_object) {
      as.numeric(ci.auc(roc_object))
    })
  })  
  
  #c66 erikseen, ei tuoda kapselista raakoja rocceja.
  c66_auc_ci <- c66_results$auc_ci %>% 
    mutate(longtag = gsub("_auc_ci", "", longtag),
           tunniste =str_split_i(longtag, "_",1) %>% gsub("round","",.) %>% as.numeric) %>%
    filter(grepl(model,longtag)) %>%
    mutate(model = model) %>%
    select(ci_lower, auc, ci_upper, model, tunniste) %>%
    as.data.frame
  rownames(c66_auc_ci) <- paste0("C66",1:nrow(c66_auc_ci))
    
    
  ci_df_list <- lapply(auc_cilist, function(sublist) data.frame(do.call(rbind, sublist), stringsAsFactors = FALSE))
  ci_df_list <- lapply(1:nrow(indeksit), function(i) cbind(ci_df_list[[i]], model=model, tunniste=indeksit$tunniste[i]))
  ci_table <- ci_df_list %>%
    bind_rows %>%
    setNames(c("ci_lower", "auc", "ci_upper", "model", "tunniste")) %>%
    rbind(.,c66_auc_ci) %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    mutate(dataset = rownames(.)) %>%
    {rownames(.) <- NULL; .} %>% 
    mutate(
      dataset = case_when(
        grepl("97",dataset) ~ "FR97",
        grepl("07",dataset) ~ "FR07",
        grepl("12",dataset) ~ "FR12",
        grepl("02",dataset) ~ "FR02",
        grepl("C66",dataset) ~ "C66"
      ),
      auc_ci = glue("[{format3d(ci_lower,3)}, {format3d(ci_upper,3)}]")
    ) %>%
    left_join(indeksit, by="tunniste") %>%
    filter(!(model=="logreg" & model_vars=="fix_none")) %>%
    mutate(model_vars=gsub("_","",model_vars))
  
  
  betacount <- get_betacounts(model, reslist, indeksit) %>% 
    mutate(tunniste=as.integer(tunniste)) %>% 
    unique
  
  auc_table <- ci_table %>% 
    left_join(betacount, by=c("model","tunniste")) %>%
    select(dataset,model_vars, followup, model,auc,auc_ci, beta) 
  
  return(auc_table)
}


kokoa_valmis_AUC_taulu <- function(reslist, c66_results) {
  indeksit <- as.data.frame(reslist[[2]])
  models <- c("ridge", "ela05", "lasso", "stepwise", "logreg")
  
  auctab_long <- lapply(models, function(model) kokoa_model_aucit(model, reslist, indeksit, c66_results)) %>% bind_rows %>% distinct
  
  return(auctab_long)
}


# Kalibraatio -----

get_Y_values <- function(tagi, calib_ylist) {
  
  tagparts <- tagi %>%
    gsub("FU","",.) %>%
    str_split(.,"_") %>%
    lapply(., function(i) gsub("fix","fix_",i)) %>%
    unlist %>%
    setNames(c("set","fu","data","mod"))
  
  
  aikapiste_indices <- grep(tagparts[["set"]], names(calib_ylist), ignore.case=T)
  
  if (length(aikapiste_indices) == 0) {
    stop("No Y values found for the specified aikapiste.")
  }
  
  if(tagparts[["fu"]] == "5") sub_index <- 1
  if(tagparts[["fu"]] == "8") sub_index <- 2
  if(tagparts[["fu"]] == "10") sub_index <- 3
  
  y_values <- calib_ylist[[aikapiste_indices]][[sub_index]]
  return(y_values)
}

modeltag_to_probvec <- function(tagi, reslist){
  tagparts <- tagi %>%
    gsub("FU","",.) %>%
    str_split(.,"_") %>%
    lapply(., function(i) gsub("fix","fix_",i)) %>%
    unlist %>%
    setNames(c("set","fu","data","mod"))
  tagparts <- gsub("onlyclinicals", "fix_clinical", tagparts)
  
  reslist_level1 <- reslist[[1]][grepl(tagparts[["data"]],names(reslist[[1]])) & 
                                   grepl(paste0("_", tagparts[["fu"]]),names(reslist[[1]]))] #haluttu datasetup ja followup
  reslist_level2 <- reslist_level1[[1]][grepl(tagparts[["mod"]],names(reslist_level1[[1]]))] #haluttu malli
  reslist_level3 <- reslist_level2[[1]]$probs #probs-lista
  
  #probs-listojen nimet on train02 tyylisiä, erilailla tehtyjä kun tägit. matchataan ne haluttuun tagiin.
  set_name_pairs <- cbind(c("FR02","FR97","FR07","FR12","C66"),c("TRAIN02","TEST97","TEST07","TEST12", ""))
  settag <- set_name_pairs[set_name_pairs[,1] == tagparts[["set"]], 2]
  settag_index <- which(names(reslist_level3)%>%gsub("TRAINDAT","TRAIN02",.) == settag)
  
  probvec <- reslist_level3[[settag_index]] %>% as.numeric
  
  return(probvec)
}




### Kalibraatiokäyrän piirtäminen ----
get_segment_dfs <- function(preds, y, length.seg, line.bins){
  x <- preds
  bins <- seq(0, min(1, max(x)), length = 101)
  x <- x[x >= 0 & x <= 1]
  f0 <- table(cut(x[y == 0], bins))
  f1 <- table(cut(x[y == 1], bins))
  bins0 <- (bins[-101])
  bins1 <- (bins[-101])
  maxf <- max(f0, f1)
  f0 <- (0.1 * f0) / max(f0)
  f1 <- (0.1 * f1) / max(f1)
  
  line_bins <- -0.05
  length_seg <- 1
  segment_df.pos <- data.frame(ystart=-0.05, yend=as.numeric(length.seg * f1 + line.bins), xstart=bins1, xend=bins1)
  segment_df.neg <- data.frame(ystart=-0.05, yend=as.numeric(length.seg * -f0 + line.bins), xstart=bins1, xend=bins1)
  
  return(list(pos=segment_df.pos, neg=segment_df.neg))
}

ggplot_calcurve <- function(calstat, preds, response, main_title, xlims=c(-0.02, 1), ylims=c(-0.115, 1), dist.label=0.04, dist.label2=0.03){
  vpxy <- calstat$CalibrationCurves$FlexibleCalibration

  # rug_cases <- data.frame(
  #   x = preds[response == 1],
  #   ystart = -0.005,  # Y start position for cases
  #   yend = -0.055,    # Y end position for cases (same as start for a vertical line)
  #   xstart = preds[response == 1],  # X start position for cases
  #   xend = preds[response == 1],    # X end position for cases (same as start for a vertical line)
  #   response = "cases"
  # )
  # 
  # rug_controls <- data.frame(
  #   x = preds[response == 0],
  #   ystart = -0.06,  # Y start position for controls
  #   yend = -0.11,    # Y end position for controls (same as start for a vertical line)
  #   xstart = preds[response == 0],  # X start position for controls
  #   xend = preds[response == 0],    # X end position for controls (same as start for a vertical line)
  #   response = "controls"
  # )
  
  gg_calcurve <- ggplot(vpxy, aes(x = x, y = y)) +
    annotate("segment", x = 0, y = 0, xend = 1, yend = 1, color = "red") +  
    geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "gray", alpha = 0.2) +
    geom_line(color = "black") +
    # geom_segment(data = rug_cases, aes(x = xstart, y = ystart, xend = xend, yend = yend)) +
    # geom_segment(data = rug_controls, aes(x = xstart, y = ystart, xend = xend, yend = yend)) +
    geom_segment(x=0, y=-0.0575, xend=1, yend=-0.0575, color="black") +
    annotate("text", x = 0 - dist.label, y = -0.0575 + dist.label2, label = "1") +
    annotate("text", x = 0 - dist.label, y = -0.0575 - dist.label2, label = "0") +
    labs(title = main_title, x = "Predicted probability", y = "Observed proportion") +
    theme_minimal() +
    coord_cartesian(xlim = xlims, ylim = ylims)
  
  
  return(gg_calcurve)
}



### Kalibraatiotaulun kokoaminen --------

#pyörittää loopilla kalibraariota kaikille annettujen aikapisteiden malleille

mle_logit <- function(y, logit) {
  loglik <- function(par, X, y) {
    eta <- X %*% par
    p <- 1 / (1 + exp(-eta))
    -sum(y * log(p) + (1 - y) * log(1 - p))
  }
  
  X <- cbind(1, logit)
  start <- rep(0, ncol(X))
  fit <- optim(start, loglik, X = X, y = y, method = "BFGS", hessian = TRUE)
  
  coef_mle <- fit$par[2] # slope = 2
  se <- sqrt(diag(solve(fit$hessian)))[2]
  ci_lower <- coef_mle - qnorm(0.975) * se
  ci_upper <- coef_mle + qnorm(0.975) * se
  
  data.frame(
    slope = coef_mle,
    ci_lower = ci_lower,
    ci_upper = ci_upper
  )
}

mle_intercept <- function(y, logit) {
  loglik <- function(intercept, offset, y) {
    eta <- intercept + offset
    p <- 1 / (1 + exp(-eta))
    -sum(y * log(p) + (1 - y) * log(1 - p))
  }
  
  # Fit model with offset and estimate intercept
  start <- 0
  fit <- optim(start, loglik, offset = logit, y = y, method = "BFGS", hessian = TRUE)
  
  intercept <- fit$par
  se <- sqrt(1 / fit$hessian)
  ci_lower <- intercept - qnorm(0.975) * se
  ci_upper <- intercept + qnorm(0.975) * se
  
  data.frame(
    intercept = intercept,
    ci_lower = ci_lower,
    ci_upper = ci_upper
  )
}

val.prob.simple <- function(p, y, cl.level = 0.95) {
  # Input validation
  if (!all(y %in% 0:1)) {
    stop("The binary outcome can only contain values 0 and 1")
  }
  if (length(p) != length(y)) {
    stop("lengths of p and y do not agree")
  }
  if (any(p > 1 | p < 0)) {
    stop("Probabilities must be between 0 and 1")
  }
  
  # Logistic calibration
  logit <- log(p/(1 - p))
  y <- y[order(p)]
  n <- length(p)
  
  #sort in order of preds
  logit <- logit[order(p)]
  p <- p[order(p)]
  i <- !is.infinite(logit)
  
  #calibration slope ja calibration intercept
  slope <- mle_logit(y,logit)
  inter <- mle_intercept(y[i], logit[i]) 
  
  # Loess smoothing
  smooth <- loess(y ~ p, degree = 2)
  smooth.pred <- predict(smooth, se = TRUE)
  
  # Calculate confidence intervals for smooth curve
  a <- (1 - cl.level)
  smooth.lower <- smooth.pred$fit - qnorm(1 - a/2) * smooth.pred$se.fit
  smooth.upper <- smooth.pred$fit + qnorm(1 - a/2) * smooth.pred$se.fit
  
  # Ensure bounds are within [0,1]
  smooth.lower <- pmax(0, pmin(1, smooth.lower))
  smooth.upper <- pmax(0, pmin(1, smooth.upper))
  
  # Calculate metrics
  # Brier score
  B <- sum((p - y)^2)/n
  Bmax <- mean(y) * (1 - mean(y))^2 + (1 - mean(y)) * mean(y)^2
  Bscaled <- 1 - B/Bmax
  
  # C-statistic (ROC AUC)
  n.pos <- sum(y)
  n.neg <- length(y) - n.pos
  concordant <- 0
  for(i in 1:length(y)) {
    if(y[i] == 1) {
      concordant <- concordant + sum(p[i] > p[y == 0])
    }
  }
  c.stat <- concordant / (n.pos * n.neg)
  
  # C-statistic confidence intervals
  # Using variance estimator from Hanley & McNeil
  q1 <- c.stat / (2 - c.stat)
  q2 <- 2 * c.stat^2 / (1 + c.stat)
  se <- sqrt((c.stat * (1 - c.stat) + (n.pos - 1) * (q1 - c.stat^2) + 
                (n.neg - 1) * (q2 - c.stat^2)) / (n.pos * n.neg))
  c.lower <- c.stat - qnorm(1 - a/2) * se
  c.upper <- c.stat + qnorm(1 - a/2) * se
  
  # Prepare results
  results <- list(
    stats = list(
      "C (ROC)" = c.stat,
      "Intercept" = inter,
      "Slope" = slope,
      "Brier" = B,
      "Brier scaled" = Bscaled
    ),
    CalibrationCurves = list(
      FlexibleCalibration = data.frame(
        x = p,
        y = smooth.pred$fit,
        ymin = smooth.lower,
        ymax = smooth.upper
      )
    )
  )
  
  #käyrän pisteet lasketaan seuraavasti:
  # y <- y[order(p)]
  # p <- p[order(p)]
  # smooth <- loess(y ~ p, degree = 2)
  # smooth.pred <- predict(smooth, se = TRUE)
  # smooth.lower <- smooth.pred$fit - qnorm(1 - a/2) * smooth.pred$se.fit
  # smooth.upper <- smooth.pred$fit + qnorm(1 - a/2) * smooth.pred$se.fit
  # smooth.lower <- pmax(0, pmin(1, smooth.lower))
  # smooth.upper <- pmax(0, pmin(1, smooth.upper))
  # FlexibleCalibration = data.frame(
  #   x = p,
  #   y = smooth.pred$fit,
  #   ymin = smooth.lower,
  #   ymax = smooth.upper
  # )
  #missä input y on todellinen vaste, ja input p on ennustemallin henkilölle laskema riski.
  
  return(results)
}



loop_calibration <- function(calib_ylist, plot_tags, all_tags, reslist){
  calstat.list<-list()
  calplots.list<-list()
  
  cat("\n|",paste0(rep("-", length(all_tags)),collapse=""),"\n")
  cat("| ")
  
  for(i in 1:length(all_tags)){
    #edistymisprinttejä
    tagi <- all_tags[i]
    if(tagi %in% plot_tags) progressmark <- "p"
    if(!tagi %in% plot_tags) progressmark <- "x"
    cat(progressmark) 
    
    preds <- modeltag_to_probvec(tagi, reslist)
    yvalue <- get_Y_values(tagi, calib_ylist)
    
    #lasketaan kalibraatio. plotit vaan jos niitä halutaan.
    # suppressMessages(suppressWarnings({
    #   calibrated <- val.prob.ci.2(preds, yvalue, smooth="loess")
    # }))
    
    calibrated <- val.prob.simple(preds, yvalue)
    calplot<-NULL
    if(tagi %in% plot_tags) calplot <- ggplot_calcurve(calstat=calibrated, preds=preds, response=yvalue, main_title=tagi)
    
    calstats <- calibrated$stats %>% 
      data.frame %>%
      setNames(c("C (ROC)", 
                 "intercept.Point estimate", "intercept.Lower confidence limit", "intercept.Upper confidence limit", 
                 "slope.Point estimate", "slope.Lower confidence limit.2.5 %", "slope.Upper confidence limit.97.5 %", 
                 "Brier","Brier scaled"))
    
    #tulokset listalla
    calstat.list[[i]] <- data.frame(tag=tagi, calstats)
    calplots.list[[i]] <- calplot
    
  }
  return(list(calstat.list=calstat.list, calplots.list=calplots.list))
}

kokoa_valmis_calibration_taulu <- function(ylist, c66_results, all_tags, reslist) {
  
  c66results_actual_contents <- ifelse(any(grepl("FRcomb", c66_results$calstat.tab$tag)), "FRcomb","C66")
  
  plot_tags <- all_tags %>% grep("FU8|FU5",.,value=T) %>% grep("ridge|logreg",., value=T)
  plot_non66tags <- plot_tags %>% gsub("C66", c66results_actual_contents,. ) %>% discard(grepl(c66results_actual_contents,.))
  all_non66tags <- all_tags %>%  gsub("C66", c66results_actual_contents,. ) %>% discard(grepl(c66results_actual_contents,.))
  
  cals <- loop_calibration(calib_ylist=ylist, plot_tags=plot_non66tags, all_tags=all_non66tags, reslist) #ajetaan tässä kalibraatiot muille kuin c66-malleille.
  calstat.list.fr <- cals$calstat.list
  calplots.list.fr <- cals$calplots.list
  calstat.tab.fr <- bind_rows(calstat.list.fr) %>% data.frame()
  rownames(calstat.tab.fr) <- NULL
  
  if(!is.null(c66_results)){
    calstat.tab.c66 <- c66_results$calstat.tab %>%
      mutate(tag = gsub("fixclinical_logreg", "onlyclinicals_logreg", tag)) %>%
      filter(tag %in% (all_tags%>%gsub("C66", c66results_actual_contents,. ))) %>% 
      select(-any_of(c('slope.Lower.confidence.limit',  'slope.Upper.confidence.limit')))
    calstat.tab <- rbind(calstat.tab.fr, calstat.tab.c66)
    calplots.list <- calplots.list.fr
  }else{
    calstat.tab <- calstat.tab.fr
    calplots.list <- calplots.list.fr
  }
  
  return(list(caltab=calstat.tab, calplot=calplots.list.fr))
}




### Kalibraatiotaulun siistimistä -------------

format3d <- function(x,d){
  formatted <- format(round(x,d),nsmall=d, trim=T, scientific=F)
  if(any(formatted < 0)){
    ifelse(round(x,d) >= 0, paste0(" ", formatted), formatted)
  }else{
    formatted
  }
}


format_ciprint <- function(x,lower,upper,d=2){
  paste0("[", format3d(lower,d), ", ", format3d(upper,d), "]")
}

kokoa_calstats_tagged <- function(calstat.tab, all_tags) {

  # c66results_actual_contents <- ifelse(any(grepl("FRcomb", c66_results$calstat.tab$tag)), "FRcomb","C66")
  
  best_mod_options <- all_tags %>%
    # gsub("C66", c66results_actual_contents,. ) %>% #temp, käytössä kun ajettiin frcomb c66 setin tilalla.
    str_split("_", simplify = TRUE) %>%
    as_tibble(.name_repair = ~ c("FR", "Followup", "kiinteet", "modeltype")) %>%
    mutate(
      Followup = str_remove(Followup, "FU"),
      kiinteet = case_when(
        kiinteet == "fixnone" ~ "none",
        kiinteet == "fixclinical" ~ "clinical",
        TRUE ~ kiinteet)
      ) %>%
    # add_column(tag = all_tags%>%gsub("C66", c66results_actual_contents,. )) %>%
    add_column(tag = all_tags) %>%
    arrange(str_sort(tag, numeric=T))
  
  
  calstat.tab.cleaned <- calstat.tab %>%
    mutate(across(-1, as.numeric)) %>%
    mutate(
      index=1:nrow(.),
      Intercept = format_ciprint(intercept.Point.estimate, intercept.Lower.confidence.limit, intercept.Upper.confidence.limit),
      Slope = format_ciprint(slope.Point.estimate, slope.Lower.confidence.limit.2.5.., slope.Upper.confidence.limit.97.5..),
      Brier = Brier,
      Brier_scaled = Brier.scaled
    ) %>%
    left_join(best_mod_options, ., by="tag") %>%
    arrange(factor(Followup, levels = c("5", "8", "10")), factor(FR, levels = c("FR02", "FR97", "FR07", "FR12", "C66"))) %>%
    select(Followup, FR, intercept.Point.estimate, Intercept, slope.Point.estimate, Slope, Brier, Brier_scaled, index,tag) %>%
    as_tibble()
  

  return(calstat.tab.cleaned)
}



# roc.test taulut -------------

kokoa_roctest_taulu <- function(problist, c66_results, dummydata) {
  
  #vastemuuttujat
  ylist <- dummydata$ylist
  #käytettävä problist
  problist_t <- problist %>% map(., ~transpose(.x))

  big_roctest <- tibble()
  for(i.dataset in 1:length(problist)){
    dataset_probs <- problist_t[[i.dataset]]
    for(j.fuset in 1:length(dataset_probs)){
      fuset_probs <- dataset_probs[[j.fuset]]
      
      fulens <- names(dataset_probs) %>% str_split(.,"FU",simp=T)%>%as.data.frame%>%pull(2)
      fulen <- fulens[j.fuset]
      ylist_colnum  <- match(fulen, c("5","8","10"))
      
      response <- ylist[[i.dataset]][,ylist_colnum]
      
      t_ridge <- roc.test(response, predictor1=fuset_probs$ridge, predictor2=fuset_probs$logreg, quiet=T, method="delong", alternative="two.sided")
      t_ela05 <- roc.test(response, predictor1=fuset_probs$ela05, predictor2=fuset_probs$logreg, quiet=T, method="delong", alternative="two.sided")
      t_lasso <- roc.test(response, predictor1=fuset_probs$lasso, predictor2=fuset_probs$logreg, quiet=T, method="delong", alternative="two.sided")
      t_stepw <- roc.test(response, predictor1=fuset_probs$stepwise, predictor2=fuset_probs$logreg, quiet=T, method="delong", alternative="two.sided")
      
      
      roctest_list <- lapply(list(t_ridge, t_ela05, t_lasso, t_stepw), FUN=function(tlist){
        diff <- sum(tlist$estimate*c(1,-1))
        auc_new <- tlist$estimate[1]
        auc_old <- tlist$estimate[2]
        zstat <- tlist$statistic
        pval <- tlist$p.value
        ci_upper <- max(tlist$conf.int)
        ci_lower <- min(tlist$conf.int)
        ttab <- tibble(diff, auc_new, auc_old, zstat, pval, ci_upper, ci_lower)
      }) %>%
        setNames(c("ridge","ela05","lasso","stepwise")) %>%
        bind_rows(.id="model") %>%
        mutate(dataset= names(problist_t)[i.dataset],
               fuset  = names(dataset_probs)[j.fuset])
      
      big_roctest <- bind_rows(big_roctest, roctest_list)
    }
  }
  
  big_roctest <- big_roctest %>%
    mutate(tagpart1 =   gsub("fix_(.*)_FU(\\d+)", "FU\\2_fix\\1_", fuset),
           tag = paste0(tagpart1, model),
           prettyCI = glue("[{sprintf('%.3f', ci_lower)}, {sprintf('%.3f', ci_upper)}]"),
           prettyDiff = glue("{sprintf('%.3f', diff)}"),
           prettyVal = glue("{prettyDiff} {prettyCI}"))
  
  if(!is.null(c66_results)){
    c66_roctest <- c66_results$roctest %>%
      mutate(dataset="c66",
             model = str_split_i(tag, "_", 3),
             prettyCI = glue("[{sprintf('%.3f', ci_lower)}, {sprintf('%.3f', ci_upper)}]"),
             prettyDiff = glue("{sprintf('%.3f', diff)}"),
             prettyVal = glue("{prettyDiff} {prettyCI}"),
             fuset = str_split_i(tag, "_",1),
             tagpart1 = str_split_i(tag, "_",2),
             tagpart1 = paste0(fuset, "_", tagpart1,"_"))
    big_roctest <- bind_rows(big_roctest, c66_roctest)
  }
  
  return(big_roctest)
}


# betataulut -----

setup_easy_betatab <- function(x){
  x %>% 
    do.call(cbind,.) %>% 
    setNames(paste0("tun",1:length(x))) %>% 
    mutate(var=rownames(.)) %>% 
    select(var,any_of(c("tun1","tun2","tun3","tun4","tun5","tun6")))
}


extract_betas <- function(method, reslist) {
  map(reslist[[1]], ~ .x[[method]][[4]] %>% select(beta))
}


kokoa_betalist <- function(reslist) {
  modnames <- c("ridge","ela05","lasso","logreg","stepwise")
  betalist_names <- c("var", "none_FU5", "clinical_FU5", "none_FU8", "clinical_FU8")
  
  beta_tables <- map(modnames, extract_betas, reslist=reslist)
  processed_betas <- map(beta_tables, setup_easy_betatab)
  betalist <- map(processed_betas, ~ .x %>% setNames(betalist_names))
  names(betalist) <- paste0(modnames, "_betas")
  
  fi_names <- c("IKA", "SUKUP2", "ALKI1", "TUPAKOINTI1", "HDL", "KOL", 
                "TRIG", "SYSm", "DIASm", "VYOTARO", "BMI", "GLC", "BP_TREAT1", 
                "LIPID_TREAT1", "DIAB_FAMILYHIST1")
  
  en_names <- c("Age","Female sex","Alcohol","Current smoking","HDL-C","Total cholesterol",
                "Triglycerides", "Systolic BP","Diastolic BP","Waist circumference","BMI","Glucose","BP medication",
                "Lipid medication","Diabetes family history")
  
  
  betalist$ridge_betas$var <- recode(betalist$ridge_betas$var, !!!setNames(en_names, fi_names))
  betalist$ela05_betas$var <- recode(betalist$ela05_betas$var, !!!setNames(en_names, fi_names))
  betalist$lasso_betas$var <- recode(betalist$lasso_betas$var, !!!setNames(en_names, fi_names))
  betalist$logreg_betas$var <- recode(betalist$logreg_betas$var, !!!setNames(en_names, fi_names))
  betalist$stepwise_betas$var <- recode(betalist$stepwise_betas$var, !!!setNames(en_names, fi_names))
  
  
  return(betalist)
}


make_sparse_betatab <- function(betalist_item){
  betalist_item %>%
    filter(rowSums(select(., 2:7) != 0) != 0) %>%
    tibble::remove_rownames() %>%
    mutate_at(2:7, ~round(as.numeric(.), 4)) %>%
    mutate_all(as.character) %>%
    mutate_all(~ifelse(. == "0", "", .))
}


kokoa_valmis_betaval_taulu <- function(betalist, scale_params) {
  
  ridge.beta <- betalist[[1]]
  ela05.beta <- betalist[[2]]
  lasso.beta <-  betalist[[3]]
  logr.beta <- betalist[[4]]
  stp.beta <- betalist[[5]]
  
  scaleinfo <- data.frame(var=names(scale_params$scale), 
                          center=format(scale_params$center, digits=3,scientific=T), 
                          scale=format(scale_params$scale,digits=3,scientific=T))
  
  trueorder <- ridge.beta$var
  
  stepwise.betatab.f <- stp.beta %>% make_sparse_betatab %>% right_join(scaleinfo, ., by="var") %>% .[order(match(.$var, trueorder)),] %>% replace_na(list(center="",scale=""))
  lasso.betatab.f <- lasso.beta %>% make_sparse_betatab %>% right_join(scaleinfo, ., by="var") %>% .[order(match(.$var, trueorder)),]%>% replace_na(list(center="",scale=""))
  ela05.betatab.f <- ela05.beta %>% make_sparse_betatab %>% right_join(scaleinfo, ., by="var") %>% .[order(match(.$var, trueorder)),]%>% replace_na(list(center="",scale=""))
  ridge.betatab.f <- ridge.beta %>% make_sparse_betatab %>% right_join(scaleinfo, ., by="var") %>% .[order(match(.$var, trueorder)),]%>% replace_na(list(center="",scale=""))
  logreg.betatab.f <- logr.beta %>% make_sparse_betatab %>% right_join(scaleinfo, ., by="var") %>% .[order(match(.$var, trueorder)),]%>% replace_na(list(center="",scale=""))
  
  return(list(stepwise.betatab.f=stepwise.betatab.f, 
              lasso.betatab.f=lasso.betatab.f, 
              ela05.betatb.f=ela05.betatab.f, 
              ridge.betatab.f=ridge.betatab.f, 
              logreg.betatab.f=logreg.betatab.f))
  
}


kokoa_betavals_by_fu <- function(betalist, clinical_vars) {
  
  ridge.beta <- betalist[[1]]
  ela05.beta <- betalist[[2]]
  lasso.beta <-  betalist[[3]]
  logr.beta <- betalist[[4]]
  stp.beta <- betalist[[5]]
  
  
  getj <- function(x,j) x%>%select(all_of(c(1,j)))
  fuset_betat <- list()
  for(j in 2:ncol(ridge.beta)){
    setti <- getj(ridge.beta,j) %>%
      left_join(., getj(ela05.beta,j), by="var") %>%
      left_join(.,getj(lasso.beta,j), by="var") %>%
      left_join(., getj(stp.beta,j),by="var") %>%
      left_join(., getj(logr.beta,j), by="var") %>%
      setNames(c("var","ridge","ela05","lasso","step","logreg")) %>%
      mutate_all(~ ifelse(is.na(.), "", ifelse(. == 0, "", as.character(.)))) %>%
      mutate_at(-c(1), ~ ifelse(. == "", "", format(round(as.numeric(.), digits=3),nsmall=3)))
    
    fixed_vars <- setdiff(clinical_vars, c("TUPAKOINTI"))
    non_fixed <- setdiff(setti$var,fixed_vars)
    setti <- setti[match(intersect(c(fixed_vars, non_fixed), setti$var), setti$var), ,drop=F]
    
    fuset_betat[[j-1]]<-setti
    
  }
  
  names(fuset_betat) <- c("fu5_fulldata", "fu5_fixdata", "fu8_fulldata", "fu8_fixdata")
  betatfu5 <- left_join(fuset_betat$fu5_fixdata, fuset_betat$fu5_fulldata, by="var", suffix=c("_fix","_full"))
  betatfu8 <- left_join(fuset_betat$fu8_fixdata, fuset_betat$fu8_fulldata, by="var", suffix=c("_fix","_full"))
  
  fu_betat <- list(betatfu5, betatfu8)
  names(fu_betat) <- c("fu5mods", "fu8mods")
  fu_betat <- lapply(fu_betat, FUN=function(x){
    x %>% 
      select(var, logreg_fix, everything()) %>%
      select(-logreg_full) %>%
      rename(logreg=logreg_fix)
  })

  return(fu_betat)
}


# Beta plotit ----

make_beta_lineplot <- function(btab,title){
  beta_order <- setdiff(btab$var,"(Intercept)")
  btab_long <- btab %>%
    filter(var != "(Intercept)") %>%
    mutate(var = factor(var, levels=beta_order)) %>%
    pivot_longer(.,cols=-1) %>%
    mutate(value = as.numeric(value))
  
  betaplot <- ggplot(btab_long, aes(x=var, y=value, color=name, group=name)) +
    geom_line(linewidth=1, alpha=0.5) + 
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=5), axis.title.x = element_blank()) +
    labs(title = title)  
  
  return(betaplot)
}

kokoa_valmiit_betaval_kuvat <- function(betalist) {
  ridge <- betalist[[1]]
  ela05 <- betalist[[2]]
  lasso <-  betalist[[3]]
  logr.beta <- betalist[[4]]
  stp.beta <- betalist[[5]]
  
  return(list(fulldat.ridgeplot = make_beta_lineplot(ridge[,c(1,grep("none",colnames(ridge)))], "Ridge, fulldat") ,
              fulldat.ela05plot = make_beta_lineplot(ela05[,c(1,grep("none",colnames(ela05)))], "Ela0.5, fulldat") ,
              fulldat.lassoplot = make_beta_lineplot(lasso[,c(1,grep("none",colnames(lasso)))], "LASSO, fulldat") ,
              fulldat.stepwplot = make_beta_lineplot(stp.beta[,c(1,grep("none",colnames(stp.beta)))], "Stepwise, fulldat"),
              fulldat.logreplot = make_beta_lineplot(logr.beta[,c(1,grep("none",colnames(logr.beta)))], "Logreg, fulldat"),
              fixdat.ridgeplot = make_beta_lineplot(ridge[,c(1,grep("clinical",colnames(ridge)))], "Ridge, fixed clinicals"),
              fixdat.ela05plot = make_beta_lineplot(ela05[,c(1,grep("clinical",colnames(ela05)))], "Ela05, fixed clinica"),
              fixdat.lassoplot = make_beta_lineplot(lasso[,c(1,grep("clinical",colnames(lasso)))], "LASSO, fixed clinicals"),
              fixdat.stepwplot = make_beta_lineplot(stp.beta[,c(1,grep("clinical",colnames(stp.beta)))], "Stepwise, fixed clinicals"),
              fixdat.logreplot = make_beta_lineplot(logr.beta[,c(1,grep("none",colnames(logr.beta)))], "Logreg, fixed clinicals")
  ))
}

# Kaplan-meierit ----
fix_breaks <- function(plot.obj) {
  orig_breaks <- ggplot_build(plot.obj)[["layout"]][["panel_params"]][[1]][["y"]][["breaks"]]
  orig_minorbreaks <- ggplot_build(plot.obj)[["layout"]][["panel_params"]][[1]][["y"]][["minor_breaks"]]
  
  new_breaks <- orig_breaks
  new_breaks[which.min(new_breaks)] <- min(plot.obj$data$estimate)
  new_breaks <- round(new_breaks, 3)
  
  new_minorbreaks <- orig_minorbreaks
  new_minorbreaks <- new_minorbreaks[-which(new_minorbreaks == min(orig_breaks))] 
  
  plot.obj <- plot.obj + scale_y_continuous(breaks = new_breaks, minor_breaks = new_minorbreaks)
  return(plot.obj)
}

make_top_bot_x_percent_km_plots <- function(suri, xperc,title, force_y_axis_range = NULL) {
  #kaplan-meier korkeimmista 10% probseista
  su.topx <- suri %>% 
    mutate(topxprobs = ifelse(suri$probs > quantile(suri$probs, 1-xperc), glue("top {xperc*100}%"), glue("bottom {(1-xperc)*100}%"))) %>%
    survfit2(Surv(time_to_event, E4_DIABETES)~topxprobs, data=.) %>% 
    ggsurvfit(type="risk") +
    add_confidence_interval() +
    labs(x = "years",y = "Risk of event") + 
    ggtitle(title) +
    scale_color_manual(values = c("tomato3", "steelblue")) + # Change line colors
    scale_fill_manual(values = c("tomato3", "steelblue"))  # Change confidence interval band colors
  
  # su.topx <- fix_breaks(su.topx)
  
  if(!is.null(force_y_axis_range) && force_y_axis_range) su.topx <- su.topx + ylim(0, 1)
  
  
  return(su.topx)
}

build_caltab_yhteensopiva_nimi <- function(ap, modname, vecname){
  fulen <- gsub("[^0-9]","", vecname)
  fixstatus <- str_split(vecname, "_", simplify=T)[,1:2] %>% paste0(.,collapse="")
  modname <-   modname %>% 
    gsub("stepwise","stepwise_",.) %>%
    gsub("logreg", "logreg_",.)
  if(!grepl("66",ap)) title <- glue("FR{ap}_FU{fulen}_{fixstatus}_{modname}")
  if(grepl("66",ap)) title <- glue("{ap}_FU{fulen}_{fixstatus}_{modname}")
  
  return(title)
}



kokoa_vuoden_km_plotit <- function(problist_yr,  yearname, survdat, followups=c(5,8,10), qgrid=F) {
  kmplotlist <- list()
  for(mod in 1:length(problist_yr)){
    modprobs <- problist_yr[[mod]] #malli-kohtaiset probsit
    modprobs <- modprobs[,grepl(paste0(followups, collapse="|"), names(modprobs))]
    modname <- names(problist_yr)[mod]
    for(st in 1:ncol(modprobs)){
      probvec <- modprobs[,st]
      vecname <- colnames(modprobs)[st]
      title <- build_caltab_yhteensopiva_nimi(yearname, modname, vecname)
      suriprob <-  cbind(survdat,probs=probvec) %>% tibble
      
      if(any(is.nan(suriprob$probs))){ #jos suri probsit on pilalla (fu8 stepwise)
        kmplotlist[[title]] <- ggplot() + ggtitle(title)
      }else{
        # kmplotlist_freeylim[[title]] <- make_top_bot_x_percent_km_plots(suriprob, 0.1, title=title)
        if(!qgrid) kmplotlist[[title]] <- make_top_bot_x_percent_km_plots(suriprob, 0.1, title=title, force_y_axis_range = T)
        if(qgrid){
          qgrid <- seq(0.05,1, by=0.05)
          sublist <- lapply(X=qgrid, FUN=function(X) list(quantile=q, make_top_bot_x_percent_km_plots(suriprob, X, title=title, force_y_axis_range = T)))
          
        }
      }
    }
  }
  return(kmplotlist)
}


get_survival_data <- function(yri,timeskip=0, censor_glucose=T) {
  fulldat <- setup_xy_DIAB(yr=yri, followup_length=5, dropID = F, fu_min=timeskip, censor_glucose = censor_glucose)$fullset 
  #kun dropID on FALSE, niin tehdään aina FU50, riippumatta siitä mitä tässä annetaan followup_lengthiksi.
  survdat  <- fulldat %>% 
    select(all_of(c("studyID", "BL_AGE", "FU_END_AGE", "E4_DIABETES_FU50", "E4_DIABETES_AGE"))) %>%
    mutate(time_to_event = E4_DIABETES_AGE - BL_AGE) %>%
    select(all_of(c("E4_DIABETES_FU50", "time_to_event")))
  names(survdat)<- gsub('E4_DIABETES_FU50', 'E4_DIABETES', names(survdat))
  return(survdat)
}


#datan kokoaminen
kokoa_valmiit_KM_kuvat <- function(iditablist, c66_results, followups, timeskip=0, censor_glucose=T) {
  problist <- iditablist$problists
  kmplotlist <- list()
  years <- c(02,97,07,12)
  for(yr in 1:length(years)){
    cat("\n\tKootaan vuotta", years[yr])
    survdat <- get_survival_data(years[yr],timeskip=timeskip[[1]], censor_glucose=censor_glucose)    
    problist_yr <- problist[[yr]]
    yearname <- years[yr]
    
    vuoden_kmplotit <- kokoa_vuoden_km_plotit(problist_yr, yearname, survdat, followups=followups, qgrid=F)
    kmplotlist <- c(kmplotlist, vuoden_kmplotit)
  }
  
  if(!is.null(c66_results)) kmplotlist <- c(kmplotlist, c66_results$kmplots66)
  
  cat("\n")
  return(kmplotlist)
}




# Tulosten yleistä siistimistä ----



paste_mean_sd <- function(x){
  means <- colMeans(x)
  sds <- apply(x, 2, sd)
  
  formated_mean <- format(means,digits=3)
  formated_sd <- format(sds,digits=3)
  formated <- glue("{formated_mean} ({formated_sd})")
  names(formated) <- names(x)
  formated
}

build_meanSD_tab <- function(dummydata, c66meansd) {
  datameansd <- 
    map(dummydata$datasets[c(1,4,5,6)], paste_mean_sd) %>% 
    bind_rows %>% 
    t %>% 
    data.frame %>% 
    mutate(var=row.names(.)) %>% 
    setnames(c("fr02","fr97","fr07","fr12","var")) %>% 
    select(var, fr97, fr07, fr02, fr12) %>% 
    as_tibble %>%
    left_join(c66meansd, by="var")
  return(datameansd)
}

#### siistitään tallennettavia taulukoita ----
beautify_clinvartab <- function(fulltbl1) {
  tab1.all <- fulltbl1
  tab1.fu5 <- tab1.all %>% remove_row_type(type="all", variables=c("FU8_cases")) %>% tbl_butcher
  tab1.fu810 <- tab1.all %>% remove_row_type(type="all", variables=!c("FU8_cases")) %>% tbl_butcher
  tab1.all <- tab1.all %>% tbl_butcher
  
  return(list(all=tab1.all,fu5=tab1.fu5,fu810=tab1.fu810))
}

beautify_auctabs <- function(auctablist) {
  widetab <- auctablist %>% 
    pivot_wider(id_cols=c(model, followup, model_vars), 
                values_from=c(auc,auc_ci,beta),
                names_from=dataset,
                names_vary = "slowest") %>%
    select(-contains("beta_FR")) %>%
    mutate(
      model_vars=case_when(
        model_vars=="fixnone" ~ "full",
        model_vars=="fixclinical" ~ "fix",
        model_vars=="onlyclinicals" ~ "clin"
      ),
      latexmod = glue("$\\textrm{<str_to_title(model)>}^{<model_vars>}$", .open="<", .close=">"),
      sortmodel = factor(model, levels=c("ridge","ela05","lasso","stepwise","logreg")),
      sortfix = factor(model_vars, levels=c("full","fix","clin"))) %>%
    arrange(followup, sortfix, sortmodel) %>%
    select(followup, latexmod, contains("auc"), contains("beta")) %>%
    rename(beta=beta_C66)
  
  wide5 <- widetab %>% filter(followup==5)
  
  
  return(list(all=widetab, fu5=wide5))
}


beautify_iditabs <- function(IDI_caltab) {
  tab3.all <- IDI_caltab
  tab3.fu5 <- tab3.all[trimws(as.character(tab3.all$Followup)) == "5", ] %>% select( -any_of(c("Brier", "Brier_scaled")))
  tab3.fu810 <- tab3.all
  return(list(fu5=tab3.fu5, fu810=tab3.fu810))
}



beautify_caltab <- function(calstats_allmods){
  prettytabs <- calstats_allmods %>% 
    mutate(tag2=tag) %>%
    separate_wider_delim(cols=tag2, delim="_", names=c("fr","fu","setup","Model")) %>%
    arrange(
      match(Followup, c("FU5","FU8","FU10")),
      match(setup, c("fixclinical","fixnone","onlyclinical")),
      match(Model, c("ridge","ela05","lasso","stepwise","logreg")),
      match(FR, c("FR02","FR97","FR07","FR12","C6646"))
      ) %>%
    select(-c(fr,fu)) %>%
    mutate(Model = str_to_title(Model),
           FR = toupper(FR),
           Followup = gsub("FU","",Followup),
           datasetup = gsub("fix","", setup) %>% gsub("none","full",.) %>% gsub("clinical","fix",.),
           LatexModel = glue("$\\textrm{[Model]}^{[datasetup]}$",.open="[",.close="]")) %>%
    select(Followup, LatexModel, FR, intercept.Point.estimate, Intercept, slope.Point.estimate, Slope, Brier, Brier_scaled, tag) %>%
    filter(Followup != 10) %>%
    as_tibble %>%
    rename(int=intercept.Point.estimate,
           int_ci = Intercept,
           slope=slope.Point.estimate,
           slope_ci=Slope)
  
  fu5tabs = prettytabs%>%filter(Followup==5)
  return(list(fu5=fu5tabs, fuall=prettytabs))
}
  
  

beautify_roctesttab <- function(roctest_tabs){
  prettytabs <- roctest_tabs %>% 
    mutate(Model = str_to_title(model),
           FR = toupper(dataset),
           Followup = str_split(fuset, "FU",simplify=T)[,2],
           datasetup = str_split(tagpart1, "_",simplify=T)[,2],
           datasetup = gsub("fix","", datasetup) %>% gsub("none","full",.) %>% gsub("clinical","fix",.),
           LatexModel = glue("$\\textrm{[Model]}^{[datasetup]}$",.open="[",.close="]")) %>%
    select(Followup, LatexModel, FR, prettyDiff, prettyCI, tag) %>%
    setNames(c("Followup", "Model", "Set", "dAUC", "CI", "tag"))
  
  
  
  mainpap_tabs <- prettytabs %>% filter(Followup == 5)
  extra_tabs <- prettytabs %>% filter(Followup != 5)
  
  return(list(fu5=mainpap_tabs, fu810=extra_tabs))
  
}

#### siistitään tallennettavia plotteja ----
overlay_clin_plots <- function(plot_list_best, plot_list_clin, 
                               fill_color_best = "orange", line_color_best = "orange", 
                               fill_color_clin = "steelblue", line_color_clin = "steelblue"){
  double_plot_list <- list()
  for (i in seq_along(plot_list_best)) {
    #kootaan uusi plotti. poistetaan 
    newplot <- plot_list_best[[i]]
    newplot[["layers"]][[2]][["aes_params"]][["fill"]] <- NA
    newplot[["layers"]][[3]][["aes_params"]][["colour"]] <- NA
    
    modname_best <- plot_list_best[[i]][["labels"]][["title"]]
    modname_clin <- plot_list_clin[[i]][["labels"]][["title"]]
    datasetname <- str_split_fixed(modname_best, "_", 2)[1]
    
    doubleplot <- newplot + 
      # Best model
      geom_line(data = plot_list_best[[i]]$data, aes(x = x, y = y), color = line_color_best) + 
      geom_ribbon(data = plot_list_best[[i]]$data, aes(ymin = ymin, ymax = ymax, fill = "best_model"), alpha = 0.2) +
      # Clinical model
      geom_line(data = plot_list_clin[[i]]$data, aes(x = x, y = y), color = line_color_clin) + 
      geom_ribbon(data = plot_list_clin[[i]]$data, aes(ymin = ymin, ymax = ymax, fill = "clinical_model"), alpha = 0.2) +
      #legend building
      scale_fill_manual(values = c("best_model" = fill_color_best, "clinical_model" = fill_color_clin), labels = c(modname_best, modname_clin)) +
      labs(fill = "", title=datasetname)+
      theme(legend.position = "top", legend.justification = "left")
    
    
    double_plot_list[[i]] <- doubleplot
  }
  return(double_plot_list)
}



clip_plot <- function(plot,xlim=5,ylim=0.3) {
  suppressMessages(plot + scale_x_continuous(limits = c(0, xlim)) + scale_y_continuous(limits = c(0, ylim)))
}


get_best_calib_titles <- function(calibrations, followup="FU5", dropLogreg=T) {
  plotnames <- calibrations[["caltab"]][["tag"]]
  
  if(dropLogreg) chosen_indices <- which(str_detect(plotnames, followup) & !str_detect(plotnames, "logreg")) #näiden järjestys on väärä
  if(!dropLogreg) chosen_indices <- which(str_detect(plotnames, followup) & str_detect(plotnames, "logreg")) 
  
  chosen_fus <- plotnames[chosen_indices] %>% str_split(., "_", simp=T) %>% gsub("FR","",.) %>% data.frame %>% pull(1)
  fu_order <- c("02","97","07","12","C66")
  
  ordered_indices<-chosen_indices[match(fu_order,chosen_fus)]
  
  return(ordered_indices)
}

get_bestplots <- function(calibrations, best_calib_titles) {
  calibrations$calplot[best_calib_titles[c(1,2,3,4,5)]]
}


get_bestmods_kmplots <- function(KM_plots, bestFUtitles) {
  bestmods_kmplots <- KM_plots[which(
    str_replace(names(KM_plots), "fixclinical_logreg_", "onlyclinicals_logreg") %in% bestFUtitles
  )]
  return(bestmods_kmplots)
}


wrap_clipped_kmplots_byFU <- function(best_tags, KM_plots, FU, xlim,ylim,logregs=F) {
  bestfutitles <- grep(FU, best_tags, value = TRUE) %>% 
    gsub("02", "2", .) %>% 
    gsub("07", "7", .) %>%
    gsub("stepwise", "stepwise_",.)
  
  if(!logregs) bestfutitles <- bestfutitles %>% discard(grepl("logreg",.))
  if(logregs) bestfutitles <- bestfutitles %>% keep(grepl("logreg",.))
  
  best_fumods_kmplots <- get_bestmods_kmplots(KM_plots, bestfutitles)
  bestkm_limrange <- best_fumods_kmplots %>% map(~clip_plot(.,xlim,ylim)) %>% wrap_plots
  bestkm_limrange
}

prep_kmplots <- function(fu_values, limits, best_tags_clins, KM_plots) {
  best_km <- map2(fu_values, limits, ~wrap_clipped_kmplots_byFU(best_tags_clins, KM_plots, glue("FU{.x}"), .x, .y))
  logreg_km <- map2(fu_values, limits, ~wrap_clipped_kmplots_byFU(best_tags_clins, KM_plots, glue("FU{.x}"), .x, .y, logregs = TRUE))
  return(list(best=best_km, logs=logreg_km))
}

# prep_calibplots <- function(fu_values, calibrations) {
#   bestplots <- map(fu_values, ~get_best_calib_titles(calibrations, glue("FU{.x}")) %>% get_bestplots(calibrations, .) %>% discard(is.null))
#   calplot_logreg <- map(fu_values, ~get_best_calib_titles(calibrations, glue("FU{.x}"), dropLogreg = F) %>% get_bestplots(calibrations, .) %>% discard(is.null))
#   doubleplots <- map2(bestplots, calplot_logreg, overlay_clin_plots)
#   return(list(best=bestplots, dups=doubleplots))
# }

#### siistitään kalibraatio plotteja ----

combine_calibration_curves <- function(plot_list_best, plot_list_clin, legend_labs, ylab="Observed Proportion", line_color_best, line_color_clin, fill_color_best, fill_color_clin) {
  plotlist <- list()
  for (i in seq_along(plot_list_best)) {
    
    #kootaan uusi plotti. vaihdetaan "parhaaseen" värit ja laitetaan "kliinisen" käyrä siihen päälle 
    plotname <- str_split_i(plot_list_best[[i]]$labels$title, "_",1)
    plotname <- gsub("FR", "FINRISK ", plotname)
    plotname <- gsub("C66", "NFBC1966", plotname)
    plotname <- gsub("FINRISK 02", "FINRISK 2002", plotname)
    plotname <- gsub("FINRISK 07", "FINRISK 2007", plotname)
    plotname <- gsub("FINRISK 12", "FINRISK 2012", plotname)
    plotname <- gsub("FINRISK 97", "FINRISK 1997", plotname)
    
    best_sans_rug <- plot_list_best[[i]]
    best_sans_rug$layers[2:8]<-NULL #rugplot pois
    best_sans_rug$layers[[1]]$aes_params$size <- 1.4
    best_sans_rug$layers[[1]]$aes_params$alpha <- 0.2
    best_sans_rug$layers[[1]]$aes_params$colour <- "grey10"
    best_sans_rug[["labels"]][["y"]] <- ylab
    
    doubleplot <- best_sans_rug +
      geom_line(data = plot_list_clin[[i]]$data, aes(x = x, y = y), color = line_color_clin, size=2, alpha=0.6) +
      geom_ribbon(data = plot_list_clin[[i]]$data, aes(ymin = ymin, ymax = ymax, fill = "clinical_model"), alpha = 0.4) +
      geom_line(data = plot_list_best[[i]]$data, aes(x = x, y = y), color = line_color_best, size=2, alpha=0.6) +
      geom_ribbon(data = plot_list_best[[i]]$data, aes(ymin = ymin, ymax = ymax, fill = "best_model"), alpha = 0.4) +
      scale_fill_manual(
        values = c("best_model" = fill_color_best, "clinical_model" = fill_color_clin),
        labels = legend_labs
      )+
      coord_cartesian(ylim=c(0,1)) +
      labs(fill = "", title=plotname)+
      theme(legend.position = "top", 
            legend.justification = "left", 
            legend.text = element_markdown(size = 12),  
            axis.title.y = element_markdown(size=12),
            axis.title.x = element_markdown(size=12),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.margin=margin(r=-8),
            legend.margin=margin(b=-10),
            plot.title = element_text(hjust=0.5))
    
    plotlist[[i]] <- doubleplot
  }
  return(plotlist)
}


prep_calibplots <- function(fu_values, calibrations) {
  
  #plotit 
  calplots <- calibrations$calplot %>% compact %>% setNames(map(.,"labels") %>% map("title") %>% unlist)
  #fu8
  logreg8_plots <- names(calplots) %>% grep("FU8.*logreg",.,value=T)
  ridgefix8_plots <- names(calplots) %>% grep("FU8.*fixclinical.ridge",.,value=T) 
  ridgefull8_plots <- names(calplots) %>% grep("FU8.*fixnone.ridge",.,value=T)
  #fu5
  logreg5_plots <- names(calplots) %>% grep("FU5.*logreg",.,value=T)
  ridgefix5_plots <- names(calplots) %>% grep("FU5.*fixclinical.ridge",.,value=T) 
  ridgefull5_plots <- names(calplots) %>% grep("FU5.*fixnone.ridge",.,value=T)
  
  #värit
  fix8col  = "#4895ef"
  full8col = "#574c89"
  fix5col  = "#e56b6f"
  full5col = "#fb8b24"
  clincol = "#adb5bd"
  

  
  double_fix_ridge_logreg8 <- suppressMessages(
    combine_calibration_curves(plot_list_best = calplots[ridgefix8_plots], 
                               plot_list_clin=calplots[logreg8_plots], 
                               legend_labs = c(
                                 "<span style='font-size:16px;'>Partially penalized ridge regression, 8-year follow-up</span>",
                                 "<span style='font-size:16px;'>Baseline, 8-year follow-up</span>"
                               ),
                               line_color_best=fix8col, line_color_clin=clincol, 
                               fill_color_best=fix8col, fill_color_clin=clincol)
  )
  
  double_full_ridge_logreg8 <- suppressMessages(
    combine_calibration_curves(plot_list_best = calplots[ridgefull8_plots], 
                               plot_list_clin=calplots[logreg8_plots], 
                               legend_labs = c(
                                 "<span style='font-size:16px;'>Fully penalized ridge regression, 8-year follow-up</span>",
                                 "<span style='font-size:16px;'>Baseline, 8-year follow-up</span>"
                               ),  
                               line_color_best=full8col, line_color_clin=clincol, 
                               fill_color_best=full8col, fill_color_clin=clincol)
  )
  
  double_fix_ridge_logreg5 <- suppressMessages(
    combine_calibration_curves(plot_list_best = calplots[ridgefix5_plots], 
                               plot_list_clin=calplots[logreg5_plots], 
                               legend_labs = c(
                                 "<span style='font-size:16px;'>Partially penalized ridge regression, 5-year follow-up</span>",
                                 "<span style='font-size:16px;'>Baseline, 5-year follow-up</span>"
                               ),
                               line_color_best=fix5col, line_color_clin=clincol, 
                               fill_color_best=fix5col, fill_color_clin=clincol)
  )
  
  double_full_ridge_logreg5 <- suppressMessages(
    combine_calibration_curves(plot_list_best = calplots[ridgefull5_plots], 
                               plot_list_clin=calplots[logreg5_plots], 
                               legend_labs = c(
                                 "<span style='font-size:16px;'>Fully penalized ridge regression, 5-year follow-up</span>",
                                 "<span style='font-size:16px;'>Baseline, 5-year follow-up</span>"
                               ),
                               line_color_best=full5col, line_color_clin=clincol, 
                               fill_color_best=full5col, fill_color_clin=clincol)
  )
  
  return(list(
    ridgefull5=double_full_ridge_logreg5,
    ridgefix5=double_fix_ridge_logreg5,
    ridgefull8=double_full_ridge_logreg8, 
    ridgefix8=double_fix_ridge_logreg8
  ))
}



#### pudotetaan plottilistoilta tavaraa ----


drop_legend_from_rest <- function(plot_list) {
  plot_list[-1] <- map(plot_list[-1], ~ .x + theme(legend.position = "none"))
  return(plot_list)
}

drop_labs_x <- function(plot_list, keep) {
  plot_list <- map(seq_along(plot_list), function(i) {
    if (i %in% keep) {
      plot_list[[i]]
    } else {
      plot_list[[i]] + labs(x=NULL)
    }
  })
  return(plot_list)
}

drop_labs_y <- function(plot_list, keep) {
  plot_list <- map(seq_along(plot_list), function(i) {
    if (i %in% keep) {
      plot_list[[i]]
    } else {
      plot_list[[i]] + labs(y = NULL)
    }
  })
  return(plot_list)
}

drop_legends <- function(plot_list, keep) {
  plot_list <- map(seq_along(plot_list), function(i) {
    if (i %in% keep) {
      plot_list[[i]]
    } else {
      plot_list[[i]] + theme(legend.position = "none")
    }
  })
  return(plot_list)
}

drop_titles <- function(plot_list, keep) {
  plot_list <- map(seq_along(plot_list), function(i) {
    if (i %in% keep) {
      plot_list[[i]]
    } else {
      plot_list[[i]] + labs(title = NULL)
    }
  })
  return(plot_list)
}


drop_axis_text_y <- function(plot_list, keep) {
  plot_list <- map(seq_along(plot_list), function(i) {
    if (i %in% keep) {
      plot_list[[i]]
    } else {
      plot_list[[i]] + theme(axis.text.y = element_blank())
    }
  })
  return(plot_list)
}

drop_axis_text_x <- function(plot_list, keep) {
  plot_list <- map(seq_along(plot_list), function(i) {
    if (i %in% keep) {
      plot_list[[i]]
    } else {
      plot_list[[i]] + theme(axis.text.x = element_blank())
    }
  })
  return(plot_list)
}


#### patchwork layoutteja ----


layout_5by1plots <- function(p){
  layout5 <- c("
            aa
            bb
            cc
            dd
            ee
            ")
  wrap_plots(p) + plot_layout(design=layout5)
}


layout_3plus2plots <- function(p){
  layout32 <- 
  c("aabbcc
     #ddee#")
  wrap_plots(p) + plot_layout(design=layout32)
}


layout_4by5plots <- function(p, p66){
  wrap_plots(p, nrow=4, ncol=4) + p66
}


combine_4by5_calplot <- function(prepped_calplots, c66_plot_fp, output_dir, target_filename) {
  
  #kootaan finriskplot, kirjoitetaan pdfksi, luetaan pdf, jota sitten muokataan kylkeen c66 plotit kapselista
  
  #filleriä
  prepped_calplots[[1]][[5]] <- prepped_calplots[[1]][[4]]
  prepped_calplots[[2]][[5]] <- prepped_calplots[[2]][[4]]
  prepped_calplots[[3]][[5]] <- prepped_calplots[[3]][[4]]
  prepped_calplots[[4]][[5]] <- prepped_calplots[[4]][[4]]
  
  #indeksit
  col1_items <- which(1:20 %% 5 == 1)
  row4_items <- tail(1:20, 5)
  row1_items <- head(1:20, 5)
  
  frplots_4by5 <- prepped_calplots %>%
    list_flatten %>% 
    drop_labs_x(keep=row4_items) %>%
    drop_labs_y(keep=col1_items) %>%
    drop_legends(keep=col1_items) %>%
    drop_axis_text_x(keep=row4_items) %>%
    drop_axis_text_y(keep=col1_items) %>%
    drop_titles(keep=row1_items) %>%
    wrap_plots(., nrow=4,ncol=5)
  
  ggsave(file.path(output_dir,"temp_doublecal_plots_4by5_finriskdata.pdf"), plot = frplots_4by5, device = "pdf", width = 300, height = 300, units = "mm")
  pdf1 <- image_read_pdf(file.path(output_dir,"temp_doublecal_plots_4by5_finriskdata.pdf"))
  cropped_plot1 <- image_crop(pdf1, geometry= paste0(2834, "x", 3542, "+0+0"))
  pdf2 <- image_read_pdf(c66_plot_fp)
  combined <- image_append(c(cropped_plot1, pdf2))
  image_write(combined, file.path(output_dir,target_filename), format="pdf") #toimii riittävän hyvin
}
