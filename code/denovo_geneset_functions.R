#################################
## De novo analysis functions  ##
#################################

## competitive tests ####
competitive_geneset_analysis <- function(new_denovos_in,
                                         published_denovos_in,
                                         genesets_in,
                                         new_denovos_mut_in,
                                         published_denovos_mut_in,
                                         output_style = 'standard') {
  
  ## Generate gene set results in new trios (Rees et al 2020)
  newTrios_setResults <-
    competitive.Obs.vs.Exp.geneset(
      denovos_in = new_denovos_in$denovos,
      N_trios = new_denovos_in$N_trios,
      mutRates_in = new_denovos_mut_in,
      genesets_in = genesets_in
    )
  
  ## Generate gene set results in previously published trios
  publishedTrios_setResults <-
    competitive.Obs.vs.Exp.geneset(
      denovos_in = published_denovos_in$denovos,
      N_trios = published_denovos_in$N_trios,
      mutRates_in = published_denovos_mut_in,
      genesets_in = genesets_in
    )
  
  ## Meta-analyse results from new and published trios
  
  ## Gather column names to merge
  vars_to_get <- c("Gene_set",
                   "Genes_tested",
                   "Remaining_genes_tested")
  ## append additional vars to get
  ## Set classes to test
  classes_to_test <- c("NS","Missense","LoF","MPC2","MPC2_LoF","NS_SNV","LoF_SNV")
  
  for (i in classes_to_test) {
    
    vars_to_get <- append(vars_to_get,
                          c(paste0(i,"_count"),
                            paste0(i,"_expected"),
                            paste0(i,"_count_remaining"),
                            paste0(i,"_expected_remaining"),
                            paste0(i,"_RR"),
                            paste0(i,"_RR_lower"),
                            paste0(i,"_RR_upper"),
                            paste0(i,"_P")))
  }
  
## Combine results from new and published trios
  combined_setResult <- merge(
    newTrios_setResults %>% 
      dplyr::select(vars_to_get) %>% 
      rename_at(vars_to_get[4:length(vars_to_get)],list(~paste0("New_Trios_",.))),
    publishedTrios_setResults %>%
      dplyr::select(vars_to_get[c(1,4:length(vars_to_get))]) %>%
      rename_at(vars_to_get[4:length(vars_to_get)],list(~paste0("Published_Trios_",.))),
    by = "Gene_set"
  )
  
  
  
  ## Generate all trios observed and expected counts
  for (class in classes_to_test) {
    
    combined_setResult[paste0("All_Trio_",class,"_obs")] <- combined_setResult[[paste0("New_Trios_",class,"_count")]] + combined_setResult[[paste0("Published_Trios_",class,"_count")]]
    combined_setResult[paste0("All_Trio_",class,"_exp")] <- combined_setResult[[paste0("New_Trios_",class,"_expected")]] + combined_setResult[[paste0("Published_Trios_",class,"_expected")]]
    combined_setResult[paste0("All_Trio_",class,"_obs_remaining")] <- combined_setResult[[paste0("New_Trios_",class,"_count_remaining")]] + combined_setResult[[paste0("Published_Trios_",class,"_count_remaining")]]
    combined_setResult[paste0("All_Trio_",class,"_exp_remaining")] <- combined_setResult[[paste0("New_Trios_",class,"_expected_remaining")]] + combined_setResult[[paste0("Published_Trios_",class,"_expected_remaining")]]
    
    ## Initiate poisson result variables
    combined_setResult[paste0("All_Trio_",class,"_P")] <- NA
    combined_setResult[paste0("All_Trio_",class,"_RR")] <- NA
    combined_setResult[paste0("All_Trio_",class,"_RR_lower")] <- NA
    combined_setResult[paste0("All_Trio_",class,"_RR_upper")] <- NA
  }
  
  
  ## Loop through each gene set and generate competitive poisson test for all trios
  for (n in 1:nrow(combined_setResult)) {
    ## generate test statistic for each mutation class
    for (class in classes_to_test) {
      ## Poisson test
      temp_poisson_result <-
        poisson.test(
          x = c(combined_setResult[[paste0("All_Trio_", class, "_obs")]][n],
                combined_setResult[[paste0("All_Trio_", class, "_obs_remaining")]][n]),
          T = c(combined_setResult[[paste0("All_Trio_", class, "_exp")]][n],
                combined_setResult[[paste0("All_Trio_", class, "_exp_remaining")]][n])
        )
    
      combined_setResult[[paste0("All_Trio_",class,"_P")]][n] <- temp_poisson_result$p.value
      combined_setResult[[paste0("All_Trio_",class,"_RR")]][n] <- temp_poisson_result$estimate
      combined_setResult[[paste0("All_Trio_",class,"_RR_lower")]][n] <- temp_poisson_result$conf.int[1]
      combined_setResult[[paste0("All_Trio_",class,"_RR_upper")]][n] <- temp_poisson_result$conf.int[2]
    }
    
  }
  
  ## Determine output style
  
  ## standard output (all variables)
  if (output_style == "standard") {
    
    output_file <- combined_setResult %>%
      select(c(1:3, contains("All_Trio"))) %>% 
      mutate_at(vars(contains(c('_RR','_P'))), .fun = list(~signif(., 3))) %>%
      rename_at(vars(contains('All_Trio')),list(~sub("All_Trio_",'',.)))
    
    return(output_file)
    
  }
  ## manuscript output (merge relative risk and 95% into single variable)
  else if (output_style == "manuscript") {
    output_file <- combined_setResult %>%
      select(c(1:3, contains("All_Trio"))) %>% 
      mutate_at(vars(contains(c('_RR','_P'))), .fun = list(~signif(., 3))) %>%
      rename_at(vars(contains('All_Trio')),list(~sub("All_Trio_",'',.)))
    
    for (i in classes_to_test) {
  
      output_file[[paste0(i,"_RR")]] <- paste0(output_file[[paste0(i,"_RR")]],
                                               " (",
                                               output_file[[paste0(i,"_RR_lower")]],
                                               ", ",
                                               output_file[[paste0(i,"_RR_upper")]],")")
      
      output_file <- output_file %>% select(-!!paste0(i,"_RR_lower"),
                                            -!!paste0(i,"_RR_upper"))
    }
    
    return(output_file)
    
  } 
  ## study breakdown output (output results separate for new, published and all trios)
  else if (output_style == "study_breakdown") {
    
    return(list(new_trio_results = newTrios_setResults,
                published_trio_results = publishedTrios_setResults,
                combined_trio_results = combined_setResult))
  }
}



##########################
## Obs.vs.Exp.geneset ####
##########################

competitive.Obs.vs.Exp.geneset <- function(denovos_in,
                                           N_trios,
                                           mutRates_in,
                                           genesets_in) {
  
  denovos_in$Confidence_2 <- "High"
  denovos_in$Annotation_2 <- denovos_in$Annotation
  ## Derive list of set names to test
  setNames <- as.list(as.character(unique(genesets_in$Gene_set)))
  
  ## Loop through gene sets and generate observed, expected and poisson statistics
  temp_results <-lapply(setNames, function(x) {
    
    ## Get set name and n genes
    temp_geneset <- genesets_in %>%
      filter(Gene_set == x)
    
    set_n <- nrow(temp_geneset)
    set_name <- as.character(x)
    
    ## Get gene names and n genes for those outside the test gene set
    ## These 'remaining' genes are used as a background set in the two.sample poisson rate ratio test
    temp_remaining_geneset <- mutRates_in %>%
      filter(!ENTREZ %in% temp_geneset$ENTREZ) %>%
      mutate(Gene_set = "Remaining_genes") %>%
      select(Gene_set,ENTREZ)
    
    remaining_set_n <- nrow(temp_remaining_geneset)
    remaining_set_name <- paste("Non",set_name,sep = "_")
    
    ## Summarise genes without mutation rates or ENTREZ IDs in gene set
    temp_summarise_geneset <- summarise.geneset(mutRates_in = mutRates_in,
                                                geneset_in = temp_geneset)
    
    ## Calculate gene set mutation rates for all classes of mutation, separately for autosomes and X chromosome. 
    temp_set_mutRates <- geneset.mutationRates.list(mutRates_in = mutRates_in,
                                                    geneset_in = temp_geneset)
    ## Count de novo variants in gene set
    temp_geneset_denovo_counts <- get.denovo.counts(denovos_in = denovos_in,
                                                    geneset_in = temp_geneset)
    ## count missense MPC variants in gene set
    temp_geneset_MPC_counts <- get.MPC.denovo.counts(denovos_in = denovos_in,
                                                     geneset_in = temp_geneset)
    ## bind standard and MPC counts
    temp_geneset_denovo_counts <- rbind(temp_geneset_denovo_counts,
                                        temp_geneset_MPC_counts)
    
    ## Remaining gene metrics ##
    
    temp_remaining_summarise_geneset <- summarise.geneset(mutRates_in = mutRates_in,
                                                          geneset_in = temp_remaining_geneset)
    
    temp_remaining_set_mutRates <- geneset.mutationRates.list(mutRates_in = mutRates_in,
                                                              geneset_in = temp_remaining_geneset)
    
    temp_remaining_geneset_denovo_counts <- get.denovo.counts(denovos_in = denovos_in,
                                                              geneset_in = temp_remaining_geneset)
    
    temp_remaining_geneset_MPC_counts <- get.MPC.denovo.counts(denovos_in = denovos_in,
                                                               geneset_in = temp_remaining_geneset)
    
    temp_remaining_geneset_denovo_counts <- rbind(temp_remaining_geneset_denovo_counts,
                                                  temp_remaining_geneset_MPC_counts)
    
    all_genesets_summary <- rbind(temp_summarise_geneset,
                                  temp_remaining_summarise_geneset)
    ###
    
    ## Count both high and medium confidence DNVs (only necessary for new trios)
    temp_geneset_highConf_results <- temp_geneset_denovo_counts %>% 
      mutate(Observed = High + Medium) %>% 
      dplyr::select(Annotation,Observed)
    
    ## Generate enrichment P values and effect sizes
    temp_geneset_highConf_results <- Obs.vs.Exp.geneset.list(observed_counts = temp_geneset_highConf_results,
                                                             N_trios = N_trios,
                                                             auto_rates = temp_set_mutRates[,c("Annotation","Auto_Rates")],
                                                             chrX_rates = temp_set_mutRates[,c("Annotation","chrX_Rates")])
    
    ## Note: since switching to a competitive test, the P values and effect size are now redundent, so removed below. Perhaps remove from function?
    temp_geneset_highConf_results <- temp_geneset_highConf_results %>%
      select(Annotation,Observed,Expected)
    
    ## Add LoF and MPC2 for NS damaging analysis
    temp_geneset_highConf_results <-
      rbind(
        temp_geneset_highConf_results,
        temp_MPC2_LoF <- data.frame(
          Annotation = "MPC2_LoF",
          Observed = temp_geneset_highConf_results$Observed[temp_geneset_highConf_results$Annotation == "MPC2"] + temp_geneset_highConf_results$Observed[temp_geneset_highConf_results$Annotation == "LoF"],
          Expected = temp_geneset_highConf_results$Expected[temp_geneset_highConf_results$Annotation == "MPC2"] + temp_geneset_highConf_results$Expected[temp_geneset_highConf_results$Annotation == "LoF"]
        )
      )
    
    
    ## Generate observed and expected DNVs for remaining genes outside test set
    temp_remaining_geneset_highConf_results <- temp_remaining_geneset_denovo_counts %>% 
      mutate(Observed = High + Medium) %>% 
      dplyr::select(Annotation,Observed)
    
    temp_remaining_geneset_highConf_results <- Obs.vs.Exp.geneset.list(observed_counts = temp_remaining_geneset_highConf_results,
                                                                       N_trios = N_trios,
                                                                       auto_rates = temp_remaining_set_mutRates[,c("Annotation","Auto_Rates")],
                                                                       chrX_rates = temp_remaining_set_mutRates[,c("Annotation","chrX_Rates")])
    
    temp_remaining_geneset_highConf_results <- temp_remaining_geneset_highConf_results %>% 
      rename(Observed_remaining = Observed,
             Expected_remaining = Expected) %>%
      select(Annotation,Observed_remaining,Expected_remaining)
    
    temp_remaining_geneset_highConf_results <-
      rbind(
        temp_remaining_geneset_highConf_results,
        temp_MPC2_LoF <- data.frame(
          Annotation = "MPC2_LoF",
          Observed_remaining = temp_remaining_geneset_highConf_results$Observed_remaining[temp_remaining_geneset_highConf_results$Annotation == "MPC2"] + temp_remaining_geneset_highConf_results$Observed_remaining[temp_remaining_geneset_highConf_results$Annotation == "LoF"],
          Expected_remaining = temp_remaining_geneset_highConf_results$Expected_remaining[temp_remaining_geneset_highConf_results$Annotation == "MPC2"] + temp_remaining_geneset_highConf_results$Expected_remaining[temp_remaining_geneset_highConf_results$Annotation == "LoF"]
        )
      )
    
    ###
    
    ## Merge gene set and remaining gene set results
    merged_results <- merge(temp_geneset_highConf_results,temp_remaining_geneset_highConf_results,by = "Annotation")
    
    ## Loop through annotations and generate competitive DNV enrichment results
    for (i in 1:nrow(merged_results)) {
      
      temp_result <- poisson.test(x = c(merged_results$Observed[i],
                                        merged_results$Observed_remaining[i]),
                                  T = c(merged_results$Expected[i],
                                        merged_results$Expected_remaining[i]),
                                  r = 1)  
      temp_P <- temp_result$p.value
      temp_RR <- temp_result$estimate
      temp_lower_RR <- temp_result$conf.int[1]
      temp_upper_RR <- temp_result$conf.int[2]
      
      merged_results$P[i] <- temp_P
      merged_results$RR[i] <- temp_RR   
      merged_results$RR_lower[i] <- temp_lower_RR   
      merged_results$RR_upper[i] <- temp_upper_RR   
    }
    
    annotations_list <- list(
      NS = "NS",
      LoF = "LoF",
      Synonymous = "Synonymous",
      Missense = "Missense",
      MPC1 = "MPC1",
      MPC2 = "MPC2",
      MPC2_LoF = "MPC2_LoF",
      NS_SNV = "NS_SNV",
      LoF_SNV = "LoF_SNV"
    )
  
  
    test_loop <- lapply(annotations_list, function(class){
  
      merged_results %>%
        filter(Annotation == class) %>%
        mutate(Gene_set = set_name,
               !!paste(class,"_count",sep = '') := Observed,
               !!paste(class,"_expected",sep = '') := Expected,
               !!paste(class,"_count_remaining",sep = '') := Observed_remaining, 
               !!paste(class,"_expected_remaining",sep = '') := Expected_remaining,
               !!paste(class,"_P",sep = '') := P,
               !!paste(class,"_RR",sep = '') := RR,
               !!paste(class,"_RR_lower",sep = '') := RR_lower,
               !!paste(class,"_RR_upper",sep = '') := RR_upper) %>%
        dplyr::select(Gene_set,
                      !!paste(class,"_count",sep = ''),
                      !!paste(class,"_expected",sep = ''),
                      !!paste(class,"_count_remaining",sep = ''),
                      !!paste(class,"_expected_remaining",sep = ''),
                      !!paste(class,"_RR",sep = ''),
                      !!paste(class,"_RR_lower",sep = ''),
                      !!paste(class,"_RR_upper",sep = ''),
                      !!paste(class,"_P",sep = ''))
      
    })
    temp_DN_results <- test_loop %>% reduce(left_join, by = "Gene_set")
  
    temp_summarise_geneset$Genes_tested <- temp_summarise_geneset$Genes_in_set - (temp_summarise_geneset$Genes_without_mutRates + temp_summarise_geneset$Genes_without_Entrez_IDs)
    temp_remaining_summarise_geneset$Remaining_genes_tested <- temp_remaining_summarise_geneset$Genes_in_set - (temp_remaining_summarise_geneset$Genes_without_mutRates + temp_remaining_summarise_geneset$Genes_without_Entrez_IDs)
    temp_summarise_geneset$Remaining_genes_tested <- temp_remaining_summarise_geneset$Remaining_genes_tested
    
    temp_summarise_geneset <- temp_summarise_geneset %>% select(Gene_set,Genes_tested,Remaining_genes_tested)
    
    
    temp_DN_results <- merge(temp_summarise_geneset,temp_DN_results,by = "Gene_set")
    temp_DN_results
  })
  
  final_result = do.call(rbind,temp_results)
  
  return(final_result)
}
##########################

###################################################
## Summarise mutation rates for genes in geneSet ##
###################################################

summarise.geneset <- function(mutRates_in,
                              geneset_in) {
  
  ## total genes, genes without mut rates, genes without extrez IDs
  geneset_name <- unique(as.character(geneset_in$Gene_set))
  
  ## remove genes which don't have ENTREZ ID ##
  N_genes_in_set <- geneset_in %>% 
    dplyr::count()
  
  N_genes_in_set_without_mutRates <- geneset_in %>% 
    filter(!ENTREZ %in% mutRates_in$ENTREZ) %>%
    dplyr::count()
  
  N_genes_no_ENTREZ <- geneset_in %>% 
    filter(is.na(ENTREZ)) %>%
    dplyr::count()
  
  summarise_geneset <- data.frame(geneset_name,
                                  N_genes_in_set,
                                  N_genes_in_set_without_mutRates,
                                  N_genes_no_ENTREZ,stringsAsFactors = F)
  colnames(summarise_geneset) <- c("Gene_set","Genes_in_set","Genes_without_mutRates","Genes_without_Entrez_IDs")
  
  return(summarise_geneset)
}

###################################################

##################################################################
## Estimate gene set mutation rates for all classes of mutation ##
##################################################################

geneset.mutationRates.list <- function(mutRates_in,geneset_in) {
  ## remove genes which don't have ENTREZ ID ##
  geneset_in <- geneset_in %>% 
    filter(!is.na(ENTREZ))
  ## get mutation rates for genes in set
  set_mutRates <- mutRates_in %>% 
    filter(ENTREZ %in% geneset_in$ENTREZ)
  
  ## produce autosome and X chrome mutation rates ##
  autosome_mutRates <- subset(set_mutRates,set_mutRates$Chr != "X" & set_mutRates$Chr != "Y")
  xChr_mutRates <- subset(set_mutRates,set_mutRates$Chr == "X")
  
  auto_syn_sum <- sum(autosome_mutRates$syn,na.rm = T)
  auto_non_sum <- sum(autosome_mutRates$non,na.rm = T)
  auto_splice_sum <- sum(autosome_mutRates$splice,na.rm = T)
  auto_frame_sum <- sum(autosome_mutRates$frameshift,na.rm = T)
  auto_inframe_sum <- sum(autosome_mutRates$inframe,na.rm = T)
  auto_mis_sum <- sum(autosome_mutRates$mis,na.rm = T)
  
  auto_MPC1_sum <- sum(autosome_mutRates$MPC1_mutRate,na.rm = T)
  auto_MPC2_sum <- sum(autosome_mutRates$MPC2_mutRate,na.rm = T)
  
  auto_LoF_sum <- auto_non_sum + auto_splice_sum + auto_frame_sum
  auto_NS_sum <- auto_non_sum + auto_splice_sum + auto_frame_sum + auto_inframe_sum + auto_mis_sum
  
  ## no indels
  auto_LoF_SNV_sum <- auto_non_sum + auto_splice_sum
  auto_NS_SNV_sum <- auto_non_sum + auto_splice_sum + auto_mis_sum
  
  auto_rates <- data.frame(row.names = c("Synonymous",
                                         "Missense",
                                         "Splice",
                                         "Nonsense",
                                         "FrameShift",
                                         "Inframe",
                                         "MPC1",
                                         "MPC2",
                                         "LoF",
                                         "NS",
                                         "LoF_SNV",
                                         "NS_SNV"), 
                          Auto_Rates = c(auto_syn_sum,
                                         auto_mis_sum,
                                         auto_splice_sum,
                                         auto_non_sum,
                                         auto_frame_sum,
                                         auto_inframe_sum,
                                         auto_MPC1_sum,
                                         auto_MPC2_sum,
                                         auto_LoF_sum,
                                         auto_NS_sum,
                                         auto_LoF_SNV_sum,
                                         auto_NS_SNV_sum)) 
  
  ## X Chromosome 
  chrx_syn_sum <- sum(xChr_mutRates$syn,na.rm = T)
  chrx_non_sum <- sum(xChr_mutRates$non,na.rm = T)
  chrx_splice_sum <- sum(xChr_mutRates$splice,na.rm = T)
  chrx_frame_sum <- sum(xChr_mutRates$frameshift,na.rm = T)
  chrx_inframe_sum <- sum(xChr_mutRates$inframe,na.rm = T)
  chrx_mis_sum <- sum(xChr_mutRates$mis,na.rm = T)
  
  chrx_MPC1_sum <- sum(xChr_mutRates$MPC1_mutRate,na.rm = T)
  chrx_MPC2_sum <- sum(xChr_mutRates$MPC2_mutRate,na.rm = T)
  
  chrx_LoF_sum <- chrx_non_sum + chrx_splice_sum + chrx_frame_sum
  chrx_NS_sum <- chrx_non_sum + chrx_splice_sum + chrx_frame_sum + chrx_inframe_sum + chrx_mis_sum
  
  ## no indels
  chrx_LoF_SNV_sum <- chrx_non_sum + chrx_splice_sum
  chrx_NS_SNV_sum <- chrx_non_sum + chrx_splice_sum + chrx_mis_sum
  
  chrx_rates <- data.frame(row.names = c("Synonymous",
                                         "Missense",
                                         "Splice",
                                         "Nonsense",
                                         "FrameShift",
                                         "Inframe",
                                         "MPC1",
                                         "MPC2",
                                         "LoF",
                                         "NS",
                                         "LoF_SNV",
                                         "NS_SNV"), 
                          chrX_Rates = c(chrx_syn_sum,
                                         chrx_mis_sum,
                                         chrx_splice_sum,
                                         chrx_non_sum,
                                         chrx_frame_sum,
                                         chrx_inframe_sum,
                                         chrx_MPC1_sum,
                                         chrx_MPC2_sum,
                                         chrx_LoF_sum,
                                         chrx_NS_sum,
                                         chrx_LoF_SNV_sum,
                                         chrx_NS_SNV_sum)) 
  
  return_rates <- merge(auto_rates,chrx_rates,by = 0)
  colnames(return_rates)[1] <- "Annotation"
  
  return(return_rates)
}
##################################################################

############################################
#### Count de novo variants in gene set ####
############################################

get.denovo.counts <- function(denovos_in,
                              geneset_in) {
  
  ## Analyse high confidence autosomal de novos (IGV exclude) from high QC trios ##
  AnnotationsToGet <- c("Missense_SNV","Frameshift","Inframe","Splice_SNV","stop_gained","Synonymous_SNV")
  LoF_annotations <- c("Frameshift","Splice_SNV","stop_gained")
  NS_annotations <- c("Frameshift","Splice_SNV","stop_gained","Inframe","Missense_SNV")
  NS_MPC1_annotations <- c("Frameshift","Splice_SNV","stop_gained","Missense_MPC1")
  NS_MPC2_annotations <- c("Frameshift","Splice_SNV","stop_gained","Missense_MPC2")
  LoF_SNV_annotations <- c("Splice_SNV","stop_gained")
  NS_SNV_annotations <- c("Splice_SNV","stop_gained","Missense_SNV")

  
  
  if (nrow(denovos_in %>% 
           filter(Entrez %in% geneset_in$ENTREZ &
                  Annotation_2 %in% AnnotationsToGet)) > 0) {
    
    denovo_count_in_set <- denovos_in %>% 
      filter(Entrez %in% geneset_in$ENTREZ) %>% ## only look at de novos in genes that have mutation rates
      dplyr::group_by(Annotation_2,Confidence_2) %>% 
      dplyr::summarise(Count=n())
    
    ## Retrieve important annotations ##
    denovo_count_in_set <- denovo_count_in_set %>% 
      filter(Annotation_2 %in% AnnotationsToGet)
    
    
    DN_counts <- dcast(data = setDT(denovo_count_in_set),formula = Annotation_2~Confidence_2,value.var = "Count")
    colnames(DN_counts)[1] <- c("Annotation")
    DN_counts[is.na(DN_counts)] <- 0
    
    if (!"High" %in% colnames(DN_counts)) {
      DN_counts$High <- 0
    }
    if (!"Low" %in% colnames(DN_counts)) {
      DN_counts$Low <- 0
    }
    if (!"Medium" %in% colnames(DN_counts)) {
      DN_counts$Medium <- 0
    }
    
    #DN_counts <- DN_counts[,c(1,2,4,3)]
    
    ## Sum NS and LoF counts ##
    DN_counts <- rbind(DN_counts, data.frame(Annotation="LoF",t(colSums(DN_counts[DN_counts$Annotation %in% LoF_annotations,-1],na.rm = T))))
    DN_counts <- rbind(DN_counts, data.frame(Annotation="NS",t(colSums(DN_counts[DN_counts$Annotation %in% NS_annotations,-1],na.rm = T)))) 
    
    ## Sum NS and LoF SNV counts ##
    DN_counts <- rbind(DN_counts, data.frame(Annotation="LoF_SNV",t(colSums(DN_counts[DN_counts$Annotation %in% LoF_SNV_annotations,-1],na.rm = T))))
    DN_counts <- rbind(DN_counts, data.frame(Annotation="NS_SNV",t(colSums(DN_counts[DN_counts$Annotation %in% NS_SNV_annotations,-1],na.rm = T))))
    
  }
  else {
    DN_counts <- data.frame(Annotation = c("Frameshift","Missense_SNV","Synonymous_SNV","stop_gained","LoF","NS","LoF_SNV","NS_SNV"),High = 0, Low = 0,Medium = 0)
  }
  
  if (length(AnnotationsToGet[!AnnotationsToGet %in% DN_counts$Annotation]) > 0) {
    temp_missing <- data.frame(Annotation = c(AnnotationsToGet[!AnnotationsToGet %in% DN_counts$Annotation]),High = 0, Low = 0,Medium = 0)
    DN_counts <- rbind(DN_counts,temp_missing)
  }
  return(DN_counts)
}
############################################

################################################
#### Count MPC de novo variants in gene set ####
################################################

get.MPC.denovo.counts <- function(denovos_in,
                                  geneset_in){
  
  ## Analyse high confidence autosomal de novos (IGV exclude) from high QC trios ##
  AnnotationsToGet <- c("Missense_MPC1","Missense_MPC2")
  
  denovos_in$MPC[denovos_in$MPC == "."] <- NA
  if (!is.numeric(denovos_in$MPC)) {
    denovos_in$MPC <- as.numeric(levels(denovos_in$MPC))[denovos_in$MPC]
  }
  denovos_in$Annotation_2 <- as.character(denovos_in$Annotation_2)
  for (i in 1:nrow(denovos_in)) {
    
    if (!is.na(denovos_in$MPC[i])) {
      if (denovos_in$MPC[i] >= 2) {
        denovos_in$Annotation_2[i] <- "Missense_MPC2"
      }
      else if (denovos_in$MPC[i] >= 1) {
        denovos_in$Annotation_2[i] <- "Missense_MPC1"
      }
      
    }
  }
  
  if (nrow(denovo_count_in_set <- denovos_in %>% 
           filter(Entrez %in% geneset_in$ENTREZ &
                  Annotation_2 %in% AnnotationsToGet)) > 0) {
    
    denovo_count_in_set <- denovos_in %>% 
      filter(Entrez %in% geneset_in$ENTREZ) %>% ## only look at de novos in genes that have mutation rates
      group_by(Confidence_2,Annotation_2) %>% 
      dplyr::summarise(Count=n())
    
    ## Retrieve important annotations ##
    denovo_count_in_set <- denovo_count_in_set %>% 
      filter(Annotation_2 %in% AnnotationsToGet)
    
    DN_counts <- dcast(data = setDT(denovo_count_in_set), formula = Annotation_2~Confidence_2, value.var = "Count")
    colnames(DN_counts)[1] <- c("Annotation")
    DN_counts[is.na(DN_counts)] <- 0
    
    if (!"High" %in% colnames(DN_counts)) {
      DN_counts$High <- 0
    }
    if (!"Low" %in% colnames(DN_counts)) {
      DN_counts$Low <- 0
    }
    if (!"Medium" %in% colnames(DN_counts)) {
      DN_counts$Medium <- 0
    }
    
  }
  else {
    DN_counts <- data.frame(Annotation = c("Missense_MPC1","Missense_MPC2"),High = 0, Low = 0,Medium = 0)
  }
  
  if (length(AnnotationsToGet[!AnnotationsToGet %in% DN_counts$Annotation]) > 0) {
    temp_missing <- data.frame(Annotation = c(AnnotationsToGet[!AnnotationsToGet %in% DN_counts$Annotation]),High = 0, Low = 0,Medium = 0)
    DN_counts <- rbind(DN_counts,temp_missing)
  }
  
  MPC1_temp <- data.frame(Annotation="Missense_MPC1",High = sum(DN_counts$High,na.rm = T),Low = sum(DN_counts$Low,na.rm = T),Medium = sum(DN_counts$Medium,na.rm = T))
  DN_counts <- rbind(DN_counts %>% filter(Annotation == "Missense_MPC2"),MPC1_temp)
  return(DN_counts)
}
################################################


## Obs.vs.Exp.geneset.list ####
Obs.vs.Exp.geneset.list <- function(observed_counts,
                                    auto_rates,
                                    chrX_rates,
                                    N_trios) {
  
  observed_counts$Annotation <-
    mapvalues(
      observed_counts$Annotation,
      from = c(
        "Frameshift",
        "Splice_SNV",
        "stop_gained",
        "Synonymous_SNV",
        "Missense_SNV",
        "Missense_MPC1",
        "Missense_MPC2"
      )
      ,
      to = c(
        "FrameShift",
        "Splice",
        "Nonsense",
        "Synonymous",
        "Missense",
        "MPC1",
        "MPC2"
      )
    )
  
  ## Generate expected counts ##
  auto_expected_counts <- data.frame(Annotation = auto_rates$Annotation, 
                                     Exp = auto_rates$Auto_Rates * (N_trios[["total_probands"]] * 2))
  
  chrX_expected_counts <- data.frame(Annotation = chrX_rates$Annotation,
                                     Exp = ((chrX_rates$chrX_Rates * (N_trios[["female_probands"]] * 2)) +
                                            (chrX_rates$chrX_Rates * N_trios[["male_probands"]])
                                            )
                                     )
  
  expected_counts <- merge(auto_expected_counts,chrX_expected_counts,by = "Annotation") %>% 
    mutate(Expected = Exp.x + Exp.y) %>%
    dplyr::select(Annotation,Expected)
  ##############################
  obs_exp <- merge(observed_counts,expected_counts,by = "Annotation")
  obs_exp[is.na(obs_exp)] <- 0
  obs_exp$P <- NA
  obs_exp$OR <- NA
  obs_exp$lower <- NA
  obs_exp$higher <- NA
  
  ## Observed vs expected tests. Note, this does not control for background rates
  for (i in 1:nrow(obs_exp)) {
    obs_exp$P[i] <- ppois(q = obs_exp$Observed[i] - 1,
                          lambda = obs_exp$Expected[i],
                          lower.tail = F)
    obs_exp$OR[i] <- poisson.test(obs_exp$Observed[i],
                                  obs_exp$Expected[i])$estimate
    obs_exp$lower[i] <- poisson.test(obs_exp$Observed[i],
                                     obs_exp$Expected[i])$conf.int[1]
    obs_exp$higher[i] <- poisson.test(obs_exp$Observed[i],
                                      obs_exp$Expected[i])$conf.int[2]
  }
  return(obs_exp)
}








