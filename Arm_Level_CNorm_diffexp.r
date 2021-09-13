#requires("edgeR")
#requires("limma")
#requires("tidyverse")
#requires("data.table")
#requires("ggpubr")
#requires("ggplot2")
#
#uses a list of dataframes generated from the ChrArm_CNA_freq_TCGA_puritycorrection.R function, "ChrArm_CNA_summary"
#
#uses segment CNA files and RSEM gene expression files from TCGA database
#
#uses purity files downloaded from https://zenodo.org/record/253193#.YT-EfDZKi-w 
#
#adjustable parameters:
#RSEM_min_count = minimum average RSEM value
#cohort_subsample_minFrac = minimum remaining fraction of cohort after subsampling (suggested value: 0.7), set to 1 for no CNorm
#arm = chromosome arm of interest 
#gain_loss = CNA direction
#acceptable_difference_percent = difference in cohort incidence to aim for (suggested value: 5)
#iter = maximum number of iterations of subsampling (suggested value: 30)
#
#vignette: 
#diffexp_gain8q_CNorm <- Arm_Level_CNorm_diffexp(BRCA_RSEM, BRCA_CNA, BRCA_purity, BRCA_CNA_summary, 10, 0.7, "1q", "gain", 5, 30)

Arm_Level_CNorm_diffexp <- function(TCGA_RSEM, TCGA_purity, TCGA_CNA, ChrArm_CNA_summary, RSEM_min_count, cohort_subsample_minFrac, arm, gain_loss, acceptable_difference_percent, iter) {
    
    #define function
    CNA_distributions_plot_func <- function(x) {
    CNA_distributions_ChrArm2 <- x[,c(1,4)]
    CNA_distributions_ChrArm3 <- x[,c(1,5)]
    CNA_distributions_ChrArm2$status <- "yes"
    CNA_distributions_ChrArm3$status <- "no"
    colnames(CNA_distributions_ChrArm3)[2] <- "yes"
    CNA_distributions_ChrArm4 <- rbind(CNA_distributions_ChrArm2,CNA_distributions_ChrArm3)
    return(CNA_distributions_ChrArm4) }

    #preparing files
    print("... preparing RSEM file ... ")
    TCGA_RSEM <- TCGA_RSEM[-1,]
    TCGA_RSEM <- separate(data = TCGA_RSEM, col = Hybridization.REF, into = c('Gene', 'x'), sep = "\\|")
    TCGA_RSEM$Gene <- ifelse(TCGA_RSEM$Gene =="?", TCGA_RSEM$x, TCGA_RSEM$Gene)
    TCGA_RSEM <- TCGA_RSEM[,-2]
    TCGA_RSEM <- data.frame(TCGA_RSEM)
    rownames(TCGA_RSEM) <- make.names(TCGA_RSEM$Gene, unique=TRUE)
    TCGA_RSEM <- TCGA_RSEM[,-1]
    as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
    for(i in 1:ncol(TCGA_RSEM)){
        TCGA_RSEM[i] <- as.numeric.factor(unlist(TCGA_RSEM[i]))
    }
    TCGA_RSEM_ids <- colnames(TCGA_RSEM)
    TCGA_RSEM_ids <- data.frame(TCGA_RSEM_ids)
    colnames(TCGA_RSEM_ids)[1] <- c("id")
    TCGA_RSEM_ids <- separate(data = TCGA_RSEM_ids, col = id, into = c('TCGA', 'TSS', 'Participant', 'SampleVial', 'portion', 'plate', 'center'), sep = "\\.")
    TCGA_RSEM_ids$id <- paste(TCGA_RSEM_ids$TCGA, TCGA_RSEM_ids$TSS, TCGA_RSEM_ids$Participant, TCGA_RSEM_ids$SampleVial, sep='-')
    TCGA_RSEM_ids$SampleVial <- gsub('.{1}$', '', TCGA_RSEM_ids$SampleVial)
    TCGA_RSEM_ids$id <- gsub('.{1}$', '', TCGA_RSEM_ids$id)
    TCGA_RSEM_ids$Sample <- colnames(TCGA_RSEM)
    TCGA_RSEM_ids$SampleVial <- as.numeric(as.character(TCGA_RSEM_ids$SampleVial))
    TCGA_RSEM_ids <- subset(TCGA_RSEM_ids, SampleVial < 6)
    TCGA_RSEM2 <- TCGA_RSEM[, names(TCGA_RSEM) %in% TCGA_RSEM_ids$Sample]
    print("success")
    #remove samples with no CNA or purity data
    print("... preparing merge of CNA, RSEM, purity data ... ")
    TCGA_purity_samples <- TCGA_purity
    TCGA_purity_samples <- separate(data = TCGA_purity_samples, col = SampleName, into = c('TCGA', 'TSS', 'Participant', 'SampleVial', 'portion', 'plate', 'center'), sep = "\\-")
    TCGA_purity_samples$id <- paste(TCGA_purity_samples$TCGA, TCGA_purity_samples$TSS, TCGA_purity_samples$Participant, TCGA_purity_samples$Sample, sep='-')
    TCGA_purity_samples$id <- gsub('.{1}$', '', TCGA_purity_samples$id)
    TCGA_purity_samples = TCGA_purity_samples[!duplicated(TCGA_purity_samples$id),]
    TCGA_CNA <- TCGA_CNA[,c(1,2)]
    TCGA_CNA <- separate(data = TCGA_CNA, col = Sample, into = c('TCGA', 'TSS', 'Participant', 'SampleVial', 'portion', 'plate', 'center'), sep = "\\-")
    TCGA_CNA$id <- paste(TCGA_CNA$TCGA, TCGA_CNA$TSS, TCGA_CNA$Participant, TCGA_CNA$SampleVial, sep='-')
    TCGA_CNA$id <- gsub('.{1}$', '', TCGA_CNA$id)
    TCGA_CNA$SampleVial <- gsub('.{1}$', '', TCGA_CNA$SampleVial)
    TCGA_CNA <- subset(TCGA_CNA, SampleVial < 6)
    TCGA_CNA_samples <- unique(TCGA_CNA$id)
    TCGA_CNA_samples <- data.frame(TCGA_CNA_samples)
    colnames(TCGA_CNA_samples) <- c("id")
    TCGA_CNA_samples$CNAcheck <- 1
    TCGA_RSEM_ids <- merge(TCGA_RSEM_ids, TCGA_purity_samples, by="id")
    TCGA_RSEM_ids <- merge(TCGA_RSEM_ids, TCGA_CNA_samples, by="id")
    TCGA_RSEM_ids <- distinct(TCGA_RSEM_ids,Sample, .keep_all= TRUE)
    TCGA_RSEM2 <- TCGA_RSEM2[, names(TCGA_RSEM2) %in% TCGA_RSEM_ids$Sample]
    print("success")
    w <- ncol(TCGA_RSEM2)/2
    #remove gene rows with low counts
    TCGA_RSEM3 <- TCGA_RSEM2[rowSums(TCGA_RSEM2 < RSEM_min_count) <=  w, , drop = FALSE]
    #grouping based on CNA status
    print("... grouping based on CNA status ...")
    `%ni%` <- Negate(`%in%`)
    ChrArm_summary_gain <- ChrArm_CNA_summary[[3]]
    ChrArm_summary_loss <- ChrArm_CNA_summary[[1]]
    ChrArm_summary_gain$id <- gsub('.{1}$', '', ChrArm_summary_gain$id)
    ChrArm_summary_loss$id <- gsub('.{1}$', '', ChrArm_summary_loss$id)
    ChrArm_summary_gain$status <- "gain"
    ChrArm_summary_loss$status <- "loss"
    ChrArm_summary_gain <- ChrArm_summary_gain[,c(1,4:5)]
    ChrArm_summary_loss <- ChrArm_summary_loss[,c(1,4:5)]
    ChrArm_summary <- rbind(ChrArm_summary_gain, ChrArm_summary_loss)
    ChrArm_summary$ChrArm <- as.character(ChrArm_summary$ChrArm)
    ChrArm_summary$ChrArm2 <- paste(ChrArm_summary$ChrArm, ChrArm_summary$status, sep="_")
    ChrArm_summary <- subset(ChrArm_summary, id %in% TCGA_RSEM_ids$id)
    ChrArmY <- paste(arm, gain_loss, sep="_")
    ids_to_grab1 <- subset(ChrArm_summary2, ChrArm == arm & status == gain_loss)$id
    ChrArm_summary_subset1 <- subset(ChrArm_summary2, id %in% ids_to_grab1)
    ChrArm_summary_subset_ELSE1 <- subset(ChrArm_summary2, id %ni% ids_to_grab1)
    min3subset <- round(length(unique(ChrArm_summary_subset1$id))*cohort_subsample_minFrac)
    min3ELSE <- round(length(unique(ChrArm_summary_subset_ELSE1$id))*cohort_subsample_minFrac)
    print("success")
        ChrArm_summary2 <- ChrArm_summary
        ChrArm_list <- unique(ChrArm_summary$ChrArm2)
        #my.list2 <- list()
    print("... performing CNorm to generate custom cohort ...")
        for(i in 1:iter) {
            
            ids_to_grab <- subset(ChrArm_summary2, ChrArm == arm & status == gain_loss)$id
            ChrArm_summary_subset <- subset(ChrArm_summary2, id %in% ids_to_grab)
            ChrArm_summary_subset_ELSE <- subset(ChrArm_summary2, id %ni% ids_to_grab)
            CNA_distributions_ChrArm <- matrix(ncol=5, nrow=length(unique(ChrArm_summary$ChrArm2)))
            colnames(CNA_distributions_ChrArm) <- c("ChrArm_status","yes_count", "no_count", "yes", "no")
            
            for(i in 1:length(ChrArm_list)){
                ChrArm_X <- ChrArm_list[i]
                CNA_distributions_ChrArm[i,1] <- ChrArm_X
                if (ChrArm_X == ChrArmY)
                    next
                test_ids <- subset(ChrArm_summary_subset_ELSE, ChrArm2 == ChrArm_X)
                test_ids2 <- subset(ChrArm_summary_subset, ChrArm2 == ChrArm_X)
                CNA_distributions_ChrArm[i,2] <- nrow(test_ids2)
                CNA_distributions_ChrArm[i,3] <- nrow(test_ids)
                CNA_distributions_ChrArm[i,4] <- 100*nrow(test_ids2)/length(unique(ChrArm_summary_subset$id))
                CNA_distributions_ChrArm[i,5] <- 100*nrow(test_ids)/length(unique(ChrArm_summary_subset_ELSE$id))
            }
            CNA_distributions_ChrArm <- data.frame(CNA_distributions_ChrArm)
            CNA_distributions_ChrArm$yes <- as.numeric(as.character(CNA_distributions_ChrArm$yes))
            CNA_distributions_ChrArm$no <- as.numeric(as.character(CNA_distributions_ChrArm$no))
            CNA_distributions_ChrArm$yes_count <- as.numeric(as.character(CNA_distributions_ChrArm$yes_count))
            CNA_distributions_ChrArm$no_count <- as.numeric(as.character(CNA_distributions_ChrArm$no_count))
            CNA_distributions_ChrArm$diff <- CNA_distributions_ChrArm$yes - CNA_distributions_ChrArm$no
            CNA_distributions_ChrArm$yes_to_remove <- CNA_distributions_ChrArm$yes_count - (CNA_distributions_ChrArm$no/100)*(length(unique(ChrArm_summary_subset$id))) 
            CNA_distributions_ChrArm$no_to_remove <- CNA_distributions_ChrArm$no_count - (CNA_distributions_ChrArm$yes/100)*(length(unique(ChrArm_summary_subset_ELSE$id)))
            CNA_distributions_ChrArm$yes_to_remove <- round(CNA_distributions_ChrArm$yes_to_remove)
            CNA_distributions_ChrArm$no_to_remove <- round(CNA_distributions_ChrArm$no_to_remove)
            CNA_distributions_ChrArm2 <- subset(CNA_distributions_ChrArm, abs(diff) > acceptable_difference_percent)
            CNA_distributions_ChrArm2 <- subset(CNA_distributions_ChrArm2, yes_count != 0 & no_count !=0)
            if (cohort_subsample_minFrac == 1)
                break
            if (nrow(CNA_distributions_ChrArm2) == 0)
                break
            if (length(unique(ChrArm_summary_subset$id)) < min3subset)
                break
            if (length(unique(ChrArm_summary_subset_ELSE$id)) < min3ELSE)
                break
            CNA_distributions_ChrArm2$same_or_else <- ifelse(CNA_distributions_ChrArm2$no_to_remove > 0, "same", "else")
            CNA_distributions_ChrArm2$ChrArm_status <- as.character(CNA_distributions_ChrArm2$ChrArm_status)
            CNA_distributions_ChrArm2 <- CNA_distributions_ChrArm2[sample(nrow(CNA_distributions_ChrArm2)),]
            ChrArm_summary_subset_ELSE <- ChrArm_summary_subset_ELSE[sample(nrow(ChrArm_summary_subset_ELSE)),]
            ChrArm_summary_subset <- ChrArm_summary_subset[sample(nrow(ChrArm_summary_subset)),]
            ELSE_ids_removable <- sort(table(subset(ChrArm_summary_subset_ELSE, ChrArm2 %ni% subset(CNA_distributions_ChrArm2, same_or_else == "else")$ChrArm_status )$id), decreasing = TRUE)
            ELSE_ids_removable <- data.frame(ELSE_ids_removable)
            ELSE_ids_removable$Var1 <- as.character(ELSE_ids_removable$Var1)
            ELSE_ids_removable <- ELSE_ids_removable$Var1
            ELSE_ids_removable2 <- sort(table(subset(ChrArm_summary_subset_ELSE, ChrArm2 %in% subset(CNA_distributions_ChrArm2, same_or_else == "same")$ChrArm_status )$id), decreasing = TRUE)
            ELSE_ids_removable2 <- data.frame(ELSE_ids_removable2)
            ELSE_ids_removable2$Var1 <- as.character(ELSE_ids_removable2$Var1)
            ELSE_ids_removable2 <- ELSE_ids_removable2$Var1
            
            subset_ids_removable <- sort(table(subset(ChrArm_summary_subset, ChrArm2 %in% subset(CNA_distributions_ChrArm2, yes_to_remove > 0)$ChrArm_status )$id), decreasing = TRUE)
            subset_ids_removable <- data.frame(subset_ids_removable)
            subset_ids_removable$Var1 <- as.character(subset_ids_removable$Var1)
            subset_ids_removable <- subset_ids_removable$Var1
            
            if (CNA_distributions_ChrArm2$same_or_else[1] == "same") {
                ids_to_remove <- ELSE_ids_removable2[1:round(abs(CNA_distributions_ChrArm2$no_to_remove[1])/3)]
                ChrArm_summary2 <- subset(ChrArm_summary2, id %ni%  ids_to_remove) 
            } else {
                ids_to_remove <- subset_ids_removable[1:round(abs(CNA_distributions_ChrArm2$yes_to_remove[1])/3)]
                ChrArm_summary2 <- subset(ChrArm_summary2, id %ni%  ids_to_remove)
            }
        }
        print("success")
        print("... performing edgeR analysis ...")
        CNA_distributions_ChrArm5 <-  CNA_distributions_plot_func(CNA_distributions_ChrArm)
        after_plot <- ggpar(ggbarplot(CNA_distributions_ChrArm5, x="ChrArm_status", y ="yes", fill="status", size = 0, palette = c("#00AFBB", "#FC4E07"), position = position_dodge(0.9)), font.tickslab = c(6), xtickslab.rt = 45)
        ChrArm_summary_subset <- subset(ChrArm_summary2, ChrArm == arm & status == gain_loss)
        CNA_ids_list <- unique(ChrArm_summary2$id)
        TCGA_RSEM_ids$CNA_status <- ifelse(TCGA_RSEM_ids$id %in% ChrArm_summary_subset$id, 2, 1)
        TCGA_RSEM_ids$remove_status <- ifelse(TCGA_RSEM_ids$id %in% ChrArm_summary2$id, 1, 0)
        to.remove <- subset(TCGA_RSEM_ids, remove_status==0)$Sample
        TCGA_RSEM3 <- subset(TCGA_RSEM3,select = names(TCGA_RSEM3) %ni% to.remove)
        TCGA_RSEM_ids <- subset(TCGA_RSEM_ids, remove_status ==1 )
        yEndtoEnd2 <- DGEList(counts=TCGA_RSEM3, genes=rownames(TCGA_RSEM3))
        labels <- data.frame("Name"=colnames(yEndtoEnd2), "TimePoint"=TCGA_RSEM_ids$CNA_status)
        timepoint<-factor(labels$TimePoint)
        design<-model.matrix(~0+timepoint)
        xglm2 <-  estimateDisp(yEndtoEnd2, design)
        fit <- glmFit(xglm2, design)
        lrt <- glmLRT(fit, contrast=c(-1,1))
        lrtres <- data.frame(lrt$genes,lrt$table)
        lrtres$neglogP <- (abs(lrtres$logFC)/lrtres$logFC)*(-1)*log10(lrtres$PValue)
        print("success")
        print("... preparing plot ...")
        ChrArm_summary2 <- ChrArm_summary
        ChrArm_list <- unique(ChrArm_summary$ChrArm2)
        ids_to_grab <- subset(ChrArm_summary2, ChrArm == arm & status == gain_loss)$id
        ChrArm_summary_subset <- subset(ChrArm_summary2, id %in% ids_to_grab)
        ChrArm_summary_subset_ELSE <- subset(ChrArm_summary2, id %ni% ids_to_grab)
        CNA_distributions_ChrArm <- matrix(ncol=5, nrow=length(unique(ChrArm_summary$ChrArm2)))
        colnames(CNA_distributions_ChrArm) <- c("ChrArm_status","yes_count", "no_count", "yes", "no")
        
        for(i in 1:length(ChrArm_list)){
            ChrArm_X <- ChrArm_list[i]
            CNA_distributions_ChrArm[i,1] <- ChrArm_X
            if (ChrArm_X == ChrArmY)
                next
            test_ids <- subset(ChrArm_summary_subset_ELSE, ChrArm2 == ChrArm_X)
            test_ids2 <- subset(ChrArm_summary_subset, ChrArm2 == ChrArm_X)
            CNA_distributions_ChrArm[i,2] <- nrow(test_ids2)
            CNA_distributions_ChrArm[i,3] <- nrow(test_ids)
            CNA_distributions_ChrArm[i,4] <- 100*nrow(test_ids2)/length(unique(ChrArm_summary_subset$id))
            CNA_distributions_ChrArm[i,5] <- 100*nrow(test_ids)/length(unique(ChrArm_summary_subset_ELSE$id))
        }
        CNA_distributions_ChrArm <- data.frame(CNA_distributions_ChrArm)
        CNA_distributions_ChrArm$yes <- as.numeric(as.character(CNA_distributions_ChrArm$yes))
        CNA_distributions_ChrArm$no <- as.numeric(as.character(CNA_distributions_ChrArm$no))
        CNA_distributions_ChrArm$yes_count <- as.numeric(as.character(CNA_distributions_ChrArm$yes_count))
        CNA_distributions_ChrArm$no_count <- as.numeric(as.character(CNA_distributions_ChrArm$no_count))
        CNA_distributions_ChrArm$diff <- CNA_distributions_ChrArm$yes - CNA_distributions_ChrArm$no
        CNA_distributions_ChrArm$yes_to_remove <- CNA_distributions_ChrArm$yes_count - (CNA_distributions_ChrArm$no/100)*(length(unique(ChrArm_summary_subset$id))) 
        CNA_distributions_ChrArm$no_to_remove <- CNA_distributions_ChrArm$no_count - (CNA_distributions_ChrArm$yes/100)*(length(unique(ChrArm_summary_subset_ELSE$id)))
        CNA_distributions_ChrArm$yes_to_remove <- round(CNA_distributions_ChrArm$yes_to_remove)
        CNA_distributions_ChrArm$no_to_remove <- round(CNA_distributions_ChrArm$no_to_remove)
        CNA_distributions_ChrArm4 <-  CNA_distributions_plot_func(CNA_distributions_ChrArm)
        before_plot <- ggpar(ggbarplot(CNA_distributions_ChrArm4, x="ChrArm_status", y ="yes", fill="status", size = 0, palette = c("#00AFBB", "#FC4E07"), position = position_dodge(0.9)), font.tickslab = c(6), xtickslab.rt = 45)
    
    theplot <- ggarrange(before_plot, after_plot, 
                         labels = c("before normalization", "after normalization"),
                         ncol = 1, nrow = 2)
    print(theplot)
    return(lrtres)
}
