#requires(tidyverse)
#requires(data.table)
#
#uses a chromosome information file such as the file chromdata.txt included in this package (i.e. "chromdata")
#
#uses segment CNA files from TCGA database, purity files downloaded from https://zenodo.org/record/253193#.YT-EfDZKi-w 
#
#generates list of four dataframes called final_output: final_output[[1]] is the Deletion Tumor Sample summary, final_output[[2]] is the Deletion Frequency summary, 
#final_output[[3]] is the Amplification Tumor Sample summary, and final_output[[4]] is the Amplification Frequency summary
#
#also prints out files with lists of samples called for each potential arm-level CNA
#
#adjustable parameters: y = name on printed text files, z = minimal fraction of chromosome arm included in CNA, w = maximal Segment_Mean value to call losses, g = minimal Segment_Mean value to call gains

ChrArm_CNA_freq_TCGA_puritycorrection <- function(TCGA_CNA, TCGA_purity, y, z, w, g, chromdata) { 
    #prep files
    print("... preparing data...")
    TCGA_purity <- separate(data = TCGA_purity, col = SampleName, into = c('TCGA', 'TSS', 'Participant', 'SampleVial', 'portion', 'plate', 'center'), sep = "\\-")
    TCGA_purity$id <- paste(TCGA_purity$TCGA, TCGA_purity$TSS, TCGA_purity$Participant, TCGA_purity$SampleVial, sep='-')
    TCGA_purity <- TCGA_purity[,-c(1:7)]
    TCGA_purity = TCGA_purity[!duplicated(TCGA_purity$id),]
    TCGA_CNA <- separate(data = TCGA_CNA, col = Sample, into = c('TCGA', 'TSS', 'Participant', 'SampleVial', 'portion', 'plate', 'center'), sep = "\\-")
    TCGA_CNA$id <- paste(TCGA_CNA$TCGA, TCGA_CNA$TSS, TCGA_CNA$Participant, TCGA_CNA$SampleVial, sep='-')
    TCGA_CNA <- subset(TCGA_CNA, TCGA_CNA$SampleVial == "01A" | TCGA_CNA$SampleVial == "01B" )
    TCGA_CNA <- TCGA_CNA[,-c(1:7)]
    print("success")
    #correct CNAs based on purity
    print("... correcting CNAs based on purity ...")
    TCGA_CNA_purity <- merge(TCGA_CNA, TCGA_purity, by = "id")
    TCGA_CNA_purity$Segment_Mean <- log2(((2^(TCGA_CNA_purity$Segment_Mean) )/TCGA_CNA_purity$Purity_InfiniumPurify) - ((1-TCGA_CNA_purity$Purity_InfiniumPurify)/(TCGA_CNA_purity$Purity_InfiniumPurify)))
    TCGA_CNA_purity$Segment_Mean[is.na(TCGA_CNA_purity$Segment_Mean)] <- 0
    TCGA_CNA_purity <- TCGA_CNA_purity[,-7]
    print("success")
    #assign segments to arms, deal with centromere-spanning segments
    print("... assigning segments to arms ...")
    TCGA_CNA_purity <- merge(TCGA_CNA_purity, chromdata, by = "Chromosome")
    TCGA_CNA_purity$length <- TCGA_CNA_purity$End - TCGA_CNA_purity$Start
    TCGA_CNA_purity$length <- as.numeric(as.character(TCGA_CNA_purity$length))
    TCGA_CNA_purity$Mid_Segment <- TCGA_CNA_purity$End - (TCGA_CNA_purity$End - TCGA_CNA_purity$Start)/2
    TCGA_CNA_purity$Mid_Segment <- as.numeric(as.character(TCGA_CNA_purity$Mid_Segment))
    TCGA_CNA_purity$Segment_Mean <- as.numeric(as.character(TCGA_CNA_purity$Segment_Mean))
    TCGA_CNA_purity$centromere_start <- as.numeric(as.character(TCGA_CNA_purity$centromere_start))
    TCGA_CNA_purity$Plength <- as.numeric(as.character(TCGA_CNA_purity$Plength))
    TCGA_CNA_purity$Qlength <- as.numeric(as.character(TCGA_CNA_purity$Qlength))
    TCGA_CNA_purity$WC <- ifelse(TCGA_CNA_purity$Start < TCGA_CNA_purity$centromere_start & TCGA_CNA_purity$End > TCGA_CNA_purity$centromere_end, 1, 0)
    TCGA_CNA_purity <- rbind(TCGA_CNA_purity, TCGA_CNA_purity %>% filter(WC == 1) %>% mutate(WC = 2))
    TCGA_CNA_purity$ChromArm <- ifelse(TCGA_CNA_purity$WC == 2, paste(TCGA_CNA_purity$Chromosome, "q", sep ="", collapse=NULL), ifelse(TCGA_CNA_purity$WC == 1, paste(TCGA_CNA_purity$Chromosome, "p", sep ="", collapse=NULL), ifelse(TCGA_CNA_purity$Mid_Segment < TCGA_CNA_purity$centromere_start, paste(TCGA_CNA_purity$Chromosome, "p", sep ="", collapse=NULL), paste(TCGA_CNA_purity$Chromosome, "q", sep ="", collapse=NULL))))
    TCGA_CNA_purity$ChromArm_length <- ifelse(TCGA_CNA_purity$ChromArm %like% "p", TCGA_CNA_purity$Plength, TCGA_CNA_purity$Qlength)
    TCGA_CNA_purity$ChromArm_length <- as.numeric(as.character(TCGA_CNA_purity$ChromArm_length))
    TCGA_CNA_purity$length <- ifelse(TCGA_CNA_purity$WC == 0, TCGA_CNA_purity$length, ifelse(TCGA_CNA_purity$WC == 2, (TCGA_CNA_purity$End - TCGA_CNA_purity$centromere_end), (TCGA_CNA_purity$centromere_start - TCGA_CNA_purity$Start) ))
    TCGA_CNA_purity$length <- as.numeric(as.character(TCGA_CNA_purity$length))
    TCGA_CNA_purity <- TCGA_CNA_purity[!(TCGA_CNA_purity$ChromArm=="13p" | TCGA_CNA_purity$ChromArm=="14p" | TCGA_CNA_purity$ChromArm=="15p" | TCGA_CNA_purity$ChromArm=="21p" | TCGA_CNA_purity$ChromArm=="22p"),]
    ids <- split(TCGA_CNA_purity, TCGA_CNA_purity$id)
    print("success")
    #arm-level loss summary
    print("... preparing chromosome loss summary ...")
    {  
        TCGA_CNA_purity_losses <- subset(TCGA_CNA_purity, TCGA_CNA_purity$Segment_Mean < w)
        TCGA_CNA_purity_losses_chrX <- split(TCGA_CNA_purity_losses, TCGA_CNA_purity_losses$ChromArm)
        TCGA_CNA_purity_losses_chrX <- TCGA_CNA_purity_losses_chrX[sapply(TCGA_CNA_purity_losses_chrX, function(x) dim(x)[1]) > 0]
        TCGA_deletion_frequencies_ALL <- matrix(ncol = 3, nrow = 1)
        TCGA_deletion_frequencies_ALL <- data.frame(TCGA_deletion_frequencies_ALL)
        colnames(TCGA_deletion_frequencies_ALL) <- c('Chromosome', 'del_frequencies', 'avg_segMean_deleted')
        Deletion_summary <- matrix(ncol=4, nrow=1)
        colnames(Deletion_summary) <- c('id', 'avg_segMean_deleted', 'Frac_chrom_deleted', 'ChrArm')
        for (j in 1:length(TCGA_CNA_purity_losses_chrX)) {
            TCGA_deletion_frequencies <- matrix(ncol = 3, nrow = 1)
            colnames(TCGA_deletion_frequencies) <- c('Chromosome', 'del_frequencies', 'avg_segMean_deleted')
            TCGA_CNA_purity_losses_chrJ <- TCGA_CNA_purity_losses_chrX[[j]]
            TCGA_CNA_purity_losses_chrJ_sample <- split(TCGA_CNA_purity_losses_chrJ, TCGA_CNA_purity_losses_chrJ$id)
            xout <- matrix(ncol = 4, nrow = length(TCGA_CNA_purity_losses_chrJ_sample))
            for (i in 1:length(TCGA_CNA_purity_losses_chrJ_sample)) {
                TCGA_CNA_purity_losses_chrJ_sample_patX <- TCGA_CNA_purity_losses_chrJ_sample[[i]]
                TCGA_CNA_purity_losses_chrJ_sample_patX$weightedSegMean <- TCGA_CNA_purity_losses_chrJ_sample_patX$length*TCGA_CNA_purity_losses_chrJ_sample_patX$Segment_Mean
                xout[i, 1] <- TCGA_CNA_purity_losses_chrJ_sample_patX$id[1]
                xout[i, 2] <- sum(TCGA_CNA_purity_losses_chrJ_sample_patX$weightedSegMean)/sum(TCGA_CNA_purity_losses_chrJ_sample_patX$length)
                xout[i, 3] <- sum(TCGA_CNA_purity_losses_chrJ_sample_patX$length)/TCGA_CNA_purity_losses_chrJ_sample_patX$ChromArm_length[1]
                xout <- data.frame(xout)
                xout$X1 <- as.character(xout$X1)
                xout$X2 <- as.numeric(as.character(xout$X2))
                xout$X3 <- as.numeric(as.character(xout$X3))
            }
            xout_sig <- subset(xout, xout$X2 < w & xout$X3 > z)
            xout_sig$X4 <- as.character(xout_sig$X4)
            if(nrow(xout_sig) == 0) xout_sig[1,] <- c(0,0,0,TCGA_CNA_purity_losses_chrJ_sample_patX$ChromArm[1])
            xout_sig$X4 <- TCGA_CNA_purity_losses_chrJ_sample_patX$ChromArm[1]
            colnames(xout_sig) <- c('id', 'avg_segMean_deleted', 'Frac_chrom_deleted', 'ChrArm')
            write.table(xout_sig, file = paste(y, "Chromosome_", TCGA_CNA_purity_losses_chrJ_sample_patX$ChromArm[1], "_deletion_samples", ".txt"), sep="\t", row.names=FALSE, quote = FALSE)
            Deletion_summary <- rbind(Deletion_summary, xout_sig)
            TCGA_deletion_frequencies[1, 1] <- TCGA_CNA_purity_losses_chrJ_sample_patX$ChromArm[1]
            TCGA_deletion_frequencies[1, 2] <- nrow(subset(xout_sig, id != 0))/length(ids)
            TCGA_deletion_frequencies[1, 3] <- mean(xout_sig$avg_segMean_deleted)
            TCGA_deletion_frequencies <- data.frame(TCGA_deletion_frequencies)
            TCGA_deletion_frequencies_ALL <- rbind(TCGA_deletion_frequencies_ALL, TCGA_deletion_frequencies)
        }
        Deletion_summary <- Deletion_summary[-1,]
        TCGA_deletion_frequencies_ALL <- TCGA_deletion_frequencies_ALL[-1,]
        write.table(TCGA_deletion_frequencies_ALL, file = paste(y, "_deletion_frequencies", ".txt",sep=""), sep="\t", row.names=FALSE, quote = FALSE)
        write.table(Deletion_summary, file = paste(y, "_Deletion_summary", ".txt",sep=""), sep="\t", row.names=FALSE,  quote = FALSE)
        
    }
    print("success")
    #arm-level gain summary
    print("... preparing chromosome gain summary ...")
    {
        TCGA_CNA_purity_gains <- subset(TCGA_CNA_purity, TCGA_CNA_purity$Segment_Mean > g)
        TCGA_CNA_purity_gains_chrX <- split(TCGA_CNA_purity_gains, TCGA_CNA_purity_gains$ChromArm)
        TCGA_CNA_purity_gains_chrX <- TCGA_CNA_purity_gains_chrX[sapply(TCGA_CNA_purity_gains_chrX, function(x) dim(x)[1]) > 0]
        
        TCGA_amplification_frequencies_ALL <- matrix(ncol = 3, nrow = 1)
        TCGA_amplification_frequencies_ALL <- data.frame(TCGA_amplification_frequencies_ALL)
        colnames(TCGA_amplification_frequencies_ALL) <- c('Chromosome', 'amp_frequencies', 'avg_segMean_amplified')
        Amplification_summary <- matrix(ncol = 4, nrow = 1)
        colnames(Amplification_summary) <- c('id', 'avg_segMean_amplified', 'Frac_chrom_amplified', 'ChrArm')
        for (j in 1:length(TCGA_CNA_purity_gains_chrX)) {
            TCGA_amplification_frequencies <- matrix(ncol = 3, nrow = 1)
            colnames(TCGA_amplification_frequencies) <- c('Chromosome', 'amp_frequencies', 'avg_segMean_amplified')
            TCGA_CNA_purity_gains_chrJ <- TCGA_CNA_purity_gains_chrX[[j]]
            TCGA_CNA_purity_gains_chrJ_sample <- split(TCGA_CNA_purity_gains_chrJ, TCGA_CNA_purity_gains_chrJ$id)
            xout2 <- matrix(ncol = 4, nrow = length(TCGA_CNA_purity_gains_chrJ_sample))
            
            for (i in 1:length(TCGA_CNA_purity_gains_chrJ_sample)) {
                TCGA_CNA_purity_gains_chrJ_sample_patX <- TCGA_CNA_purity_gains_chrJ_sample[[i]]
                TCGA_CNA_purity_gains_chrJ_sample_patX$weightedSegMean <- TCGA_CNA_purity_gains_chrJ_sample_patX$length*TCGA_CNA_purity_gains_chrJ_sample_patX$Segment_Mean
                xout2[i, 1] <- TCGA_CNA_purity_gains_chrJ_sample_patX$id[1]
                xout2[i, 2] <- sum(TCGA_CNA_purity_gains_chrJ_sample_patX$weightedSegMean)/sum(TCGA_CNA_purity_gains_chrJ_sample_patX$length)
                xout2[i, 3] <- sum(TCGA_CNA_purity_gains_chrJ_sample_patX$length)/TCGA_CNA_purity_gains_chrJ_sample_patX$ChromArm_length[1]
                xout2 <- data.frame(xout2)
                xout2$X1 <- as.character(xout2$X1)
                xout2$X2 <- as.numeric(as.character(xout2$X2))
                xout2$X3 <- as.numeric(as.character(xout2$X3))
            }
            
            xout_sig2 <- subset(xout2, xout2$X2 > g & xout2$X3 > z)
            xout_sig2$X4 <- as.character(xout_sig2$X4)
            if(nrow(xout_sig2) == 0) xout_sig2[1,] <- c(0,0,0,TCGA_CNA_purity_gains_chrJ_sample_patX$ChromArm[1])
            xout_sig2$X4 <- TCGA_CNA_purity_gains_chrJ_sample_patX$ChromArm[1]
            colnames(xout_sig2) <- c('id', 'avg_segMean_amplified', 'Frac_chrom_amplified', 'ChrArm')
            write.table(xout_sig2, file = paste(y, "Chromosome_", TCGA_CNA_purity_gains_chrJ_sample_patX$ChromArm[1], "_amplification_samples", ".txt",sep=""), sep="\t", row.names=FALSE, quote = FALSE)
            Amplification_summary <- rbind(Amplification_summary, xout_sig2)
            TCGA_amplification_frequencies[1, 1] <- TCGA_CNA_purity_gains_chrJ_sample_patX$ChromArm[1]
            TCGA_amplification_frequencies[1, 2] <- nrow(subset(xout_sig2, id != 0))/length(ids)
            TCGA_amplification_frequencies[1, 3] <- mean(xout_sig2$avg_segMean_amplified)
            TCGA_amplification_frequencies <- data.frame(TCGA_amplification_frequencies)
            TCGA_amplification_frequencies_ALL <- rbind(TCGA_amplification_frequencies_ALL, TCGA_amplification_frequencies)
        }
        Amplification_summary <- Amplification_summary[-1,]
        TCGA_amplification_frequencies_ALL <- TCGA_amplification_frequencies_ALL[-1,]
        write.table(TCGA_amplification_frequencies_ALL, file = paste(y, "_amplification_frequencies", ".txt",sep=""), sep="\t", row.names=FALSE, quote = FALSE)
        write.table(Amplification_summary, file = paste(y, "_Amplification_summary", ".txt",sep=""), sep="\t", row.names=FALSE, quote = FALSE)
        
    }
    final_output <- list(Deletion_summary, TCGA_deletion_frequencies_ALL, Amplification_summary, TCGA_amplification_frequencies_ALL)
    final_output[[1]] <- final_output[[1]][final_output[[1]]$id != 0, ]
    final_output[[2]][is.na(final_output[[2]])] <- 0
    final_output[[3]] <- final_output[[3]][final_output[[3]]$id != 0, ]
    final_output[[4]][is.na(final_output[[4]])] <- 0
    return(final_output)
    print("success")
}
