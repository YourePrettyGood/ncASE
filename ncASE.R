#usage: replicate_1_par_1 replicate_1_par_2 ... filter read_length output

library(tidyverse)

args<-commandArgs(TRUE)

#mydist <- function(data, snp1, snp2){
#  coord1 <- data$position[data$variant == snp1]
#  coord2 <- data$position[data$variant == snp2]
#  return(object)
#}

#myweight <- function(data, snp){
#  gene <- data$transcript[data$variant == snp]
#  snp_pos <- data$position[data$variant == snp]
#  all_pos <- data$position[data$transcript == gene & data$position - snp_pos < read_length & data$position - snp_pos > -read_length]
#  s <- sum(1-abs(all_pos-snp_pos)/read_length)
#  return(1/s)
#}

#We assume the input will be in sorted order (of both transcript and position)
snp_weight <- function(x, read_length) {
   dist_mat <- 1 - abs(outer(x$position, x$position, "-"))/read_length
   dist_mat[dist_mat < 0] <- 0
   weights <- 1/colSums(dist_mat)
   return(data.frame(transcript=x$transcript, position=x$position, weight=weights, stringsAsFactors=FALSE))
}

if(length(args) < 5)
   stop("Need at least one pair of replicate files and output filename")
if((length(args)-1) %% 2 != 0)
   stop("Need two files per replicate")

read_length <- as.numeric(args[length(args)-1]) 

max_replicate <- (length(args)-3)/2

transcript_count_list <- vector("list", max_replicate)

for(replicate in 1:max_replicate) {

   cat("Processing replicate", replicate, "\n")

   par1_cols <- c("transcript", "position", "par1ref", "par1alt", "par1ref_count", "par1alt_count")
   par2_cols <- c("transcript", "position", "par2ref", "par2alt", "par2ref_count", "par2alt_count")
   par_col_classes <- c("character", "integer", "character", "character", "integer", "integer")

   par1 <- read.table(args[2*replicate-1], quote="\"", col.names=par1_cols, colClasses=par_col_classes)
   cat(paste(nrow(par1), "par1 SNPs from input file\n"))
   par2 <- read.table(args[2*replicate], quote="\"", col.names=par2_cols, colClasses=par_col_classes)
   cat(paste(nrow(par2), "par2 SNPs from input file\n"))

   #Filter for par[12]ref_count > 0 && par[12]alt_count > 0, and
   # set up a column as a join key, eliminating the redundant transcript
   # and position columns:
   par1 <- par1 %>% filter(par1ref_count * par1alt_count > 0) %>%
      mutate(variant=paste(transcript, position, sep="_")) %>%
      select(-c(transcript, position))
   cat(paste(nrow(par1), "par1 SNPs after non-zero ref*alt count filter\n"))
   #species_1 <- species_1[species_1$V5 * species_1$V6 > 0,]
   par2 <- par2 %>% filter(par2ref_count * par2alt_count > 0) %>%
      mutate(variant=paste(transcript, position, sep="_")) %>%
      select(-c(transcript, position))
   cat(paste(nrow(par2), "par2 SNPs after non-zero ref*alt count filter\n"))
   #species_2 <- species_2[species_2$V5 * species_2$V6 > 0,]

   #Perform an inner join on the sets of SNPs identified by mapping to each
   # parent:
   snps <- inner_join(par1, par2, by="variant")
   cat(paste(nrow(snps), "SNPs after inner join of the two parent mappings\n"))
   #snps <- inner_join(species_1 %>% mutate(match_col=paste(V1,V2,sep="_")), species_2 %>% mutate(match_col=paste(V1,V2,sep="_")),by="match_col")
   #colnames(snps) <- c("transcript","position","species_1_ref","species_1_alt","species_1_species_1","species_1_species_2","variant","V1","V2","species_2_ref","species_2_alt","species_2_species_2","species_2_species_1")

   #Filter for SNPs with reciprocal ref and alt on the two parental
   # transcriptomes (e.g. par1ref == par2alt && par1alt == par2ref:
   #Note: This will filter out sites with 0 support for one of the two
   # alleles, which might be the strongest signals of ASE.
   snps <- snps %>% filter(par1ref == par2alt & par1alt == par2ref) %>%
      separate(variant, into=c("transcript", "position"), sep="_", remove=FALSE) %>%
      mutate(position=as.integer(position)) %>%
      arrange(transcript, position)
   cat(paste(nrow(snps), "SNPs after reciprocal parental allele filter\n"))
   #snps <- snps %>% filter((species_1_ref == species_2_alt | species_2_alt == ".") & (species_1_alt == species_2_ref | species_1_alt == ".") & (species_1_alt != "." | species_2_alt != "."))
   #snps <- snps %>% filter(species_1_ref == species_2_alt & species_1_alt == species_2_ref)
   #snps <- snps %>% select("variant","transcript", "species_1_species_1","species_1_species_2","species_2_species_2","species_2_species_1","position")

   #Set thresholds for proximity of reciprocal mapping counts:
   filter_threshold <- as.numeric(args[length(args)-2])
   #The old threshold+filter combination is equivalent to the following:
   #low_threshold <- (1 - filter_threshold)/(1 + filter_threshold)
   #high_threshold <- (1 + filter_threshold)/(1 - filter_threshold)
   #New thresholds, representing % absolute deviance in either direction:
   low_threshold <- 1 - filter_threshold
   high_threshold <- 1 + filter_threshold
   #Old thresholds for old filters:
   #low_threshold <- 0.5 - filter_threshold/2
   #high_threshold <- 0.5 + filter_threshold/2

   #
   snps <- snps %>% filter(par1ref_count > low_threshold * par2alt_count)
   cat(paste(nrow(snps), "SNPs after high-pass par1ref filter\n"))
   snps <- snps %>% filter(par1ref_count < high_threshold * par2alt_count)
   cat(paste(nrow(snps), "SNPs after adding low-pass par1ref filter\n"))
   snps <- snps %>% filter(par2ref_count > low_threshold * par1alt_count)
   cat(paste(nrow(snps), "SNPs after adding high-pass par1alt filter\n"))
   snps <- snps %>% filter(par2ref_count < high_threshold * par1alt_count)
   cat(paste(nrow(snps), "SNPs after all reciprocal mapping counts filters\n"))
   #snps <- snps[snps$species_1_species_1/(snps$species_1_species_1 + snps$species_2_species_1) > low_threshold
   #                & snps$species_1_species_1/(snps$species_1_species_1 + snps$species_2_species_1) < high_threshold
   #                & snps$species_2_species_2/(snps$species_1_species_2 + snps$species_2_species_2) > low_threshold
   #                & snps$species_2_species_2/(snps$species_1_species_2 + snps$species_2_species_2) < high_threshold,]

   #Calculate a per-SNP weight based on distance in the transcript:
   transcripts <- unique(snps$transcript)
   weight_df <- data.frame(transcript=c(), position=c(), weight=c())
   for (tx in transcripts) {
      weight_df <- rbind(weight_df, snp_weight(snps[snps$transcript == tx,], read_length))
   }
   cat(paste("Calculated SNP weights for", length(transcripts), "transcripts\n"))
   cat(paste("snps$transcript class", class(snps$transcript), "snps$position class", class(snps$position), "weight_df$transcript class", class(weight_df$transcript), "weight_df$position class", class(weight_df$position), "\n"))
   snps <- full_join(snps, weight_df, by=c("transcript", "position"))
   cat(paste("Joined SNP weights into main data.frame for replicate", replicate, "\n"))
   #snps$weight <- sapply(snps$variant,myweight,data=snps)

   #Calculate weighted counts for each transcript:
   transcript_counts <- snps %>% group_by(transcript) %>%
      summarize(par1ref_counts=sum(par1ref_count*weight),
         par1alt_counts=sum(par1alt_count*weight),
         par2ref_counts=sum(par2ref_count*weight),
         par2alt_counts=sum(par2alt_count*weight)) %>%
      transmute(transcript=transcript,
         par1_counts=(par1ref_counts+par2alt_counts)/2,
         par2_counts=(par1alt_counts+par2ref_counts)/2)
   cat(paste("Summarized per-parent transcript counts for", nrow(transcript_counts), "transcripts for replicate", replicate, "\n"))
   #transcript_counts <- transcript_counts %>%
   #   rename_if(is_numeric, function(x, replicate) {paste0("rep", replicate, sep="_")}, replicate)
   colnames(transcript_counts)[2:3] <- c(paste("rep", replicate, "par1", "counts", sep="_"), paste("rep", replicate, "par2", "counts", sep="_"))
   cat("Renamed count columns to include replicate\n")
   transcript_count_list[[replicate]] <- transcript_counts
   #x <- aggregate(snps$species_1_species_1*snps$weight, by=list(snps$transcript), "sum")
   #y <- aggregate(snps$species_2_species_1*snps$weight, by=list(snps$transcript), "sum")
   #w <- aggregate(snps$species_1_species_2*snps$weight, by=list(snps$transcript), "sum")
   #z <- aggregate(snps$species_2_species_2*snps$weight, by=list(snps$transcript), "sum")

   #combined <- data.frame(x$Group.1, x$x, y$x, w$x, z$x)

   #if (replicate == 1){
   #   replicate_counts <- list(data.frame(x$Group.1,x.x=(x$x+y$x)/2,z.x=(w$x + z$x)/2))
   #} else {
   #   replicate_counts <- c(replicate_counts, list(data.frame(x$Group.1,x.x=(x$x+y$x)/2,z.x=(w$x + z$x)/2)))
   #}
}


#temp <- replicate_counts[[1]]$x.Group.1

#if (max_replicate > 1)
#   for(replicate in 2:max_replicate) {
#      temp <- replicate_counts[[replicate]]$x.Group.1[replicate_counts[[replicate]]$x.Group.1 %in% temp]
#   }

#for(replicate in 1:max_replicate){
#  i <- replicate_counts[[replicate]][replicate_counts[[replicate]]$x.Group.1 %in% temp,]

#  i <- i[order(as.character(i$x.Group.1)),]
#  if(replicate==1){ 
#    combined_i <- data.frame(i$x.Group.1, i$x.x, i$z.x)
#  }else{
#  combined_i <- data.frame(combined_i, i$x.x, i$z.x)
#  }
#}


#labels <- c("transcript")
#for(replicate in 1:max_replicate){
#  labels <- c(labels, paste("replicate",replicate,"species_1",sep="_"),paste("replicate",replicate,"species_2",sep="_"))
#}

#names(combined_i) <- labels

#write.table(combined_i, file=args[length(args)], quote=FALSE, sep="\t", row.names=FALSE) 

#Join all the transcript count tables together and output:
i <- 1
final_transcript_counts <- transcript_count_list[[i]]
cat(paste(nrow(final_transcript_counts), "transcripts from replicate", i, "\n"))
i <- i + 1
while (i <= max_replicate) {
   final_transcript_counts <- final_transcript_counts %>% full_join(transcript_count_list[[i]], by="transcript")
   cat(paste(nrow(final_transcript_counts), "transcripts after joining with replicate", i, "\n"))
   i <- i + 1
}

write.table(final_transcript_counts, file=args[length(args)], quote=FALSE, sep="\t", row.names=FALSE)
