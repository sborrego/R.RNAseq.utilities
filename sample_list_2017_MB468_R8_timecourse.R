# Making a sample_list.txt
# see README for more information

# UNIX
# cat *.txt | sed 's/[#;,]//g' > combo_stats.txt

working.dir <- "~/Data/RNASeq"
sample.list.dir <- "~/Data/RNASeq/sample_list"

setwd(working.dir)

combo_stats6 <- read.delim("~/Data/RNASeq/4R044-L6-Statistics/combo_stats6.txt",
                          header=FALSE)
combo_stats7 <- read.delim("~/Data/RNASeq/4R044-L7-Statistics/combo_stats7.txt",
                          header=FALSE)

combo_stats6$V2 <- sub("(\\w+).*", "\\1", combo_stats6$V1)
combo_stats6$V0 <- sub("\\w+", "\\1", combo_stats6$V1)
combo_stats7$V2 <- sub("(\\w+).*", "\\1", combo_stats6$V1)
combo_stats7$V0 <- sub("\\w+", "\\1", combo_stats7$V1)

df6 <- data.frame(cat=combo_stats6$V2,
                  data=combo_stats6$V0)
df7 <- data.frame(cat=combo_stats7$V2,
                  data=combo_stats7$V0)
df <- rbind(df6, df7)

library(dplyr)
files <- filter(df, cat=="Files") %>% select(data)
barcode <- filter(df, cat=="Barcode") %>% select(data)
reads <- filter(df, cat=="Reads") %>% select(data)
lib <- filter(df, cat=="Library") %>% select(data)
cycles <- filter(df, cat=="Cycles") %>% select(data)

original_sample_list <- data.frame(files, barcode, reads, lib, cycles)
colnames(original_sample_list) <- c("File", "Barcode",
                                    "Reads", "Library", "Cycles")

setwd(sample.list.dir)
write.csv(as.data.frame(original_sample_list),
          file="original_sample_list_2017_MB468_R8_timecourse.csv")

combo_list <- read.delim("sample_list.txt", header=TRUE)

filename <- combo_list %>% filter(grepl("combo*", combo_list[,1]))
total_reads1 <- combo_list %>% filter(!grepl("combo*", combo_list[,1]))
total_reads2 <- as.numeric(as.vector(unlist(total_reads1)))
total_reads <- total_reads2 / 4

l6_reads <- as.numeric(as.vector(unlist(original_sample_list[1:24, "Reads"])))
l7_reads <- as.numeric(as.vector(unlist(original_sample_list[26:49,"Reads"])))
lane_read_sum <-  l6_reads + l7_reads

sample_list <- data.frame(Filename=filename,
                          Total_Reads=total_reads,
                          Lane_Reads_Sums=lane_read_sum,
                          L6_sample_name=original_sample_list[1:24,
                                                              "File"],
                          L6_reads=l6_reads,
                          L7_sample_name=original_sample_list[26:49,
                                                              "File"],
                          L7_reads=l7_reads,
                          Barcode=original_sample_list[1:24, "Barcode"])

cell <- c("MB468", "R8")
minute <- c("000", "030", "120", "720")
rep <- c("01", "02", "03")

sample_names <- list()
sample_cell <- list()
sample_time <- list()
sample_rep <- list()

for (MIN in minute){
    for (CELL in cell){
        for (REP in rep){
            samples <- paste(CELL, MIN, REP, sep="_")
            sample_names[[samples]] <- samples
        }
    }
}

sample_list$Sample.name <- unlist(sample_names)
sample_list$Cell.line <- rep(cell, each=3)
sample_list$Time <- rep(minute, each=6)
sample_list$Replicate <- rep(rep, each=1)

colnames(sample_list)[1] <- "Filename"
refcols <- c("Filename", "Sample.name", "Cell.line", "Time", "Replicate")
sample_list <- sample_list[, c(refcols, setdiff(names(sample_list), refcols))]

write.csv(as.data.frame(sample_list),
          file="sample_list_2017_MB468_R8_timecourse.csv")
