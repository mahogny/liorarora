library(sqldf)

dat <- read.csv("tss.csv",sep="\t",stringsAsFactors = FALSE)
colnames(dat)<-c("gene","start","end","chromosome")
dat <- dat[dat$chromosome %in% c(sprintf("%s",1:19),"X","Y"),]
dat$chromosome <- sprintf("chr%s",dat$chromosome)

dat <- sqldf("select gene, min(start) as start, max(end) as end, chromosome from dat group by gene")

out <- data.frame(chrom=dat$chromosome, start=dat$start-5000, end=dat$start, name=dat$gene)
write.table(out, "out/upstream.bed",sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

out <- data.frame(chrom=dat$chromosome, start=dat$end, end=dat$end+5000, name=dat$gene)
write.table(out, "out/downstream.bed",sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

#out$end-out$start
