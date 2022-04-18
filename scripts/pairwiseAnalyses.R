# 8/31/17

# TO DO:
# in SNP file, there are still too many SNPs, e.g. at chrI 2146
# need to add the regional block for ends of chromosomes
# in SNP file, the SNPs are out of order, which makes me worried about the filters. Some indexing problem?
# SNP file also has chromosomes out of order. Not a big deal, but will need to take into account when plotting
# FIX multipeaks dump problem
# FIX mclapply


# can then call this from external
i = 1
alignmentDir <- "~/Desktop/data/illumina/np2us"
SNPs <- read.table("~/Desktop/data/illumina/np2us/SNPs_Maggie_170809_BY_positions.txt", stringsAsFactors=FALSE, head=FALSE)
# see comments above. As of 8/31/17, the SNPs seem not to be fully filtered, and are out of sorting order
for (thisChr in unique(SNPs[,1])){SNPs[SNPs[,1] == thisChr, 2] <- sort(SNPs[SNPs[,1] == thisChr, 2])}
SNPs <- rbind(SNPs[SNPs[,1] == "chrI",], SNPs[SNPs[,1] == "chrII",], SNPs[SNPs[,1] == "chrIII",], SNPs[SNPs[,1] == "chrIV",], SNPs[SNPs[,1] == "chrV",], SNPs[SNPs[,1] == "chrVI",], SNPs[SNPs[,1] == "chrVII",], SNPs[SNPs[,1] == "chrVIII",], SNPs[SNPs[,1] == "chrIX",], SNPs[SNPs[,1] == "chrX",], SNPs[SNPs[,1] == "chrXI",], SNPs[SNPs[,1] == "chrXII",], SNPs[SNPs[,1] == "chrXIII",], SNPs[SNPs[,1] == "chrXIV",], SNPs[SNPs[,1] == "chrXV",], SNPs[SNPs[,1] == "chrXVI",])

experimentFile <- read.table("~/Desktop/data/illumina/np2us/comparison-table-short3.txt", stringsAsFactors=FALSE, head=TRUE)
resultsFolder <- "~/Desktop/data/illumina/np2us/results"

withMultipool <- TRUE


# common annotations, functions, etc

library("VariantAnnotation")
source("~/Desktop/data/illumina/np2us/gTest.R")
source("~/Desktop/data/illumina/np2us/x_qtl_seq_functions_170831.R")
source("~/Desktop/data/illumina/np2us/mp_JB_170901.R")
source("~/Desktop/data/illumina/np2us/peaksFromVector.R")
geneInfo = read.table("~/Desktop/data/illumina/np2us/ensemblGenes_ensembl83_160307_MOD.txt", stringsAsFactors=FALSE, sep="\t", header=TRUE)
rownames(geneInfo) <- geneInfo[,"geneID"]
allNames <- geneInfo[, "geneName"]
names(allNames) <- geneInfo[,1]
allNames[which(allNames == "")] <- names(allNames)[which(allNames == "")]
allNamesInv <- names(allNames)
names(allNamesInv) <- allNames

sepBetweenChr <- 1e5
trimFromEnd = 15e3
obsMin <- 10
LoessSpan = 0.1
AFThres = 0.09653124 # same as in Albert 2014
multiThres = 4.5 # LOD threshold ofr multipool, if run with N=1000



# done reading annotation files, now get to work

for (i in (1:nrow(experimentFile))){   

plotName <- paste0(experimentFile[i, ], collapse="_")
mainPlotName <- paste0(experimentFile[i, ], collapse=" ")
thisGene <- experimentFile[i, "Gene"]


highFile <- dir(alignmentDir, pattern=paste0(experimentFile[i, "high"], ".*vcf$"), full.names=TRUE)
lowFile <- dir(alignmentDir, pattern=paste0(experimentFile[i, "low"], ".*vcf$"), full.names=TRUE)
highVCF <- readVcf(highFile)
lowVCF <- readVcf(lowFile)
highVCFCounts <- t(sapply(info(highVCF)$AD, function(x){c(x[[1]], x[[2]])}))
lowVCFCounts <- t(sapply(info(lowVCF)$AD, function(x){c(x[[1]], x[[2]])}))
rownames(highVCFCounts) <- sapply(names(rowRanges(highVCF)), function(x){strsplit(x, "_")[[1]][1]})
rownames(lowVCFCounts) <- sapply(names(rowRanges(lowVCF)), function(x){strsplit(x, "_")[[1]][1]})
# these counts were read from the vcf, so they don't contain counts for sites without reads
# need to build a common SNP table
theseCounts <- data.frame(SNPs, matrix(0, nrow=nrow(SNPs), ncol=4))
rownames(theseCounts) <- paste(theseCounts[,1], theseCounts[,2], sep=":")
colnames(theseCounts) <- c("chr", "pos", "high_ref", "high_alt", "low_ref", "low_alt")
theseCounts[rownames(highVCFCounts), c("high_ref", "high_alt")] <- highVCFCounts
theseCounts[rownames(lowVCFCounts), c("low_ref", "low_alt")] <- lowVCFCounts

gcoords= getGcoords(theseCounts$chr, theseCounts$pos, sepBetweenChr)
names(gcoords) = rownames(theseCounts)

# get plot coordinates for the chromosome line dividers
chrCutoffs <- sapply(unique(theseCounts$chr), function(x){gcoords[theseCounts$chr == x][1] - sepBetweenChr/2})
names(chrCutoffs) <- unique(theseCounts$chr)
chrLabels = sapply(1:(length(chrCutoffs)-1), function(i)(chrCutoffs[i] + chrCutoffs[i+1])/2)
# add half the length of chrXVI
chrLabels = c(chrLabels, chrCutoffs[16] + sepBetweenChr + 948066/2)
names(chrLabels)[16] = "chrXVI"


# write out the counts to be analyzed
save(theseCounts, file=paste(paste0(experimentFile[i, c(1, 3, 4, 5, 6)], collapse="_"), ".RData", sep=""))

# compute BY allele frequencies
theseCounts$highBYAF <- theseCounts$high_ref / (theseCounts$high_ref + theseCounts$high_alt)
theseCounts$lowBYAF <- theseCounts$low_ref / (theseCounts$low_ref + theseCounts$low_alt)
theseCounts$highCoverage <- theseCounts$high_ref + theseCounts$high_alt
theseCounts$lowCoverage <- theseCounts$low_ref + theseCounts$low_alt

#theseCounts <- cbind(theseCounts,
#    t(apply(theseCounts[,c("high_ref", "high_alt", "low_ref", "low_alt")], 1, function(x){
#        ret = c(NA, NA)
#        testMatrix = cbind(x[1:2], x[3:4])
#        if (sum(testMatrix) > 0){
#            gT = g.test(testMatrix)
#            lgP = -log10(gT$p.value)
#            lgP[lgP == Inf] <- max(lgP[lgP != Inf], na.rm=T)
#            lgP[lgP == -Inf] <- min(lgP[lgP != -Inf], na.rm=T)
#            ret = c(lgP, gT$statistic)
#        }
#        names(ret) = c("minusLog10PValue", "G")
#        ret
#    })))

medianCov = apply(theseCounts[,c("highCoverage", "lowCoverage")], 2, median)
SNPsAtMinObs = length(theseCounts[,"highCoverage"] >= obsMin & theseCounts[,"lowCoverage"] >= obsMin)



##############
# now we're ready to plot
pdf(file = paste(resultsFolder, plotName, "_all_raw.pdf", sep=""), width=11, height=8)
for(thisPop in c("high", "low")){
    plot(gcoords, rep(0.5, length(gcoords)), ylim=c(0,1), main = paste(mainPlotName, thisPop, collapse=" "), xaxt='n', xlab="chromosome", ylab="BY allele frequency", type="n")
    
    # CAN1 is on chr05, pos 33466
    # MATALPHA is on chr03, 198671
    abline(v = getGcoords("chrV", 33466, sepBetweenChr), lwd=2, col="green")
    abline(v = getGcoords("chrIII", 198671, sepBetweenChr), lwd=2, col="green")
    abline(v = getGcoords(geneInfo[thisGene, "chr"], mean(as.numeric(geneInfo[thisGene, c("start", "end")])), sepBetweenChr), lwd=2, col="purple")

    points(gcoords, theseCounts[,paste0(thisPop, "BYAF", collapse="")], cex=.2, col="#00000022")
    
    roll = rollLoessByChrWithWeights(data.frame(theseCounts[,"chr"], theseCounts[,paste0(thisPop, "BYAF", collapse="")], gcoords, median(theseCounts[,paste0(thisPop, "Coverage", collapse="")]), stringsAsFactors=FALSE), LoessSpan)
    names(roll) = rownames(theseCounts)
    
    for (j in unique(theseCounts[,"chr"])){
        points(gcoords[theseCounts$chr == j], roll[theseCounts$chr == j], type="l", lwd=2)
    }
    
    for (j in chrCutoffs){
        abline(v = j, lty=2, col="light blue")
    }
    abline(h = 0.5, lty=2, col="light blue")
    
    legend("topright", legend = c(paste0("median coverage: ", median(theseCounts[,paste0(thisPop, "Coverage", collapse="")]))), box.lty=0)
    axis(side=1, at = chrLabels, labels = as.roman(1:16), tick=F)
}
dev.off()




# plot the BY-RM lines for the experiments in one plot
pdf(file = paste(resultsFolder, plotName, "_all_combined.pdf", sep=""), width=11, height=8)
plot(gcoords, rep(0.5, length(gcoords)), ylim=c(0,1), cex=.2, col="grey", main = mainPlotName, type="n", xaxt='n', xlab="chromosome", ylab="BY allele freq")

abline(v = getGcoords(geneInfo[thisGene, "chr"], mean(as.numeric(geneInfo[thisGene, c("start", "end")])), sepBetweenChr), lwd=2, col="purple")

cols <- cbind(c("red", "blue"), c("#FF000022", "#0000FF22"))
rownames(cols) <- c("high", "low")
colnames(cols) <- c("lines", "points")

for (thisPop in c("high", "low")){
    points(gcoords, theseCounts[,paste0(thisPop, "BYAF", collapse="")], cex=.2, col=cols[thisPop, "points"])
}

for(thisPop in c("high", "low")){
    roll = rollLoessByChrWithWeights(data.frame(theseCounts[,"chr"], theseCounts[,paste0(thisPop, "BYAF", collapse="")], gcoords, median(theseCounts[,paste0(thisPop, "Coverage", collapse="")]), stringsAsFactors=FALSE), LoessSpan)
    for (j in unique(theseCounts[,"chr"])){
        points(gcoords[theseCounts[,"chr"] == j], roll[theseCounts[,"chr"] == j], type="l", lwd=2, col=cols[thisPop, "lines"])
    }
}
for (j in chrCutoffs){
    abline(v = j, lty=2, col="light blue")
}
abline(h = 0.5, lty=2, col="light blue")

#legend("topleft", legend=c("high tail", "low tail"), col=cols, lty=1, lwd=2, box.lty=0)    legend("topright", legend = c(paste0("median coverage: ", medianCov[covDict[thisPop]]), box.lty=0))
axis(side=1, at = chrLabels, labels = as.roman(1:16), tick=F)

dev.off()


# difference plot
pdf(file = paste(resultsFolder, plotName, "_difference.pdf", sep=""), width=11, height=8)

rollHigh = rollLoessByChrWithWeights(data.frame(theseCounts[,"chr"], theseCounts[,"highBYAF"], gcoords, median(theseCounts[,"highCoverage"]), stringsAsFactors=FALSE), LoessSpan)
rollLow = rollLoessByChrWithWeights(data.frame(theseCounts[,"chr"], theseCounts[,"lowBYAF"], gcoords, median(theseCounts[,"lowCoverage"]), stringsAsFactors=FALSE), LoessSpan)


plot(gcoords, rollHigh - rollLow, cex=.2, col="grey", main = mainPlotName, type="n", xaxt='n', xlab="chromosome", ylab="High - Low Population", ylim=c(-1,1))
abline(v = getGcoords(geneInfo[thisGene, "chr"], mean(as.numeric(geneInfo[thisGene, c("start", "end")])), sepBetweenChr), lwd=2, col="purple")

points(gcoords, theseCounts[,"highBYAF"] - theseCounts[,"lowBYAF"], cex=.2, col="#00000022")

# these difference thresholds from the nullSorts:
abline(h = AFThres, lty = 2, lwd=2, col = "red")
abline(h = -AFThres, lty = 2, lwd=2, col = "red")

for (j in unique(SNPs[,1])){
    points(gcoords[theseCounts[,"chr"] == j], (rollHigh - rollLow)[theseCounts[,"chr"] == j], type="l", lwd=2, col="black")
}

for (j in chrCutoffs){
    abline(v = j, lty=2, col="light blue")
}
abline(h = 0, lty=2, col="light blue")

axis(side=1, at = chrLabels, labels = as.roman(1:16), tick=F)

dev.off()





#####################
# multipool from here on

# make a switch to kill plotting etc here if no multipool desired
# if(!withMultipool){break()} # not actually sure this would work - need to test

# run multipool
# per chromosome!

# TURNS OUT THIS IS NOT A GOOD IDEA ON OUR SERVER USING MCLAPPLY: IT DOES NTO TERMINATE THE CHILD PROCESSES AND THEY SIT THERE AND TAKE UP MEMORY!
# were able to kill like this: for pid in $(ps -ef | grep "R_PATH_HERE" | awk '{print $2}'); do kill -9 $pid; done
# further debugging shows that the problem is not with screen, and not with the R function we call. looks like a mclapply problem on our server?
#multipoolOutput <- mclapply(unique(theseCounts[,"chr"]), function(j){
#    doMultiPoolFromWithinR(theseCounts[theseCounts$chr == j, c("pos", "high_ref", "high_alt")], theseCounts[theseCounts$chr == j, c("pos", "low_ref", "low_alt")])
#}, mc.cores=2)
multipoolOutput <- lapply(unique(theseCounts[,"chr"]), function(j){
    doMultiPoolFromWithinR(theseCounts[theseCounts$chr == j, c("pos", "high_ref", "high_alt")], theseCounts[theseCounts$chr == j, c("pos", "low_ref", "low_alt")])
})


# call peaks from multipool LODs
multiPeaks <- lapply(multipoolOutput, function(x){
    thisLODTrace = x[[2]]
    theseChrPeaks = callPeaks(thisLODTrace[,2], 4.5, 2)
    theseChrPeaks
})

save(multiPeaks, multipoolOutput, file=paste(paste0(experimentFile[i, c(1, 3, 4, 5, 6)], collapse="_"), "_multipoolResults.RData", sep=""))
zz <- file(paste(paste0(experimentFile[i, c(1, 3, 4, 5, 6)], collapse="_"), "_multipoolPeaks.txt", sep=""), open="wt")
sink(zz)
multiPeaks
sink()


# plot difference over MulitpoolLOD
pdf(paste(resultsFolder, plotName, "_difference_withMultipool.pdf", sep=""), width=11, height=11)
par(mfrow=c(2,1))

plot(gcoords, rep(0, length(gcoords)), cex=.2, col="grey", main = mainPlotName, type="n", xaxt='n', xlab="chromosome", ylab="High - Low Population", ylim=c(-1,1))
points(gcoords, theseCounts[,"highBYAF"] - theseCounts[,"lowBYAF"], cex=.2, col="#00000022")

abline(v = getGcoords(geneInfo[thisGene, "chr"], mean(as.numeric(geneInfo[thisGene, c("start", "end")])), sepBetweenChr), lwd=2, col="purple")

# difference thresholds from the Albert 2014 nullSorts:
abline(h = AFThres, lty = 2, lwd=2, col = "red")
abline(h = -AFThres, lty = 2, lwd=2, col = "red")

for (j in unique(SNPs[,1])){
    points(gcoords[theseCounts[,"chr"] == j], (rollHigh - rollLow)[theseCounts[,"chr"] == j], type="l", lwd=2, col="black")
}

for (j in chrCutoffs){
    abline(v = j, lty=2, col="light blue")
}
abline(h = 0, lty=2, col="light blue")

axis(side=1, at = chrLabels, labels = as.roman(1:16), tick=F)


ylimMax = max(c(multiThres, sapply(multipoolOutput, function(x){max(x[[2]][,2])}))) + 1
plot(gcoords, rep(0, length(gcoords)), main = mainPlotName, type="n", xaxt='n', xlab="chromosome", ylab="Multipool LOD", ylim=c(0, ylimMax))

abline(v = getGcoords(geneInfo[thisGene, "chr"], mean(as.numeric(geneInfo[thisGene, c("start", "end")])), sepBetweenChr), lwd=2, col="purple")

for (j in 1:16){
    points(getGcoords(paste0("chr", as.roman(j)), multipoolOutput[[j]][[2]][,1], sepBetweenChr), multipoolOutput[[j]][[2]][,2], type="l", lwd=2, col="black")
    # add stars at peaks
    if(!is.null(multiPeaks[[j]])){
        for (thisPeak in 1:nrow(multiPeaks[[j]])){
            thisPeakPlotPos <- getGcoords(j, multiPeaks[[j]][thisPeak, "maxIndex"], sepBetweenChr)
            text(thisPeakPlotPos, multiPeaks[[j]][thisPeak, "maxValue"] + 0.5, labels="*", col="red", cex=5)
        }
    }

}


abline(h = 0, lty=2, col="light blue")
for (j in chrCutoffs){
    abline(v = j, lty=2, col="light blue")
}

axis(side=1, at = chrLabels, labels = as.roman(1:16), tick=F)
abline(h = multiThres, lty=2, col="red", lwd=2)

dev.off()

}


#######
# ADD:
# overplot eQTL & earlier XQTL with these results
# load big table with all those results; if gene is present, fish out and plot






