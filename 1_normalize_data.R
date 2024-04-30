###
#
# Original script: https://github.com/cbmporter/CandidaKOanalysis
# This is an R script to process, analyze, and visualize genetic knockouts in C. albicans experiment
# Used in the paper "A CRISPR Cas9-based gene drive platform for genetic interaction analysis in Candida albicans? by Shapiro et al.
# Edited by: Myra Paz Masinas, Boone Lab, University of Toronto, November 2023
#
###


# load libraries
rm(list = ls())
library(amap)
library(gplots)
library("extrafont")
library("optparse")

# script arguments
option_list = list(
	make_option(c("-d", "--data"), type="character", default=NULL, help="Directory where the data is located", metavar="character"),
    make_option(c("-f", "--files"), type="character", default=NULL, help="Directory where the supplemental files are located", metavar="character"),
	make_option(c("-o", "--output"), type="character", default=NULL, help="Directory where to save output files", metavar="character"),
	make_option(c("-r", "--replicate"), type="integer", default=3, help="Number of replicates. Default is 3.", metavar="number"),
	make_option(c("-w", "--wildtype"), type="integer", default=2, help="Number of wildtype matrix. Default is 2.", metavar="number")
);

# parse arguments 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

data_dir <- opt$data
files_dir <- opt$files
output_dir <- opt$output
num_rep <- opt$replicate
num_WT <- opt$wildtype

# load and process index for data - see example index file 
setwd(files_dir)

geneIndex = read.table('geneIndex.csv', sep=",")
lowerIndex = as.matrix(geneIndex[lower.tri(geneIndex, diag=TRUE)])
upperIndex = as.matrix(geneIndex[upper.tri(geneIndex, diag=TRUE)])
lowerIndex_i = lowerIndex
upperIndex_i = upperIndex

# load index to gene name mapping
geneNames <- read.table('geneNames.csv', sep=",")

# reverse the string order of the name for upper genes
strReverse <- function(x)
        sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")
upperIndex = strReverse(upperIndex)

# keep string order for lower genes
strCopy <- function(x)
        sapply(lapply(strsplit(x, NULL), paste), paste, collapse="")
lowerIndex = strCopy(lowerIndex)

# get index of aligned gene combinations 
sortedLower = sort(lowerIndex, index.return=TRUE)
sortedUpper = sort(upperIndex, index.return=TRUE)

# load and process data: should be a folder with .csv files of 
# OD measurements for each condition; see example file
setwd(data_dir)

conditions <- list.files()

# define variables 
noReps <- num_rep * 2 
noGenes <- nrow(geneNames)
noCond <- length(conditions)
noStrains <- length(lowerIndex)

# initialize matrices
allData <- matrix(data=NA, nrow=noStrains, ncol=noReps*noCond)
names <- matrix(data=NA,nrow=1,ncol=noReps*noCond)
avgData <- matrix(data=NA, nrow=noStrains, ncol=noCond)
avgLower <- matrix(data=NA, nrow=noStrains, ncol=noCond)
avgUpper <- matrix(data=NA, nrow=noStrains, ncol=noCond)

# assign data to appropratiate matrices 
# colums will be conditions, rows will be strains
noWT <- num_WT
for (i in 1:noCond){
        # pull out data for current condition 
        currentData <- read.table(conditions[i], sep=",")
        # separate replicates 
        controls <- currentData[,(noGenes+1):ncol(currentData)]
        knockouts <- currentData[,1:noGenes]

        # save only rows from file that contain data 
        #knockouts <- knockouts[complete.cases(knockouts),]
		is_finite = lapply(knockouts, is.finite)
		is_nan = lapply(knockouts, is.nan)
		is_finite_or_nan <- append(is_finite, is_nan)
		knockouts_is_finite_or_nan <- Reduce("|", is_finite_or_nan)
		knockouts <- knockouts[knockouts_is_finite_or_nan,]
        controls <- controls[complete.cases(controls),]
        
        # initialize temporary matrices to store current data 
        allCurrentData <- matrix(data=NA, nrow=noStrains, ncol=noReps)
        currentNames <- matrix(data=NA, nrow=1, ncol=noReps)
		currentCondition <- gsub(".csv", "", conditions[i])
        for (j in 1:(noReps/2)){
                # pull out replicate data 
                repData <- knockouts[(j*noGenes-(noGenes-1)):(j*noGenes),]
				# pull out control data and take average
                ctrlData <- controls[(j*noWT-noWT+1):(j*noWT),]
                ctrlAvg <- mean(as.matrix(ctrlData), na.rm=TRUE)
				# pull out data for forward and reverse KO cases 
                lower <- as.matrix(repData[lower.tri(repData, diag = TRUE)])
                upper <- as.matrix(repData[upper.tri(repData, diag = TRUE)])
				# average strain matches
                orderedLower <- as.matrix(lower[sortedLower$ix])
                orderedUpper <- as.matrix(upper[sortedUpper$ix])
                allCurrentData[,(j*2-1)] <- orderedLower/ctrlAvg
                allCurrentData[,(j*2)] <- orderedUpper/ctrlAvg
				# set column names with condition, replicate and lower/upper information
				currentNames[(j*2-1)] <- paste0(currentCondition, '_R', j, '_lower')
				currentNames[(j*2)] <- paste0(currentCondition, '_R', j, '_upper')
				
        }
        # fill final matrices and build average data matrices
        allData[,((i*noReps-(noReps-1)):(i*noReps))] <- allCurrentData
        avgData[,i] <- rowMeans(allCurrentData, na.rm=TRUE)
        names[,((i*noReps-(noReps-1)):(i*noReps))] <- currentNames
        avgLower[,i] <- rowMeans(orderedLower/ctrlAvg, na.rm=TRUE)
        avgUpper[,i] <- rowMeans(orderedUpper/ctrlAvg, na.rm=TRUE)
}


# name rows and columns 
geneList <- sortedLower$x
geneList_i <- geneList
singles <- diag(as.matrix(geneIndex))
singlesIndex <- match(singles, geneList)
geneList <- paste0(paste0('-', geneList),'-')
for (i in 1:nrow(geneNames)){
        geneList <- gsub(paste0('^-', geneNames[i,1]), paste0(geneNames[i,2], '-'), geneList)
        geneList <- gsub(paste0(geneNames[i,1], '-$'), geneNames[i,2], geneList)
}
geneList <- gsub('-', '', geneList)

# rename columns and rows
avgData <- as.matrix(avgData)
avgLower <- as.matrix(avgLower)
avgUpper <- as.matrix(avgUpper)
rownames(allData)<-geneList
colnames(allData)<-names
rownames(avgData)<-geneList
colnames(avgData)<-gsub(".csv", "", conditions)
rownames(avgLower)<-sortedLower$x
colnames(avgLower)<-gsub(".csv", "", conditions)
rownames(avgUpper)<-upperIndex_i[sortedUpper$ix]
colnames(avgUpper)<-gsub(".csv", "", conditions)

# remove empty columns (for conditions that had fewer reps)
completeConditions <- conditions

# save data to csv files
setwd(output_dir)
write.csv(allData, "allData.csv")
write.csv(avgData, "avgData.csv")

# make heatmaps of "raw" data in plate shape
triangleHM <- function(triangleData, triangleType, colors, extension){
        tDataCol <- ncol(triangleData)
        for (i in 1:tDataCol){
                triangleGeneIndex <- geneIndex
                rownames(triangleGeneIndex) <- geneNames[,2]
                colnames(triangleGeneIndex) <- geneNames[,2]
                
                # reverse the gene index place holder for half triangle before sorting
                if (triangleType=="half"){}
                
                # sort the data
                sortRow <- sort(rownames(triangleGeneIndex), index.return=TRUE)
                sortCol <- sort(colnames(triangleGeneIndex), index.return=TRUE)
                triangleGeneIndex <- triangleGeneIndex[sortRow$ix, sortCol$ix]
                
                # replace upper trinagle with NaN
                if (triangleType=="half"){
                        triangleGeneIndex[upper.tri(triangleGeneIndex)]<-NaN
                }
                triangleGeneIndex <- as.matrix(triangleGeneIndex)
                
                
                # replace letters with data values
                tDataCol <-ncol(triangleData)
                for (j in 1:nrow(triangleData)){
                        triangleGeneIndex <- gsub(paste0(paste0('^', rownames(triangleData)[j]),'$'),triangleData[j,i],triangleGeneIndex)
                }
                
                if (triangleType=="half"){
                        tDataCol <-ncol(triangleData)
                        for (j in 1:nrow(triangleData)){
                                triangleGeneIndex <- gsub(paste0(paste0('^', rownames(avgUpper)[j]),'$'),triangleData[j,i],triangleGeneIndex)
                        }
                }
                
                
                # convert from character to double for heatmap 
                triangleGeneIndex <- apply(triangleGeneIndex, c(1,2), as.numeric)
                
                # generate pdf of heatmap (no clustering)
                pdf(paste0(paste0(paste0(colnames(triangleData)[i], extension), triangleType), '.pdf'),
                    width = 6, height = 6, family="ArialMT", useDingbats=FALSE)
                heatmap.2(triangleGeneIndex, Rowv=FALSE, Colv=FALSE, trace="none", 
                          dendrogram="none", keysize=1.1, key.title='', key.xlab='',
					      cexRow=0.6, cexCol=0.6,
                          main=colnames(triangleData)[i], col=my_palette, breaks=colors,
                          density.info='none', na.color='gray80')
                dev.off()
                        
        }
}

avgData <- as.matrix(avgData)

# triangle HM of OD data, half
tData <- avgData
rownames(tData) <- geneList_i

# remove outliers if they exist
setwd(files_dir)
nan_lower_file <- read.csv("nan_lower_gene_index.csv", header = FALSE)
nan_upper_file <- read.csv("nan_upper_gene_index.csv", header = FALSE)
nan_lower <- as.list(nan_lower_file[[1]])
nan_upper <- as.list(nan_upper_file[[1]])
outlier_upper <- c(which(geneList_i %in% nan_upper))
outlier_lower <- c(which(geneList_i %in% nan_lower))
avgUpper[outlier_upper,]<-NaN
avgLower[outlier_lower,]<-NaN

# triangle HM of OD data, full
setwd(output_dir)
tData <- rbind(avgLower, avgUpper)
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 48)
colors = unique(c(seq(min(tData, na.rm=TRUE),1,length=25),seq(1,max(tData, na.rm=TRUE),length=25)))
triangleHM(tData, "full", colors, "_OD")
print(paste("Done saving data and plotting OD heatmap"))

