##### Illustration with NVPA

library(xcms)
library(ggplot2)
#library(mixOmics)
devtools::install_github("Beirnaert/MetaboMeeseeks" )
#library(MetaboMeeseeks)
library(purrr)
library(RColorBrewer)

rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)

## selection of datadirectories (change accordingly)
datadir_lipid_neg_acute = "/Users/Charlie/Dropbox/PhD/Projects/TOXICO_matthias/NVPA_Toxico_biomarker/data/lipidomics/negative/acute/"
datadir_lipid_neg_chronic = "/Users/Charlie/Dropbox/PhD/Projects/TOXICO_matthias/NVPA_Toxico_biomarker/data/lipidomics/negative/chronic/"

## merging datadirectories in list (change accordingly)
DataDirectories = list(acute = datadir_lipid_neg_acute, 
                       chronic = datadir_lipid_neg_chronic)

## getting individual datafiles in the directories 
Datafiles <- map(.x = DataDirectories,
                 .f = ~ list.files(.x, 
                                   recursive = TRUE, 
                                   full.names = TRUE)
)





## XCMS peak picking
ppm_threshold=20
xcmsSets = map(.x = Datafiles[1:2], 
                   .f = ~ xcmsSet(files=.x, 
                                  method="centWave", 
                                  ppm = ppm_threshold, 
                                  peakwidth = c(5,30), 
                                  mzdiff = -0.01, 
                                  prefilter = c(3,5000),
                                  integrate = 2, 
                                  snthresh = 6, 
                                  noise = 3000, 
                                  fitgauss = TRUE, 
                                  verbose.columns = T)
)


## setting batch and class variable names (just to get some info on the xcms names)
walk(.x = xcmsSets,
     .f = ~ print(head(.x@phenoData))
)

identical.phenoData = readline("Are all column names in phonoData identical?")

if( identical.phenoData %in% c(TRUE, T, "yes", "y")){
    class <- readline(paste("Which of the following classes describes the individual groups: ",paste(colnames(xcmsSets[[1]]@phenoData), collapse=" or ")))
    batch <- readline(paste("Which of the following classes describes the batches: ",paste(colnames(xcmsSets[[1]]@phenoData), collapse=" or ")))
} else{
    warning("the classlabels are not identical over all xcmsSet object. The purrr::walk and purrr::map statements will have to be changed accordingly")
    class <- walk( .x = xcmsSets,
                   .f = ~ readline(paste("Which of the following classes describes the individual groups: ",paste(colnames(.x@phenoData), collapse=" or ")))
    )
    batch <- walk(.x = xcmsSets,
                  .f = ~ readline(paste("Which of the following classes describes the batches: ",paste(colnames(.x@phenoData), collapse=" or ")))
    )
}


## QC plots for the peak detection
walk2(.x = xcmsSets,
      .y = names(xcmsSets),
      .f = ~ QC_plots_peaks(.x, className = class, plottitle = .y)
)


map(.x = xcmsSets,
    .f = ~ ggplot(data.frame(.x@peaks), aes(log(into),egauss)) + 
        stat_bin2d(bins=25) + 
        scale_fill_gradientn(colours=r)
)

map(.x = xcmsSets,
    .f = ~ ggplot(data.frame(.x@peaks), aes(rtmax-rtmin)) + 
        geom_histogram()
)

## Peak filtering part 1 
xcmsSets_peakfiltered = map(.x = xcmsSets,
                            .f = ~ PeakFiltR(.x, 
                                             min_RT_width = 5, 
                                             max_RT_width = 150,  
                                             gauss_fail_threshold = 0.4)
)
## updated QC plots after the peak filtering
map2(.x = xcmsSets_peakfiltered,
     .y = names(xcmsSets_peakfiltered),
     .f = ~ ggplot(data.frame(.x@peaks), aes(rtmax-rtmin)) + 
         geom_histogram() + 
         ggtitle(.y)
)

walk2(.x = xcmsSets_peakfiltered,
      .y = names(xcmsSets_peakfiltered),
      .f = ~ MetaboMeeseeks::QC_plots_peaks(.x, className = class, plottitle = .y)
)


# first grouping 
xcmsSets_grouped = map(.x = xcmsSets_peakfiltered, #xcmsSets,#xcmsSets_peakfiltered,
                       .f = ~ group(.x, 
                                    method="density", 
                                    bw=5, 
                                    mzwid = 0.015, 
                                    minfrac=0.50, 
                                    max=100, 
                                    minsamp = 3)
)


## splitting of batches for single batch analysis

# xcmsSets_batchSplit = map(.x = xcmsSets_grouped,
#                           .f = ~ MetaboMeeseeks::BatchSplitR(xcmsObject = .x ,
#                                                              BatchClassLabel = batch)
# )
# 
# 
# BatchProcessingStructure = MetaboMeeseeks::getBatchStructure(xcmsSets_batchSplit)
# 
# xcmsSets_batchSplit = MetaboMeeseeks::unStructureBatches(xcmsSets_batchSplit, BatchProcessingStructure)

##

# batchSets_grouped = map(.x = xcmsSets_batchSplit,#xcmsSets_peakfiltered,
#                         .f = ~ group(.x,
#                                      method="density",
#                                      bw=5,
#                                      mzwid = 0.015,
#                                      minfrac=0.50,
#                                      max=100,
#                                      minsamp = 3)
# )
# 
# batchSets_grouped_filled = map2(.x = batchSets_grouped,
#                                           .y = xcmsSets_batchSplit,
#                                           .f = ~ batchPeakFilling(xcmsObject = .x,
#                                                                   unStructureBatchesObject = .y)
# )
# 
# batchSets_grouped_filled_filtered = map(.x = batchSets_grouped_filled,
#                                           .f = ~ GroupFiltR(.x,
#                                                             max_RT_width = 100,
#                                                             max_ppm_diff = 50,
#                                                             max_dupl_prop = 0.3,
#                                                             max_missing_prop = 0.3,
#                                                             thresholds_to_pass = 4)
# 
# )

###
xcmsSets_grouped_filled = map2(.x = xcmsSets_grouped,  # the grouped ones
                                .y = xcmsSets_peakfiltered, # the ungrouped ones 
                                .f = ~ batchPeakFilling(xcmsObject = .x, 
                                                        unStructureBatchesObject = .y)
)

# batchSets_grouped_filled_filtered = map(.x = batchSets_grouped_filled, 
#                                         .f = ~ GroupFiltR(.x, 
#                                                           max_RT_width = 100, 
#                                                           max_ppm_diff = 50, 
#                                                           max_dupl_prop = 0.3, 
#                                                           max_missing_prop = 0.3, 
#                                                           thresholds_to_pass = 4)
#                                         
# )

xcmsSets_grouped_filled_filtered = map(.x = xcmsSets_grouped_filled, 
                                        .f = ~ GroupFiltR(.x, 
                                                          max_RT_width = 100, 
                                                          max_ppm_diff = 25, 
                                                          max_dupl_prop = 0.3, 
                                                          max_missing_prop = 0.3, 
                                                          thresholds_to_pass = 4)
                                        
)


# GroupData =  map(.x = batchSets_grouped_filled_filtered, 
#                  .f = ~ .x@groups)
 GroupData =  map(.x = xcmsSets_grouped_filled_filtered, 
                  .f = ~ .x@groups)


# FilledData_Matrices = map(.x = batchSets_grouped_filled_filtered, 
#                           .f = ~ t(groupval(.x, value = "into", method = "medret"))
# )
FilledData_Matrices = map(.x = xcmsSets_grouped_filled_filtered, 
                          .f = ~ t(groupval(.x, value = "into", method = "medret"))
)



# pmap(.l = list(x = FilledData_Matrices, 
#                y = batchSets_grouped_filled_filtered, 
#                z = names(batchSets_grouped_filled_filtered))  , 
#      .f = function(x,y,z)  QC_plots_features(FeatureMatrix = x, 
#                                              XCMSobject = y, 
#                                              className = class, 
#                                              plottitle = z)
# )

pmap(.l = list(x = FilledData_Matrices, 
               y = xcmsSets_grouped_filled_filtered, 
               z = names(xcmsSets_grouped_filled_filtered))  , 
     .f = function(x,y,z)  QC_plots_features(FeatureMatrix = x, 
                                             classVector = sampclass(y), 
                                             plottitle = z)
)




FilledData_Matrices_log2 = map(.x = FilledData_Matrices,
                               .f = ~ log2(.x +1)
)



# classlabels = map(.x = batchSets_grouped_filled_filtered,
#                   .f = ~ as.character(.x@phenoData[class][,1])
# )
classlabels = map(.x = xcmsSets_grouped_filled_filtered,
                  .f = ~ as.character(.x@phenoData[class][,1])
)


original_column_indices =  map(.x = FilledData_Matrices_log2,
                               .f = ~ seq(1,ncol(.x))
)


subclasses = list("ctrl", "diluted", "IC10", "QC") # noem het COI




################################### missing values filter ########################################

missing_pass = map2(.x = FilledData_Matrices_log2,
                    .y = classlabels,
                    .f = ~ MissingValuesCheckR(DataMatrix = .x, 
                                                subclasses_list = subclasses, 
                                                prop_missing_threshold = 0.2, 
                                                NA_or_numeric_limit = 1, # peak filling imputeert een nul als hij niks vindt. In dit geval zet je deze limiet op 1, en alles kleiner dan of gelijk aan 1 wordt aanzien al missing
                                                labels = .y, 
                                                feature_orientation = "columns", 
                                                groups_ok_threshold = 1)
)

original_column_indices = map2(.x = original_column_indices,
                               .y = missing_pass,
                               .f = ~ .x[.y]
)

FilledData_Matrices_log2 = map2(.x = FilledData_Matrices_log2,
                                .y = missing_pass,
                                .f = ~ .x[,.y]
)


walk2(.x = FilledData_Matrices_log2,
      .y = xcmsSets_grouped_filled_filtered,
      .f = ~ MetaboMeeseeks::QC_plots_features(FeatureMatrix = .x, 
                                               classVector = sampclass(.y))
)

walk2(.x = FilledData_Matrices_log2,
      .y = classlabels,
      .f = ~ MetaboMeeseeks::QC_plots_features(FeatureMatrix = .x, 
                                               classVector = .y)
)


################################### median RSD filter ########################################
# this one is not working
RSD_pass = map2(.x = FilledData_Matrices_log2,
                .y = classlabels,
                .f = ~ mRSD_check(DataMatrix = .x, 
                                  subclasses = unlist(subclasses), 
                                  RSD_threshold = 0.05, 
                                  ignore_missing_values = TRUE, 
                                  NA_or_numeric_limit = 1, 
                                  labels = .y,
                                  include_RSD_matrix = FALSE)
)

original_column_indices = map2(.x = original_column_indices,
                               .y = RSD_pass,
                               .f = ~ .x[.y]
)

FilledData_Matrices_log2 = map2(.x = FilledData_Matrices_log2,
                                .y = RSD_pass,
                                .f = ~ .x[,.y]
)
walk2(.x = FilledData_Matrices_log2,
      .y = xcmsSets_grouped_filled_filtered,
      .f = ~ MetaboMeeseeks::QC_plots_features(FeatureMatrix = .x, 
                                               classVector = sampclass(.y))
)
walk2(.x = FilledData_Matrices_log2,
      .y = classlabels,
      .f = ~ MetaboMeeseeks::QC_plots_features(FeatureMatrix = .x, 
                                               classVector = .y)
)



FilledData_Matrices_log2_scaled = map2(.x = FilledData_Matrices_log2,
                                       .y = classlabels,
                                       .f = ~ speaq::SCANT(.x[.y %in% unlist(subclasses),], type = c( "center"), what = "columns")
)




##### PCA
part = 1

dim(FilledData_Matrices_log2[[part]])

PCADAMATRIX = FilledData_Matrices_log2_scaled[[part]]
pca_colnames = colnames(PCADAMATRIX)
classlabels.subset = classlabels[[part]][ classlabels[[part]] %in% unlist(subclasses)]
plot.marks = paste(xcmsSets_grouped_filled_filtered[[part]]@phenoData[classlabels[[part]] %in% unlist(subclasses),1],xcmsSets_grouped_filled_filtered[[part]]@phenoData[classlabels[[part]] %in% unlist(subclasses),2], sep =".")

colnames(PCADAMATRIX) = paste("feat",as.character(seq(1,ncol(PCADAMATRIX))), sep = "" )
colnames(PCADAMATRIX) = as.character(seq(1,ncol(PCADAMATRIX)))


PCA.plotcol = classlabels[[part]][classlabels[[part]] %in% unlist(subclasses) ]
PCA.plotshape = as.character(xcmsSets_grouped_filled_filtered[[1]]@phenoData["V14"][,1])[classlabels[[part]] %in% unlist(subclasses) ]
#PCA.plotcol = classlabels[[part]][classlabels[[part]] %in% unlist(subclasses)[1:3] ]
#PCA.plotshape = as.character(xcmsSets_grouped_filled_filtered[[1]]@phenoData["V14"][,1])[classlabels[[part]] %in% unlist(subclasses)[1:3] ]


common.pca <- prcomp(PCADAMATRIX)


loadings <- common.pca$rotation
scores <- common.pca$x
varExplained <- common.pca$sdev^2

barplot(varExplained/sum(varExplained), 
        main="Scree Plot",
        ylab="Proportion of variance explained", 
        xlab = "Principal comonent", 
        names.arg = as.character(seq(1,length(varExplained))) )

PCA_labelsnum = rep(NA, length(PCA.plotcol))
PCA_labelsnum[PCA.plotcol == "QC"] = 3
PCA_labelsnum[PCA.plotcol == "ctrl"] = 1
PCA_labelsnum[PCA.plotcol == "diluted"] = 8
PCA_labelsnum[PCA.plotcol == "IC10"] = 2
cp1 <- 1
cp2 <- 2 
pdf(file = "/Users/Charlie/Dropbox/PhD/Thesis/Chapter3_meeseeks/images/NVPA_PCA_1_2.pdf", width = 7, height = 4.5)
par(mar = c(4.5, 4.5, 2, 2)) # bottom, left, top, right
plot(scores[,cp1]/max(scores[,cp1]), scores[,cp2]/max(scores[,cp2]),
     main=paste("Score plot, PC",cp1," vs. PC",cp2,sep=""),
     xlab=paste("PC",cp1,"  (",round(varExplained[cp1]/sum(varExplained),digits=2),")",sep = ""),
     ylab=paste("PC",cp2,"  (",round(varExplained[cp2]/sum(varExplained),digits=2),")",sep = ""),
     col = as.numeric(as.factor(PCA.plotshape)),
     pch = PCA_labelsnum)
lines(x = c(-100,100), y = c(0,0))
lines(x = c(0,0), y = c(-100,100))
legend("topright", legend=c("Batch 1", "Batch 2", "","QC", "control", "diluted", "IC10"),
       col=c("black", "red", "black","black", "black", "black", "black"), 
       text.col = c("black", "red", rep("black", 5)),
       pch = c(NA,NA,NA, 3, 1, 8 , 2), cex=0.8)
dev.off()

cp1 <- 3
cp2 <- 4 
pdf(file = "/Users/Charlie/Dropbox/PhD/Thesis/Chapter3_meeseeks/images/NVPA_PCA_3_4.pdf", width = 7, height = 4.5)
par(mar = c(4.5, 4.5, 2, 2)) # bottom, left, top, right
plot(scores[,cp1]/max(scores[,cp1]), scores[,cp2]/max(scores[,cp2]),
     main=paste("Score plot, PC",cp1," vs. PC",cp2,sep=""),
     xlab=paste("PC",cp1,"  (",round(varExplained[cp1]/sum(varExplained),digits=2),")",sep = ""),
     ylab=paste("PC",cp2,"  (",round(varExplained[cp2]/sum(varExplained),digits=2),")",sep = ""),
     col = as.numeric(as.factor(PCA.plotshape)),
     pch = PCA_labelsnum)
lines(x = c(-100,100), y = c(0,0))
lines(x = c(0,0), y = c(-100,100))
legend("topright", legend=c("Batch 1", "Batch 2", "","QC", "control", "diluted", "IC10"),
       col=c("black", "red", "black","black", "black", "black", "black"), 
       text.col = c("black", "red", rep("black", 5)),
       pch = c(NA,NA,NA, 3, 1, 8 , 2), cex=0.8)
dev.off()

par(mar = c(5.1, 4.1, 4.1, 2.1)) # bottom, left, top, right

important_loadings = as.numeric(which((loadings[,4] > 0.1 | loadings[,4] < -0.1)))
imfeats_pca = pca_colnames[important_loadings]


tet = MetaboMeeseeks::Meeseeks.RF(PCADAMATRIX[classlabels[[part]] %in% unlist(subclasses)[c(1,3)],], classlabels[[part]][ classlabels[[part]] %in% unlist(subclasses)[c(1,3)]])

varimpPlot_df = data.frame(var = seq(1,ncol(PCADAMATRIX)), impmean = matrixStats::colMeans2(tet$varImportance), impsds = matrixStats::colSds(tet$varImportance))
plot(varimpPlot_df$impmean, varimpPlot_df$impsds)
impvars13 = which(varimpPlot_df$impmean > 1)
imfeats_rf = pca_colnames[impvars13 ]
plot(as.factor(classlabels[[part]]), PCADAMATRIX[,impvars[1]])


















common.pca <- prcomp(t(testing), center = T, scale. = T)


loadings <- common.pca$rotation
scores <- common.pca$x
varExplained <- common.pca$sdev^2

barplot(varExplained/sum(varExplained), 
        main="Scree Plot",
        ylab="Proportion of variance explained", 
        xlab = "Principal comonent", 
        names.arg = as.character(seq(1,length(varExplained))) )

PCA_labelsnum = rep(NA, length(PCA.plotcol))
PCA_labelsnum[PCA.plotcol == "QC"] = 3
PCA_labelsnum[PCA.plotcol == "ctrl"] = 1
PCA_labelsnum[PCA.plotcol == "diluted"] = 8
PCA_labelsnum[PCA.plotcol == "IC10"] = 2
cp1 <- 1
cp2 <- 2 
plot(scores[,cp1]/max(scores[,cp1]), scores[,cp2]/max(scores[,cp2]),
     main=paste("Score plot, PC",cp1," vs. PC",cp2,sep=""),
     xlab=paste("PC",cp1,"  (",round(varExplained[cp1]/sum(varExplained),digits=2),")",sep = ""),
     ylab=paste("PC",cp2,"  (",round(varExplained[cp2]/sum(varExplained),digits=2),")",sep = ""),
     col = as.numeric(as.factor(PCA.plotshape)),
     pch = PCA_labelsnum)
lines(x = c(-100,100), y = c(0,0))
lines(x = c(0,0), y = c(-100,100))
legend("topright", legend=c("Batch 1", "Batch 2", "","QC", "control", "diluted", "IC10"),
       col=c("black", "red", "black","black", "black", "black", "black"), 
       text.col = c("black", "red", rep("black", 5)),
       pch = c(NA,NA,NA, 3, 1, 8 , 2), cex=0.8)


cp1 <- 5
cp2 <- 6 
plot(scores[,cp1]/max(scores[,cp1]), scores[,cp2]/max(scores[,cp2]),
     main=paste("Score plot, PC",cp1," vs. PC",cp2,sep=""),
     xlab=paste("PC",cp1,"  (",round(varExplained[cp1]/sum(varExplained),digits=2),")",sep = ""),
     ylab=paste("PC",cp2,"  (",round(varExplained[cp2]/sum(varExplained),digits=2),")",sep = ""),
     col = as.numeric(as.factor(PCA.plotshape)),
     pch = PCA_labelsnum)
lines(x = c(-100,100), y = c(0,0))
lines(x = c(0,0), y = c(-100,100))
legend("topright", legend=c("Batch 1", "Batch 2", "","QC", "control", "diluted", "IC10"),
       col=c("black", "red", "black","black", "black", "black", "black"), 
       text.col = c("black", "red", rep("black", 5)),
       pch = c(NA,NA,NA, 3, 1, 8 , 2), cex=0.8)


important_loadings = as.numeric(which((loadings[,4] > 0.1 | loadings[,4] < -0.1)))



tet = MetaboMeeseeks::Meeseeks.RF(t(testing)[classlabels[[part]] %in% unlist(subclasses)[c(1,3)],], classlabels[[part]][ classlabels[[part]] %in% unlist(subclasses)[c(1,3)]])

varimpPlot_df = data.frame(var = seq(1,ncol(t(testing))), impmean = matrixStats::colMeans2(tet$varImportance), impsds = matrixStats::colSds(tet$varImportance))
plot(varimpPlot_df$impmean, varimpPlot_df$impsds)
impvars_testing = rownames(testing)[which(varimpPlot_df$impmean > 0.5)]


