require(cluster)
require(doParallel)
require(foreach)
require(e1071)

library(rlang)
library(tictoc)
library(RWeka)
library(xgboost)
library(stringr)
library(optparse)

source('spn.value.r')
source('spn.value.int.r')
source('spn.evidence.r')
source('spn.evidence.int.r')
source('spn.learn.r')
source('spn.print.r')
source('spn.test.r')
source('spn.utils.r')
source('spn.clara2.r')
source('spn.clustering.r')

getDummies <- function(data) {
    for(var in colnames(data)){
        data[[var]] <- as.integer(as.factor(data[[var]]))
    }
    return(data)
}

datasets <- c('zoo',
              'bridges_version1',
              'lymph',
              'flags',
              'autos',
              'breast-cancer',
              'heart-h',
              'ecoli',
              'liver-disorders',
              'dermatology',
              'colic',
              'balance-scale',
              'soybean',
              'diabetes',
              'vehicle',
              'tic-tac-toe',
              'vowel',
              'solar-flare_2',
              'cmc',
              'car',
              'segment',
              'sick',
              'hypothyroid',
              'spambase',
              'nursery')

option_list = list(
  make_option(c("-d", "--dataset"), type="character", default=NULL,
              help="Dataset file name.", metavar="character"),
	make_option(c("-e", "--evaluation"), type="character", default="acc",
              help="Type of experiment", metavar="character"),
  make_option(c("-o", "--output_dir"), type="character", default=NULL,
              help="Output dir", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

nruns = 1
ncores = 1
files <- list.files(path="../data/arff", full.names = TRUE)
if (is.null(opt$output_dir)) {
    output_dir <- opt$evaluation
} else output_dir <- opt$output_dir
if (!is.null(opt$dataset)) {
    datasets = c(opt$dataset)
}

if(is.element('mem', opt$evaluation)) {
    cat("Running memory comparison on the following datasets: ", datasets, '\n')
} else if(is.element('acc', opt$evaluation)) {
    cat("Running accuracy comparison on the following datasets: ", datasets, '\n')
} else {
    stop('Unknown evaluation method. Please set -e to either acc or mem.')
}

if (!dir.exists(output_dir)) {dir.create(output_dir)}
seed=42
for(file in files){
    name <- tools::file_path_sans_ext(basename(file))
    if(is.element(name, datasets)) {
      cat('Next dataset: ', name, '\n')
      data <- read.arff(file)
      data <- getDummies(data)
      if(is.element('mem', opt$evaluation)) {
          res <- spn.cv(data, nfolds=5, classcol=ncol(data), nruns=nruns,
                        ncores=ncores, use_memory=FALSE, seed=seed, run.rob=TRUE,
                        path=paste(output_dir, '/', name, '_nomem.csv', sep=""))
          res_mem <- spn.cv(data, nfolds=5, classcol=ncol(data), nruns=nruns,
                            ncores=ncores, use_memory=TRUE, seed=seed, run.rob=TRUE,
                            path=paste(output_dir, '/', name, '_mem.csv', sep=""))
      } else {
          spn.cv_acc(data, nfolds=5, classcol=ncol(data), nruns=nruns,
                            ncores=ncores, use_memory=TRUE, run.rob=TRUE,
                            path=paste(output_dir, '/', name, sep=""))
      }
  }
}
