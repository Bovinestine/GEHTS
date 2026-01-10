# GE-HTS project pipeline for in situ sequencing data analysis
# author: Nathan Wooseok Lee

# saving data
pub_path <- './Rdata/pub'

save(mac.pub, file = paste0(pub_path, '/seurat_mac.pub.Rdata'))
save(sin.pub, file = paste0(pub_path, '/seurat_sin.pub.Rdata'))
save(cmb.pub, file = paste0(pub_path, '/seurat_cmb.pub.Rdata'))
save(rna.pub, file = paste0(pub_path, '/seurat_rna.pub.Rdata'))

save(mac.sct, file = paste0(pub_path, '/seurat_mac.sct.Rdata'))
save(sin.sct, file = paste0(pub_path, '/seurat_sin.sct.Rdata'))
save(cmb.sct, file = paste0(pub_path, '/seurat_cmb.sct.Rdata'))
save(rna.sct, file = paste0(pub_path, '/seurat_rna.sct.Rdata'))


# Load necessary data by activating the following code:
# Actual variable name in the middle of file name: [prefix]_[variable name].Rdata
dir_path <- './Rdata/pub'
# Raw data before normalization
load(paste0(dir_path, '/seurat_mac.pub.Rdata')) # variable name: mac # refreshed on 251230 about 230508
load(paste0(dir_path, '/seurat_sin.pub.Rdata')) # variable name: singleAll
load(paste0(dir_path, '/seurat_cmb.pub.Rdata'))
load(paste0(dir_path, '/seurat_rna.pub.Rdata')) # RNA seq data

# SCT normalized for application of gene counts
load(paste0(dir_path, '/seurat_mac.sct.Rdata'))
load(paste0(dir_path, '/seurat_sin.sct.Rdata'))
load(paste0(dir_path, '/seurat_cmb.sct.Rdata')) 
load(paste0(dir_path, '/seurat_rna.sct.Rdata')) 

library(SeuratDisk)
macSebastian <- readRDS('./RNA-seq/macCtrl_sebastian.rds')# healthy mac RNA-seq data (filtered and SCTransformed) from Sebastian 2021

# probability of efficacy of drugs
# random forest model probability
# the loaded variables do not have the subfix ("_rf" things)
load(paste0(dir_path, '/prediction/prob_sng10_rf.Rdata'))
load(paste0(dir_path, '/prediction/prob_sng0.1_rf.Rdata'))
load(paste0(dir_path, '/prediction/prob_cmb0.1_rf.Rdata'))

# Prediction models 
"Elastic Net", "Random Forest", "XGBoost", "SVR", "GPR"

en_model <- readRDS("./Rdata/prediction/Elastic Net_model.rds")
rf_model <- readRDS("./Rdata/prediction/Random Forest_model.rds")
xgb_model <- readRDS("./Rdata/prediction/XGBoost_model.rds")
svr_model <- readRDS("./Rdata/prediction/SVR_model.rds")
gpr_model <- readRDS("./Rdata/prediction/GPR_model.rds")

# universal scale for radar chart
load(paste0(dir_path, '/uni_scale.Rdata'))
