#Step 1 

install.packages("devtools")

library(devtools)


devtools::install_github("https://github.com/lr98769/doppelgangerIdentifier")

#install sva package due to error prompt for sva dependnecy above

if(!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("sva")

# check if DI package is installed 

library(doppelgangerIdentifier)

#viewing all the functions / objects in the package 
lsf.str("package:doppelgangerIdentifier")


# 4 dataset (cited)  - ("Paper: Doppelganger Spotting in Biomedical Gene Expression Data- Table 2")
#- rc- renal cell carcinoma proteomics data set
#- dmd - duchenne muscular dystrophy (dmd) micro array data set
#lk -leukemia data set
#all- acute lymphoblastic leukemia miroarray dataset

#Step 2 : import renal carcinoma data & metadata to identify dd and ascertain if they have an inflationary effect on model accuracy

#importing microarray rc  dataset
data(rc)
#import metadata for rc dataset
data(rc_metadata)

start_time = Sys.time()


#getPPCCDoppelgangers detects ppcc dd btw 2 batches of data or within a batch
#possible doppelgangers are samples from the same class but different patients. samples pairs of same class and patient are
#indicative of leakage while sample pairs of different class and patient are indicative of negative controls
ppccDoppelgangerResults_rc = getPPCCDoppelgangers(rc, rc_metadata)



#Step 3 : 

#dataset before ppcc calcualtion for min-max normalization - step 1 of identification method as explained in paper
View(ppccDoppelgangerResults_rc$Processed_data)


#dataset of ppcc matrix btw samples pairs of different batch (seperate datasets). PPCC is a metric of similarity- step 2 of identification method
View(ppccDoppelgangerResults_rc$PPCC_matrix)

# group the sample pairs - step 3 of identification method as explained in paper
View(ppccDoppelgangerResults_rc$PPCC_df)

# calcualting the cutoff for the identification of dd which are defined as sample pairs from same class and different patients with ppcc > cut-off point
paste("PPCC cut-off:", ppccDoppelgangerResults_rc$cut_off)

# Step 5 - visualize the output via a univariate scatterplot of dd by passing the output from the function getPPCCDopplegangers
visualisePPCCDoppelgangers(ppccDoppelgangerResults_rc)




ppccDoppelgangerResults_leuk = getPPCCDoppelgangers(leuk, leuk_metadata)

visualisePPCCDoppelgangers(ppccDoppelgangerResults_leuk)


end_time = Sys.time()

end_time - start_time

# Additional - create a list of datasets and their metadata

ds_all <- c(rc, leuk, dmd, all)

ds_metadata_all <- c(rc_metadata, leuk_metadata, dmd_metadata, all_metadata)






#--------------------------------------------------------------------------------------------------------
#Part 2 ; Verify if PPCC dd are also functional doppelgangers

# Step 1 of identificaton method 


# train random knn models per user-defined experiment plan to verify the confounding effects of ppcc dd identifed by getPPCCDoppleganger function above. 
verificationResults_rc = verifyDoppelgangers("https://github.com/lr98769/doppelgangerIdentifier/blob/main/tutorial/experiment_plans/rc_ex_plan.csv", rc, rc_metadata)



#preprocessing stage to get the correct the batch and perform min-max normalization, just before model training
View(verificationResults_rc$combat_minmax)

#random feature selection to get 10 randomly generated feature sets and 2 feature sets with features of lowest and highest variance
#each feature set is a vector of gene labels
View(verificationResults_rc$feature_sets)

#seperate the samples into training and validation datasets and view the accuracies of training-validation in matrix format
#It is the column which shows model accuracies for training-validation set
View(verificationResults_rc$accuracy_mat)

#now view the acuracies of training-validation dataset in dataframe form
#each row shows the accuracy of model trained and evaluated on specific training-validation set with specific feature set.
View(verificationResults_rc$accuracy_df)


ori_train_valid_names = c("Doppel_0","Doppel_2", "Doppel_4", "Doppel_6", "Doppel_8", "Neg_Con", "Pos_Con")

new_train_valid_names = c("0 Doppel", "2 Doppel", "4 Doppel", "6 Doppel", "8 Doppel", "Binomial", "Perfect Leakage")

visualiseVerificationResults(verificationResults_rc, ori_train_valid_names, new_train_valid_names)


verificationResults_leuk = verifyDoppelgangers("https://github.com/lr98769/doppelgangerSpotting/blob/master/experiment_plans/leuk_experiment_plan.csv", leuk, leuk_metadata)

ori_train_valid_names_2 = c("Doppel_0", "Doppel_2", "Doppel_4", "Doppel_6")

new_train_valid_names_2 = c("0 Doppel", "2 Doppel", "4 Doppel", "6 Doppel" )

visualiseVerificationResults(verificationResults_leuk, ori_train_valid_names_2, new_train_valid_names_2)



---------------------------------------------------------------------------------------------------------------
#check working directory and import the bc datasets as copying the entire url gives compressed error
getwd()
  
  
#importing bc RNA sequence dataset to your local
  
bc = readRDS('bc_her2_tut.rds')

bc_meta = readRDS('bc_her2_meta_tut.rds')


#identify dd in the bc RNA-Seq data set (similar steps to the rc dataset)
bc_dd = getPPCCDoppelgangers(raw_data = bc, meta_data = bc_meta, do_batch_corr = TRUE, do_min_max = TRUE, batch_corr_method = "ComBat_seq")

visualisePPCCDoppelgangers(bc_dd)

#verify if detected ppcc dd is also fd. Experiment plan csv file incrementally increases the number of ppcc dd samples
#in the validation dataset


verify_bc = verifyDoppelgangers("https://github.com/lr98769/doppelgangerIdentifier/blob/main/tutorial/experiment_plans/bc_ex_plan.csv", raw_data = bc,
                                meta_data = bc_meta, batch_corr_method = "ComBat_seq", k=9, size_of_val_set = 48, feature_set_portion = 0.01)


#increasing number of dd samples in validation data set is used to illustrate if identified ppcc dd are also functional dopelgangers

ori_train_valid_names = c("Doppel_0","Doppel_6", "Doppel_12", "Doppel_18", "Doppel_24", "Neg_Con", "Pos_Con_24")

new_train_valid_names = c("0 Doppel", "6 Doppel", "12 Doppel", "18 Doppel", "24 Doppel", "Binomial", "24 Perfect Leakage")

visualiseVerificationResults(verify_bc, original_train_valid_names = ori_train_valid_names, new_train_valid_names = new_train_valid_names)


