cf <- c(
  0.2904047, 0.2904047, 0.2070468, 0.2070468, 0.2676134,
  0.2676134, 0.3284704, 0.3284704, 0.3979398, 0.3979398,
  0.4863355, 0.4863355, 0.6306967, 0.6306967, 0.6906705, 0.6906705
)

###### susceptibility ####
# (based on Davies et al, Nature paper)
sus <- c(
  0.3956736, 0.3956736, 0.3815349, 0.3815349, 0.7859512,
  0.7859512, 0.8585759, 0.8585759, 0.7981468, 0.7981468,
  0.8166960, 0.8166960, 0.8784811, 0.8784811, 0.7383189, 0.7383189
)

#### load contact matrices ####
for(x in 2:5){
  load(paste0("data/", list.files("data", pattern = "contact"))[x])
}

# countrycode::countrycode("Thailand", "country.name", "iso3c")
cm_matrices[["Thailand"]]$home <- as.matrix(contact_home$THA) 
dimnames(cm_matrices[["Thailand"]]$home) <- 
  list(rownames(cm_matrices[["Brazil"]]$home),
       colnames(cm_matrices[["Brazil"]]$home))

cm_matrices[["Thailand"]]$work <- as.matrix(contact_work$THA) 
dimnames(cm_matrices[["Thailand"]]$work) <- 
  list(rownames(cm_matrices[["Brazil"]]$work),
       colnames(cm_matrices[["Brazil"]]$work))

cm_matrices[["Thailand"]]$school <- as.matrix(contact_school$THA) 
dimnames(cm_matrices[["Thailand"]]$school) <- 
  list(rownames(cm_matrices[["Brazil"]]$school),
       colnames(cm_matrices[["Brazil"]]$school))

cm_matrices[["Thailand"]]$other <- as.matrix(contact_others$THA) 
dimnames(cm_matrices[["Thailand"]]$other) <- 
  list(rownames(cm_matrices[["Brazil"]]$other),
       colnames(cm_matrices[["Brazil"]]$other))

