#' Get FishLife life history traits
#'
#' Retrieves life history trait predictions from FishLife (Thorson et al. 2018). This
#' is a wrapper for the Plot_taxa() function in the FishLife R package.
#'
#' @param species A character vector of species scientific names
#' @return A dataframe with life history trait predictions from FishLife for each species
#' @examples
#' # Look up life history traits
#' species <- c("Gadus morhua", "Centropristis striata", "Paralichthys dentatus")
#' fishlife(species)
#' @export
#' 
 install.packages("remotes")
 remotes::install_github("James-Thorson/FishLife")
 
library(FishLife)
library(dplyr)


#### example of species vector:
species = c('Cynoscion acoupa',
           'Mugil rubrioculus')

save_it = FALSE


fishlife <- function(species){

  # Setup container
  spp <- sort(unique(species))
  fl <- data.frame(species=spp, linf_cm=NA, k=NA, winf_g=NA, tmax_yr=NA, tmat_yr=NA,
                   m=NA, lmat_cm=NA, temp_c=NA, stringsAsFactors=F)

  # Loop through species
  for(i in 1:nrow(fl)){

    # Get spp info
    sciname <- fl$species[i]
    genus <- stringr::word(sciname, 1)
    nwords_in_spp <- length(strsplit(sciname, " ")[[1]])
    species <- stringr::word(sciname, start=2, end=nwords_in_spp)
    species <- ifelse(species=="spp", "predictive", species)

    # Try looking up in FishLife
    spp_info <- try(FishLife::Plot_taxa(FishLife::Search_species(Genus=genus, Species=species)$match_taxonomy))
    if(inherits(spp_info, "try-error")){
      # Record blanks
#      fl[i,2:ncol(fl)] <- rep(NA, ncol(fl)-1)
    }else{
      # Values are in log-scale except temperature
      spp_lh_vals_log <- spp_info[[1]]$Mean_pred
      spp_lh_vals <- c(exp(spp_lh_vals_log[1:7]), spp_lh_vals_log[8],spp_lh_vals_log[9:20])
 #     fl[i,2:ncol(fl)] <- spp_lh_vals
    }

  }

  # Return
#  return(fl)
  return(spp_lh_vals)
}

### Update 2

##### modified function to loop through species and compile traits in a single data frame,
##### with species as rows and trait values as columns
#### check species input 

fishlife_for_many = function(species){
  
  # Setup container
  spp = sort(unique(species))
  fl = data.frame(species=spp, linf_cm=NA, k=NA, winf_g=NA, tmax_yr=NA, tmat_yr=NA,
                  m=NA, lmat_cm=NA, temp_c=NA, stringsAsFactors=FALSE)
  
  # Column names for the life history values
  col_names = c("Loo", "K", "Winfinity", "tmax", "tm", "M", "Lm", "Temperature",
                "ln_var", "rho", "ln_MASPS", "ln_margsd", "h", "logitbound_h",
                "ln_Fmsy_over_M", "ln_Fmsy", "ln_r", "r", "ln_G", "G")
  
  # Initialize the final data frame
  final_data = data.frame(matrix(ncol = length(col_names) + 1, nrow = 0))
  colnames(final_data) = c("species", col_names)
  
  # Loop through species
  for(i in 1:nrow(fl)){
    
    # Get spp info
    sciname = fl$species[i]
    genus = stringr::word(sciname, 1)
    nwords_in_spp = length(strsplit(sciname, " ")[[1]])
    species = stringr::word(sciname, start=2, end=nwords_in_spp)
    species = ifelse(species=="spp", "predictive", species)
    
    # Try looking up in FishLife
    spp_info = try(FishLife::Plot_taxa(FishLife::Search_species(Genus=genus, Species=species)$match_taxonomy), silent=TRUE)
    if(inherits(spp_info, "try-error")){
      # Record blanks
      species_data = data.frame(matrix(NA, nrow=1, ncol=length(col_names)))
      colnames(species_data) = col_names
    } else {
      # Values are in log-scale except temperature
      spp_lh_vals_log = spp_info[[1]]$Mean_pred
      spp_lh_vals = c(exp(spp_lh_vals_log[1:7]), spp_lh_vals_log[8], spp_lh_vals_log[9:20])
      # Create a data frame with the life history values
      species_data = data.frame(matrix(spp_lh_vals, nrow=1))
    }
    
    colnames(species_data) = col_names
    species_data = cbind(species=sciname, species_data)
    
    # Final data frame
    final_data = rbind(final_data, species_data)
  }
  return(final_data)
}

if (save_it == TRUE) { 
  
  fishlife_for_many(species) %>%
    write.csv('fishlife_results.csv') 
  
  fishlife_for_many(species) %>%
    View()
  
  print('Done!')
} else {
  fishlife_for_many(species) %>%
    View()
  
  print('Done!')}


