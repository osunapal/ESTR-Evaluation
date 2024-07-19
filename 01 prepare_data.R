library(dplyr)
library(rjson)
############################################################################
#If there are duplicate family IDs, they should be made unique:
#to use it: fam$FamilyID <- deduplicate_family_ids(fam$FamilyID)
deduplicate_family_ids <- function(family_ids) {
  # Count occurrences of each family ID
  id_counts <- table(family_ids)
  
  # Create a named vector for suffixes
  suffixes <- setNames(rep("", length(family_ids)), family_ids)
  
  # Iterate over each family ID
  for (id in unique(family_ids)) {
    if (id_counts[id] > 1) {
      # Assign suffixes for duplicates
      duplicates <- which(family_ids == id)
      suffixes[duplicates] <- paste0("_", seq_along(duplicates) - 1)
    }
  }
  # Generate unique IDs by appending suffixes
  unique_family_ids <- paste0(family_ids, suffixes)
  return(unique_family_ids)
}

############################################################################
#Load the data:
d <- fromJSON(file="merged_full.multisample_profile.json")
IDs	<- read.table( "familyIDs.txt",h=T,as.is=T )
#Deduplicate the FamilyIDs
IDs$FamilyID <- deduplicate_family_ids(IDs$FamilyID)
uniqueIDs <- unique(c(IDs$FatherID,IDs$MotherID,IDs$SampleID))
results <- data.frame( matrix(nrow=0,ncol=length(uniqueIDs)) )
colnames( results ) <- uniqueIDs

locus <- ""
chr <- ""
pos <- 0
motifs <- names( d$Counts )
total_iterations <- length(motifs)
pb <- txtProgressBar(min = 0, max = total_iterations, style = 3)
counter <- 0
for( motif in motifs )
{
  counter <- counter + 1
	whichrow <- nrow(results) + 1
	regions <- names( d$Counts[[motif]]$RegionsWithIrrAnchors )
	for( region in regions )
	{
		result <- d$Counts[[motif]]$RegionsWithIrrAnchors[[region]]
		for( person in names(result) )
		{
			results[whichrow,person]<-result[[person]]
		}
		locus[whichrow] <- region
		locusinfo <- unlist(strsplit(region,":"))
		chr[whichrow] <- locusinfo[1]
		pos[whichrow] <- mean(as.numeric(unlist(strsplit(locusinfo[2],"-"))))
	}
	setTxtProgressBar(pb, counter)
}
close(pb)

############################################################################
#Split the data into father, mother and child

father<-data.frame(matrix(NA,nrow=nrow(results),ncol=nrow(IDs)))
names(father)<-IDs$FamilyID
mother<-father
child<-father

for(i in 1:nrow(IDs))
{
  familyID <- IDs$FamilyID[i]
  fatherID <- IDs$FatherID[i]
  motherID <- IDs$MotherID[i]
  childID  <- IDs$SampleID[i]
  
  if( fatherID %in% names(results) )
  {
    result   <- results[,names(results)==fatherID]
    father[,i] <- result
  }
  if( motherID %in% names(results) )
  {
    result   <- results[,names(results)==motherID]
    mother[,i] <- result
  }
  if( childID %in% names(results) )
  {
    result   <- results[,names(results)==childID]
    child[,i] <- result
  }
}

