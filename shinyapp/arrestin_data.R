cons <- readr::read_delim(I(gsub("\t+", "\t", # Substitute multiple tabs by single tab for reading
                                 readLines("4JQIA_consurf_summary.txt"), 
                                 perl = T)),
                          delim = "\t", skip = 13, # First 13 lines contain column descriptions
                          skip_empty_rows = T) %>%
  filter(SCORE != "(normalized)") %>% # Filter out second row of column names/info
  filter(!grepl("*", COLOR, fixed=T)) %>% # Filter out non-confident colors
  type_convert() %>%
  # Fixing import problems
  rename(ID = `    3LATOM`, POS = ` POS`) %>%
  filter(ID != "-") %>%
  # Creating Residue variable
  separate(ID, into=c("Res", NA), sep=":") %>%
  separate(Res, into=c("resname", "resnum"), sep=3) %>%
  unite("Residue", resname, resnum, sep=":", remove=F) %>%
  # Filling in residues of the protein not in ConSurf
  bind_rows(map_dfr(c("GLY:309", "ALA:310", "ASN:311"), ~tibble_row(Residue=., 
                                                                    COLOR=NA, 
                                                                    POS=as.numeric( str_split(., ":")[[1]][2] ) + 8
                                                                    ))) %>%
  arrange(POS)


resvec <- map_int(filter(cons, COLOR == 9 | COLOR == 8)$Residue, ~as.integer(str_split(., ":")[[1]][2]))-7



jumps <- which(diff(resvec) != 1)

test <- data.frame(
  ini=resvec[c(1, jumps + 1)],
  end=resvec[c(jumps, length(resvec))],
  cons=1
)

for (row in 1:nrow(test)) {
  if (row != 1) {
    end_of_previous = test[row-1,]$end
    start_of_present = test[row,]$ini
    
    if ((start_of_present - end_of_previous) > 1) {
      test <- rbind(test,
                    data.frame(ini=end_of_previous+1, end=start_of_present-1, cons=0))
    } else {
      print(c(row, start_of_present - end_of_previous))
    }
  }
}

top_cons <- arrange(test, ini)

if (top_cons[nrow(test),]$end < length(pull(cons, Residue))) {
  top_cons[nrow(test)+1,] <- data.frame(ini=top_cons[nrow(test),]$end+1, end=length(pull(cons, Residue)), cons=0)
}

all_cons <- pull(cons, COLOR)







uniprot <- jsonlite::read_json("P29066.json")

sites <- c() # Binding sites single-position
regions <- NULL # Interaction with proteins and disordered regions
ss <- NULL # Secondary structure


for (elem in uniprot$features) {
  
  if (elem$type == "Binding site") {
    sites <- c(sites, elem$location$start$value)
    
  } else if (elem$type == "Region") {
    regions <- rbind(regions, tibble_row(name = gsub("Interaction with ", "Int. with ", elem$description),
                                         start = elem$location$start$value, end = elem$location$end$value))
    
  } else if (elem$type %in% c("Beta strand", "Helix", "Turn")) {
    ss <- rbind(ss, tibble_row(name = elem$type,
                               start = elem$location$start$value, end = elem$location$end$value))
    
  }
}



regions <- regions %>%
  # Regions have to be in different annotations because they overlap
  mutate(group = case_when(name %in% c("Int. with SRC", "Int. with TRAF6") ~ 1,
                           name %in% c("Int. with CHRM2", "Disordered") ~ 2)) %>%
  # Add the NAs already for group 2 so that the vector is correctly generated later
  add_row(name=NA, group=2, start=1, end=44) %>%
  arrange(start)


regions1 <- filter(regions, group == 1) %>% select(!group)
regions2 <- filter(regions, group == 2) %>% select(!group)



get_vector <- function(reg) {
  library(magrittr)
  
  for (row in 1:nrow(reg)) {
    if (row != 1) {
      end_of_previous = reg[row-1,]$end
      start_of_present = reg[row,]$start
      
      if ((start_of_present - end_of_previous) > 1) {
        # Add a range of NAs for positions without annotation
        reg <- rbind(reg,
                     data.frame(start=end_of_previous+1, end=start_of_present-1, name=NA))
        
      } else if ((start_of_present - end_of_previous) < 1) {
        # Sanity check
        print(c(row, start_of_present - end_of_previous))
      }
    }
  }
  return(reg %>% 
           arrange(start) %$% # magrittr exposition pipe
           rep(name, end-start+1)) # This produces directly the vector
}

regions1v <- get_vector(regions1)
regions2v <- get_vector(regions2)

regions2v[sites] = rep("Binding site", length(sites)) # Put the binding site sites in group 2, in whose locations there are NAs



ss <- ss %>%
  add_row(name=NA, start=1, end=6, .before=1) %>%
  add_row(name=NA, start=393, end=418)

ssv <- get_vector(ss)


# The UniProt sequence for which we have data for starts at position 6 and ends before this large snippet
# 6, 418-nchar("HREVPESETPVDTNLIELDTNDDDIVFEDFARQRLKGMKDDKDEEDDGTGSPHLNNR")=361

regions1v <- regions1v[6:361]
regions2v <- regions2v[6:361]
ssv <- ssv[6:361]