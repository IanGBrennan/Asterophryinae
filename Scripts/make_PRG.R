require(ips)

# these functions assume alignments are in FASTA format

# if the path directory has alignments that only include a single taxon/sample
make_PRG_single <- function(taxon=NULL, path=NULL, file.type=c(".fasta",".fas")){
  # change the wd
  setwd(path)
  
  # grab the locus names
  file.names <- dir(".", pattern=file.type)
  loci.names <- sapply(file.names, function(x) strsplit(x, file.type)[[1]][1])
  names(loci.names) <- NULL
  
  # concatenate all the sequences together
  new.file <- paste0(taxon,"_PRG.fasta")
  concat.call <- paste0("cat *",file.type, " >> ", new.file)
  system(concat.call)
  
  # read in the new file
  prg <- read.FASTA(file=new.file, type="DNA")
  names(prg) <- loci.names
  write.FASTA(prg, file=new.file)
}

# if the path directory has alignments that include multiple species/samples
make_PRG <- function(taxon, grep=F,
                     filetype=c(".fasta",".phylip"), replace.missing=T,
                     path=NULL, comment.missing=F,
                     same.length=T, sequence=c("DNA","AA")){
  # make a vector of the alignment files by name
  file.names <- dir(path, pattern = filetype)
  # keep the original taxon name
  taxon.og <- taxon
  
  # make a directory to hold the PRG(s)
  PRG.dir <- paste0(path, "/PRG")
  dir.create(PRG.dir)
  
  # make an empty file
  cat(NULL,file=paste0(PRG.dir, "/", taxon.og, "_All_Loci.fasta"))
  
  for (p in 1:length(file.names)) {
    print(paste("working on alignment", file.names[[p]]))
    # first remove the first line with the number of taxa and characters
    file <- paste(path, "/", file.names[p], sep="")
    short.file <- file.names[p]
    shortie <- strsplit(short.file, filetype)[[1]]
    #removed <- paste("sed -i '1d'", file)
    #system(removed)
    
    # read in the alignment and grab the sample names
    if(filetype==".fasta" && same.length==T){
      full_alignment <- read.dna(file, format="fasta", as.matrix=TRUE)
      alignment_names <- row.names(full_alignment)}
    if(filetype==".fasta" && same.length==F){
      if(sequence=="DNA"){}
      full_alignment <- read.FASTA(file, type=sequence)
      alignment_names <- names(full_alignment)}
    else if(filetype==".phylip"){#taxon <- paste0(taxon.og,"\t"); 
      full_alignment <- read.dna(file, format="sequential")
      alignment_names <- row.names(full_alignment)}
    
    if(!taxon == "longest" && grep == T){taxon <- alignment_names[grep(taxon.og,alignment_names)[1]]}
    
    if(taxon == "longest" && same.length == T){
      # get the number of missing bases per sample
      missing.sum <- apply(full_alignment, MARGIN = 1, FUN = function(x) length(which(as.numeric(x) %in% c(4, 240, 2))))
      new.taxon <- names(missing.sum[which(missing.sum == min(missing.sum))])
      if(length(missing.sum > 1)){new.taxon <- new.taxon[[1]]}
      
      # pull out just the new target taxon
      target.alignment <- full_alignment[new.taxon,]
      
      # rename the alignment with the locus name
      row.names(target.alignment) <- shortie
      
      # write the file (appending each new sequence)
      if(filetype==".fasta") {write.FASTA(target.alignment, paste0(PRG.dir, "/", taxon, "_All_Loci.fasta"), append=T)}
      else if(filetype==".phylip") {ips::write.fas(target.alignment, paste0(PRG.dir, "/", taxon.og, "_All_Loci.fasta"), append=T)}
    }
    
    if(taxon == "longest" && same.length == F){
      # figure out which sequence is longest
      seq.lengths <- sapply(full_alignment, length)
      longest.seq <- which(seq.lengths==max(seq.lengths))
      if(length(longest.seq)>1){longest.seq <- longest.seq[[1]]}
      
      # pull out just the new target sequences
      target.alignment <- full_alignment[longest.seq]
      
      # rename the alignment with the locus name
      names(target.alignment) <- shortie
      
      # write the file (appending each new sequence)
      if(filetype==".fasta") {write.FASTA(target.alignment, paste0(PRG.dir, "/", taxon, "_All_Loci.fasta"), append=T)}
      else if(filetype==".phylip") {ips::write.fas(target.alignment, paste0(PRG.dir, "/", taxon.og, "_All_Loci.fasta"), append=T)}
    }
    
    else if (!taxon == "longest"){
      if(!taxon %in% alignment_names){
        if(replace.missing==F && comment.missing==T){
          # tell us that the taxon is missing from the current alignment
          missing <- paste(taxon, "is not in alignment", shortie)
          print(missing)
        }
        else if(replace.missing==T){
          new.taxon <- alignment_names[[1]]
          
          # tell us that the taxon is missing from the current alignment, and what we'll use instead
          missing <- paste(taxon, "is not in alignment", shortie, "replacing with sequence from", new.taxon)
          print(missing)
          
          # pull out just the new target taxon
          target.alignment <- full_alignment[new.taxon,]
          
          # rename the alignment with the locus name
          row.names(target.alignment) <- shortie
          
          # write the file (appending each new sequence)
          if(filetype==".fasta") {write.FASTA(target.alignment, paste0(PRG.dir, "/", taxon, "_All_Loci.fasta"), append=T)}
          else if(filetype==".phylip") {ips::write.fas(target.alignment, paste0(PRG.dir, "/", taxon.og, "_All_Loci.fasta"), append=T)}
        }
        if(comment.missing==T){
          # write the info to file
          write.table(missing, file=paste0(path,"/",taxon,"_MISSING_LOCI.txt"), append=T, row.names=F, col.names=F)
        }
      }
      else if(taxon %in% alignment_names){
        # pull out just the target taxon
        if(same.length==T){target.alignment <- full_alignment[taxon,]}
        if(same.length==F){target.alignment <- full_alignment[taxon]}
        
        # rename the alignment with the locus name
        if(same.length==T){row.names(target.alignment) <- shortie}
        if(same.length==F){names(target.alignment) <- shortie}
        
        # write the file (appending each new sequence)
        if(filetype==".fasta") {write.FASTA(target.alignment, paste0(PRG.dir, "/", taxon.og, "_All_Loci.fasta"), append=T)}
        else if(filetype==".phylip") {ips::write.fas(target.alignment, paste0(PRG.dir, "/", taxon.og, "_All_Loci.fasta"), append=T)}
      }
    }
  }
}

# make_PRG(taxon="Austrochaperina_fryi_1", filetype=".fasta", 
#          replace.missing=F, path=getwd(), sequence="DNA", same.length=T)


# parallelized version for making lots of PRGs
PRGs_aln <- function(name.vec, aln.path=getwd(), filetype=c(".fasta", ".phylip"), 
                     same.length=T, sequence=c("DNA","AA")){
  #mcsapply(name.vec, function(x) make_PRG(taxon=x, replace.missing=F,
  #                                        path=aln.path, filetype = filetype,
  #                                        comment.missing=F), mc.cores=8)
  sapply(name.vec, function(x) make_PRG(taxon=x, replace.missing=F,
                                          path=aln.path, filetype = filetype,
                                          comment.missing=F, 
                                        same.length = same.length,
                                        sequence = sequence))
}
# make sure you supply name.vec as a vector of sample names!


# expand PRG to separate sequence files
expand_PRG <- function(prg, taxon=NULL){
  if(is.null(taxon)){stop("must provide a taxon name")}
  dir.create(paste0(getwd(),"/EXPANDED")) # make a directory for the new files
  for(k in 1:length(prg)){
    new.seq <- prg[k]
    names(new.seq) <- taxon
    write.FASTA(new.seq, file=file.path(getwd(),"EXPANDED",paste0(names(prg)[[k]],".fasta")))
  }
}

# An mc-version of the sapply function.
mcsapply <- function (X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE) {
  FUN <- match.fun(FUN)
  answer <- parallel::mclapply(X = X, FUN = FUN, ...)
  if (USE.NAMES && is.character(X) && is.null(names(answer))) 
    names(answer) <- X
  if (!isFALSE(simplify) && length(answer)) 
    simplify2array(answer, higher = (simplify == "array"))
  else answer
}

# Remove missing bases from a DNAbin sequence
deleteMissing <- function(DNAbin, nset = c("-","n","?")){
  iupac <- c(n = 240, `?` = 2, `-` = 4, r = 192, y = 48, s = 96, 
             w = 144, k = 80, m = 160, b = 112, d = 208, h = 176, 
             v = 224)
  nset <- iupac[nset]
  nset <- as.raw(nset)
  isNotEmpty <- function(x, nset) {
    ifelse(all(unique(x) %in% nset), FALSE, TRUE)
  }
  DNAbin[[1]] <- DNAbin[[1]][which(!DNAbin[[1]] %in% nset)]
  return(DNAbin)
}

# Remove missing data from a DNAbin object (many DNAbin sequences)
deleteEmptyCellsPRG <- function (PRGDNAbin, nset = c("-", "n", "?"), 
          quiet = FALSE) 
{
  if (!inherits(PRGDNAbin, "DNAbin")) 
    stop("'DNAbin' is not of class 'DNAbin'")

  # a simple for-loop to correct each sequence
  for(k in 1:length(PRGDNAbin)){PRGDNAbin[k] <- deleteMissing(PRGDNAbin[k], nset=c("-","n","?"))}
  # because I couldn't make the lapply work
  # testy <- lapply(1:length(PRGDNAbin), function(x) deleteMissing(PRGDNAbin[x], nset=nset))
  PRGDNAbin
}

 