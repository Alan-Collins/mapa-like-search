require(RCurl)
require(dplyr)

########LIST OF FUNCTIONS###########


pNCBISummary=function(input.file){ #this partitions the NCBI summary file into organisms and protein accesion. I compare the organism lists to find LapGD organisms. Just selecting organisms probably would have worked just fine in retrospect.
  data.input<-readLines(input.file)
  input.lines<-which(grepl("[", data.input, fixed=TRUE)==TRUE)
  access.lines<-input.lines+2
  access.lines<-strsplit(data.input[access.lines]," ", fixed=TRUE)
  seperate<-strsplit(data.input[input.lines], "[", fixed=TRUE)
  
  ORG=NULL
  ACCESS=NULL
  
  for (i in 1:length(access.lines)){
    org<-substr(seperate[[i]][2],0,nchar(seperate[[i]][2])-1)
    ORG<-c(ORG,org)
    access<-access.lines[[i]][1]
    ACCESS<-c(ACCESS,access)
  }
  out<-data.frame(ORG,ACCESS,stringsAsFactors = FALSE)
  
  return(out)
}


searchable.faa=function(faa.file){ #converts genome faa file into a table with the accession, AA squence, putative function, and AA size. can pull info from this table about RTX status, cleavage site, etc
  test1<-readLines(faa.file)
  inputs<-which(substr(test1,0,1)==">") 
  
  ACC<-NULL
  SEQQ<-NULL
  pFUN<-NULL
  SIZE<-NULL
  for (i in 1:length(inputs)){
    slp.inputs<-gregexpr(" ", test1[inputs[i]], fixed=TRUE)
    species.split<-gregexpr("[", test1[inputs[i]], fixed=TRUE)
    acc<-substr(test1[inputs[i]], 2, slp.inputs[[1]][1]-1)
    p.fun<-substr(test1[inputs[i]],slp.inputs[[1]][1]+1, species.split[[1]]-1)
    first<-inputs[i]+1
    
    if (i+1>length(inputs)){
      last<-length(test1)
    } else { last<-inputs[i+1]-1 }
    
    seq=NULL
    
    for(a in first:last){
      seq<-paste(seq, test1[a], sep="")
      size<-nchar(seq)
    }
    
    ACC<-c(ACC,acc)
    pFUN<-c(pFUN,p.fun)
    SEQQ<-c(SEQQ,seq)
    SIZE<-c(SIZE,size)
  }
  
  tbl<-data.frame(ACC,pFUN, SEQQ, SIZE, stringsAsFactors = FALSE)
  names(tbl)<-c("accession", "putative.function","protein.sequence", "protein.length")
  return(tbl)
}

verify.files=function(input.url){
  require("RCurl")
  temp.url<-paste(input.url, "/", sep="")
  getfiles<-getURL(temp.url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
  getfiles<-strsplit(getfiles,"\n")
  sp<-strsplit(input.url,"/", fixed=TRUE)
  get.l<-length(sp[[1]])
  temp.protein<-paste(sp[[1]][get.l], "_protein.faa.gz", sep="")
  temp.table<-paste(sp[[1]][get.l], "_feature_table.txt.gz", sep="")
  find.table<-any(getfiles[[1]]==temp.table)
  find.protein<-any(getfiles[[1]]==temp.protein)
  if (find.table==TRUE & find.protein==TRUE){
    return(TRUE)
  } else { return(FALSE)}
}

s.url=function(url){
  #  verify.files(url)
  sp<-strsplit(url,"/", fixed=TRUE)
  get.l<-length(sp[[1]])
  new.url<-paste(url, "/", sp[[1]][get.l], sep="")
  return(new.url)
}

any.RTX=function(y){ #Retursn proteins with RTX motifs
  #D-x-[LI]-x(4)-G-x-D-x-[LI]-x-G-G-x(3)-D
  RTX<-c("D[[:alpha:]]L[[:alpha:]][[:alpha:]][[:alpha:]][[:alpha:]]G[[:alpha:]]D[[:alpha:]]L[[:alpha:]]GG[[:alpha:]][[:alpha:]][[:alpha:]]D",
         "D[[:alpha:]]I[[:alpha:]][[:alpha:]][[:alpha:]][[:alpha:]]G[[:alpha:]]D[[:alpha:]]L[[:alpha:]]GG[[:alpha:]][[:alpha:]][[:alpha:]]D",
         "D[[:alpha:]]L[[:alpha:]][[:alpha:]][[:alpha:]][[:alpha:]]G[[:alpha:]]D[[:alpha:]]I[[:alpha:]]GG[[:alpha:]][[:alpha:]][[:alpha:]]D",
         "D[[:alpha:]]I[[:alpha:]][[:alpha:]][[:alpha:]][[:alpha:]]G[[:alpha:]]D[[:alpha:]]I[[:alpha:]]GG[[:alpha:]][[:alpha:]][[:alpha:]]D")

  checker<-lapply(RTX, function(x) grepl(x,y))
  results<-any(checker==TRUE)
  return(results)
}

find.cleavage=function(anchor){ ##returns T or F to whateve AA string you feed it. In the script its residues 80-150 of proteins that have RTX motifs and are at least 1500 AA long
  
  TAAG<-grepl("TAAG", anchor)
  PAAG<-grepl("PAAG", anchor)
  AAAG<-grepl("AAAG", anchor)
  TAAV<-grepl("TAAV", anchor)
  
  
  if (TAAG==TRUE|PAAG==TRUE|AAAG==TRUE|TAAV==TRUE){
    return(TRUE)} else { return(FALSE) }
}


##################################
##################################
##################################
##################################

#######STEP 1 PROCESS NCBI SUMMARY FILES####### 
  #Download LapD_MoxY_N domain containing strain summary list: https://www.ncbi.nlm.nih.gov/sites/entrez?LinkName=cdd_protein&db=cdd&cmd=Link&from_uid=318615
    #8700 hits as of 10-20-2018
  
  #Download LapG comtaining strain summary list: https://www.ncbi.nlm.nih.gov/sites/entrez?LinkName=cdd_protein&db=cdd&cmd=Link&from_uid=310549
    #20685 hits as of 10-20-2018

  ###if you download another format the pNCBI summary function won't work###

####CONVERT SUMMARY FILES#####
pSUM.lapD<-pNCBISummary("Data/LapD.txt")
pSUM.lapG<-pNCBISummary("Data/LapG.txt")

# 
##############################

#######STEP 2 - Identify LapGD-encoding bugs, do a little analysis########
lapD.organism.tally<- length(unique(pSUM.lapD$ORG)) #how many unique genomes have LapD-like protein? 4024
lapG.organism.tally<- length(unique(pSUM.lapG$ORG)) #how many unique genomes have LapG-like protein? 5545
lapGD.organisms<-unique(intersect(pSUM.lapG$ORG,pSUM.lapD$ORG)) #select organisms found on both lists. 3592 automatically removes replcates so no need for unique function.
length(lapGD.organisms)
lapD.percent.overlap<-(length(lapGD.organisms)/lapD.organism.tally)*100 #How many LapD-encoding genomes also encode LapG? 96%
#############################################################


#######STEP 3 - Download Genome file#######S
all.genomes<-read.csv("Data/ALL_genomes_proks.csv", stringsAsFactors = FALSE, header=TRUE) #file contains all bacterial genomes as of 10/30/2018


#######STEP 4 #######

###################GET FTP INFO FOR LAPGD BUGS###################

#############create dataframe containing all the lapD/G encoding species with subgroup and sequence info################

target.rows<-lapply(lapGD.organisms, function(x) which(all.genomes$X.Organism.Name==x))
get.rows<-sapply(target.rows, function(x) x[1]) #returns only one value per 
get.rows<-na.omit(get.rows)

lapd.g.species.info<-data.frame(all.genomes[get.rows,c(1,7,18,21)]) #create dataset of only organisms and info of interest
names(lapd.g.species.info)<-c("organism", "subgroup", "level", "RefSeq.FTP") #name cols

write.csv(lapd.g.species.info, "Output/lapDG_Species_info.csv")



#######STEP 5#######

#############create data frame containing species names and protein and table source files######################

# https://www.ncbi.nlm.nih.gov/assembly/GCA_000828895.1/GCA_000828895.1_protein.faa.gz

# GCA_000466945.1


tbl<-sapply(lapd.g.species.info$RefSeq.FTP, function(x) paste(s.url(x), "_feature_table.txt.gz", sep=""))
prt<-sapply(lapd.g.species.info$RefSeq.FTP, function(x) paste(s.url(x), "_protein.faa.gz", sep=""))
new.d<-data.frame(lapd.g.species.info,tbl,prt)
write.csv(new.d,"Output/dg.genome.url.files.csv")






#####download FAA files######
  ##THE NCBI DB DOESNT LIKE IT WHEN YOU TRY TO DOWNLOAD THOUSANDS OF GEMONES. THERES PROBABLY A BATCH WAY TO DO THIS, BUT THIS IS HOW THE SCRIPT WORKS.
  ##THE DOWNLOAD WILL PROBABLY TIME OUT. I JUST RESTART THE FOR LOOP BELOW AT WHERE EVER IT TIMED OUT.
errors <- list()

for (i in 1:length(new.d$organism)){ #
  filename = paste("./Output/Genomes/", new.d$organism[i], "_protein.faa.gz", sep="")
  if(file.exists(filename)){
    next
  }
  tryCatch({
    message("Downloading: ", new.d$prt[i],"\n")
    download.file(as.character(new.d$prt[i]), destfile = filename)   #download genome faa
 #   Sys.sleep(10)
    message("Download Complete.\n\n")
  },error = function(e){ 
    errors <<- c(errors, new.d$organism[i])
    })
}



#I REMOVED THE SECTION THAT ALSO DOWNLOADED THE FEATURE TABLE.

#######STEP 6 #######

######ONCE YOU HAVE DOWNLOADED THE GENOME FAA FILES YOU CAN LOOK THROUGH THEM. THIS WILL SELECT LARGE PROTEINS WITH RTX MOTIFS AND A LAPG CLEAVAGE SITE###########
#######GENERATE LIST OF LARGE RTX PROTEINS. THEN LOOK FOR CLEAVAGE. SEPERATE BY CLEAVAGE STATUS###############
setwd("/Users/AlanCollins/Github/mapa-like-search/Output/Genomes/")
faa.files<-list.files() 
my.table<-data.frame(ACC=character(), pFUN=character(), SEQQ=character(), SIZE=numeric(), stringsAsFactors=FALSE) 

ptm <- proc.time() #for fun I timed this
for (i in 1:length(faa.files)){
#  limit.size = NULL
  closeAllConnections()
  species.name<-strsplit(faa.files[i], "_protein.faa.gz") #get species name
  species.name<-species.name[[1]][1]
  
  message("processing....", species.name, " faa file\t", i, " of ", length(faa.files))
  
  f.proteins<-searchable.faa(gzfile(faa.files[i])) #convert faa file to searchable format.
  
  RTX.status<-lapply(f.proteins$protein.sequence, function(x) any.RTX(x)) #look through protein sequence column for proteins with RTX
  
  if (any(RTX.status==TRUE)==TRUE) { ##only interogate genomes with large RTX adhesins.
    putative.RTX<-which(RTX.status==TRUE)
    limit.size<-subset(f.proteins[putative.RTX,], protein.length>=1000)
    org.rep<-rep(species.name, length(limit.size$protein.length))
    limit.size<-data.frame(organism=org.rep,limit.size)
  }
  
  my.table<-rbind(my.table, limit.size) 
  message("complete\n\n")
  
}
motif_status<-lapply(my.table$protein.sequence, function(x) find.cleavage(substr(x, 80,150))) ##search list of large RTX proteins for LapG cleavage site. Return T or F for putative cleavage site.
new_my.table<-mutate(my.table, lapG=motif_status, loose.lapG =loose_motif_status, loosest.lapG = loosest_motif_status)
#write.xlsx2(new_my.table, "LapA-Like_V3.xlsx", append=TRUE, sheetName="All Large RTX")
closeAllConnections()
proc.time()-ptm #it takes a pretty long time
substrates2<-subset(new_my.table[which(new_my.table$lapG==TRUE),1:5], stringsAsFactors=FALSE)
no.substrates2<-subset(new_my.table[which(new_my.table$lapG==FALSE),])
#write.xlsx2(substrates2, "LapA-Like_V3.xlsx", append=TRUE, sheetName="Substrates")
#write.xlsx2(no.substrates2, "LapA-Like_V3.xlsx", append=TRUE, sheetName="Not_Substrates")

###########################

#write.csv(substrates2, file = "all_lapG_substrates.csv")

substrates.unique <- substrates2[!duplicated(substrates2$accession),]
substrates.unique$organism <- as.character(substrates.unique$organism)
row.names(substrates.unique) <- seq(length(row.names(substrates.unique)))
protein.count <- data.frame(table(unlist(substrates.unique$organism)))


multiples <- protein.count[protein.count$Var1[which(protein.count$Freq >1)],]

doubles <-protein.count[protein.count$Var1[which(protein.count$Freq == 2)],]

singles <- protein.count[protein.count$Var1[which(protein.count$Freq == 1)],]

multiple.rtx <- substrates.unique[which(substrates.unique$organism == multiples$Var1[1]),]

double.rtx <- substrates.unique[which(substrates.unique$organism == doubles$Var1[1]),]

single.rtx <- substrates.unique[which(substrates.unique$organism == singles$Var1[1]),]


for ( i in 2: length(multiples$Var1)){
  multiple.rtx <- rbind(multiple.rtx, substrates.unique[which(substrates.unique$organism == multiples$Var1[i]),])
  
}

for ( i in 2: length(doubles$Var1)){
  double.rtx <- rbind(double.rtx, substrates.unique[which(substrates.unique$organism == doubles$Var1[i]),])
  
}

for ( i in 2: length(singles$Var1)){
  single.rtx <- rbind(single.rtx, substrates.unique[which(substrates.unique$organism == singles$Var1[i]),])
  
}

pseudomonas.multi <- multiple.rtx[ grepl( "Pseudomonas" , multiple.rtx$organism ), ]
pseudomonas.single <- single.rtx[ grepl( "Pseudomonas" , single.rtx$organism ), ]

multiple.rtx$Adhesin.num <- "z"

multi.organisms <- as.character(levels(as.factor(multiple.rtx$organism)))

length(multi.organisms)

for (i in 1:length(multi.organisms)){
  org <- multi.organisms[i]
  sizes <- NULL
  sizes <- multiple.rtx$protein.length[which(multiple.rtx$organism == org)]
  sizes <- sizes[order(sizes, decreasing = TRUE)]
  for ( j in 1:length(sizes)){
    multiple.rtx$Adhesin.num[which(multiple.rtx$organism == org & multiple.rtx$protein.length == sizes[j])] <- j
  }

}

biggest.rtx <- multiple.rtx[which(multiple.rtx$Adhesin.num == 1),]


doubles$ratio <- 0
doubles$biggest <- 0
doubles$smallest <- 0
for(i in 1: length(doubles$Var1)){
  sizes <- double.rtx$protein.length[which(double.rtx$organism == doubles$Var1[i])]
  sizes <- sizes[order(sizes, decreasing = TRUE)]
  doubles$ratio[i] <- sizes[1]/sizes[2]
  doubles$biggest[i] <- sizes[1]
  doubles$smallest[i] <- sizes[2]
  
}

Shewanella.multi <- multiple.rtx[ grepl( "Shewanella" , multiple.rtx$organism ), ]
Shewanella.single <- single.rtx[ grepl( "Shewanella" , single.rtx$organism ), ]

Vibrio.multi <- multiple.rtx[ grepl( "Vibrio" , multiple.rtx$organism ), ]
Vibrio.single <- single.rtx[ grepl( "Vibrio" , single.rtx$organism ), ]

###
###
###
### Plots
###
###
###
###


library(ggplot2);library(reshape2)

ggplot(data = multiple.rtx, aes(protein.length, fill = 'multiple RTX')) + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +  
  geom_density(alpha=0.25) + scale_x_continuous(limits = c(0,16500), expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  geom_density(data = single.rtx, aes(protein.length, fill = 'single RTX'), alpha=0.25) + 
  labs(title = "All substrate size distributions", x = "Protein Length", y = "Density", fill = "") 
  ggsave("overall substrate sizes.png", width = 20, height = 10, units = "cm", dpi = 300)

ggplot(data = multiple.rtx, aes(Adhesin.num, protein.length, fill = Adhesin.num))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), legend.position="none",
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_boxplot() +
  geom_jitter(aes(alpha = 0.1), color = "lightblue")+
  ggsave("All RTX proteins boxplot by size rank in org.png")


ggplot(data = doubles, aes(biggest, ratio))+
  geom_jitter()

ggplot(data = doubles, aes(biggest, smallest)) +
  geom_jitter()

ggplot(data = doubles, aes("a", ratio)) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  geom_boxplot()+
  geom_jitter(aes(alpha = 0.1), color = "lightblue")+
  ggsave("Ratios in double RTX organisms.png")

ggplot(data = pseudomonas.multi, aes(protein.length, fill = 'multiple RTX')) + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +  
  geom_density(alpha=0.25) + scale_x_continuous(limits = c(0,12000), expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +  
  geom_density(alpha=0.25)+
  geom_density(data = pseudomonas.single, aes(protein.length, fill = 'single RTX'), alpha=0.25) + 
  labs(title = "Pseudomonas substrate size distributions", x = "Protein Length", y = "Density", fill = "") +
  ggsave("Pseudomonas substrate sizes.png", width = 20, height = 10, units = "cm", dpi = 300)

ggplot(data = Shewanella.multi, aes(protein.length, fill = 'multiple RTX')) + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +  
  geom_density(alpha=0.25) + scale_x_continuous(limits = c(0,16500), expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +  
  geom_density(alpha=0.25)+
  geom_density(data = Shewanella.single, aes(protein.length, fill = 'single RTX'), alpha=0.25) + 
  labs(title = "Shewanella substrate size distributions", x = "Protein Length", y = "Density", fill = "")+ 
  ggsave("Shewanella substrate sizes.png", width = 20, height = 10, units = "cm", dpi = 300)

ggplot(data = Vibrio.multi, aes(protein.length, fill = 'multiple RTX')) + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +  
  geom_density(alpha=0.25) + scale_x_continuous(limits = c(0,16500), expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +  
  geom_density(alpha=0.25)+
  geom_density(data = Vibrio.single, aes(protein.length, fill = 'single RTX'), alpha=0.25) + 
  labs(title = "Vibrio substrate size distributions", x = "Protein Length", y = "Density", fill = "")+ 
  ggsave("Vibrio substrate sizes.png", width = 20, height = 10, units = "cm", dpi = 300)
