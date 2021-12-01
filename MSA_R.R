
##Install LATEX
install.packages("rmarkdown", dep = TRUE)


###Install MSA
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("msa")

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("seqinr")

##Load Library
library(msa)
library(seqinr)

system.file("tex", "texshade.sty", package="msa")

#Example MSA

##WORKING MSA command
mySequenceFile <- system.file("examples", "exampleAA.fasta", package="msa")
mySequences <- readAAStringSet(mySequenceFile)
mySequences

myFirstAlignment <- msa(mySequences)
myFirstAlignment 

msaPrettyPrint(myFirstAlignment, alFile="./myFile.fasta", output="pdf", showNames="none",
               showLogo="none", askForOverwrite=FALSE, verbose=TRUE)

##READING SEQ in files?
bdir = "C:/Users/User/Desktop/GitHub/Bioinformatic/Publication"
workdir = file.path(bdir,"MSA_Test/")
outdir = file.path(bdir,"Alignements")

myFiles = list.files(path = workdir,full.names = TRUE)
for (i in 1:66){
  newdir = substr(basename(myFiles[i]), start = 1 ,stop = 6)
  output = file.path(outdir, newdir)
  dir.create(output,showWarnings = FALSE)
  current_seq <- readAAStringSet(myFiles[i])
  current_align <- msa(current_seq)
  fasta_file = paste(outdir,'/',newdir,".fasta",sep="")
  tex_file = paste(outdir,"/",newdir,".tex",sep="")
  print(pdf_file)
  msaPrettyPrint(current_align, output="tex", showConsensus = "none", askForOverwrite=FALSE, verbose=FALSE,
                 file = tex_file, alFile = fasta_file )
  texi2pdf(tex_file, clean=TRUE)
}


###ClustalW or Clustal OMEGA or MUSCLE
#ClustalW
myClustalWAlignment <- msa(mySequences, "ClustalW")
myClustalWAlignment

#Clustal Omega
myClustalOmegaAlignment <- msa(mySequences, "ClustalOmega")
myClustalOmegaAlignment

#MUSCLE
myMuscleAlignment <- msa(mySequences, "Muscle")
myMuscleAlignment


#Print
help("print,MsaDNAMultipleAlignment-method")

print(myFirstAlignment)
print(myFirstAlignment, show="complete")
print(myFirstAlignment, showConsensus=FALSE, halfNrow=3)
print(myFirstAlignment, showNames=FALSE, show="complete")


#To set row or column masks an IRanges object must be supplied:
myMaskedAlignment <- myFirstAlignment
colM <- IRanges(start=1, end=100)
colmask(myMaskedAlignment) <- colM
myMaskedAlignment 
unmasked(myMaskedAlignment)


conMat <- consensusMatrix(myFirstAlignment)
dim(conMat)
conMat[, 101:110]

conMat <- consensusMatrix(myMaskedAlignment)
conMat[, 95:104]

#Consensus Sequences and Consertavtion Scores
printSplitString <- function(x, width=getOption("width") - 1)
{
  starts <- seq(from=1, to=nchar(x), by=width)
  for (i in 1:length(starts))
    cat(substr(x, starts[i], starts[i] + width - 1), "\n")
}
printSplitString(msaConsensusSequence(myFirstAlignment))


printSplitString(msaConsensusSequence(myFirstAlignment, type="upperlower",
                                      thresh=c(40, 20)))

printSplitString(msaConsensusSequence(myMaskedAlignment, type="upperlower",
                                      thresh=c(40, 20)))
data(BLOSUM62)
msaConservationScore(myFirstAlignment, BLOSUM62)

msaConservationScore(myFirstAlignment, BLOSUM62, gapVsGap=0,
                     type="upperlower", thresh=c(40, 20))

msaConservationScore(myMaskedAlignment, BLOSUM62, gapVsGap=0,
                     type="upperlower", thresh=c(40, 20))

#Interface with others Packages
hemoSeq <- readAAStringSet(system.file("examples/HemoglobinAA.fasta",
                                       package="msa"))
hemoAln <- msa(hemoSeq)
hemoAln

hemoAln2 <- msaConvert(hemoAln, type="seqinr::alignment")
