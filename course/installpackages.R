# usually you can just do this
source("https://bioconductor.org/biocLite.R")
biocLite("KEGG.db")
library("KEGG.db")

# before you proceed further, set up your working directory
# A. check your current directory
  getwd()
  # the installation will not work if you are running it on the network drives
  # such as '\\icnas3.cc.ic.ac.uk\username' or 'H:\\username\'

# B. set your working directory
  setwd("C:\\Users\\yfwang\\Desktop\\whateverfolderyouwant")

# However, we don't have the administor permission of the lib computer
# We have to procedd this option instead

# step (1) decide where you are going to install this package
  mydir="C:\\Users\\yfwang\\Desktop\\R_packages\\"

# step (2) 
  # download the packages into the folder that you created in step (1) 
  # (download A and B first,and hopefully you don't need to download C)
# A. link for the KEGG.db
  # http://bioconductor.org/packages/3.3/data/annotation/html/KEGG.db.html
# B. link for the goseq  
  # http://bioconductor.org/packages/3.3/bioc/html/goseq.html
# C. link for the org.Mm.eg.db
  # http://bioconductor.org/packages/3.3/data/annotation/html/org.Mm.eg.db.html

# step (3) install the KEGG.db first  
  
  keggloc<-paste0(mydir,"KEGG.db_3.2.3.tar.gz")
  install.packages(goseqloc, lib=mydir, repos = NULL, type = "source")
  library(KEGG.db, lib.loc=mydir)
 
# step (4) install the goseq
  source("https://bioconductor.org/biocLite.R")
  biocLite("goseq")
  library("goseq")
  # if it works, then congrats!
  # if not, try the following 
  goseqloc<-paste0(mydir,"goseq_1.24.0.tar.gz")
  install.packages(goseqloc, lib=mydir, repos = NULL, type = "source")
  library(goseq, lib.loc=mydir)

  
  orgmmloc<-paste0(mydir,"org.Mm.eg.db_3.3.0.tar.gz")
  install.packages(orgmmloc, lib=mydir, repos = NULL, type = "source")
  library(org.Mm.eg.db, lib.loc=mydir)
  

