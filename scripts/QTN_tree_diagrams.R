## author: Christian Brion 2019-12-18

## -----
## setup and load packages
library("plyr")
library("reshape2")
library("ape")
library("tidyverse")

base_dir <- "~/emacs/N-end_Rule_QTL_paper/variant_analysis/"
plots_dir <- paste0(base_dir, "plots/") 

tree <- read.tree(paste0(base_dir, "1011_matrix.tree.newick"))

## get information about the 1011 Sc strains, provided 
Strain1002 <- read.csv(paste0(base_dir, "1002Strains.txt"),
                       header = T, 
                       sep="\t", 
                       comment.char = "", 
                       quote = "")

## rotate one branch of the tree for better presentation
tree <- rotate(tree, 1126, polytom = c(1, 2))

## this makes the basic tree plot 
plot(tree, type = "unrooted",
     use.edge.length = TRUE, #ploting the tree as control
     node.pos = NULL,
     show.tip.label = FALSE,
     show.node.label = FALSE,
     edge.color = "black",
     edge.width = 1,
     edge.lty = 1, font = 3,
     cex = par("cex"), adj = NULL,
     srt = 0, no.margin = FALSE,
     root.edge = FALSE,
     label.offset = 0, underscore = FALSE,
     x.lim = NULL, y.lim = NULL, direction = "rightwards",
     lab4ut = NULL, tip.color = "black", plot = TRUE,
     rotate.tree = 0, open.angle = 0, node.depth = 1,
     align.tip.label = FALSE)


## -----
## find RM and BY on the tree
tree2 <- tree
tree2$tip.label[!(tree2$tip.label %in% c("AAA","ADT"))] <- " "

my_red <- "#B65050"
my_blu <- "#446F9A"
my_pur <- "#7C416B"

colbranch<-c(rep("lightgrey",nrow(tree2$edge)))
for (i in 1:nrow(tree2$edge)) {
  edgeN<-c(1:nrow(tree2$edge))[tree2$edge[,2]==i]
  if (i <= length(tree2$tip.label)) {
    if (tree2$tip.label[i]=="AAA") {
      colbranch[edgeN]<- my_red
    } else if (tree2$tip.label[i]=="ADT") {
      colbranch[edgeN]<- my_blu
    } else {
      colbranch[edgeN]<-"lightgrey"
    }
  }
}

## produce the tree with the position of the parental strains
plot(tree, type = "unrooted",
     use.edge.length = TRUE,
     node.pos = NULL,
     show.tip.label = FALSE,
     show.node.label = FALSE,
     edge.color = colbranch,
     edge.width = 1,
     edge.lty = 1, font = 3,
     cex = par("cex"), adj = NULL,
     srt = 0, no.margin = FALSE,
     root.edge = FALSE,
     label.offset = 0, underscore = FALSE,
     x.lim = NULL, y.lim = NULL, direction = "rightwards",
     lab4ut = NULL, tip.color = "black", plot = TRUE,
     rotate.tree = 0, open.angle = 0, node.depth = 1,
     align.tip.label = FALSE,main="parent")

## load the genotypes at the variants of interest
tree_data <- read.table(paste0(base_dir, "QTN_tree_data.csv"),
                        header = T, sep = ",", stringsAsFactors = F)
## load(paste0(base_dir, "IRA2_varGenosStringsOnly_210716.RData"))

## provide strains information to the table above
tree_table <- merge(Strain1002, tree_data, by.x = "Standardized.name", by.y = 0, all.x=T)
head(tree_table)

dcast(table2, Clades~var3) #allelic count per cluster

## add a joined vector for the three variants
table2$all3vars <- paste(table2$var3, table2$var4, table2$var6, sep="_")

## clade coutns for these:
dcast(table2[!str_detect(table2$all3vars, "\\."),], Clades~all3vars) #allelic count per cluster


#=======================
#coloring the edge of the tree to follow variant evolution

#QTNs <- c("var3", "var4", "var6")

QTNs <- colnames(tree_table)[20:50]

for(QTN in QTNs){

if(QTN %in%  colnames(tree_table)[20:50]){
    colbranch<-c(rep("grey", nrow(tree$edge)))
    for (i in 1:nrow(tree$edge)) {
        ## coloring the last edges (top of the
        ## tree) according to the allele of the
        ## site of interest in the corresponding strains
      edgeN<-c(1:nrow(tree$edge))[tree$edge[,2]==i]
      if (i <= length(tree$tip.label)) {
        if (is.na(tree_table[,QTN][tree_table$Standardized.name==tree$tip.label[i]])) {
          colbranch[edgeN] <- "grey"
        } else if (tree_table[,QTN][tree_table$Standardized.name==tree$tip.label[i]]=="0/0") {
          colbranch[edgeN] <- my_blu
        } else if (tree_table[,QTN][tree_table$Standardized.name==tree$tip.label[i]]=="1/1") {
          colbranch[edgeN] <- my_red
        } else if (tree_table[,QTN][tree_table$Standardized.name==tree$tip.label[i]]=="0/1") {
          colbranch[edgeN] <- my_pur
        }
      }
    }
}

    ## next section colors in branches
    ## according to their tips, if unambiguous
    ## print(date())
    for (j in 1:100) {
        ##coloring the edges within the tree according to
        ##the colors of the two out-coming edges, going
        ## down at each iteration (100 iterations)
  for (i in (length(tree$tip.label)+1):nrow(tree$edge)) {
    edgeN<-c(1:nrow(tree$edge))[tree$edge[,2]==i]
    if (length(levels(as.factor(colbranch[tree$edge[,1]==i])))==1) {
        ## control same color of the two out-coming edges
        colbranch[edgeN]<-levels(as.factor(colbranch[tree$edge[,1]==i]))
        ## give this color to the inner-edge (if not, stay black)
    }
  }
}
## print(date()) runs in 10–15secs

    ## generating the tree with the edges colored
    ## by allele version in a svg format (for Inkskape)
    pdf(paste0(plots_dir, "tree_",
               strsplit(QTN, "X")[[1]][2],
               ".pdf"),
        width = 10, height = 10)
    
plot(tree, type = "unrooted", use.edge.length = TRUE,
     node.pos = NULL, show.tip.label = FALSE, show.node.label = FALSE,
     edge.color = colbranch, edge.width = 0.75, edge.lty = 1, font = 3,
     cex = par("cex"), adj = NULL, srt = 0, no.margin = FALSE,
     root.edge = FALSE, label.offset = 0, underscore = FALSE,
     x.lim = NULL, y.lim = NULL, direction = "rightwards",
     lab4ut = NULL, tip.color = "black", plot = TRUE,
     rotate.tree = 0, open.angle = 0, node.depth = 1,
     align.tip.label = FALSE,
     main = strsplit(QTN, "X")[[1]][2])
## legend("bottomright", inset=0.3,
##     legend=c("0/0_0/0_0/0", "1/1_1/1_1/1", "1/1_1/1_0/0", "0/0_1/1_1/1", "0/0_1/1_0/0", "0/1_1/1_1/1", "0/1_1/1_0/0", "0/1_1/1_0/1", "0/1_0/1_0/1", "other")
##     ,fill=c("blue", "red", "firebrick1", "royalblue", "royalblue3", "hotpink", "hotpink1", "hotpink2", "hotpink3", "grey"))
dev.off()
}
