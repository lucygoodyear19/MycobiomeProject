########################################################################################
#################################### Subset FastTree ###################################
########################################################################################


# Author: Lucy Goodyear (lucy.goodyear19@imperial.ac.uk)
# Version: 0.0.1

# clear workspace
rm(list=ls())


########################################################################################
###################################### Set up ##########################################


# load packages
library(phyloseq)
library(Biostrings)
library(phytools)
library(adephylo) # for listTips() function

# load filtered phyloseq object for both DADA2 pipeline
dada2 <- readRDS(paste0(dada2_data_path, "physeqob_DADA2.rds"))

# load FastTree tree
tree.unrooted <- read_tree(paste0(dada2_data_path, "Tree/fasttree_tree"))

longest_edge <- which.max(tree$edge.length)
long_nodes <- tree$edge[longest_edge,]     #this should usually include one tip
new_outgroup <- long_nodes[long_nodes < Ntip(tree)]
tree_rooted <- root(tree, outgroup=new_outgroup, resolve.root=TRUE)


pick_new_outgroup <- function(tree.unrooted){
  require("magrittr")
  require("data.table")
  require("ape") # ape::Ntip
  # tablify parts of tree that we need.
  treeDT <- 
    cbind(
      data.table(tree.unrooted$edge),
      data.table(length = tree.unrooted$edge.length)
    )[1:Ntip(tree.unrooted)] %>% 
    cbind(data.table(id = tree.unrooted$tip.label))
  # Take the longest terminal branch as outgroup
  new.outgroup <- treeDT[which.max(length)]$id
  return(new.outgroup)
}
new.outgroup <- pick_new_outgroup(tree)
tree_rooted <- ape::root(tree, outgroup=new.outgroup)#, resolve.root=TRUE)

# find longest branch length and use furthest ASV as outgroup for root
#root_asv <- tree$tip.label[max(tree$edge.length)]
# root tree
#tree_rooted <- root(tree, root_asv, resolve.root = TRUE)

dada2 <- merge_phyloseq(dada2, 
                        phy_tree(tree_rooted))


###################################################################################
######################### Subset phylogenetic tree ################################


# subset tree via (mostly) taxonomically independent method

# create tree subsetting by order
chytrids <- subset_taxa(dada2, Order == "o__Rhizophydiales")

# find the tip labels at each end
longest_edge_chyt <- max(phy_tree(chytrids)$edge.length)
nm <- phy_tree(chytrids)$tip.label[phy_tree(chytrids)$edge.length == longest_edge_chyt]
long_nodes_chyt <- phy_tree(chytrids)$edge[longest_edge_chyt,] # this should usually include one tip
new_outgroup_chyt <- long_nodes_chyt[long_nodes_chyt < Ntip(phy_tree(chytrids))]
rhizo_node <- phy_tree(chytrids)$tip.label

mrca_node <- getMRCA(as.phylo(phy_tree(chytrids)), c("ASV34440", two))

tips <- as.data.frame(listTips(phy_tree(chytrids))[[1]])
names(tips) <- "node_nums"
one <- rownames(tips)[tips$node_nums == 1]
two <- rownames(tips)[tips$node_nums == nrow(tips)]

# use ASV ends to find MRCA node number
mrca_node <- getMRCA(as.phylo(phy_tree(dada2)), c(one, "ASV34323"))
# extract clade from node number
clade_not_chytrids <- extract.clade(as.phylo(phy_tree(dada2)), mrca_node)
try <- tree.extract.clade(as.phylo(phy_tree(dada2)), mrca_node)
# subset dada2 by clade
chytrids_final <- subset_taxa(dada2, taxa_names(dada2) %in% clade_not_chytrids$tip.label)

# save chytrids to separate phyloseq object
saveRDS(chytrids_final, paste0(results_path, "Tree/physeqob_DADA2_chytrids.rds"))


phy <- phy_tree(chytrids)
rootID <- length(phy$tip.label) + 1;
counter <- 1;
for (i in rootID:max(phy$edge[,1])) {
  clade <- extract.clade(phy, i);
  # do something. just printing clade properties here
  print(paste0("clade ", counter, " (node ", i, ") has ",
               length(clade$tip.label), " tips."));
  counter <- counter + 1;
}


## end of script
