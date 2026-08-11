# truei_trees.r -----------------------------------------------------------
#
# MAB 9/15/25
# Code for phylogenetic tree visualization for the following manuscript:
# Becker MA, Bell KC, Lopez-Ortiz S, Maldonado JE, Edwards CW, Castañeda-Rico S. (2026) 
# Genomic detection of a novel Peromyscus truei invasion on San Clemente Island, California. 
# Therya 17:277-294. doi: 10.12933/therya.2026.6265
#

# load packages -----------------------------------------------------------
library(BiocManager)
library(treeio)
library(tidyverse)
library(ggtree)
library(ggpubr)
library(tidytree)
library(patchwork)
library(glue)

# cytb_tips -------------------------------------------------------------

cytb_tips <- read.csv("C:/Users/madel/OneDrive/Documents/Therya_truei_ms/analyses/subset_truei_tips.csv") %>% 
  select(c("Sample","Group"))

cytb_labs <- c(expression(paste(italic("P. truei")," II CA")), 
               expression(paste(italic("P. truei")," II ID")),
               expression(paste(italic("P. truei")," II NV")),
               expression(paste(italic("P. truei")," II UT")),
               "San Clemente Island", 
               "Outgroups")

cytb_tips$Group <- factor(cytb_tips$Group, levels=c("CA","ID","NV","UT","SCL","OUT"))

# MrBayes_cytb ------------------------------------------------------------

BI_cytb <- read.newick("C:/Users/madel/OneDrive/Documents/Therya_truei_ms/analyses/rooted_mrbayes_truei_cytb_unique_aln_9_3.contree",node.label = 'support') %>% 
  tree_subset(node = 82, levels_back = 0)

BI_cytb@data$support <- as.integer(BI_cytb@data$support)


ggtree(BI_cytb) %<+% cytb_tips %>% 
  
  # collapse truei I
  scaleClade(76, .4) %>%collapse(76, 'max', color = "black", fill = "white") +
  
  # add padding on the right
  #xlim(NA, 0.0001) +
  
  # tip labels
  #geom_tiplab() +
  # get node numbers
  #geom_text2(aes(subset=isTip, label=node, fontsize=.5), hjust=-.3) +

  # support
  geom_point2(aes(subset=(support > 95)), shape=20, size=2, colour="black") +
  
  geom_cladelab(node=93, label="P. truei", fontface = 3, barsize = 1, barcolour = "black", fontsize = 4.5, offset = 2.4, offset.text =0.15 , extend=0.5, align = TRUE) +
  geom_cladelab(node=93, label="              II", barsize = 1, barcolour = NA, fontsize = 4.5, offset = 2.4, align = TRUE) +
  geom_cladelab(node=76, label="P. truei", fontface = 3, barcolour = NA, fontsize = 4.5, align = TRUE) +
  geom_cladelab(node=76, label="             I", barcolour = NA, fontsize = 4.5, align = TRUE) +
  geom_cladelab(node=107, label="P. gratus", fontface = 3, barcolour = NA, fontsize = 4.5, align = TRUE) +

  geom_cladelab(node=69, label="LACM 5833", barcolour = NA, fontsize = 4.5, align = TRUE) +
  geom_cladelab(node=70, label="LACM 5834", barcolour = NA, fontsize = 4.5, align = TRUE) +
  
  
  # tip colors and shapes
  geom_tippoint(aes(color=Group,shape=Group), size=2) + 
  scale_color_manual(name = "Taxon/Locality", values = c("#66CCEE","#228833", "#AA3377","#CCBB44","#EE6677", "#BBBBBB"),
                     labels = cytb_labs) +
  scale_shape_manual(name = "Taxon/Locality", values=c(16,16,16,16,17,15), 
                     labels = cytb_labs) +
  
  xlim(0, 13) +
  theme(legend.position = c(0.2,0.78),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) 



# IQtree_cytb -------------------------------------------------------------

ML_cytb <- read.iqtree("C:/Users/madel/OneDrive/Documents/Therya_truei_ms/analyses/subset_truei_cytb_unique_aln_9_3.nex.contree") 

ggtree(ML_cytb) %<+% cytb_tips %>%  
  
  # collapse truei I
  scaleClade(117, .4) %>%collapse(117, 'max', color = "black", fill = "white") +
  
  # add padding on the right
  xlim(NA, 0.15) +
  
  # get node numbers
  #geom_text2(aes(subset=!isTip, label=node, fontsize=.5), hjust=-.3) 
  
  # get tip numbers
  #geom_text2(aes(subset=isTip, label=node, fontsize=.5), hjust=-.3) 
  
  geom_cladelab(node=76, label="P. truei", fontface = 3, barsize = 1, barcolour = "black", offset = 0.01, fontsize = 4.5) +
  geom_cladelab(node=76, label="            II", barsize = 1, barcolour = "black", offset = 0.01, fontsize = 4.5) +
  geom_cladelab(node=117, label="P. truei", fontface = 3, barcolour = NA, offset = 0.012, fontsize = 4.5) +
  geom_cladelab(node=117, label="  I", barcolour = NA, offset = 0.025, fontsize = 4.5) +
  geom_cladelab(node=144, label="P. gratus", fontface = 3, barcolour = NA, offset = 0.003, fontsize = 4.5) +
  geom_cladelab(node=111, label="LACM", barcolour = NA, offset = - 0.008, vjust = 1, hjust = 0, fontsize = 4.5) +
  geom_cladelab(node=111, label="\n\nspecimens", barcolour = NA, offset = - 0.011, hjust = 0, fontsize = 4.5) +
  
  # support
  geom_point2(aes(subset=(UFboot > 95)), shape=20, size=2, colour="black") +
  
  # tip colors and shapes
  geom_tippoint(aes(color=Group,shape=Group), size=2) + 
  scale_color_manual(name = "Taxon/Locality", values = c("#66CCEE","#228833", "#AA3377","#CCBB44","#EE6677", "#BBBBBB"),
                     labels = cytb_labs) +
  scale_shape_manual(name = "Taxon/Locality", values=c(16,16,16,16,17,15), 
                     labels = cytb_labs) +
  theme(legend.position = c(0.25,0.7),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        plot.margin = margin(l = -0.5, r = 0.5, t = 0.5, b = 0.5, unit = "cm"))   


# mitogenome_mrbayes ------------------------------------------------------


BI_MT <-read.mrbayes("C:/Users/madel/OneDrive/Documents/Therya_truei_ms/analyses/Alignment_Truei_77samples.nex.con.tre") %>% 
  tree_subset(node = 98, levels_back = 0) 

# break up labels
BI_labels<- data.frame(label = BI_MT@phylo$tip.label) %>%
  separate_wider_delim(cols = "label", names = c("g", "s", "n"), delim = "_") %>%
  transmute(new = if_else((n =="SCL"),
                          paste(g,s,n,sep = "_"),         # no italics for these
                          glue("italic({g})~italic({s})~'({n})'")))

# rename tree's labels
BI_MT@phylo$tip.label = BI_labels$new

# plot tree
ggtree(BI_MT) +
  
  # add padding on the right
  xlim(0, 0.52) +
  
  # get node numbers
  #geom_text2(aes(subset=isTip, label=node, fontsize=.5), hjust=-.3) +
  
  #italic attempt
  #geom_richtext(data=td_filter(isTip), aes(label=name), label.color=NA) +
  
  # tip labels
  geom_tiplab(parse = TRUE) +
  
  # support
  geom_point2(aes(subset=!isTip), shape=20, size=3, colour="black") +
  
  # clade labels
  geom_cladelab(node=57, label="Clade A", barsize = 1, barcolour = "black", align = TRUE, offset = 0.16, offset.text = 0.005, extend = 0.2, fontsize = 4.5) +
  geom_cladelab(node=75, label="Clade B", barsize = 1, barcolour = "black", align = TRUE, offset = 0.16, offset.text = 0.005, extend = 0.2, fontsize = 4.5) +
  geom_cladelab(node=66, label="Clade C", barsize = 1, barcolour = "black", align = TRUE, offset = 0.16, offset.text = 0.005, extend = 0.2, fontsize = 4.5) +
  geom_cladelab(node=82, label="Clade H", barsize = 1, barcolour = "black", align = TRUE, offset = 0.16, offset.text = 0.005, extend = 0.2, fontsize = 4.5) +
  geom_cladelab(node=86, label="Clade P", barsize = 1, barcolour = "black", align = TRUE, offset = 0.16, offset.text = 0.005, extend = 0.2, fontsize = 4.5) +
  
  geom_cladelab(node=79, label="Habromys", fontface = 3, barsize = 1, barcolour = "black", align = TRUE, offset = 0.16, offset.text = 0.005, extend = 0.2, fontsize = 4.5) +
  geom_cladelab(node=55, label="Neotomodon", fontface = 3, barsize = 1, barcolour = "black", align = TRUE, offset = 0.16, offset.text = 0.005, extend = 0.2, fontsize = 4.5) + 
  geom_cladelab(node=54, label="Podomys", fontface = 3, barsize = 1, barcolour = "black", align = TRUE, offset = 0.16, offset.text = 0.005, extend = 0.2, fontsize = 4.5) +
  
  # tip labels
  geom_cladelab(node=30, label="Megadontomys", fontface = 3, barsize = 1, barcolour = "black", align = TRUE, offset = 0.16, offset.text = 0.005, extend = 0.2, fontsize = 4.5) +
  geom_cladelab(node=31, label="Clade E", barsize = 1, barcolour = "black", align = TRUE, offset = 0.16, offset.text = 0.005, extend = 0.2, fontsize = 4.5) +
  geom_cladelab(node=32, label="Osgoodomys", fontface = 3, barsize = 1, barcolour = "black", align = TRUE, offset = 0.16, offset.text = 0.005, extend = 0.2, fontsize = 4.5) +
  
  # bars for my taxa
  geom_cladelab(node=69, label="P. truei", fontface = 3, barsize = 1, barcolour = "#BB5566", offset = 0.113, offset.text = 0.002, extend = 0.4, fontsize = 4) +
  geom_cladelab(node=39, label="P. gambelii", fontface = 3, barsize = 1, barcolour = "#4477AA", offset = 0.16, fontsize = 4, offset.text = 0.002, extend = 0.4) +
  
  # boxes for my clades
  geom_highlight(node = 86, fill="#4477AA", alpha=.3, extend = 0.128) +
  geom_highlight(node = 66, fill="#EE6677", alpha=.3, extend = 0.11)


# mitogenome iqtree -------------------------------------------------------

ML_MT <- read.iqtree("C:/Users/madel/OneDrive/Documents/Therya_truei_ms/analyses/Mitogenome_truei_partitioned_iqtree.contree") %>% 
  tree_subset(node = 103, levels_back = 0)

# break up labels
ML_labels<- data.frame(label = ML_MT@phylo$tip.label) %>%
  separate_wider_delim(cols = "label", names = c("g", "s", "n"), delim = "_") %>%
  transmute(new = if_else((n =="SCL"),
                          paste(g,s,n,sep = "_"),         # no italics for these
                          glue("italic({g})~italic({s})~'({n})'")))

# rename tree's labels
ML_MT@phylo$tip.label = ML_labels$new

ggtree(ML_MT) +
  
  # add padding on the right
  xlim(NA, 0.463) +
  
  # get node numbers
  #geom_text2(aes(subset=!isTip, label=node, fontsize=.5), hjust=-.3) +
  
  # get tip numbers
  #geom_text2(aes(subset=isTip, label=node, fontsize=.5), hjust=-.3) 
  
  # tip labels
  geom_tiplab(parse=TRUE) +
  
  # support
  geom_point2(aes(subset=(UFboot > 95)), shape=20, size=3, colour="black") +

  # clade labels
  geom_cladelab(node=57, label="Clade A", barsize = 1, barcolour = "black", align = TRUE, offset = 0.16, offset.text = 0.005, extend = 0.2, fontsize = 4.5) +
  geom_cladelab(node=75, label="Clade B", barsize = 1, barcolour = "black", align = TRUE, offset = 0.16, offset.text = 0.005, extend = 0.2, fontsize = 4.5) +
  geom_cladelab(node=66, label="Clade C", barsize = 1, barcolour = "black", align = TRUE, offset = 0.16, offset.text = 0.005, extend = 0.2, fontsize = 4.5) +
  geom_cladelab(node=82, label="Clade H", barsize = 1, barcolour = "black", align = TRUE, offset = 0.16, offset.text = 0.005, extend = 0.2, fontsize = 4.5) +
  geom_cladelab(node=86, label="Clade P", barsize = 1, barcolour = "black", align = TRUE, offset = 0.16, offset.text = 0.005, extend = 0.2, fontsize = 4.5) +
  
  geom_cladelab(node=79, label="Habromys", fontface = 3, barsize = 1, barcolour = "black", align = TRUE, offset = 0.16, offset.text = 0.005, extend = 0.2, fontsize = 4.5) +
  geom_cladelab(node=55, label="Neotomodon", fontface = 3, barsize = 1, barcolour = "black", align = TRUE, offset = 0.16, offset.text = 0.005, extend = 0.2, fontsize = 4.5) + 
  geom_cladelab(node=54, label="Podomys", fontface = 3, barsize = 1, barcolour = "black", align = TRUE, offset = 0.16, offset.text = 0.005, extend = 0.2, fontsize = 4.5) +
  
  # tip labels
  geom_cladelab(node=30, label="Megadontomys", fontface = 3, barsize = 1, barcolour = "black", align = TRUE, offset = 0.16, offset.text = 0.005, extend = 0.2, fontsize = 4.5) +
  geom_cladelab(node=31, label="Clade E", barsize = 1, barcolour = "black", align = TRUE, offset = 0.16, offset.text = 0.005, extend = 0.2, fontsize = 4.5) +
  geom_cladelab(node=32, label="Osgoodomys", fontface = 3, barsize = 1, barcolour = "black", align = TRUE, offset = 0.16, offset.text = 0.005, extend = 0.2, fontsize = 4.5) +
  
  # bars for my taxa
  geom_cladelab(node=69, label="P. truei", fontface = 3, barsize = 1, barcolour = "#BB5566", offset = 0.11, offset.text = 0.002, extend = 0.4, fontsize = 4) +
  geom_cladelab(node=42, label="P. gambelii", fontface = 3, barsize = 1, barcolour = "#4477AA", offset = 0.149, fontsize = 4, offset.text = 0.002, extend = 0.4) +
    
  
  # boxes for my clades
  geom_highlight(node = 86, fill="#4477AA", alpha=.3, extend = 0.128) +
  geom_highlight(node = 66, fill="#EE6677", alpha=.3, extend = 0.11)



# beast tree --------------------------------------------------------------

beast <-read.beast("C:/Users/madel/OneDrive/Documents/Therya_truei_ms/analyses/BEAST_MitoTruei_77sPosiStrict_3PointsT1R1_2.tree") %>% 
  tree_subset(node = 85, levels_back = 0)

# break up labels
beast_labels<- data.frame(label = beast@phylo$tip.label) %>%
  separate_wider_delim(cols = "label", names = c("g", "s", "n"), delim = "_") %>%
  transmute(new = if_else((n =="SCL"),
                          paste(g,s,n,sep = "_"),         # no italics for these
                          glue("italic({g})~italic({s})~'({n})'")))

# rename tree's labels
beast@phylo$tip.label = beast_labels$new

ggtree(beast)  +
  
  # add padding on the right
  xlim(0, 7) +
  
  # get node numbers
  #geom_text2(aes(subset=isTip, label=node, fontsize=.5), hjust=-.3) +
  
  # tip labels
  geom_tiplab(size=3.5, parse = TRUE) +
  
  # median ages
  geom_text2(aes(subset=!isTip, label=format(round(CAheight_median, digits = 3), nsmall = 3)), vjust = -0.5, hjust = 1, size = 3.5) +
  
  # 95% HPD bars
  geom_range(range = 'CAheight_0.95_HPD', color='darkgreen', size=2, alpha=.3, center='CAheight_median') +
  
  # bars for my taxa
  geom_cladelab(node=65, label="P. truei", fontface = 3, size = 3.5, barsize = 1, barcolour = "#BB5566", offset = 1.113, offset.text = 0.002, extend = 0.4, fontsize = 4) +
  geom_cladelab(node=5, label="P. gambelii", fontface = 3, size = 3.5, barsize = 1, barcolour = "#4477AA", offset = 1.16, fontsize = 4, offset.text = 0.002, extend = 0.4) 
  

