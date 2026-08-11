
# measurements_boxplots.R -------------------------------------------------
#
# MAB 10/31/25
# Code for morphological measurement analyses for the following manuscript:
# Becker MA, Bell KC, Lopez-Ortiz S, Maldonado JE, Edwards CW, Castañeda-Rico S. (2026) 
# Genomic detection of a novel Peromyscus truei invasion on San Clemente Island, California. 
# Therya 17:277-294. doi: 10.12933/therya.2026.6265

# load packages -----------------------------------------------------------

library(tidyverse)
library(ggpubr)

# make von Bloeker only dataset -------------------------------------------


vB <- read.csv("CI_meas_tibble_imm.csv") %>% 
  select(catalogNo,
         island = Island,
         genetic_ID,
         sex = Sex,	
         total_length,
         tail_length,	
         foot_length = hindfoot,	
         ear_length = ear_from_notch,
         immature,
         pregnant,
         preparator = MamPreparator,
         year = Year) %>% 
  mutate(island = case_when(island %in% c("AI-E", "AI-M", "AI-W") ~ "AI",TRUE ~ island)) %>% 
  mutate(tail_percent = tail_length/(total_length-tail_length)) %>% 
  filter(preparator == "von Bloeker, J C") %>% 
  mutate(island = recode(island,
                         "AI" = "Anacapa",
                         "CAT" = "Santa Catalina",
                         "SBI" = "Santa Barbara",
                         "SCL" = "San Clemente",
                         "SCZ" = "Santa Cruz",
                         "SMI" = "San Miguel",
                         "SNI" = "San Nicolas",
                         "SRI" = "Santa Rosa")) %>% 
  filter(!catalogNo %in% c(6213,6623))


# get averages ----------------------------------------------------------------

vB %>% 
  filter(island == "San Clemente") %>% 
  summarize(total_length_mean = mean(total_length, na.rm = TRUE))


# establish ggplot theming ------------------------------------------------


theming <-
  theme_classic() +
  theme(legend.text = element_text(size = 12, face = "italic"),
        legend.title = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,color = "black", size = 10),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) 


# box plots ---------------------------------------------------------------


# box/dot of ear length 
ear <- vB %>% 
  ggplot(aes(x=island,y=ear_length)) +
  geom_boxplot() +
  geom_dotplot(aes(fill = genetic_ID),
               binaxis='y', 
               stackdir='center', 
               dotsize = 0.8) +
  labs(fill = "Genetic Species ID", y = "Ear length (mm)") +
  scale_y_continuous(breaks=seq(14,24,2)) +
  scale_fill_manual(labels = c("P. gambelii", "P. truei"), values = c("#4477AA","#EE6677")) +
  theming

# box/dot of total length 
total <- vB %>% 
  ggplot(aes(x=island,y=total_length)) +
  geom_boxplot() +
  geom_dotplot(aes(fill = genetic_ID),
               binaxis='y', 
               stackdir='center', 
               dotsize = 0.8) +
  labs(fill = "Genetic Species ID", y = "Total body length (mm)") +
  #scale_y_continuous(breaks=seq(14,24,2)) +
  scale_fill_manual(labels = c("P. gambelii", "P. truei"), values = c("#4477AA","#EE6677")) +
  theming

# box/dot of tail percentage
percent <-vB %>% 
  filter(total_length > 150) %>% 
  ggplot(aes(x=island,y=tail_percent)) +
  geom_boxplot() +
  geom_dotplot(aes(fill = genetic_ID),
               binaxis='y', 
               stackdir='center', 
               dotsize = 0.8) +
  geom_hline(yintercept = 0.9, linetype = "dashed") +
  labs(fill = "Genetic Species ID", y = "Tail length/Body length ratio") +
  
  #scale_y_continuous(breaks=seq(14,24,2)) +
  scale_fill_manual(labels = c("P. gambelii", "P. truei"), values = c("#4477AA","#EE6677")) +
  theming

# violin of total length (shows SCL/SRI pattern)
violin <- vB %>% 
  ggplot(aes(x=island,y=total_length)) +
  geom_violin() +
  geom_dotplot(aes(fill = genetic_ID),
               binaxis='y', 
               stackdir='center', 
               dotsize = 0.8) +
  labs(fill = "Genetic Species ID", y = "Total body length (mm)") +
  scale_fill_manual(labels = c("P. gambelii", "P. truei"), values = c("#4477AA","#EE6677")) +
  theming

# Figure 5 plot grid
ggarrange(ear, percent, total, violin, 
          labels = "AUTO", 
          common.legend = TRUE, 
          legend = "bottom")


# supplementary figure box plots ------------------------------------------


# what about pregnancy?
# box/dot of total length 
p <-
vB %>% 
  ggplot(aes(x=island,y=total_length)) +
  geom_boxplot() +
  geom_dotplot(aes(fill = pregnant),
               binaxis='y', 
               stackdir='center', 
               dotsize = .6) +
  labs(fill = "Pregnant?", y = "Total body length (mm)") +
  scale_fill_brewer(palette="PiYG") +
  theme_classic() +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,color = "black", size = 10),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) 

# what about immaturity?
# box/dot of total length 
i <-
vB %>% 
  ggplot(aes(x=island,y=total_length)) +
  geom_boxplot() +
  geom_dotplot(aes(fill = immature),
               binaxis='y', 
               stackdir='center', 
               dotsize = .6) +
  labs(fill = "Immature?", y = "Total body length (mm)") +
  scale_fill_brewer(palette="YlGn") +
  theme_classic() +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,color = "black", size = 10),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) 


# supplementary box plots
ggarrange(p, i, 
          labels = "AUTO")


# what about sex in general?
# box/dot of total length 
vB %>% 
  filter(sex %in% c("F","M")) %>% 
  ggplot(aes(x=island,y=total_length, fill = sex)) +
  geom_boxplot() +
  labs(fill = "Sex", y = "Total body length (mm)") +
  theming


# box/dot of tail length 
vB %>% 
  ggplot(aes(x=island,y=tail_length)) +
  geom_boxplot() +
  geom_dotplot(aes(fill = genetic_ID),
               binaxis='y', 
               stackdir='center', 
               dotsize = .5) +
  scale_fill_manual(labels = c("P. gambelii", "P. truei"), values = c("#4477AA","#EE6677")) +
  theming


