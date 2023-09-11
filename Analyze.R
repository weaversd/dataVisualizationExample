#libraries
#if you're running this for the first time, you'll have to install
#these packages. You can google how to install packages in R, or you can
#click on the package tab to the right and mannually install each by name.
#once installed, you won't have to install them again
library(stringr) #working with strings (text)
library(ggplot2) #creating plots
library(dplyr) #reshaping data
library(tidyr) #reshaping data
library(seqinr) #working with biological sequences
library(ggpubr) #statistics


conflicts_prefer(dplyr::select) #if there is a conflict, this tells
#R to use 'select' from 'dplyr'


setwd("~/04_Champion_Lab/VisualizationWorkshop") #change to wherever you downloaded
#these files to

DDAproteins <- read.csv("DDAresults/db.proteins.csv") #import DDA results from CSV

#pull out columns that contain area values
areas <- select(DDAproteins, matches("Area.DDA_300Smeg_D_[1-3]"))

#loop through each row and calculate area, sd, and RSD. Also identify missing values
for (i in 1:nrow(DDAproteins)) {
  DDAproteins$AccessionID[i] <- 
    str_split(DDAproteins$Accession[i], "\\|")[[1]][1] #pull out only accession IDs
  
  areaVector <- as.numeric(areas[i,]) #pull out specific row of areas
  
  DDAproteins$missingAreaCount[i] <- sum(is.na(areaVector)) #count missing values
  
  if (DDAproteins$missingAreaCount[i] < 3) {
    #calculate values and assign to new columns
    DDAproteins$AreaAvg[i] <- mean(areaVector, na.rm = T)
    DDAproteins$AreaSD[i] <- sd(areaVector, na.rm = T)
    DDAproteins$AreaRSD[i] <- (DDAproteins$AreaSD[i] / DDAproteins$AreaAvg[i]) * 100
  } else {
    #if there are too many missing values, leave as NA
    DDAproteins$AreaAvg[i] <- NA
    DDAproteins$AreaSD[i] <- NA
    DDAproteins$AreaRSD[i] <- NA
  }
}

#pull out only the columns we're interested in
DDAproteins <- select(DDAproteins, c(1, (ncol(DDAproteins) - 4):(ncol(DDAproteins))))

#creat the initial plot
plot <- ggplot(DDAproteins) + #which data to use
  geom_histogram(aes(x = AreaRSD), binwidth = 1, #which variables to use
                 fill = "green4", color = "black", alpha = 0.5) + #color and size
  theme_bw(base_size = 25) + #thematic elements
  theme(panel.grid = element_blank()) +
  labs(y = "Protein Count", #labels
       x = "%RSD (Area, n = 3)",
       title = "DDA proteins %RSDs")
show(plot) #display the plot


DDAproteinMeanCV <- mean(DDAproteins$AreaRSD, na.rm = T) #calculate mean RSD
DDAproteinMeanCVlabel <- paste0("Mean: ",
                                round(DDAproteinMeanCV, 2), "%") #create label

plot + geom_vline(xintercept = DDAproteinMeanCV, linewidth = 2,
             linetype = 2, color = "grey50") + #add line to graph at mean
  geom_text(x = DDAproteinMeanCV + 4, y = 330, label = DDAproteinMeanCVlabel, hjust = 0,
            color = "grey50", size = 5, angle = 270) #add label to graph


#import DIA proteins
DIAproteins <- read.csv("DIAresults/report.pg_matrix.tsv", sep = "\t")

#pull out area columns
areas <- select(DIAproteins, matches("timsTOF_files"))

#similar loop as above for calculating attributes
for (i in 1:nrow(DIAproteins)) {
  
  areaVector <- as.numeric(areas[i,])
  
  DIAproteins$missingAreaCount[i] <- sum(is.na(areaVector))
  
  if (DIAproteins$missingAreaCount[i] < 3) {
    DIAproteins$AreaAvg[i] <- mean(areaVector, na.rm = T)
    DIAproteins$AreaSD[i] <- sd(areaVector, na.rm = T)
    DIAproteins$AreaRSD[i] <- (DIAproteins$AreaSD[i] / DIAproteins$AreaAvg[i]) * 100
  }
}

#pull out columns we care about
DIAproteins <- select(DIAproteins, c(1, (ncol(DIAproteins) - 3):(ncol(DIAproteins))))

#label and mean
DIAproteinMeanCV <- mean(DIAproteins$AreaRSD, na.rm = T)
DIAproteinMeanCVlabel <- paste0("Mean: ",
                                round(DIAproteinMeanCV, 2), "%")

#create plot, same as above
plot <- ggplot(DIAproteins) +
  geom_histogram(aes(x = AreaRSD), binwidth = 1,
                 fill = "green4", color = "black", alpha = 0.5) +
  theme_bw(base_size = 25) +
  theme(panel.grid = element_blank()) +
  labs(y = "Protein Count",
       x = "%RSD (Area, n = 3)",
       title = "DIA proteins %RSDs")
show(plot)

#add text and line for mean
plot + geom_vline(xintercept = DIAproteinMeanCV, linewidth = 2,
             linetype = 2, color = "grey50") +
  geom_text(x = DIAproteinMeanCV + 1.6, y = 600, label = DIAproteinMeanCVlabel, hjust = 0,
            color = "grey50", size = 5, angle = 270)



#add column clarifying which dataset each row comes from
DDAproteins$method <- "DDA"
DDAproteins$Protein.Group <- as.character(DDAproteins$Protein.Group) #turn to string for binding
DIAproteins$method <- "DIA"

#combine dataframes
combinedProteins <- bind_rows(DDAproteins, DIAproteins)

#create combined plot
plot <- ggplot(combinedProteins) + #which data
  geom_histogram(aes(x = AreaRSD, fill = method), binwidth = 1 #variables
                 , color = "black", alpha = 0.5) + #colors sizes etc
  theme_bw(base_size = 25) +
  theme(panel.grid = element_blank()) + #theme elements
  labs(y = "Protein Count", #labels
       x = "%RSD (Area, n = 3)",
       title = "proteins %RSDs") +
  facet_grid(method~.) #create a grid with plots based on which method
show(plot) #show plots


#create a unique row for each accession ID in DIA dataset
DIAproteins$AccessionID <- str_split(DIAproteins$Protein.Group, ";")
DIAproteins <- DIAproteins %>% unnest(AccessionID)

#recombine
combinedProteins <- bind_rows(DDAproteins, DIAproteins)

#plot recombined data, same as above
plot <- ggplot(combinedProteins) +
  geom_histogram(aes(x = AreaRSD, fill = method), binwidth = 1
                 , color = "black", alpha = 0.5) +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank()) +
  labs(y = "Protein Count",
       x = "%RSD (Area, n = 3)",
       title = "proteins %RSDs") +
  facet_grid(method~.)
show(plot)

#pull in protein masses from uniprot. Look up how to install "UniProt.ws" from
#BioConductor if you have not installed it before
BiocManager::install("UniProt.ws")
library(UniProt.ws)


#create column for mass
combinedProteins$mass <- NA
#loop through rows in chunks of 100 (largest import size that is reasonable in time)
for (i in seq(1,nrow(combinedProteins), 100)) {
  if (i + 99 < nrow(combinedProteins)) {
    #pull info from uniprot
    masses <- queryUniProt(query = combinedProteins$AccessionID[i:(i+99)], fields = c("mass", "accession"))
    print(nrow(masses))
    combinedProteins$mass[i:(i+99)] <- masses$Mass #save in dataframe
  } else {
    #edge case for last set that isn't a multiple of 100
    masses <- queryUniProt(query = combinedProteins$AccessionID[i:nrow(combinedProteins)], fields = "mass")
    combinedProteins$mass[i:nrow(combinedProteins)] <- masses$Mass
  }
}

#convert to kDa
combinedProteins$mass <- combinedProteins$mass / 1000
#create mass bins of 10 kDa width
combinedProteins$MassBin <- cut(combinedProteins$mass, seq(0, max(combinedProteins$mass, na.rm = T), 10), include.lowest = T)

#plot histogram by mass
ggplot(combinedProteins[combinedProteins$mass < 50,]) + #which data (only masses below 50)
  geom_histogram(aes(x = AreaRSD,
                     fill = method), binwidth = 1, color = "black", alpha = 0.5) +
  theme_bw(base_size = 20,) +
  theme(panel.grid = element_blank()) +
  labs(y = "Protein Count",
       x = "%RSD (Area, n = 3)",
       title = "Proteins",
       color = element_blank()) +
  coord_cartesian(xlim = c(0, 100)) +
facet_grid(MassBin~method) #create a grid of plots for each method/Mass bin combo

#create density plot
ggplot(combinedProteins[combinedProteins$mass < 50,]) +
  geom_density(aes(x = AreaRSD,
                     fill = method),
               alpha = 0.5,
               linewidth = 1.3) +
  theme_bw(base_size = 20,) +
  theme(panel.grid = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) +
  labs(y = "Density",
       x = "%RSD (Area, n = 3)",
       title = "Proteins",
       color = element_blank()) +
  coord_cartesian(xlim = c(0, 100)) +
  facet_grid(MassBin~.)








#create violin, box, and scatter plot (oh my)
ggplot(combinedProteins, aes(x = method, y = AreaRSD)) + #which data
  geom_jitter(aes(color = method), size = 0.5, height = 0, #add the dots (jitter adds random spread so they don't all overlap)
              alpha = 0.5) +
  geom_violin(alpha = 0, size = 0.8) + #add violin plot
  theme_bw(base_size = 25) +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  labs(x =element_blank(),
       y = "%RSD quantification (n =3)") +
  stat_compare_means(comparisons = list(c("DDA", "DIA")), #calculate significance using default comparison mean
                     label = "p.signif") +
  geom_boxplot(width = 0.03, outlier.shape = NA, #add box plot on top
               size = 0.8, alpha = 0.6)


#same as above but split by different mass bins
ggplot(combinedProteins[combinedProteins$mass <= 60,], aes(x = method, y = AreaRSD)) +
  geom_jitter(aes(color = method), size = 0.5, height = 0,
              alpha = 0.5) +
  geom_violin(alpha = 0, size = 0.8) +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  labs(x =element_blank(),
       y = "%RSD quantification (n =3)") +
  stat_compare_means(comparisons = list(c("DDA", "DIA")),
                     label = "p.signif",
                     size = 5,bracket.size = 0.8) +
  geom_boxplot(width = 0.03, outlier.shape = NA,
               size = 0.8, alpha = 0.6) +
  facet_wrap(~MassBin) + #creat the grid by mass bins
  coord_cartesian(ylim = c(0,200)) +
  theme( #additional thematic customizations
    axis.ticks = element_blank(),
    axis.line = element_line(colour = "grey50"),
    panel.grid = element_line(color = "#b4aea9"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(linetype = "dashed"),
    panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
    plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4")
  )

#save as png
ggsave("percRSD DIA vs DDA by mass.png")




#create scatter plot (RSD by mass)
ggplot(combinedProteins,
       aes(x = AreaRSD, y = mass, base = 10)) +
  geom_point(aes(color = method), size = 1, alpha = 0.5) +
  theme_bw(base_size = 25) +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  labs(y = "Protein Mass (kDa)",
       x = "%RSD quantification (n =3)") +
  coord_cartesian(xlim = c(0,75), y = c(0,100)) + #set plot boundries
  facet_grid(method~.)

# same as above but with density layer added
ggplot(combinedProteins,
       aes(x = AreaRSD, y = mass)) +
  geom_point(size = 0.5)+
  geom_density_2d_filled(contour_var = "ndensity", alpha = 0.7) +
  theme_bw(base_size = 25) +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  labs(y = "Protein Mass (kDa)",
       x = "%RSD quantification (n = 3)") +
  coord_cartesian(xlim = c(0,50), y = c(0,100)) +
  facet_grid(method~.)


