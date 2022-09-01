# packages ----------------------------------------------------------------

library(ggplot2)
library(ggpubr)

# cwm vs periods ----------------------------------------------------------

##defining fonts, faces and sizes for all the graphics
mynamestheme <- theme(strip.text = element_text(family = "serif", size = (10)),
                      axis.title = element_text(family = "serif", size = (11),
                                                face = "bold"),
                      axis.text = element_text(family = "serif", size = (11)))

###colors for periods in ridgelines
Per <- c("Dry" = "#E59866", "Dry-Rain" = "#F6DDCC", "Rain" = "#5499C7",
        "Rain-Dry" = "#D6EAF8")


##png object
png("CWMvsPeriods.png", width = 16, height = 12, units = "cm",
    pointsize = 8,
    res = 300)

###CWM phenols vs periods
CWMphenols.vs.periods <- ggplot(data = data.litter, 
                                aes(x = Period, y = CWPhenols))+
         geom_boxplot(width=0.5, color="darkgrey",
                                         alpha= 0.8)+
     geom_jitter(size = 2, position = position_jitter(width = 0.1))+
    geom_text(x = "Dry", y = 42, label = "(a)", family = "serif",  
            size = 5)+##subsection
     geom_text(x = "Dry", y = 36, label = "a", family = "serif",  
            size = 4)+
     geom_text(x = "Dry-Rain", y = 36, label = "a", family = "serif", 
            size = 4)+
   geom_text(x = "Rain", y = 36, label = "ab", family = "serif", 
            size = 4)+
   geom_text(x = "Rain-Dry", y = 36, label = "b", family = "serif",  
            size = 4)+
   theme_classic()+
   labs(x = "Periods", 
       y = "CWM Phenols")+ mynamestheme

### CWM phenols ridgeline

CWM.phe.ridge <- ggplot(data.litter, aes(x = CWPhenols, y = Period))+
  geom_density_ridges(aes(fill = Period))+
  geom_text(x = 22, y = "Rain-Dry", label = "(b)", family = "serif",  
            size = 5)+##subsection
  theme_ridges()+
  theme_classic()+ scale_fill_manual(values = Per, guide = "none")+
  scale_color_manual("black")+ mynamestheme

######CWM Leaf width vs periods
CWMwidth.vs.periods <- ggplot(data = data.litter, 
                             aes(x = Period, y = CWwidth))+
  geom_boxplot(width=0.5, color="darkgrey", alpha= 0.8)+
  geom_text(x = "Dry", y = 3.45, label = "(c)", family = "serif",  
            size = 5)+##subsection
  geom_jitter(size = 2, position = position_jitter(width = 0.1))+
  geom_text(x = "Dry", y = 3.1, label = "a", family = "serif", hjust = "left",  
            size = 4)+
  geom_text(x = "Dry-Rain", y = 3.15, label = "a", family = "serif",
            hjust = "left", size = 4)+
  geom_text(x = "Rain", y = 2.8, label = "ab", family = "serif", hjust = "left", 
            size = 4)+
  geom_text(x = "Rain-Dry", y = 2.8, label = "b", family = "serif", 
            hjust = "left",  size = 4)+
  theme_classic()+
  labs(x = "Periods", 
       y = "CWM Leaf Width")+ mynamestheme


###CWM Leaf width ridgeline
CWMwidth.ridge <- ggplot(data.litter, aes(x = CWwidth, y = Period))+
  geom_density_ridges(aes(fill = Period))+
  geom_text(x = 1.8, y = "Rain-Dry", label = "(d)", family = "serif",  
            size = 5)+##subsection
  theme_ridges()+
  theme_classic()+
  labs(y = "Periods", 
       x = "CWM Leaf Width")+
  scale_fill_manual(values = Per, guide = "none")+
  scale_color_manual("black")+ mynamestheme

ggarrange(CWMphenols.vs.periods, CWM.phe.ridge, CWMwidth.vs.periods,
          CWMwidth.ridge, ncol = 2, nrow = 2)

dev.off()


# biodiversity components vs inputs ---------------------------------------

##excluding lateral input because it was not evaluated in post hoc
data.litter.vt <- data.litter %>% filter(Input != "lateral")

##png object
png("components.vs.Inputs.png", width = 18, height = 14, units = "cm",
    pointsize = 8,
    res = 300)

###cwm nitrogen vs input
CWMnit.vs.Input <- ggplot(data = data.litter.vt, 
                          aes(x = Input, y = CWNitro))+
  geom_boxplot(width=0.4, color="darkgrey", alpha= 0.8)+
  geom_jitter(size = 2, position = position_jitter(width = 0.1))+
  #geom_text(x = "terrestrial", y = 0.89, label = "a", family = "serif", 
            #hjust = "left",size = 4)+
  geom_text(x = "vertical", y = 1.07, label = "b", family = "serif", 
            hjust = "left", size = 4)+
  theme_classic()+
  labs(x = "Input", 
       y = "CWM Nitrogen")+
  mynamestheme

##CWM leaf width vs input
CWwidth.vs.Input <- ggplot(data = data.litter.vt, 
                           aes(x = Input, y = CWwidth))+
  geom_boxplot(width=0.4, color="darkgrey", alpha= 0.8)+
  geom_jitter(size = 2, position = position_jitter(width = 0.1))+
  geom_text(x = "terrestrial", y = 3.15, label = "a", family = "serif", 
            hjust = "left",size = 4)+
  geom_text(x = "vertical", y = 2.8, label = "b", family = "serif", 
            hjust = "left", size = 4)+
  theme_classic()+
  labs(x = "Input", 
       y = "CWM Leaf Width")+
   mynamestheme

###functional diversity vs input
FDiv.vs.Input <- ggplot(data = data.litter.vt, 
                        aes(x = Input, y = FDiv))+
  geom_boxplot(width=0.4, color="darkgrey", alpha= 0.8)+
  geom_jitter(size = 2, position = position_jitter(width = 0.1))+
  #geom_text(x = "terrestrial", y = 0.92, label = "a", family = "serif", 
            #hjust = "left",
            #size = 4)+
  geom_text(x = "vertical", y = 0.72, label = "b", family = "serif", 
            hjust = "left", size = 4)+
  theme_classic()+
  labs(x = "Input", 
       y = "Functional Divergence (FDiv)")+
  mynamestheme

ggarrange(CWMnit.vs.Input,CWwidth.vs.Input, FDiv.vs.Input,
          nrow = 2, ncol =2)

dev.off()
