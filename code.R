
# Litter diversity components ---------------------------------------------

library(tidyverse)
library(vegan)
library(FD)
library(mFD)
library(ggplot2)
library(lme4)
library(readxl)
library(lmerTest)
library(broom.mixed)
library(broom)
library(emmeans)
library(ggridges)
library(viridis)
library(xlsx)
library(gridExtra)

# litter input matrix -----------------------------------------------------
litter <- read_excel("Riparian_community.xlsx", sheet = "OM_dynamic")

##separating the variables for the dbFD process
variables <- litter %>% select(Community, Month, Input, Period, Point)

##cleaning the litter matrix
litter.matrix <- litter %>% select(-c(Community, Month, Input, Period, Point))

##removing rows with 0 abundance in every species
litter.matrix.1 <- litter.matrix[which(rowSums(litter.matrix)>0),]

##removing for the variables also
variables <- variables[which(rowSums(litter.matrix)>0),]

# litter traits -----------------------------------------------------------

traits <- read_excel("Riparian_community.xlsx", sheet = "litter.traits")

traits <- as.data.frame(traits)

rownames(traits) <- traits$Species

traits <- traits[,-c(1,7)]


# Calculating functional indices ------------------------------------------

indices <- dbFD(traits, litter.matrix.1)

indices.matrix <- do.call(cbind.data.frame, indices[1:8])

indices.matrix <- cbind(indices.matrix, indices$CWM)


# joining variables and indices -------------------------------------------


data.litter <- cbind(variables, indices.matrix)%>% rename(richness = nbsp,
                                                         CWPhenols = Phenols,
                                                         CWNitro = Nitrogen,
                                                         CWCarbon = Carbon,
                                                         CWlength = Leaf_Length,
                                                         CWwidth = Leaf_Width)



####Pielou index

#shannon index
H <- diversity(litter.matrix.1, index = "shannon")

##Pielou
J <- H/log(specnumber(litter.matrix.1))

data.litter <- cbind(data.litter, J)

# modeling ----------------------------------------------------------------

##richness
mRichness <- lmer(richness ~ Period * Input + (1|Point), data = data.litter)

summary(mRichness)

arich <- anova(mRichness)

##FRic
mFric <- lmer(FRic ~ Period * Input + (1|Point), data = data.litter)

summary(mFric)

aFric <- anova(mFric)
plot(mFric)
plot(mFric, FRic ~ fitted(.))
qqnorm(residuals(mFric))


#FEve
mFeve <- lmer(FEve ~ Period * Input  + (1|Point), data = data.litter)

summary(mFeve)

afeve <- anova(mFeve)

#FDiv
mFdiv <- lmer(FDiv ~ Period * Input + (1|Point), data = data.litter)

summary(mFdiv)

afdiv <- anova(mFdiv) ##input differences

c.mFdiv <- emmeans(mFdiv, "Input")

contrast(c.mFdiv, "tukey") ##significative differences between TI and VI

###likelihood test between models
mFdiv.null <- lmer(FDiv ~ 1 + (1|Point), data = data.litter)

anova(mFdiv, mFdiv.null)


#CWPhenols

mPhe <- lmer(CWPhenols ~ Period *Input + (1|Point), data = data.litter)

summary(mPhe)

a.mphe <- anova(mPhe) ##significative in Month and Input

#determining group differences among period

c.mPhe <- emmeans(mPhe, "Period")

contrast(c.mPhe, "tukey") ##Significative differences between Rain-Dry and Dry
                          ## and Dry-Rain. Rain differences were not estimated

##contrast among Input

c.in.Phe <- emmeans(mPhe, "Input")
contrast(c.in.Phe, "tukey")

##########likelihood test
mPhe.null <- lmer(CWPhenols ~ 1 + (1|Point), data = data.litter)
anova(mPhe, mPhe.null) ##significant


#CWNitro
mNit <- lmer(CWNitro ~ Period * Input + (1|Point), data = data.litter)

summary(mNit)

a.mnit <- anova(mNit)#slighty significative in Period and significative in Input

##contrast
c.Mnit <- emmeans(mNit, "Input")
contrast(c.Mnit, "tukey") ##significative difference
                           ##between vertical and terrestrial

###likelihood test
mNit.null <- lmer(CWNitro ~ 1 + (1|Point), data = data.litter)
anova(mNit.null, mNit) ##significant

#CWCarbon
mCarb <- lmer(CWCarbon ~ Period * Input + (1|Point), data = data.litter)

summary(mCarb)

amCarb<- anova(mCarb)##input is significant

#contrast
c.mCarb <- emmeans(mCarb, "Input")
contrast(c.mCarb, "tukey")##difference between terrestrial and vertical

##likelihood test
mCarb.null <- lmer(CWCarbon ~ 1 + (1|Point), data = data.litter)
anova(mCarb, mCarb.null) ##not significant

#CWlength
mLen <- lmer(CWlength ~ Period * Input + (1|Point), data = data.litter)

summary(mLen)

a.mLen <- anova(mLen)

#CWwidth
mWidth <- lmer(CWwidth ~ Period * Input + (1|Point), data = data.litter)

a.mwidth <- anova(mWidth)##significant for both month and input

##contrast among periods
c.mWidth <- emmeans(mWidth, "Period")
contrast(c.mWidth, "tukey") #Rain-Dry significative difference between Dry and
                            #Dry-rain

##contrast between inputs
c.input.width <- emmeans(mWidth, "Input")
contrast(c.input.width, "tukey") ##terrestrial and vertical

##likelihood  test
mWidth.null <- lmer(CWwidth ~ 1 + (1|Point), data = data.litter)

anova(mWidth, mWidth.null) ##significant

#Pielou
mPie <- lmer(J ~ Period * Input + (1|Point), data = data.litter)

a.mPie <- anova(mPie)

# Graphics ----------------------------------------------------------------

##excluding lateral input because it was not evaluated in post hoc
data.litter.vt <- data.litter %>% filter(Input != "lateral")


##defining fonts, faces and sizes for all the graphics
mynamestheme <- theme(strip.text = element_text(family = "serif", size = (12)),
                      axis.title = element_text(family = "serif", size = (13),
                                                face = "bold"),
                      axis.text = element_text(family = "serif", size = (12)))

###FDiv vs Input

##png object
png("FDiv.vs.Input.png", width = 16, height = 12, units = "cm", pointsize = 8,
    res = 300)

##ggplot basic information
FDiv.vs.Input <- ggplot(data = data.litter.vt, 
                        aes(x = Input, y = FDiv))
##ggplot complement

p = FDiv.vs.Input+ geom_boxplot(width=0.4, color="darkgrey", alpha= 0.8)+
  geom_jitter(size = 2, position = position_jitter(width = 0.1))+
  geom_text(x = "terrestrial", y = 0.92, label = "a", family = "serif", 
            hjust = "left",
            size = 5)+
  geom_text(x = "vertical", y = 0.72, label = "b", family = "serif", 
            hjust = "left", size = 5)+
  theme_classic()+
  labs(x = "Input", 
       y = "Litter Functional Divergence (FDiv)")


p + mynamestheme

dev.off()


##CWM phenois vs Input 

##png object
png("CWMphenols.vs.Input.png", width = 16, height = 12, units = "cm",
    pointsize = 8,
    res = 300)

##ggplot basic information
CWMphenols.vs.Input <- ggplot(data = data.litter.vt, 
                        aes(x = Input, y = CWPhenols))
##ggplot complement

p = CWMphenols.vs.Input + geom_boxplot(width=0.4, color="darkgrey", alpha= 0.8)+
  geom_jitter(size = 2, position = position_jitter(width = 0.1))+
  geom_text(x = "terrestrial", y = 32.8, label = "a", family = "serif", 
            hjust = "left",
            size = 5)+
  geom_text(x = "vertical", y = 33.3, label = "b", family = "serif", 
           hjust = "right", size = 5)+
  theme_classic()+
  labs(x = "Input", 
       y = "CWM Phenols")


p + mynamestheme

dev.off()


##########CWM phenois vs Periods 

##png object
png("CWMphenols.vs.Periods.png", width = 16, height = 12, units = "cm",
    pointsize = 8,
    res = 300)

##ggplot basic information
CWMphenols.vs.periods <- ggplot(data = data.litter, 
                              aes(x = Period, y = CWPhenols))
##ggplot complement

p = CWMphenols.vs.periods + geom_boxplot(width=0.5, color="darkgrey",
                                         alpha= 0.8)+
  geom_jitter(size = 2, position = position_jitter(width = 0.1))+
  geom_text(x = "Dry", y = 36, label = "a", family = "serif",  
            size = 5)+
  geom_text(x = "Dry-Rain", y = 36, label = "a", family = "serif", 
            size = 5)+
  geom_text(x = "Rain", y = 36, label = "ab", family = "serif", 
            size = 5)+
  geom_text(x = "Rain-Dry", y = 36, label = "b", family = "serif",  
            size = 5)+
  theme_classic()+
  labs(x = "Periods", 
       y = "CWM Phenols")


p + mynamestheme

dev.off()


##########Periods vs CWM Phenois ridgeline

##png object
png("CWMphenols.vs.Periods.ridgeline.png", width = 16, height = 12,
    units = "cm", pointsize = 8,
    res = 300)


Per <- c("Dry" = "#E59866", "Dry-Rain" = "#F6DDCC", "Rain" = "#5499C7",
         "Rain-Dry" = "#D6EAF8")


o <- ggplot(data.litter, aes(x = CWPhenols, y = Period))+
  geom_density_ridges(aes(fill = Period))+
  theme_ridges()+
  theme_classic()

o + scale_fill_manual(values = Per, guide = "none")+
  scale_color_manual("black")+ mynamestheme

dev.off()

########CWM nitrogen vs Input boxplot

##png object
png("CWMnit.vs.Input.png", width = 16, height = 12, units = "cm", pointsize = 8,
    res = 300)

##ggplot basic information
CWMnit.vs.Input <- ggplot(data = data.litter.vt, 
                              aes(x = Input, y = CWNitro))
##ggplot complement

p = CWMnit.vs.Input + geom_boxplot(width=0.4, color="darkgrey", alpha= 0.8)+
  geom_jitter(size = 2, position = position_jitter(width = 0.1))+
  geom_text(x = "terrestrial", y = 0.89, label = "a", family = "serif", 
            hjust = "left",size = 5)+
  geom_text(x = "vertical", y = 0.98, label = "b", family = "serif", 
           hjust = "left", size = 5)+
  theme_classic()+
  labs(x = "Input", 
       y = "CWM Nitrogen")


p + mynamestheme

dev.off()


########CWM leaf width vs Input boxplot

##png object
png("CWwidth.vs.Input.png", width = 16, height = 12, units = "cm",
    pointsize = 8,
    res = 300)

##ggplot basic information
CWwidth.vs.Input <- ggplot(data = data.litter.vt, 
                          aes(x = Input, y = CWwidth))
##ggplot complement

p = CWwidth.vs.Input + geom_boxplot(width=0.4, color="darkgrey", alpha= 0.8)+
  geom_jitter(size = 2, position = position_jitter(width = 0.1))+
  geom_text(x = "terrestrial", y = 3.15, label = "a", family = "serif", 
            hjust = "left",size = 5)+
  geom_text(x = "vertical", y = 2.8, label = "b", family = "serif", 
            hjust = "left", size = 5)+
  theme_classic()+
  labs(x = "Input", 
       y = "CWM Leaf Width")


p + mynamestheme

dev.off()

##########CWM leaf width vs Periods 

##png object
png("CWwidth.vs.Periods.png", width = 16, height = 12, units = "cm", 
    pointsize = 8,
    res = 300)

##ggplot basic information
CWwidth.vs.periods <- ggplot(data = data.litter, 
                                aes(x = Period, y = CWwidth))
##ggplot complement

p = CWwidth.vs.periods + geom_boxplot(width=0.5, color="darkgrey", alpha= 0.8)+
  geom_jitter(size = 2, position = position_jitter(width = 0.1))+
  geom_text(x = "Dry", y = 3.05, label = "a", family = "serif", hjust = "left",  
          size = 5)+
  geom_text(x = "Dry-Rain", y = 3.12, label = "a", family = "serif",
            hjust = "left", size = 5)+
  geom_text(x = "Rain", y = 2.8, label = "ab", family = "serif", hjust = "left", 
          size = 5)+
  geom_text(x = "Rain-Dry", y = 2.8, label = "b", family = "serif", 
            hjust = "left",  size = 5)+
  theme_classic()+
  labs(x = "Periods", 
       y = "CWM Leaf Width")


p + mynamestheme

dev.off()



###periods vs CWM leaf width ridgeline and CWM phenols



##png object
png("CWPhenols.CWwidth.periods.ridgeline.png", width = 16, height = 12,
    units = "cm", pointsize = 8,
    res = 300)


Per <- c("Dry" = "#E59866", "Dry-Rain" = "#F6DDCC", "Rain" = "#5499C7",
         "Rain-Dry" = "#D6EAF8")

CWph <- ggplot(data.litter, aes(x = CWPhenols, y = Period))+
  geom_density_ridges(aes(fill = Period))+
  theme_ridges()+
  theme_classic()+
             labs(y = "Periods", 
              x = "CWM Phenols")

phe <- CWph + scale_fill_manual(values = Per, guide = "none")+
          scale_color_manual("black")+ mynamestheme

CWwidth <- ggplot(data.litter, aes(x = CWwidth, y = Period))+
  geom_density_ridges(aes(fill = Period))+
  theme_ridges()+
  theme_classic()+
  labs(y = "Periods", 
       x = "CWM Leaf Width")

width <- CWwidth + scale_fill_manual(values = Per, guide = "none")+
            scale_color_manual("black")+ mynamestheme

grid.arrange(phe, width, ncol=2)

dev.off()
