
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
litter <- read_table("litter_dynamic.txt")

##separating the variables for the dbFD process
variables <- litter %>% select(Community, Month, Input, Period, Point)

##cleaning the litter matrix
litter.matrix <- litter %>% select(-c(Community, Month, Input, Period, Point))

##removing rows with 0 abundance in every species
litter.matrix.1 <- litter.matrix[which(rowSums(litter.matrix)>0),]

##removing for the variables also
variables <- variables[which(rowSums(litter.matrix)>0),]

# litter traits -----------------------------------------------------------

traits <- read_table("litter_traits.txt")

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

