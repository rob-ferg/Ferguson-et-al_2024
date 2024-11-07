# Set a workspace
setwd('D:/R Code')

install.packages("ggplot2")
install.packages("ggpubr")
install.packages("dplyr")
install.packages("car")
install.packages("lmtest")
install.packages("vegan")
install.packages("DescTools")
install.packages("PMCMRplus")
install.packages("boot")
install.packages("stats")
install.packages("coin")
install.packages("dunn.test")
install.packages("cowplot")
install.packages("mblm")
install.packages("lme4")
install.packages("lmerTest")
install.packages("mgcv")

library(ggplot2)
library(ggpubr)
library(dplyr)
library(car)
library(pracma)
library(lmtest)
library(vegan)
library(DescTools)
library(PMCMRplus)
library(boot)
library(stats)
library(coin)
library(dunn.test)
library(cowplot)
library(mblm)
library(lme4)
library(lmerTest)
library(mgcv)

## CREATE CUSTOM DATAFRAMES ## ------------------------------------------------
# Set working directory
c.bio.df <- read.csv("D:/Ch1_Data")
c.bio.df$Strain <- factor(c.bio.df$Strain, levels = c("Control", "A4", "A5", "G1", "SL1", "Mix", "330", "414", "B3", "C2"))
c.bio.df$Inoculation <- factor(c.bio.df$Inoculation, levels = c("Control", "AMF"))
# Create a dataframe that uses only the strains we can use to compare clades
clade.df <- c.bio.df[c.bio.df$Clade %in% c("I", "II", "VI"),]
# Create a dataframe without Mix and Control, order by organization
c.strain.df <- c.bio.df[c.bio.df$Organization %in% c("Homokaryon", "Heterokaryon"),]
c.strain.df$Strain <- factor(c.strain.df$Strain, levels = c("A4", "A5", "G1", "SL1", "330", "414", "B3", "C2"))
# Create a dataframe without Mix, order by organization
c.nomix.df <- c.bio.df[c.bio.df$Strain %in% c("330", "414", "B3", "C2", "Control", "A4", "A5", "G1", "SL1"),]
c.nomix.df$Strain <- factor(c.nomix.df$Strain, levels = c("Control", "A4", "A5", "G1", "SL1", "330", "414", "B3", "C2"))
# Create a dataframe without Control, order by organization
c.nocon.df <- c.bio.df[c.bio.df$Strain %in% c("330", "414", "B3", "C2", "Mix", "A4", "A5", "G1", "SL1"),]
c.nocon.df$Strain <- factor(c.nocon.df$Strain, levels = c("A4", "A5", "G1", "SL1", "Mix", "330", "414", "B3", "C2"))


## TOTAL DM VS ORGANIZATION ## ------------------------------------------------

# Create plot, choose df, choose variables
total.org.gg <- ggplot(c.nomix.df, aes(x=Organization, y=Total_DM, fill=Organization)) + 
  # Add whiskers
  stat_boxplot(geom = "errorbar", width = 0.15, color = "grey1") + 
  # Add colours
  geom_boxplot(fill = c("#bcbcbc", "firebrick3", "royalblue3"),
               color = "grey1", outlier.colour = "grey1") + 
  # Add axis titles and limits
  labs(x = "AMF Nuclear Organization", y = "Total Dry Mass (g)", base_size = 40, base_family = "") + 
  # Add a custom legend
  scale_fill_manual(values = c("#bcbcbc", "firebrick3", "royalblue3"), 
                    labels = c("Control", "Heterokaryon", "Homokaryon")) +
  # Choose pubr theme
  theme_pubr(base_size = 16) 

# Call the plot
total.org.gg

## Stats: 

# Check homoscedasticity assumption
leveneTest(Total_DM ~ Organization, c.nomix.df)
# Data is homoscedastic

# Create ANOVA model
total.org.anova <- aov((Total_DM) ~ Organization, data = c.nomix.df)
# Check normality of residuals
shapiro.test(total.org.anova$residuals)
hist(total.org.anova$residuals)
qqnorm(total.org.anova$residuals)
# Failed normality test, use nonparametric test
kruskal.test(Total_DM ~ Organization, c.nomix.df)
  # P-value = 0.0246, SIGNIFICANT EFFECT OF ORG ON TOTAL DM 
dunn.test(c.nomix.df$Total_DM, g=c.nomix.df$Organization, method="bonferroni")
  # Con-Het (p=0.0103), Con-Hom (p=0.0274)
DunnettTest(Total_DM ~ Organization, c.nomix.df, conf.level = 0.95)

aggregate(Root_DM ~ Organization, data = c.nomix.df, mean)

## ROOT DM VS ORGANIZATION ## ------------------------------------------------

# Create plot, choose df, choose variables
root.org.gg <- ggplot(c.nomix.df, aes(x=Organization, y= (Root_DM), fill=Organization)) + 
  # Add whiskers
  stat_boxplot(geom = "errorbar", width = 0.15, color = "grey1") + 
  # Add colours
  geom_boxplot(fill = c("#bcbcbc", "firebrick3", "royalblue3"), 
               color = "grey1", outlier.colour = "grey1") + 
  # Add axis titles and limits
  labs(x = "AMF Nuclear Organization", y = "Root Dry Mass (g)", base_size = 40, base_family = "") + 
  # Add a custom legend
  scale_fill_manual(values = c("#bcbcbc", "firebrick3", "royalblue3"), 
                    labels = c("Control", "Heterokaryon", "Homokaryon")) +
  theme(legend.position = "top", legend.key = "rectangular") +
  # Choose pubr theme
  theme_pubr(base_size = 16)
# Call the plot
root.org.gg

## Stats: 

# Check homoscedasticity assumption
leveneTest(Root_DM ~ Organization, c.nomix.df)
# Data is homoscedastic

# Create ANOVA model
root.org.anova <- aov(log(Root_DM) ~ Organization, data = c.nomix.df)
# Check normality of residuals
shapiro.test(root.org.anova$residuals)
hist(root.org.anova$residuals)
qqnorm(root.org.anova$residuals)
# Check summary of model
summary(root.org.anova)
  # P-value = 0.096, MARGINALLY SIGNIFICANT EFFECT OF ORGANIZATION ON ROOT DM #############
TukeyHSD(root.org.anova, conf.level = 0.95)
as.factor(c.nomix.df$Organization)
dunn.test(c.nomix.df$Root_DM, g=c.nomix.df$Organization, method="bonferroni")
DunnettTest(log(Root_DM) ~ Organization, c.nomix.df, control = "Control", alternative = "greater", conf.level = 0.95)
  # P-value = 0.0538: HETEROKARYON MARGINALLY GREATER ROOT DM THAN CONTROL ##############

## ROOT DM VS STRAINS ## -----------------------------------------------

# Create plot, choose df, choose variables
root.strain.gg <- ggplot(c.nomix.df, aes(x=Strain, y= (Root_DM), fill=Strain)) + 
  # Add whiskers
  stat_boxplot(geom = "errorbar", width = 0.15, color = "grey1") + 
  # Add colours
  geom_boxplot(fill = c("#bcbcbc", "firebrick3", "firebrick3", "firebrick3", "firebrick3", "royalblue3", "royalblue3", "royalblue3", "royalblue3"), 
               color = "grey1", outlier.colour = "grey1") + 
  # Add axis titles and limits
  labs(x = "AMF Strain", y = "Root Dry Mass (g)", base_size = 40, base_family = "") + 
  ggtitle("Title") +
  # Choose pubr theme
  theme_pubr(base_size = 16)
# Call the plot
root.strain.gg

## Stats:

# Check homoscedasticity assumption
leveneTest(Root_DM ~ Strain, c.nomix.df)
# Data is homoscedastic

# Create ANOVA model
root.strain.anova <- aov(log(Root_DM) ~ Strain, data = c.nomix.df)
# Check normality of residuals
shapiro.test(root.strain.anova$residuals)
hist(root.strain.anova$residuals)
qqnorm(root.strain.anova$residuals)
# Check summary of model
summary(root.strain.anova)
# P-value = 0.0331: SIGNFICANT EFFECT OF STRAIN ON ROOT DM #######################
# Check as a linear mixed effects model
root.strain.lmer <- lmer(log(Root_DM) ~ (1 | Strain), data = c.nomix.df)
ranova(root.strain.lmer)
# P-value = 0.084: MARGINALLY SIGNIFICANT EFFECT ON ROOT DM ####################
# Use model for multiple comparisons
TukeyHSD(root.strain.anova, conf.level = 0.95)
# No significant differences in multiple comparisons
dunn.test(c.nomix.df$Root_DM, g=c.nomix.df$Strain, method="bonferroni")
DunnettTest(Root_DM ~ Strain, c.nomix.df, control = "Control", alternative = "less", conf.level = 0.95)
# P-value = 0.0378: SL1 SIGNIFICANTLY GREATER THAN CONTROL #######################

## SHOOT P VS STRAIN ## --------------------------------------------------------

# Create plot, choose df, choose variables
p.strain.gg <- ggplot(c.nomix.df, aes(x=Strain, y= P_Conc, fill=Strain)) + 
  # Add whiskers
  stat_boxplot(geom = "errorbar", width = 0.15, color = "grey1") + 
  # Add colours
  geom_boxplot(fill = c("#bcbcbc", "firebrick3", "firebrick3", "firebrick3", "firebrick3", "royalblue3", "royalblue3", "royalblue3", "royalblue3"), 
               color = "grey1", outlier.colour = "grey1") + 
  # Add axis titles and limits
  labs(x = "AMF Strain", y = "[Shoot Phosphorous] (mg/g shoot tissue)", base_size = 40, base_family = "") + 
  # Add a custom title
  ggtitle("Title") +
  # Choose pubr theme
  theme_pubr(base_size = 16)
# Call the plot
p.strain.gg

## Stats:

# Check homoscedasticity assumption
leveneTest(P_Conc ~ Strain, c.bio.df)
# Data is homoscedastic

# Create ANOVA model
p.strain.anova <- aov((P_Conc) ~ Organization, data = c.nomix.df)
# Check normality of residuals
shapiro.test(p.strain.anova$residuals)
hist(p.strain.anova$residuals)
qqnorm(p.strain.anova$residuals)
# Check summary of model
summary(p.strain.anova)
# P-value = 0.000602: SIGNIFICANT EFFECT OF STRAIN ON P LEVELS #########################
P.strain.lmer <- lmer(P_Conc ~ (1 | Strain), data = c.nomix.df)
ranova(P.strain.lmer)
# P-value = 3.28x10^-6: SIGNIFICANT EFFECT OF STRAIN ON P LEVELS #########################
# Conduct a dunnet test against the control
DunnettTest(P_Conc ~ Organization, c.nomix.df, control = "Control", conf.level = 0.95)
# ALL SIGNIFICANT DUNNET P VALUES, SEE RESULTS #############################################
TukeyHSD(p.strain.anova)

aggregate(P_Conc ~ Inoculation, data = c.nomix.df, mean)

## Colonization VS P ## ----------------------------------

# Create plot, add aesthetics, transform variable
col.P.gg <- ggplot(c.bio.df, aes(x=X..Hyphae, y=P_Conc )) +
  # Add point and regression line to plot
  geom_point() + stat_smooth(method = lm, colour = "royalblue3") +
  # Add axis titles to plot
  labs(x = "% AMF Colonization", y = "[Shoot Phosphorous] (mg/g tissue)") +
  # Set theme to pubr
  theme_pubr(base_size = 16)

# Call the ggplot
col.P.gg

## Stats: 

# Create regression model
col.P.lm <- lm((P_Conc) ~ X..Hyphae, data = c.bio.df)
# Check homoscedasticity assumption
bptest(col.P.lm)
# Data is homoscedastic
# Check autocorrelation
dwtest(col.P.lm)
# Data is not auto-correlated
# Check normality of residuals
shapiro.test(col.P.lm$residuals)
hist(col.P.lm$residuals)
qqnorm(col.P.lm$residuals)
# Data is normally distributed
# Check summary of model
summary(col.P.lm)
# P-value = 4.44e^-6; SIGNIFICANT EFFECT OF COLONIZATION ON P -------------------

## SHOOT N VS STRAIN ## --------------------------------------------------------

# Create plot, choose df, choose variables
n.strain.gg <- ggplot(c.nomix.df, aes(x=Strain, y= (N_Conc), fill=Strain)) + 
  # Add whiskers
  stat_boxplot(geom = "errorbar", width = 0.15, color = "grey1") + 
  # Add colours
  geom_boxplot(fill = c("#bcbcbc", "firebrick3", "firebrick3", "firebrick3", "firebrick3", "royalblue3", "royalblue3", "royalblue3", "royalblue3"), 
               color = "grey1", outlier.colour = "grey1") + 
  # Add axis titles and limits
  labs(x = "AMF Strain", y = "[Shoot Nitrogen] (mg/g tissue)", base_size = 40, base_family = "") + 
  # Add a custom legend
  ggtitle("Title") +
  # Choose pubr theme
  theme_pubr(base_size = 16)
# Call the plot
n.strain.gg

## Stats:

# Check homoscedasticity assumption
leveneTest(N_Conc ~ Strain, c.nomix.df)
# Data is homoscedastic

# Create ANOVA model
n.strain.anova <- aov(log(N_Conc) ~ Organization, data = c.nomix.df)
# Check normality of residuals
shapiro.test(n.strain.anova$residuals)
hist(n.strain.anova$residuals)
qqnorm(n.strain.anova$residuals)
# Check summary of model
summary(n.strain.anova)
# P-value = 0.00312: SIGNIFICANT EFFECT OF STRAIN ON N LEVELS #########################
N.strain.lmer <- lmer(log(N_Conc) ~ (1 | Strain), data = c.nomix.df)
ranova(N.strain.lmer)
# P-value = 0.00178: SIGNIFICANT EFFECT OF STRAIN ON N LEVELS #########################
TukeyHSD(n.strain.anova, conf.level = 0.95)
DunnettTest(N_Conc ~ Strain, c.nomix.df, control = "Control", conf.level = 0.95)

## Colonization VS N ## ----------------------------------

# Create plot, add aesthetics, transform variable
col.N.gg <- ggplot(c.bio.df, aes(x=X..Hyphae, y=(N_Conc) )) +
  # Add point and regression line to plot
  geom_point() + stat_smooth(method = lm, colour = "royalblue3") +
  # Add axis titles to plot
  labs(x = "% AMF Colonization", y = "[Shoot Nitrogen] (mg/g tissue)") +
  # Set theme to pubr
  theme_pubr(base_size = 16)

# Call the ggplot
col.N.gg

## Stats: 

# Create regression model
col.N.lm <- lm(log(N_Conc) ~ X..Hyphae, data = c.bio.df)
# Check homoscedasticity assumption
bptest(col.N.lm)
# Data is homoscedastic
# Check autocorrelation
dwtest(col.N.lm)
# Data is not auto-correlated
# Check normality of residuals
shapiro.test(col.N.lm$residuals)
hist(col.N.lm$residuals)
qqnorm(col.N.lm$residuals)
# Data is normally distributed
# Check summary of model
summary(col.N.lm)
# P-value = 0.00578, SIGNIFICANT EFFECT OF % COL ON [N] --------------------------

## SHOOT Mn VS STRAIN ## -------------------------------------------------------

# Create plot, choose df, choose variables
mn.strain.gg <- ggplot(c.nomix.df, aes(x=Strain, y= (Mn_Conc), fill=Strain)) + 
  # Add whiskers
  stat_boxplot(geom = "errorbar", width = 0.15, color = "grey1") + 
  # Add colours
  geom_boxplot(fill = c("#bcbcbc", "firebrick3", "firebrick3", "firebrick3", "firebrick3", "royalblue3", "royalblue3", "royalblue3", "royalblue3"), 
               color = "grey1", outlier.colour = "grey1") + 
  # Add axis titles and limits
  labs(x = "AMF Strain", y = "[Shoot Manganese] (µg/g tissue)", base_size = 40, base_family = "") + 
  # Add a custom legend
  ggtitle("Title") +
  # Choose pubr theme
  theme_pubr(base_size = 16)
# Call the plot
mn.strain.gg

## Stats:

# Check homoscedasticity assumption
leveneTest(Mn_Conc ~ Strain, c.nomix.df)
# Data is homoscedastic

# Create ANOVA model
mn.strain.anova <- aov(log(Mn_Conc) ~ Strain, data = c.nomix.df)
# Check normality of residuals
shapiro.test(mn.strain.anova$residuals)
hist(mn.strain.anova$residuals)
qqnorm(mn.strain.anova$residuals)
# Check summary of model
summary(mn.strain.anova)
# P-value = 0.00356: SIGNIFICANT EFFECT OF STRAIN ON MN CONC ------------------------
Mn.strain.lmer <- lmer(log(Mn_Conc) ~ (1 | Strain), data = c.nomix.df)
ranova(Mn.strain.lmer)
# P-value = 0.009589: SIGNIFICANT EFFECT OF STRAIN ON MN CONC ------------------------
DunnettTest(log(Mn_Conc) ~ Strain, c.nomix.df)
# P-values < 0.05: A4, A5, SL1, 414 AND C2 GREATER THAN CONTROL --------------------
TukeyHSD(mn.strain.anova, conf.level = 0.95)

## Colonization VS Mn ## ----------------------------------

# Create plot, add aesthetics, transform variable
col.Mn.gg <- ggplot(c.bio.df, aes(x=X..Hyphae, y=(Mn_Conc) )) +
  # Add point and regression line to plot
  geom_point() + stat_smooth(method = lm, , colour = "royalblue3") +
  # Add axis titles to plot
  labs(x = "% AMF Colonization", y = "[Shoot Manganese] (µg/g tissue)") +
  # Set theme to pubr
  theme_pubr(base_size = 16)

# Call the ggplot
col.Mn.gg

## Stats: 

# Create regression model
col.Mn.lm <- lm(log(Mn_Conc) ~ X..Hyphae, data = c.bio.df)
# Check homoscedasticity assumption
bptest(col.Mn.lm)
# Data is homoscedastic
# Check autocorrelation
dwtest(col.Mn.lm)
# Data is not auto-correlated
# Check normality of residuals
shapiro.test(col.Mn.lm$residuals)
hist(col.Mn.lm$residuals)
qqnorm(col.Mn.lm$residuals)
# Data is normally distributed
# Check summary of model
summary(col.Mn.lm)
# P-value = 0.0095: SIGNIFICANT EFFECT OF COL ON MANGANESE ----------------------

## SHOOT Fe VS STRAIN ## -------------------------------------------------------

# Create plot, choose df, choose variables
fe.strain.gg <- ggplot(c.nomix.df, aes(x=Strain, y= (Fe_Conc), fill=Strain)) + 
  # Add whiskers
  stat_boxplot(geom = "errorbar", width = 0.15, color = "grey1") + 
  # Add colours
  geom_boxplot(fill = c("#bcbcbc", "firebrick3", "firebrick3", "firebrick3", "firebrick3", "royalblue3", "royalblue3", "royalblue3", "royalblue3"), 
               color = "grey1", outlier.colour = "grey1") + 
  # Add axis titles and limits
  labs(x = "AMF Strain", y = "[Shoot Iron] (µg/g tissue)", base_size = 40, base_family = "") + 
  # Add a custom legend
  ggtitle("Title") +
  # Choose pubr theme
  theme_pubr(base_size = 16)
# Call the plot
fe.strain.gg

## Stats:

# Check homoscedasticity assumption
leveneTest(Fe_Conc ~ Strain, c.nomix.df)
# Data is homoscedastic

# Create ANOVA model
fe.strain.anova <- aov(log(Fe_Conc) ~ Strain, data = c.nomix.df)
# Check normality of residuals
shapiro.test(fe.strain.anova$residuals)
hist(fe.strain.anova$residuals)
qqnorm(fe.strain.anova$residuals)
# Check summary of model
summary(fe.strain.anova)
# P-value = 0.0436: SIGNIFICANT EFFECT OF STRAIN ON FE CONC -------------------
Fe.strain.lmer <- lmer((Fe_Conc) ~ (1 | Strain), data = c.nomix.df)
ranova(Fe.strain.lmer)
TukeyHSD(fe.strain.anova, conf.level = 0.95)
DunnettTest(Fe_Conc ~ Strain, c.nomix.df)
# P-value > 0.05: G1 AND 414 GREATER THAN CONTROL ------------------------------

## Colonization VS Fe ## ----------------------------------

# Create plot, add aesthetics, transform variable
col.Fe.gg <- ggplot(c.bio.df, aes(x=X..Hyphae, y=(Fe_Conc) )) +
  # Add point and regression line to plot
  geom_point() + stat_smooth(method = lm, , colour = "royalblue3") +
  # Add axis titles to plot
  labs(x = "% AMF Colonization", y = "[Shoot Iron] (µg/g tissue)") +
  # Set theme to pubr
  theme_pubr(base_size = 16)

# Call the ggplot
col.Fe.gg

## Stats: 

# Create regression model
col.Fe.lm <- lm(log(Fe_Conc) ~ X..Hyphae, data = c.bio.df)
# Check homoscedasticity assumption
bptest(col.Fe.lm)
# Data is homoscedastic
# Check autocorrelation
dwtest(col.Fe.lm)
# Data is not auto-correlated
# Check normality of residuals
shapiro.test(col.Fe.lm$residuals)
hist(col.Fe.lm$residuals)
qqnorm(col.Fe.lm$residuals)
# Data is normally distributed
# Check summary of model
summary(col.Fe.lm)
# P-value = 0.0023: SIGNIFICANT EFFECT OF COL ON FE

## SHOOT Mg VS STRAIN ## -------------------------------------------------------

# Create plot, choose df, choose variables
mg.strain.gg <- ggplot(c.nomix.df, aes(x=Strain, y= (Mg_Conc), fill=Strain)) + 
  # Add whiskers
  stat_boxplot(geom = "errorbar", width = 0.15, color = "grey1") + 
  # Add colours
  geom_boxplot(fill = c("#bcbcbc", "firebrick3", "firebrick3", "firebrick3", "firebrick3", "royalblue3", "royalblue3", "royalblue3", "royalblue3"), 
               color = "grey1", outlier.colour = "grey1") + 
  # Add axis titles and limits
  labs(x = "AMF Strain", y = "[Shoot Magnesium] (µg/g tissue)", base_size = 40, base_family = "") + 
  # Add a custom legend
  ggtitle("Title") +
  # Choose pubr theme
  theme_pubr(base_size = 16)
# Call the plot
mg.strain.gg

## Stats:

# Check homoscedasticity assumption
leveneTest(Mg_Conc ~ Strain, c.nomix.df)
# Data is homoscedastic

# Create ANOVA model
mg.strain.anova <- aov(log(Mg_Conc) ~ Strain, data = c.nomix.df)
# Check normality of residuals
shapiro.test(mg.strain.anova$residuals)
hist(mg.strain.anova$residuals)
qqnorm(mg.strain.anova$residuals)
# Check summary of model
summary(mg.strain.anova)
# P-value = 0.0436: SIGNIFICANT EFFECT OF STRAIN ON FE CONC -------------------
Mg.strain.lmer <- lmer((Mg_Conc) ~ (1 | Strain), data = c.nomix.df)
ranova(Mg.strain.lmer)
TukeyHSD(mg.strain.anova, conf.level = 0.95)
DunnettTest(Mg_Conc ~ Strain, c.nomix.df)
# P-value > 0.05: G1 AND 414 GREATER THAN CONTROL ------------------------------

## Colonization VS Mg ## ----------------------------------

# Create plot, add aesthetics, transform variable
col.Mg.gg <- ggplot(c.bio.df, aes(x=X..Hyphae, y=(Mg_Conc) )) +
  # Add point and regression line to plot
  geom_point() + stat_smooth(method = lm, , colour = "royalblue3") +
  # Add axis titles to plot
  labs(x = "% AMF Colonization", y = "[Shoot Magnesium] (µg/g tissue)") +
  # Set theme to pubr
  theme_pubr(base_size = 16)

# Call the ggplot
col.Mg.gg

## Stats: 

# Create regression model
col.Mg.lm <- lm(log(Mg_Conc) ~ X..Hyphae, data = c.bio.df)
# Check homoscedasticity assumption
bptest(col.Mg.lm)
# Data is homoscedastic
# Check autocorrelation
dwtest(col.Mg.lm)
# Data is not auto-correlated
# Check normality of residuals
shapiro.test(col.Mg.lm$residuals)
hist(col.Mg.lm$residuals)
qqnorm(col.Mg.lm$residuals)
# Data is normally distributed
# Check summary of model
summary(col.Mg.lm)
# P-value = 0.0895: MARGINAL EFFECT OF COLONIZATION ON MAGNESIUM ----------------

## N VS Fe ## ----------------------------------

# Create plot, add aesthetics, transform variable
N.Fe.gg <- ggplot(c.nomix.df, aes(x=N_Conc, y=(Fe_Conc))) +
  # Add point and regression line to plot
  geom_point() + stat_smooth(method = lm, , colour = "royalblue3") +
  # Add axis titles to plot
  labs(x = "[Shoot Nitrogen] (mg/g tissue)", y = "[Shoot Iron] (µg/g tissue)") +
  # Set theme to pubr
  theme_pubr(base_size = 16)

# Call the ggplot
N.Fe.gg

## Stats: 

cor(c.nomix.df$N_Conc, c.nomix.df$Fe_Conc)

# Create regression model
N.Fe.lm <- lm(log(Fe_Conc) ~ N_Conc, data = c.nomix.df)
# Check homoscedasticity assumption
bptest(N.Fe.lm)
# Data is homoscedastic
# Check autocorrelation
dwtest(N.Fe.lm)
# Data is not auto-correlated
# Check normality of residuals
shapiro.test(N.Fe.lm$residuals)
hist(N.Fe.lm$residuals)
qqnorm(N.Fe.lm$residuals)
# Data is normally distributed
# Check summary of model
summary(N.Fe.lm)
  # P = 0.000477, significant effect of N on Fe ---------------------------------

## BULK MYCORRHIZAL C VS ORGANIZATION ## --------------------------------------------------------

# Create plot, choose df, choose variables
bulk.C1.org.gg <- ggplot(c.nomix.df, aes(x=Organization, y= Bulk_Mycorrhizal_C1, fill=Organization)) + 
  # Add whiskers
  stat_boxplot(geom = "errorbar", width = 0.15, color = "grey1") + 
  # Add colours
  geom_boxplot(fill = c("#bcbcbc", "firebrick3", "royalblue3"), 
               color = "grey1", outlier.colour = "grey1") + 
  # Add axis titles and limits
  labs(x = "AMF Nuclear Organization", y = "Mycorrhizal C (mg/g soil)", base_size = 40, base_family = "") + 
  # Add a custom legend
  # Choose pubr theme
  theme_pubr(base_size = 16)
# Call the plot
bulk.C1.org.gg

kruskal.test(Bulk_Mycorrhizal_C1 ~ Strain, data = c.strain.df)
# P-value > 0.05: NO EFFECT OF ORGANIZATION

## INCUBATED MYCORRHIZAL MAOC VS STRAINS (NO CONTROL/MIX) ## --------------------

# Create plot, choose df, choose variables
incMAOM.C1.strain.gg <- ggplot(c.strain.df, aes(x=Strain, y= (IncMAOM_Mycorrhizal_C1), fill=Strain)) + 
  # Add whiskers
  stat_boxplot(geom = "errorbar", width = 0.15, color = "grey1") + 
  # Add colours
  geom_boxplot(fill = c("firebrick3", "firebrick3", "firebrick3", "firebrick3", "royalblue3", "royalblue3", "royalblue3", "royalblue3"), 
               color = "grey1", outlier.colour = "grey1") + 
  # Add axis titles and limits
  labs(x = "AMF Strain", y = "Mycorrhizal MAOC (mg/g MAOM soil)", base_size = 40, base_family = "") + 
  # Add a custom legend
  ggtitle("Title") +
  # Choose pubr theme
  theme_pubr(base_size = 16)
# Call the plot
incMAOM.C1.strain.gg

## Stats:

# Check homoscedasticity assumption
leveneTest(IncMAOM_Mycorrhizal_C1 ~ Strain, c.strain.df)
# Data is homoscedastic

# Create ANOVA model
incMAOM.C1.strain.anova <- aov(log(IncMAOM_Mycorrhizal_C1) ~ Organization, data = c.strain.df)
# Check normality of residuals
shapiro.test(incMAOM.C1.strain.anova$residuals)
hist(incMAOM.C1.strain.anova$residuals)
qqnorm(incMAOM.C1.strain.anova$residuals)
# Check summary of model
summary(incMAOM.C1.strain.anova)
# P-value = 0.01, SIGNIFICANT EFFECT OF STRAIN ################################################

MAOC.strain.lmer <- lmer(log(IncMAOM_Mycorrhizal_C1) ~ (1 | Strain), data = c.strain.df)
ranova(MAOC.strain.lmer)
# P-value = 0.04, SIGNIFICANT EFFECT OF STRAIN ################################################

# Use model for multiple comparisons
TukeyHSD(incMAOM.C1.strain.anova, conf.level = 0.95)
# G1 > B3, 330 > B3 #################################################################

## BULK MYCORRHIZAL C VS TOTAL DM ## ----------------------------------

# Create plot, add aesthetics, transform variable
bulk.c1.total.gg <- ggplot(c.bio.df, aes(x=Total_DM, y=(Bulk_Mycorrhizal_C1))) +
  # Add point and regression line to plot
  geom_point() + stat_smooth(method = lm, , colour = "royalblue3") +
  # Add axis titles to plot
  labs(x = "Total Dry Mass (g)", y = "Mycorrhizal C (mg/g soil)") +
  # Set theme to pubr
  theme_pubr(base_size = 16)

# Call the ggplot
bulk.c1.total.gg

# Stats

bulk.c1.total.lm <- lm(log(Bulk_Mycorrhizal_C1) ~ Total_DM, data = c.bio.df)
# Check homoscedasticity assumption
bptest(bulk.c1.total.lm)
# Data is homoscedastic
# Check autocorrelation
dwtest(bulk.c1.total.lm)
# Data is not auto-correlated
# Check normality of residuals
shapiro.test(bulk.c1.total.lm$residuals)
hist(bulk.c1.total.lm$residuals)
qqnorm(bulk.c1.total.lm$residuals)
par(mfrow = c(2, 2), las = 1)
plot(bulk.c1.total.lm)
# Data not entirely normal, but should be okay
# Check summary of model
summary(bulk.c1.total.lm)
# P-value = 0.06, MARGINAL POSITIVE EFFECT OF TOTAL DM ON C INPUTS ##########################

## HETEROKARYON HAPLOTYPE RATIOS ## --------------------------------------------

# Create a dataframe for the ratio data
ratios.df <- read.csv("D:/Ch1_Ratios_R.csv")

# Create plot, choose df, choose variables
all.ratios.gg2 <- ggplot(ratios.df, aes(x=Treatment, y= Ratio, fill=Treatment)) + 
  stat_boxplot(geom = "errorbar", width = 0.15, color = "grey1") +
  geom_segment(aes(x = 0, xend = 5, y = 1, yend = 1), color = "grey1") +
  geom_boxplot(color = "grey1", outlier.colour = "grey1") + 
  labs(x = "AMF Strain", y = "MAT-A:MAT-B", fill = NULL) + ylim(c(0, 1.5)) + 
  scale_fill_manual(values = c("firebrick3", "royalblue3", "sienna2", "sienna2"),
                    labels = c("MAT1:MAT2", "MAT6:MAT3", "MAT1:MAT5", "MAT1:MAT5")) + 
  theme(legend.position = "right", legend.location = "plot", 
        legend.key = element_rect(), 
        legend.key.size = unit(4, "lines")) +  # Adjust the size here
  # Add a custom key
  guides(fill = guide_legend(title.position = "top",
                             title.hjust = 0,
                             direction = "horizontal",
                             keyheight = unit(2, "lines"),
                             keywidth = unit(2, "lines"))) +
  theme_pubr(base_size = 16)
# Call the plot
all.ratios.gg2

## CREATE CUSTOM LEGENDS FOR FIGURES ## ------------------------------------------

## CUSTOM LEGEND: ORGANIZATION ## ------------------------------------------------

# Create plot, choose df, choose variables
leg.org.gg <- ggplot(c.nomix.df, aes(x=Organization, y=Total_DM, fill=Organization)) + 
  # Add whiskers
  stat_boxplot(geom = "errorbar", width = 0.15, color = "grey1") + 
  # Add colours
  geom_boxplot(color = "grey1", outlier.colour = "grey1") + 
  # Add axis titles and limits
  labs(x = "AMF Nuclear Organization", y = "Total Dry Mass (g)", base_size = 40, base_family = "", fill = NULL) + 
  # Add a custom legend
  scale_fill_manual(values = c("#bcbcbc", "firebrick3", "royalblue3"), 
                    labels = c("Control", "Heterokaryon", "Homokaryon")) +
  theme(legend.position = "top", 
        legend.key = element_rect(), 
        legend.key.size = unit(4, "lines")) +  # Adjust the size here
  # Add a custom key
  guides(fill = guide_legend(title.position = "top",
                             title.hjust = 0,
                             direction = "horizontal",
                             keyheight = unit(2, "lines"),
                             keywidth = unit(2, "lines"))) +
  # Choose pubr theme
  theme_pubr(base_size = 20)

# Call the plot
leg.org.gg

## CUSTOM LEGEND: ORGANIZATION (NO CONTROL) ## ------------------------------------------------

# Create plot, choose df, choose variables
leg.org2.gg <- ggplot(c.strain.df, aes(x=Organization, y=Total_DM, fill=Organization)) + 
  # Add whiskers
  stat_boxplot(geom = "errorbar", width = 0.15, color = "grey1") + 
  # Add colours
  geom_boxplot(color = "grey1", outlier.colour = "grey1") + 
  # Add axis titles and limits
  labs(x = "AMF Nuclear Organization", y = "Total Dry Mass (g)", base_size = 40, base_family = "", fill = NULL) + 
  # Add a custom legend
  scale_fill_manual(values = c("firebrick3", "royalblue3"), 
                    labels = c("Heterokaryon", "Homokaryon")) +
  theme(legend.position = "top", 
        legend.key = element_rect(), 
        legend.key.size = unit(4, "lines")) +  # Adjust the size here
  # Add a custom key
  guides(fill = guide_legend(title.position = "top",
                             title.hjust = 0,
                             direction = "horizontal",
                             keyheight = unit(2, "lines"),
                             keywidth = unit(2, "lines"))) +
  # Choose pubr theme
  theme_pubr(base_size = 25)

# Call the plot
leg.org2.gg

