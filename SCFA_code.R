# Code used in the statistical analyses of the paper:
# Jacobo de la Cuesta-Zuluaga, Noel T. Mueller, Rafael Álvarez-Quintero, Eliana P. Velásquez-Mejía, Jelver A. Sierra, Vanessa Corrales-Agudelo, Jenny A. Carmona, José M. Abad, Juan S. Escobar.
# 2018
# Higher fecal short-chain fatty acid levels are associated with gut microbiome dysbiosis, obesity, hypertension and cardiometabolic disease risk factors 
# Nutrients. doi: .
# www.
# Vidarium (c) 2018.


# Initial useful commands ----
# Clean the workspace
#rm(list = ls())

# Seed for random generation
set.seed(5600)

# Color-blind-friendly palette
cbPalette = c("#999999","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7")


# Some useful functions ----
# Copy-paste R tables into Excel
write.excel <- function(x,row.names=TRUE,col.names=TRUE,...) {
  write.table(x,"clipboard",sep="\t",row.names=row.names,col.names=col.names,...)
}

# Calculate standard error of the mean
sem = function(x) sd(x)/sqrt(length(x))


# Libraries ----
library(GUniFrac) # To calculate UniFrac distances
library(phytools) # To load phylogenetic tree
library(reshape2) # Function melt
library(dplyr) # Function ddply
library(yarrr) # To make pirateplots
library(rms) # To perform restricted cubic splines
library(Hmisc) # To plot restricted cubic splines
library(cowplot) # To put togeher different plots (plot_grid)
library(psych) # Correlation tests with FDR adjustment (corr.test)
library(car) # Type-II anova
library(qvalue) # To calculate FDR-adjusted p-values
library(NMF) # To make heatmaps
library(geepack) # To calculate prevalence risks
library(nnet) # To estimate parameters of multinomial regressions
library(sandwich) # To calculate robust variance estimates for generalized linear models

# Loading data ----
# Change directory name to where data are saved in your computer
setwd(dir = "XXXX")

# Metadata
all.meta = read.table(file = "D:/Vidarium/GitHub/scfa/scfa.meta", header = T, sep = "\t", dec = ".", row.names = 1)

  # Calculate mean arterial pressure
mean_bp = ((2*all.meta$diastolic_bp)+all.meta$systolic_bp)/3
all.meta = cbind(all.meta, mean_bp)

# Microbiota data
  # OTU table with raw read counts
microbio.otus = read.table(file = "D:/Vidarium/Publicaciones/CAGs/reproducibility/analisis/microbio_selected.otus", header = T, sep = "\t", row.names = 1)

  # OTU table with relative abundances
microbio.relative = t(microbio.otus/rowSums(microbio.otus))

  # Rarefied OTU table
  # Be sure to set the seed if you want reproducible data (see line 15 above)
microbio.rare = Rarefy(microbio.otus)$otu.tab.rff

  # Distance tree (fast neighbor-joining tree)
microbio.tree = read.newick(file = "D:/Vidarium/Publicaciones/CAGs/reproducibility/analisis/microbio_selected.tre")

  # Greengenes 13_8_99 taxonomy for each OTU
microbio.taxonomy = read.table(file = "D:/Vidarium/Publicaciones/CAGs/reproducibility/analisis/microbio_selected.taxonomy", sep = "\t", row.names = 1, header = T)

  # In the files, there are five samples that were sequenced twice
replicate_samples = c("MI_008_H2", "MI_093_H12", "MI_130_H2", "MI_198_H2", "MI_458_H2")
replicate_positions = c(9, 95, 132, 201, 445)

# Remove replicates for the final datasets
microbio.otus = microbio.otus[-replicate_positions,]
microbio.rare = microbio.rare[-replicate_positions,]
microbio.relative = microbio.relative[,-replicate_positions]
all.meta = all.meta[-replicate_positions,]

# Calculate tertiles for all SCFAs and microbiota alpha diversity (OTU richness)
all.meta = within(all.meta, SCFA_tertile <- as.integer(cut(all.meta$SCFA, quantile(all.meta$SCFA, probs = 0:3/3), include.lowest = TRUE)))
all.meta = within(all.meta, acetate_tertile <- as.integer(cut(all.meta$acetate, quantile(all.meta$acetate, probs = 0:3/3), include.lowest = TRUE)))
all.meta = within(all.meta, propionate_tertile <- as.integer(cut(all.meta$propionate, quantile(all.meta$propionate, probs = 0:3/3), include.lowest = TRUE)))
all.meta = within(all.meta, butyrate_tertile <- as.integer(cut(all.meta$butyrate, quantile(all.meta$butyrate, probs = 0:3/3), include.lowest = TRUE)))
all.meta = within(all.meta, isobutyrate_tertile <- as.integer(cut(all.meta$isobutyrate, quantile(all.meta$isobutyrate, probs = 0:3/3), include.lowest = TRUE)))
all.meta = within(all.meta, richness_tertile <- as.integer(cut(all.meta$richness, quantile(all.meta$richness, probs = 0:3/3), include.lowest = TRUE)))

# Subset of participants that did not take medication and did not smoke (sensitivity analyses)
clean.meta = all.meta[all.meta$smoker =="No" & all.meta$medicament == "No",]
clean.otus = read.table(file = "D:/Vidarium/Publicaciones/SCFA/clean.otus", header = T, sep = "\t", row.names = 1)
clean.rare = Rarefy(clean.otus)$otu.tab.rff

clean.meta = clean.meta[complete.cases(clean.meta[,146]),]


# Analysis of raw data (no adjustment for confounders) ----
non_numeric = c(1,3,24:27) # non-numeric variables
quantitative = all.meta[,-non_numeric] # subset of quantitative data only

# transformed data for further ANOVA
tr_quantitative = data.frame(age = quantitative$age, calories = log(quantitative$calories),
                             fiber = log(quantitative$fiber), physical_activity = log(1+quantitative$physical_activity),
                             bmi = log(quantitative$bmi), body_fat = sqrt(asin(quantitative$body_fat/100)),
                             waist = log(quantitative$waist), HDL = log(quantitative$HDL),
                             LDL = quantitative$LDL, VLDL = log(quantitative$VLDL), triglycerides = log(quantitative$triglycerides),
                             hsCRP = log(quantitative$hsCRP), glucose = log(quantitative$glucose),
                             HbA1c = sqrt(asin(quantitative$HbA1c/100)), insulin = log(quantitative$insulin), 
                             HOMA_IR = log(quantitative$HOMA_IR), leptin = log(quantitative$leptin),
                             adiponectin = log(0.5+quantitative$adiponectin), LBP = quantitative$LBP,
                             systolic_bp = quantitative$systolic_bp, diastolic_bp = quantitative$diastolic_bp, mean_bp = quantitative$mean_bp, 
                             richness = quantitative$richness, SCFA = log(0.1+quantitative$SCFA),
                             acetate = log(0.1+quantitative$acetate), propionate = log(0.05+quantitative$propionate),
                             butyrate = log(0.01+quantitative$butyrate), isobutyrate = log(0.001+quantitative$isobutyrate),
                             SCFA_tertile = quantitative$SCFA_tertile, acetate_tertile = quantitative$acetate_tertile,
                             propionate_tertile = quantitative$propionate_tertile, butyrate_tertile = quantitative$butyrate_tertile,
                             isobutyrate_tertile = quantitative$isobutyrate_tertile, richness_tertile = quantitative$richness_tertile)

# table variable names for further use
table1_headings = names(quantitative)
tr_table1_headings = names(tr_quantitative)

  # Table 1 ----
  # Means and SEM overall data
apply(quantitative, MARGIN = 2, FUN = mean, na.rm = T)
apply(quantitative, MARGIN = 2, FUN = sem)
all.meta %>% 
     group_by(sex) %>%
     summarise(N = length(sex)/nrow(all.meta))

  # Means and SEM per butyrate tertiles
table1A_mean = lapply(table1_headings, function(x) {
  aggregate(quantitative[x], by = c(quantitative['butyrate_tertile']),
            FUN = mean, na.rm = T)
})
write.excel(table1A_mean)


table1A_sem = lapply(table1_headings, function(x) {
  aggregate(quantitative[x], by = c(quantitative['butyrate_tertile']),
            FUN = sem)
})
write.excel(table1A_sem)

  # To calculate sex proportions
all.meta %>% 
  group_by(butyrate_tertile, sex) %>%
  summarise(n = n()/147) # 147 is the number of observations in each tertile

    # Testing for differences by butyrate tertile
table1A_aov = lapply(tr_table1_headings, function(x) {
  mF = formula(paste(x, "~ as.factor(quantitative$butyrate_tertile)"))
  summary(aov(mF, data = tr_quantitative))
  })

chisq.test(all.meta$sex, as.factor(all.meta$butyrate_tertile)) # for sex


  # Means and SEM per OTU richness tertiles
table1B_mean = lapply(table1_headings, function(x) {
  aggregate(quantitative[x], by = c(quantitative['richness_tertile']),
            FUN = mean, na.rm = T)
})
write.excel(table1B_mean)

table1B_sem = lapply(table1_headings, function(x) {
  aggregate(quantitative[x], by = c(quantitative['richness_tertile']),
            FUN = sem)
})
write.excel(table1B_sem)

    # Testing for differences by OTU richness tertile
table1B_aov = lapply(tr_table1_headings, function(x) {
  mF = formula(paste(x, "~ as.factor(quantitative$richness_tertile)"))
  summary(aov(mF, data = tr_quantitative))
})

chisq.test(all.meta$sex, as.factor(all.meta$richness_tertile)) # for sex


  # Supplementary Table S2 ----
  # Means and SEM per total SCFA tertiles
tableS2A_mean = lapply(table1_headings, function(x) {
  aggregate(quantitative[x], by = c(quantitative['SCFA_tertile']),
            FUN = mean, na.rm = T)
})
write.excel(tableS2A_mean)

tableS2A_sem = lapply(table1_headings, function(x) {
  aggregate(quantitative[x], by = c(quantitative['SCFA_tertile']),
            FUN = sem)
})
write.excel(tableS2A_sem)

    # Testing for differences by total SCFA tertile
tableS2A_aov = lapply(tr_table1_headings, function(x) {
  mF = formula(paste(x, "~ as.factor(quantitative$SCFA_tertile)"))
  summary(aov(mF, data = tr_quantitative))
})

chisq.test(all.meta$sex, as.factor(all.meta$SCFA_tertile)) # for sex


  # Means and SEM per acetate tertiles
tableS2B_mean = lapply(table1_headings, function(x) {
  aggregate(quantitative[x], by = c(quantitative['acetate_tertile']),
            FUN = mean, na.rm = T)
})
write.excel(tableS2B_mean)

tableS2B_sem = lapply(table1_headings, function(x) {
  aggregate(quantitative[x], by = c(quantitative['acetate_tertile']),
            FUN = sem)
})
write.excel(tableS2B_sem)

    # Testing for differences by acetate tertile
tableS2B_aov = lapply(tr_table1_headings, function(x) {
  mF = formula(paste(x, "~ as.factor(quantitative$acetate_tertile)"))
  summary(aov(mF, data = tr_quantitative))
})

chisq.test(all.meta$sex, as.factor(all.meta$acetate_tertile)) # for sex


  # Means and SEM per propionate tertiles
tableS2C_mean = lapply(table1_headings, function(x) {
  aggregate(quantitative[x], by = c(quantitative['propionate_tertile']),
            FUN = mean, na.rm = T)
})
write.excel(tableS2C_mean)

tableS2C_sem = lapply(table1_headings, function(x) {
  aggregate(quantitative[x], by = c(quantitative['propionate_tertile']),
            FUN = sem)
})
write.excel(tableS2C_sem)

    # Testing for differences by propionate tertile
tableS2C_aov = lapply(tr_table1_headings, function(x) {
  mF = formula(paste(x, "~ as.factor(quantitative$propionate_tertile)"))
  summary(aov(mF, data = tr_quantitative))
})

chisq.test(all.meta$sex, as.factor(all.meta$propionate_tertile)) # for sex


  # Means and SEM per isobutyrate tertiles
tableS2D_mean = lapply(table1_headings, function(x) {
  aggregate(quantitative[x], by = c(quantitative['isobutyrate_tertile']),
            FUN = mean, na.rm = T)
})
write.excel(tableS2D_mean)

tableS2D_sem = lapply(table1_headings, function(x) {
  aggregate(quantitative[x], by = c(quantitative['isobutyrate_tertile']),
            FUN = sem)
})
write.excel(tableS2D_sem)

    # Testing for differences by isobutyrate tertile
tableS2C_aov = lapply(tr_table1_headings, function(x) {
  mF = formula(paste(x, "~ as.factor(quantitative$isobutyrate_tertile)"))
  summary(aov(mF, data = tr_quantitative))
})

chisq.test(all.meta$sex, as.factor(all.meta$isobutyrate_tertile)) # for sex


# Analysis of adjusted variables ----
  # Adjustment by age, city of origin, physical activity, caloric intake and fiber intake

# Calculating residuals from adjusted models
adjusted_SCFA = residuals(lm(log(0.1+SCFA)~age+city+calories+physical_activity+fiber, data = all.meta))
adjusted_acetate = residuals(lm(log(0.1+acetate)~age+city+calories+physical_activity+fiber, data = all.meta))
adjusted_propionate = residuals(lm(log(0.05+propionate)~age+city+calories+physical_activity+fiber, data = all.meta))
adjusted_butyrate = residuals(lm(log(0.01+butyrate)~age+city+calories+physical_activity+fiber, data = all.meta))
adjusted_isobutyrate = residuals(lm(log(0.001+isobutyrate)~age+city+calories+physical_activity+fiber, data = all.meta))
adjusted_richness = residuals(lm(richness~age+city+calories+physical_activity+fiber, data = all.meta))
adjusted_bmi = residuals(lm(bmi~age+city+calories+physical_activity+fiber, data = all.meta))
adjusted_body_fat = residuals(lm(sqrt(asin(body_fat/100))~age+city+calories+physical_activity+fiber, data = all.meta))
adjusted_waist = residuals(lm(log(waist)~age+city+calories+physical_activity+fiber, data = all.meta))
adjusted_HDL = residuals(lm(log(HDL)~age+city+calories+physical_activity+fiber, data = all.meta))
adjusted_LDL = residuals(lm(LDL~age+city+calories+physical_activity+fiber, data = all.meta))
adjusted_VLDL = residuals(lm(log(VLDL)~age+city+calories+physical_activity+fiber, data = all.meta))
adjusted_triglycerides = residuals(lm(log(triglycerides)~age+city+calories+physical_activity+fiber, data = all.meta))
adjusted_hsCRP = residuals(lm(log(hsCRP)~age+city+calories+physical_activity+fiber, data = all.meta))
adjusted_glucose = residuals(lm(log(glucose)~age+city+calories+physical_activity+fiber, data = all.meta))
adjusted_HbA1c = residuals(lm(sqrt(asin(HbA1c/100))~age+city+calories+physical_activity+fiber, data = all.meta))
adjusted_insulin = residuals(lm(log(insulin)~age+city+calories+physical_activity+fiber, data = all.meta))
adjusted_HOMA_IR = residuals(lm(log(HOMA_IR)~age+city+calories+physical_activity+fiber, data = all.meta))
adjusted_leptin = residuals(lm(log(leptin)~age+city+calories+physical_activity+fiber, data = all.meta))
adjusted_adiponectin = residuals(lm(log(0.5+adiponectin)~age+city+calories+physical_activity+fiber, data = all.meta))
adjusted_LBP = residuals(lm(LBP~age+city+calories+physical_activity+fiber, data = all.meta))
adjusted_systolic_bp = residuals(lm(systolic_bp~age+city+calories+physical_activity+fiber, data = all.meta))
adjusted_diastolic_bp = residuals(lm(diastolic_bp~age+city+calories+physical_activity+fiber, data = all.meta))
adjusted_mean_bp = residuals(lm(mean_bp~age+city+calories+physical_activity+fiber, data = all.meta))

all.meta = cbind(all.meta, adjusted_SCFA, adjusted_acetate, adjusted_propionate, adjusted_butyrate, adjusted_isobutyrate, adjusted_richness)
all.meta = within(all.meta, adj_SCFA_tertile <- as.integer(cut(adjusted_SCFA, quantile(adjusted_SCFA, probs=0:3/3), include.lowest=TRUE)))
all.meta = within(all.meta, adj_acetate_tertile <- as.integer(cut(adjusted_acetate, quantile(adjusted_acetate, probs=0:3/3), include.lowest=TRUE)))
all.meta = within(all.meta, adj_propionate_tertile <- as.integer(cut(adjusted_propionate, quantile(adjusted_propionate, probs=0:3/3), include.lowest=TRUE)))
all.meta = within(all.meta, adj_butyrate_tertile <- as.integer(cut(adjusted_butyrate, quantile(adjusted_butyrate, probs=0:3/3), include.lowest=TRUE)))
all.meta = within(all.meta, adj_isobutyrate_tertile <- as.integer(cut(adjusted_isobutyrate, quantile(adjusted_isobutyrate, probs=0:3/3), include.lowest=TRUE)))
all.meta = within(all.meta, adj_richness_tertile <- as.integer(cut(adjusted_richness, quantile(adjusted_richness, probs=0:3/3), include.lowest=TRUE)))

  # Figure 1 ----
  # Panel A
pirateplot(formula = richness ~ adj_SCFA_tertile,
           data = all.meta,
           #  pal = "google",
           pal = c("#999999","#E69F00","#56B4E9"),
           gl = 0, point.o = 0.8,
           inf.method = "ci",
           xlab = "Adjusted total SCFA tertiles",
           ylab = "OTU richness")
summary(aov(richness ~ adj_SCFA_tertile, data=all.meta))

  # Panel B
pirateplot(formula = richness ~ adj_acetate_tertile,
           data = all.meta,
           #  pal = "google",
           pal = c("#999999","#E69F00","#56B4E9"),
           gl = 0, point.o = 0.8,
           inf.method = "ci",
           xlab = "Adjusted acetate tertiles",
           ylab = "OTU richness")
summary(aov(richness ~ adj_acetate_tertile, data=all.meta))

  # Panel C
pirateplot(formula = richness ~ adj_propionate_tertile,
           data = all.meta,
           #  pal = "google",
           pal = c("#999999","#E69F00","#56B4E9"),
           gl = 0, point.o = 0.8,
           inf.method = "ci",
           xlab = "Adjusted propionate tertiles",
           ylab = "OTU richness")
summary(aov(richness ~ adj_propionate_tertile, data=all.meta))

  # Panel D
pirateplot(formula = richness ~ adj_butyrate_tertile,
           data = all.meta,
           #  pal = "google",
           pal = c("#999999","#E69F00","#56B4E9"),
           gl = 0, point.o = 0.8,
           inf.method = "ci",
           xlab = "Adjusted butyrate tertiles",
           ylab = "OTU richness")
summary(aov(richness ~ adj_butyrate_tertile, data=all.meta))


  # Figure 2 ----
  # Panel A
pirateplot(formula = bmi ~ adj_butyrate_tertile,
           data = all.meta,
           pal = c("#999999","#E69F00","#56B4E9"),
           gl = 0, point.o = 0.8,
           inf.method = "ci",
           ylab = "bmi (kg/m2)",
           xlim = c(0.5,6.5))
pirateplot(formula = bmi ~ adj_richness_tertile,
           data = all.meta,
           pal = c("#999999","#E69F00","#56B4E9"),
           gl = 0, point.o = 0.8,
           inf.method = "ci",
           at = c(4,5,6),
           add = TRUE)
summary(aov(bmi ~ adj_butyrate_tertile, data=all.meta))
summary(aov(bmi ~ adj_richness_tertile, data=all.meta))

  # Panel B
pirateplot(formula = waist ~ adj_butyrate_tertile,
           data = all.meta,
           pal = c("#999999","#E69F00","#56B4E9"),
           gl = 0, point.o = 0.8,
           inf.method = "ci",
           ylab = "Waist circumference (cm)",
           xlim = c(0.5,6.5))
pirateplot(formula = waist ~ adj_richness_tertile,
           data = all.meta,
           pal = c("#999999","#E69F00","#56B4E9"),
           gl = 0, point.o = 0.8,
           inf.method = "ci",
           at = c(4,5,6),
           add = TRUE)
summary(aov(waist ~ adj_butyrate_tertile, data=all.meta))
summary(aov(waist ~ adj_richness_tertile, data=all.meta))

  # Panel C
pirateplot(formula = triglycerides ~ adj_butyrate_tertile,
           data = all.meta,
           pal = c("#999999","#E69F00","#56B4E9"),
           gl = 0, point.o = 0.8,
           inf.method = "ci",
           ylim = c(0,400),
           ylab = "Triglycerides (mg/dL)",
           xlim = c(0.5,6.5))
pirateplot(formula = triglycerides ~ adj_richness_tertile,
           data = all.meta,
           pal = c("#999999","#E69F00","#56B4E9"),
           gl = 0, point.o = 0.8,
           inf.method = "ci",
           ylim = c(0,400),
           at = c(4,5,6),
           add = TRUE)
summary(aov(triglycerides ~ adj_butyrate_tertile, data=all.meta))
summary(aov(triglycerides ~ adj_richness_tertile, data=all.meta))

  # Panel D
pirateplot(formula = insulin ~ adj_butyrate_tertile,
           data = all.meta,
           pal = c("#999999","#E69F00","#56B4E9"),
           gl = 0, point.o = 0.8,
           inf.method = "ci",
           ylab = "Insulin (uU/mL)",
           xlim = c(0.5,6.5))
pirateplot(formula = insulin ~ adj_richness_tertile,
           data = all.meta,
           pal = c("#999999","#E69F00","#56B4E9"),
           gl = 0, point.o = 0.8,
           inf.method = "ci",
           at = c(4,5,6),
           add = TRUE)
summary(aov(insulin ~ adj_butyrate_tertile, data=all.meta))
summary(aov(insulin ~ adj_richness_tertile, data=all.meta))

  # Panel E
pirateplot(formula = mean_bp ~ adj_butyrate_tertile,
           data = all.meta,
           pal = c("#999999","#E69F00","#56B4E9"),
           gl = 0, point.o = 0.8,
           inf.method = "ci",
           ylab = "Mean blood pressure (mm Hg)",
           xlim = c(0.5,6.5))
pirateplot(formula = mean_bp ~ adj_richness_tertile,
           data = all.meta,
           pal = c("#999999","#E69F00","#56B4E9"),
           gl = 0, point.o = 0.8,
           inf.method = "ci",
           at = c(4,5,6),
           add = TRUE)
summary(aov(mean_bp ~ adj_butyrate_tertile, data=all.meta))
summary(aov(mean_bp ~ adj_richness_tertile, data=all.meta))

  # Panel F
pirateplot(formula = hsCRP ~ adj_butyrate_tertile,
           data = all.meta,
           pal = c("#999999","#E69F00","#56B4E9"),
           gl = 0, point.o = 0.8,
           inf.method = "ci",
           ylim = c(0,17),
           ylab = "hs-CRP (mg/L)",
           xlim = c(0.5,6.5))
pirateplot(formula = hsCRP ~ adj_richness_tertile,
           data = all.meta,
           pal = c("#999999","#E69F00","#56B4E9"),
           gl = 0, point.o = 0.8,
           inf.method = "ci",
           ylim = c(0,17),
           at = c(4,5,6),
           add = TRUE)
summary(aov(hsCRP ~ adj_butyrate_tertile, data=all.meta))
summary(aov(hsCRP ~ adj_richness_tertile, data=all.meta))


  # Restricted cubic splines ----
ddist <- datadist(all.meta)
options(datadist='ddist')

    # Supplementary Figure S1 ----
fit <- ols(richness ~ rcs(adjusted_SCFA,3), data=all.meta, x=TRUE, y=TRUE)
anova(fit)
rich_SCFA = ggplot(Predict(fit)) +
  geom_point(data=all.meta, aes(x=adjusted_SCFA, y=richness)) +
  coord_cartesian(ylim = c(25,250)) +
  labs(x = "Adjusted total fecal SCFAs", y = "OTU richness")

fit <- ols(richness ~ rcs(adjusted_acetate,3), data=all.meta, x=TRUE, y=TRUE)
anova(fit)
rich_acet = ggplot(Predict(fit)) +
  geom_point(data=all.meta, aes(x=adjusted_acetate, y=richness)) +
  coord_cartesian(ylim = c(25,250)) +
  labs(x = "Adjusted fecal acetate", y = "OTU richness")

fit <- ols(richness ~ rcs(adjusted_propionate,3), data=all.meta, x=TRUE, y=TRUE)
anova(fit)
rich_prop = ggplot(Predict(fit)) +
  geom_point(data=all.meta, aes(x=adjusted_propionate, y=richness)) +
  coord_cartesian(ylim = c(25,250)) +
  labs(x = "Adjusted fecal propionate", y = "OTU richness")

fit <- ols(richness ~ rcs(adjusted_butyrate,3), data=all.meta, x=TRUE, y=TRUE)
anova(fit)
rich_buty = ggplot(Predict(fit)) +
  geom_point(data=all.meta, aes(x=adjusted_butyrate, y=richness)) +
  coord_cartesian(ylim = c(25,250)) +
  labs(x = "Adjusted fecal butyrate", y = "OTU richness")

plot_grid(rich_SCFA,rich_acet,rich_prop,rich_buty, nrow = 2, ncol = 2, labels="AUTO")


    # Supplementary Figure S2 ----
fit <- ols(bmi ~ rcs(adjusted_butyrate,3), data=all.meta, x=TRUE, y=TRUE)
anova(fit)
but_bmi = ggplot(Predict(fit)) +
  geom_point(data=all.meta, aes(x=adjusted_butyrate, y=bmi)) +
  coord_cartesian(ylim = c(18,48)) +
  labs(x = "Adjusted fecal butyrate", y = "BMI (kg/m2)")

fit <- ols(waist ~ rcs(adjusted_butyrate,3), data=all.meta, x=TRUE, y=TRUE)
anova(fit)
but_waist = ggplot(Predict(fit)) +
  geom_point(data=all.meta, aes(x=adjusted_butyrate, y=waist)) +
  coord_cartesian(ylim = c(65,140)) +
  labs(x = "Adjusted fecal butyrate", y = "Waist circumference (cm)")

fit <- ols(triglycerides ~ rcs(adjusted_butyrate,3), data=all.meta, x=TRUE, y=TRUE)
anova(fit)
but_triglycerides = ggplot(Predict(fit)) +
  geom_point(data=all.meta, aes(x=adjusted_butyrate, y=triglycerides)) +
  coord_cartesian(ylim = c(25,400)) +
  labs(x = "Adjusted fecal butyrate", y = "Triglycerides (mg/dL)")

fit <- ols(insulin ~ rcs(adjusted_butyrate,3), data=all.meta, x=TRUE, y=TRUE)
anova(fit)
but_insulin = ggplot(Predict(fit)) +
  geom_point(data=all.meta, aes(x=adjusted_butyrate, y=insulin)) +
  coord_cartesian(ylim = c(0,60)) +
  labs(x = "Adjusted fecal butyrate", y = "Insulin (uU/mL)")

fit <- ols(mean_bp ~ rcs(adjusted_butyrate,3), data=all.meta, x=TRUE, y=TRUE)
anova(fit)
but_mean_bp = ggplot(Predict(fit)) +
  geom_point(data=all.meta, aes(x=adjusted_butyrate, y=mean_bp)) +
  coord_cartesian(ylim = c(60,150)) +
  labs(x = "Adjusted fecal butyrate", y = "Mean arterial \npressure (mm Hg)")

fit <- ols(hsCRP ~ rcs(adjusted_butyrate,3), data=all.meta, x=TRUE, y=TRUE)
anova(fit)
but_hsCRP = ggplot(Predict(fit)) +
  geom_point(data=all.meta, aes(x=adjusted_butyrate, y=hsCRP)) +
  coord_cartesian(ylim = c(0,16)) +
  labs(x = "Adjusted fecal butyrate", y = "hs-CRP (mg/L)")

plot_grid(but_bmi,but_waist,but_triglycerides,but_insulin,but_mean_bp,but_hsCRP, nrow = 3, ncol = 2, labels="AUTO")


    # Supplementary Figure S3 ----
fit <- ols(bmi ~ rcs(adjusted_richness,3), data=all.meta, x=TRUE, y=TRUE)
anova(fit)
rich_bmi = ggplot(Predict(fit)) +
  geom_point(data=all.meta, aes(x=adjusted_richness, y=bmi)) +
  coord_cartesian(ylim = c(18,48)) +
  labs(x = "Adjusted OTU richness", y = "BMI (kg/m2)")

fit <- ols(waist ~ rcs(adjusted_richness,3), data=all.meta, x=TRUE, y=TRUE)
anova(fit)
rich_waist = ggplot(Predict(fit)) +
  geom_point(data=all.meta, aes(x=adjusted_richness, y=waist)) +
  coord_cartesian(ylim = c(65,140)) +
  labs(x = "Adjusted OTU richness", y = "Waist circumference (cm)")

fit <- ols(triglycerides ~ rcs(adjusted_richness,3), data=all.meta, x=TRUE, y=TRUE)
anova(fit)
rich_triglycerides = ggplot(Predict(fit)) +
  geom_point(data=all.meta, aes(x=adjusted_richness, y=triglycerides)) +
  coord_cartesian(ylim = c(25,400)) +
  labs(x = "Adjusted OTU richness", y = "Triglycerides (mg/dL)")

fit <- ols(insulin ~ rcs(adjusted_richness,3), data=all.meta, x=TRUE, y=TRUE)
anova(fit)
rich_insulin = ggplot(Predict(fit)) +
  geom_point(data=all.meta, aes(x=adjusted_richness, y=insulin)) +
  coord_cartesian(ylim = c(0,60)) +
  labs(x = "Adjusted OTU richness", y = "Insulin (uU/mL)")

fit <- ols(mean_bp ~ rcs(adjusted_richness,3), data=all.meta, x=TRUE, y=TRUE)
anova(fit)
rich_mean_bp = ggplot(Predict(fit)) +
  geom_point(data=all.meta, aes(x=adjusted_richness, y=mean_bp)) +
  coord_cartesian(ylim = c(60,150)) +
  labs(x = "Adjusted OTU richness", y = "Mean arterial \npressure (mm Hg)")

fit <- ols(hsCRP ~ rcs(adjusted_richness,3), data=all.meta, x=TRUE, y=TRUE)
anova(fit)
rich_hsCRP = ggplot(Predict(fit)) +
  geom_point(data=all.meta, aes(x=adjusted_richness, y=hsCRP)) +
  coord_cartesian(ylim = c(0,16)) +
  labs(x = "Adjusted OTU richness", y = "hs-CRP (mg/L)")

plot_grid(rich_bmi,rich_waist,rich_triglycerides,rich_insulin,rich_mean_bp,rich_hsCRP, nrow = 3, ncol = 2, labels="AUTO")


  # Figure 3 ----
# PA: to obtain this figure for the subset of participants that did not take medication or smoke (Supplementary Figure S5), change the all.meta parameter by clean.meta

# Quasipoisson GLM for the association between microbiota and SCFAs
# Implements a quasipoisson generalized linear model to evaluate the association between OTU abundance an a given factor
microbiota_glm = function(otu_abs_freq, otu_rel_freq, variable, qvalue_cutoff, median_otu_abund, taxonomy, trans = NA){
  # Loads OTU tables and select OTUs according to median_otu_abund
  median_abund = apply(otu_rel_freq, MARGIN = 1, FUN = median)
  abundant_otus = t(otu_abs_freq)[median_abund >= (median_otu_abund/100), ]
  # Transforms the factor to log, asin_sqrt or none
  if(trans == "asin"){
    y.trans = y.trans = asin(sqrt(variable/100))
  } else if (trans == "log") {
    y.trans = log(variable)
  } else {
    y.trans = variable
  }
  # Runs glm on each OTU
  lm_meta = apply(abundant_otus, MARGIN = 1, FUN = function(x) glm_otu = glm(x ~ y.trans, family = quasipoisson, maxit = 100))
  #
  # Creates output table with taxonomy, coefficients and p- and q-values
  coef_lm = lapply(X = lm_meta, FUN = function(x) x$coefficients[2])
  anova_lm_meta = lapply(X = lm_meta, FUN = function(x) Anova(x, type = "II"))
  otus_meta = sapply(anova_lm_meta, FUN = function(x) x$`Pr(>Chisq)`[1])
  fdr_corrected = qvalue(otus_meta)
  sig_otus_meta = subset(anova_lm_meta, fdr_corrected$qvalues <= qvalue_cutoff)
  sig_otus_coef = subset(coef_lm, fdr_corrected$qvalues <= qvalue_cutoff)
  result_lm = sapply(X = sig_otus_meta, FUN = function(x) x$`Pr(>Chisq)`[1])
  result_lm = data.frame(Tax = taxonomy[row.names(taxonomy) %in% names(result_lm),2], Coef = unlist(sig_otus_coef, use.names = F), p = result_lm, q = fdr_corrected$qvalues[fdr_corrected$qvalues <= qvalue_cutoff])
  result_lm
}

otus_acetate = microbiota_glm(otu_abs_freq = microbio.rare, otu_rel_freq = microbio.relative, taxonomy = microbio.taxonomy, variable = adjusted_acetate, qvalue_cutoff = 0.1, median_otu_abund = 0.001, trans = "NA")
otus_butyrate = microbiota_glm(otu_abs_freq = microbio.rare, otu_rel_freq = microbio.relative, taxonomy = microbio.taxonomy, variable = adjusted_butyrate, qvalue_cutoff = 0.1, median_otu_abund = 0.001, trans = "NA")
otus_propionate = microbiota_glm(otu_abs_freq = microbio.rare, otu_rel_freq = microbio.relative, taxonomy = microbio.taxonomy, variable = adjusted_propionate, qvalue_cutoff = 0.1, median_otu_abund = 0.001, trans = "NA")
otus_isobutyrate = microbiota_glm(otu_abs_freq = microbio.rare, otu_rel_freq = microbio.relative, taxonomy = microbio.taxonomy, variable = adjusted_isobutyrate, qvalue_cutoff = 0.1, median_otu_abund = 0.001, trans = "NA")

otus_acetate$rho = apply(microbio.rare[,rownames(otus_acetate)], MARGIN = 2, FUN = function(x) cor(x, adjusted_acetate, method = "s"))
otus_butyrate$rho = apply(microbio.rare[,rownames(otus_butyrate)], MARGIN = 2, FUN = function(x) cor(x, adjusted_butyrate, method = "s"))
otus_propionate$rho = apply(microbio.rare[,rownames(otus_propionate)], MARGIN = 2, FUN = function(x) cor(x, adjusted_propionate, method = "s"))
otus_isobutyrate$rho = apply(microbio.rare[,rownames(otus_isobutyrate)], MARGIN = 2, FUN = function(x) cor(x, adjusted_isobutyrate, method = "s"))

otus_scfa_names = unique(c(rownames(otus_acetate[otus_acetate$q < 0.05,]), rownames(otus_butyrate[otus_butyrate$q < 0.05,]), rownames(otus_propionate[otus_propionate$q < 0.05,]), rownames(otus_isobutyrate[otus_isobutyrate$q < 0.05,])))

otus_rho_acetate = apply(microbio.rare[,otus_scfa_names], MARGIN = 2, FUN = function(x) cor(x, adjusted_acetate, method = "s"))
otus_rho_butyrate = apply(microbio.rare[,otus_scfa_names], MARGIN = 2, FUN = function(x) cor(x, adjusted_butyrate, method = "s"))
otus_rho_propionate = apply(microbio.rare[,otus_scfa_names], MARGIN = 2, FUN = function(x) cor(x, adjusted_propionate, method = "s"))
otus_rho_isobutyrate = apply(microbio.rare[,otus_scfa_names], MARGIN = 2, FUN = function(x) cor(x, adjusted_isobutyrate, method = "s"))

otus_scfa_rho = data.frame(Acetate = otus_rho_acetate, Butyrate = otus_rho_butyrate, Propionate = otus_rho_propionate, Isobutyrate = otus_rho_isobutyrate)
scfa_abs_rho = apply(abs(otus_scfa_rho), MARGIN = 1, FUN = max)
otus_strong_rho = otus_scfa_rho[scfa_abs_rho>0.2,]

otus_scfa_qvalue = cbind(otus_acetate[otus_scfa_names,]$q, otus_butyrate[otus_scfa_names,]$q, otus_propionate[otus_scfa_names,]$q, otus_isobutyrate[otus_scfa_names,]$q)

otus_scfa_qvalue = apply(otus_scfa_qvalue, MARGIN = 2, function(x) ifelse(x < 0.05, yes = "*", no = NA))
otus_strong_qvalue = otus_scfa_qvalue[scfa_abs_rho>0.2,]

# OTU taxonomy for complete dataset
otus_strong_tax = c("Otu00003 | Gemmiger formicilis", "Otu00006 | Akkermansia muciniphila", "Otu00009 | Enterobacter hormaechei", "Otu00010 | Clostridium celatum", "Otu00011 | Methanobrevibacter", "Otu00012 | Oscillospira", "Otu00025 | Bacteroides", "Otu00029 | (Clostridiaceae) SMB53", "Otu00073 | Haemophilus parainfluenzae", "Otu00080 | Alistipes finegoldii", "Otu00149 | (Clostridiaceae) 02d06", "Otu00151 | Paenibacillus ginsengarvi", "Otu00422 | Christensenellaceae", "Otu00005 | Faecalibacterium prausnitzii", "Otu00014 | Roseburia faecis", "Otu00020 | Streptococcus", "Otu00022 | Streptococcus", "Otu00033 | Oscillospira", "Otu00037 | Blautia", "Otu00043 | Oscillospira", "Otu00050 | Bacteroides uniformis", "Otu00056 | Subdoligranulum variabile", "Otu00060 | 02d06", "Otu00147 | Oscillospira", "Otu00264 | Clostridium lavalense", "Otu00308 | Bacillus solfatarensis", "Otu00002 | Prevotella copri")

# Heatmap
aheatmap(otus_strong_rho , color = "-Spectral:100", scale = "none", breaks = NA, Rowv=T, Colv=NA, width=30, height=8, hclustfun = "ward", distfun = "euclidean", txt = otus_strong_qvalue, labRow = otus_strong_tax)


  # Table 2 ----
# Remove samples with NAs (only complete cases)
all.meta = all.meta[complete.cases(all.meta[,"LBP"]),]
incomplete_samples = c("MI_017_H", "MI_078_H", "MI_354_H")
incomplete_positions = c(17, 78, 347)
all.meta = all.meta[-incomplete_positions,]

# Variables adjusted for pertinent covariates (rerun adjusted models above: lines 272-295)
all.meta_adjusted = data.frame(adjusted_bmi, adjusted_body_fat, adjusted_waist,
                               adjusted_HDL, adjusted_LDL, adjusted_VLDL, adjusted_triglycerides,
                               adjusted_hsCRP, adjusted_glucose, adjusted_HbA1c, adjusted_insulin,
                               adjusted_HOMA_IR, adjusted_leptin, adjusted_adiponectin, adjusted_LBP,
                               adjusted_systolic_bp, adjusted_diastolic_bp, adjusted_mean_bp, adjusted_SCFA, adjusted_acetate,
                               adjusted_propionate, adjusted_butyrate, adjusted_isobutyrate, adjusted_richness)

table2 = corr.test(all.meta_adjusted, all.meta_adjusted[19:24], method = "pearson", adjust = "fdr")
write.excel(table2$r)
write.excel(table2$p) # the upper diagonal corresponds to FDR-adjusted p-values (q-values)



  # Table 3: prevalence ratios ----
# This table was constructed with complete cases (no NAs: N=431)
# Multivariable-adjusted risk of prevalence (from doi: 10.1093/ije/dyv137)

# PA1: to obtain this table for the other SCFAs (Supplementary Table S5) modify the butyrate_tertile parameter accordingly (SCFA_tertile, acetate_tertile, propionate_tertile)
# PA2: to obtain this table for the subset of participants that did not take medication or smoke (Supplementary Table S6), modify the all.meta parameter by clean.meta

# Remove samples with NAs (only complete cases): WATCHOUT! Don't do this if already done for calculating Table 2 above
# Verify the number of rows with nrow(all.meta)
nrow(all.meta) # Must be 431, otherwise run the code below:
#all.meta = all.meta[complete.cases(all.meta[,"LBP"]),]
#incomplete_samples = c("MI_017_H", "MI_078_H", "MI_354_H")
#incomplete_positions = c(17, 78, 347)
#all.meta = all.meta[-incomplete_positions,]

# Recalculate tertiles for SCFAs and microbiota alpha diversity with N=431 individuals
all.meta = within(all.meta, SCFA_tertile <- as.integer(cut(all.meta$SCFA, quantile(all.meta$SCFA, probs=0:3/3), include.lowest=TRUE)))
all.meta = within(all.meta, acetate_tertile <- as.integer(cut(all.meta$acetate, quantile(all.meta$acetate, probs=0:3/3), include.lowest=TRUE)))
all.meta = within(all.meta, propionate_tertile <- as.integer(cut(all.meta$propionate, quantile(all.meta$propionate, probs=0:3/3), include.lowest=TRUE)))
all.meta = within(all.meta, butyrate_tertile <- as.integer(cut(all.meta$butyrate, quantile(all.meta$butyrate, probs=0:3/3), include.lowest=TRUE)))
all.meta = within(all.meta, isobutyrate_tertile <- as.integer(cut(all.meta$isobutyrate, quantile(all.meta$isobutyrate, probs=0:3/3), include.lowest=TRUE)))
all.meta = within(all.meta, richness_tertile <- as.integer(cut(all.meta$richness, quantile(all.meta$richness, probs=0:3/3), include.lowest=TRUE)))


    # Upper panel (obesity) ----
      # Subset data using bmi (lean vs. obese)
all.meta$bmi_clasif = ifelse(all.meta$bmi<25, "lean", ifelse(all.meta$bmi>=25 & all.meta$bmi<30, "overweight", "obese"))
lean_ob = all.meta[all.meta$bmi_clasif!="overweight",]

      # Convert bmi into a boolean (0 = lean, 1 = obese)
lean_ob$bmi_bool[lean_ob$bmi<=25] <- 0
lean_ob$bmi_bool[lean_ob$bmi>=30] <- 1
lean_ob$bmi_bool<-as.integer(lean_ob$bmi_bool)

      # Per butyrate tertiles ----
        # Count occurrences per levels of butyrate:
lean_ob %>% 
  group_by(butyrate_tertile) %>%
  summarise(no_rows = length(butyrate_tertile))

        # Fit a log binomial model (unadjusted model)
unadj_model = glm(bmi_bool~as.factor(butyrate_tertile), family=binomial(link="log"), data=lean_ob)
        # Risk ratio
exp(coef(unadj_model))
vcovHC(unadj_model, type="HC0") # Robust variance-covariance matrix (https://stat.ethz.ch/pipermail/r-help/2008-May/161610.html)

        # Calculate denominators used in inverse probability weights (PA: if the exposure variable is binary, use glm with family=binomial(link=logit); if it has more than two levels, use multinom(E~Z1+Z2) or ordered logistic regression [polr])
          # Confounder-adjusted model
E.out1 = multinom(butyrate_tertile ~ age + city + calories + physical_activity + fiber, data=lean_ob) # the exposure variable has more than 2 levels
ps1 = predict(E.out1, type = "probs")

            # Create stabilized weights, using a null model with E as the dependent variable
but_bool_ps1 = cbind(lean_ob$butyrate_tertile,ps1)
multinom_sptw <- function(x) {
  ifelse(x[,1]==1,mean(x[,1]==1,na.rm=T)/x[,2],ifelse(x[,1]==2,mean(x[,1]==2,na.rm=T)/x[,3],mean(x[,1]==3,na.rm=T)/x[,4]))
}
sptw1 = multinom_sptw(but_bool_ps1) # multinomial exposure variable

            # Fit a log binomial model (confounder-adjusted model)
adj_model_1 = glm(bmi_bool~as.factor(butyrate_tertile), family=binomial(link="log"), weight=sptw1, data=lean_ob)
            # Risk ratio
exp(coef(adj_model_1))
vcovHC(adj_model_1, type="HC0") # Robust variance-covariance matrix (https://stat.ethz.ch/pipermail/r-help/2008-May/161610.html)

          # Confounder-adjusted + LBP model
E.out2 = multinom(butyrate_tertile ~ age + city + calories + physical_activity + fiber + LBP, data=lean_ob) # the exposure variable has more than 2 levels
ps2 = predict(E.out2, type = "probs")

            # Create stabilized weights, using a null model with E as the dependent variable
but_bool_ps2 = cbind(lean_ob$butyrate_tertile,ps2)
multinom_sptw <- function(x) {
  ifelse(x[,1]==1,mean(x[,1]==1,na.rm=T)/x[,2],ifelse(x[,1]==2,mean(x[,1]==2,na.rm=T)/x[,3],mean(x[,1]==3,na.rm=T)/x[,4]))
}
sptw2 = multinom_sptw(but_bool_ps2) # multinomial exposure variable

            # Fit a log binomial model (confounder-adjusted model)
adj_model_2 = glm(bmi_bool~as.factor(butyrate_tertile), family=binomial(link="log"), weight=sptw2, data=lean_ob)
            # Risk ratio
exp(coef(adj_model_2))
vcovHC(adj_model_2, type="HC0") # Robust variance-covariance matrix (https://stat.ethz.ch/pipermail/r-help/2008-May/161610.html)


      # Per OTU richness tertiles ----
        # Count occurrences per levels of OTU richness:
lean_ob %>% 
  group_by(richness_tertile) %>%
  summarise(no_rows = length(richness_tertile))

        # Fit a log binomial model (unadjusted model)
unadj_model = glm(bmi_bool~as.factor(richness_tertile), family=binomial(link="log"), data=lean_ob)
        # Risk ratio
exp(coef(unadj_model))
vcovHC(unadj_model, type="HC0") # Robust variance-covariance matrix (https://stat.ethz.ch/pipermail/r-help/2008-May/161610.html)

        # Calculate denominators used in inverse probability weights (PA: if the exposure variable is binary, use glm with family=binomial(link=logit); if it has more than two levels, use multinom(E~Z1+Z2) or ordered logistic regression [polr])
          # Confounder-adjusted model
E.out1 = multinom(richness_tertile ~ age + city + calories + physical_activity + fiber, data=lean_ob) # the exposure variable has more than 2 levels
ps1 = predict(E.out1, type = "probs")

            # Create stabilized weights, using a null model with E as the dependent variable
but_bool_ps1 = cbind(lean_ob$richness_tertile,ps1)
multinom_sptw <- function(x) {
  ifelse(x[,1]==1,mean(x[,1]==1,na.rm=T)/x[,2],ifelse(x[,1]==2,mean(x[,1]==2,na.rm=T)/x[,3],mean(x[,1]==3,na.rm=T)/x[,4]))
}
sptw1 = multinom_sptw(but_bool_ps1) # multinomial exposure variable

            # Fit a log binomial model (confounder-adjusted model)
adj_model_1 = glm(bmi_bool~as.factor(richness_tertile), family=binomial(link="log"), weight=sptw1, data=lean_ob)
            # Risk ratio
exp(coef(adj_model_1))
vcovHC(adj_model_1, type="HC0") # Robust variance-covariance matrix (https://stat.ethz.ch/pipermail/r-help/2008-May/161610.html)

          # Confounder-adjusted + LBP model
E.out2 = multinom(richness_tertile ~ age + city + calories + physical_activity + fiber + LBP, data=lean_ob) # the exposure variable has more than 2 levels
ps2 = predict(E.out2, type = "probs")

            # Create stabilized weights, using a null model with E as the dependent variable
but_bool_ps2 = cbind(lean_ob$richness_tertile,ps2)
multinom_sptw <- function(x) {
  ifelse(x[,1]==1,mean(x[,1]==1,na.rm=T)/x[,2],ifelse(x[,1]==2,mean(x[,1]==2,na.rm=T)/x[,3],mean(x[,1]==3,na.rm=T)/x[,4]))
}
sptw2 = multinom_sptw(but_bool_ps2) # multinomial exposure variable

            # Fit a log binomial model (confounder-adjusted model)
adj_model_2 = glm(bmi_bool~as.factor(richness_tertile), family=binomial(link="log"), weight=sptw2, data=lean_ob)
            # Risk ratio
exp(coef(adj_model_2))
vcovHC(adj_model_2, type="HC0") # Robust variance-covariance matrix (https://stat.ethz.ch/pipermail/r-help/2008-May/161610.html)


    # Middle panel (central obesity) ----
    # Analysis using waist circumference
all.meta$waist_bool = ifelse(all.meta$sex=="male" & all.meta$waist<102, 0, ifelse(all.meta$sex=="female" & all.meta$waist<88, 0, 1))
lean_ob = all.meta

      # Per butyrate tertiles ----
      # Count occurrences per levels or butyrate:
lean_ob %>% 
  group_by(butyrate_tertile) %>%
  summarise(no_rows = length(butyrate_tertile))

      # Fit a log binomial model (unadjusted model)
unadj_model = glm(waist_bool~as.factor(butyrate_tertile), family=binomial(link="log"), data=lean_ob)
      # Risk ratio
exp(coef(unadj_model))
vcovHC(unadj_model, type="HC0") # Robust variance-covariance matrix (https://stat.ethz.ch/pipermail/r-help/2008-May/161610.html)

      # Calculate denominators used in inverse probability weights (PA: if the exposure variable is binary, use glm with family=binomial(link=logit); if it has more than two levels, use multinom(E~Z1+Z2) or ordered logistic regression [polr])
        # Confounder-adjusted model
E.out1 = multinom(butyrate_tertile ~ age + city + calories + physical_activity + fiber, data=lean_ob) # the exposure variable has more than 2 levels
ps1 = predict(E.out1, type = "probs")

          # Create stabilized weights, using a null model with E as the dependent variable
but_bool_ps1 = cbind(lean_ob$butyrate_tertile,ps1)
multinom_sptw <- function(x) {
  ifelse(x[,1]==1,mean(x[,1]==1,na.rm=T)/x[,2],ifelse(x[,1]==2,mean(x[,1]==2,na.rm=T)/x[,3],mean(x[,1]==3,na.rm=T)/x[,4]))
}
sptw1 = multinom_sptw(but_bool_ps1) # multinomial exposure variable

          # Fit a log binomial model (confounder-adjusted model)
adj_model_1 = glm(waist_bool~as.factor(butyrate_tertile), family=binomial(link="log"), weight=sptw1, data=lean_ob)
          # Risk ratio
exp(coef(adj_model_1))
vcovHC(adj_model_1, type="HC0") # Robust variance-covariance matrix (https://stat.ethz.ch/pipermail/r-help/2008-May/161610.html)

        # Confounder-adjusted + LBP model
E.out2 = multinom(butyrate_tertile ~ age + city + calories + physical_activity + fiber + LBP, data=lean_ob) # the exposure variable has more than 2 levels
ps2 = predict(E.out2, type = "probs")

          # Create stabilized weights, using a null model with E as the dependent variable
but_bool_ps2 = cbind(lean_ob$butyrate_tertile,ps2)
multinom_sptw <- function(x) {
  ifelse(x[,1]==1,mean(x[,1]==1,na.rm=T)/x[,2],ifelse(x[,1]==2,mean(x[,1]==2,na.rm=T)/x[,3],mean(x[,1]==3,na.rm=T)/x[,4]))
}
sptw2 = multinom_sptw(but_bool_ps2) # multinomial exposure variable

        # Fit a log binomial model (confounder-adjusted model)
adj_model_2 = glm(waist_bool~as.factor(butyrate_tertile), family=binomial(link="log"), weight=sptw2, data=lean_ob)
        # Risk ratio
exp(coef(adj_model_2))
vcovHC(adj_model_2, type="HC0") # Robust variance-covariance matrix (https://stat.ethz.ch/pipermail/r-help/2008-May/161610.html)


      # Per OTU richness tertiles ----
      # Count occurrences per levels or butyrate:
lean_ob %>% 
  group_by(richness_tertile) %>%
  summarise(no_rows = length(richness_tertile))

      # Fit a log binomial model (unadjusted model)
unadj_model = glm(waist_bool~as.factor(richness_tertile), family=binomial(link="log"), data=lean_ob)
      # Risk ratio
exp(coef(unadj_model))
vcovHC(unadj_model, type="HC0") # Robust variance-covariance matrix (https://stat.ethz.ch/pipermail/r-help/2008-May/161610.html)

      # Calculate denominators used in inverse probability weights (PA: if the exposure variable is binary, use glm with family=binomial(link=logit); if it has more than two levels, use multinom(E~Z1+Z2) or ordered logistic regression [polr])
        # Confounder-adjusted model
E.out1 = multinom(richness_tertile ~ age + city + calories + physical_activity + fiber, data=lean_ob) # the exposure variable has more than 2 levels
ps1 = predict(E.out1, type = "probs")

          # Create stabilized weights, using a null model with E as the dependent variable
but_bool_ps1 = cbind(lean_ob$richness_tertile,ps1)
multinom_sptw <- function(x) {
  ifelse(x[,1]==1,mean(x[,1]==1,na.rm=T)/x[,2],ifelse(x[,1]==2,mean(x[,1]==2,na.rm=T)/x[,3],mean(x[,1]==3,na.rm=T)/x[,4]))
}
sptw1 = multinom_sptw(but_bool_ps1) # multinomial exposure variable

          # Fit a log binomial model (confounder-adjusted model)
adj_model_1 = glm(waist_bool~as.factor(richness_tertile), family=binomial(link="log"), weight=sptw1, data=lean_ob)
          # Risk ratio
exp(coef(adj_model_1))
vcovHC(adj_model_1, type="HC0") # Robust variance-covariance matrix (https://stat.ethz.ch/pipermail/r-help/2008-May/161610.html)

        # Confounder-adjusted + LBP model
E.out2 = multinom(richness_tertile ~ age + city + calories + physical_activity + fiber + LBP, data=lean_ob) # the exposure variable has more than 2 levels
ps2 = predict(E.out2, type = "probs")

          # Create stabilized weights, using a null model with E as the dependent variable
but_bool_ps2 = cbind(lean_ob$richness_tertile,ps2)
multinom_sptw <- function(x) {
  ifelse(x[,1]==1,mean(x[,1]==1,na.rm=T)/x[,2],ifelse(x[,1]==2,mean(x[,1]==2,na.rm=T)/x[,3],mean(x[,1]==3,na.rm=T)/x[,4]))
}
sptw2 = multinom_sptw(but_bool_ps2) # multinomial exposure variable

          # Fit a log binomial model (confounder-adjusted model)
adj_model_2 = glm(waist_bool~as.factor(richness_tertile), family=binomial(link="log"), weight=sptw2, data=lean_ob)
          # Risk ratio
exp(coef(adj_model_2))
vcovHC(adj_model_2, type="HC0") # Robust variance-covariance matrix (https://stat.ethz.ch/pipermail/r-help/2008-May/161610.html)


    # Lower panel (hypertension) ----
    # Analysis using blood pressure
    # AHA 2017 Guideline for High Blood Pressure in Adults: Normal BP is defined as <120/<80 mm Hg; elevated BP 120-129/<80 mm Hg; hypertension stage 1 is 130-139 or 80-89 mm Hg, and hypertension stage 2 is >=140 or >=90 mm Hg
    # We also include in the hypertension category individuals with previous diagnosis of high blood pressure or taking antihypertensive medications
lean_ob = all.meta
lean_ob$HT_bool = ifelse(lean_ob$meds_hypertension==TRUE | lean_ob$systolic_bp>=130 | lean_ob$diastolic_bp>=80, 1, 0)

      # Per butyrate tertiles ----
      # Count occurrences per levels or butyrate:
lean_ob %>% 
  group_by(butyrate_tertile) %>%
  summarise(no_rows = length(butyrate_tertile))

      # Fit a log binomial model (unadjusted model)
unadj_model = glm(HT_bool~as.factor(butyrate_tertile), family=binomial(link="log"), data=lean_ob)
      # Risk ratio
exp(coef(unadj_model))
vcovHC(unadj_model, type="HC0") # Robust variance-covariance matrix (https://stat.ethz.ch/pipermail/r-help/2008-May/161610.html)

      # Calculate denominators used in inverse probability weights (PA: if the exposure variable is binary, use glm with family=binomial(link=logit); if it has more than two levels, use multinom(E~Z1+Z2) or ordered logistic regression [polr])
        # Confounder-adjusted model
E.out1 = multinom(butyrate_tertile ~ age + city + calories + physical_activity + fiber, data=lean_ob) # the exposure variable has more than 2 levels
ps1 = predict(E.out1, type = "probs")

        # Create stabilized weights, using a null model with E as the dependent variable
but_bool_ps1 = cbind(lean_ob$butyrate_tertile,ps1)
multinom_sptw <- function(x) {
  ifelse(x[,1]==1,mean(x[,1]==1,na.rm=T)/x[,2],ifelse(x[,1]==2,mean(x[,1]==2,na.rm=T)/x[,3],mean(x[,1]==3,na.rm=T)/x[,4]))
}
sptw1 = multinom_sptw(but_bool_ps1) # multinomial exposure variable

        # Fit a log binomial model (confounder-adjusted model)
adj_model_1 = glm(HT_bool~as.factor(butyrate_tertile), family=binomial(link="log"), weight=sptw1, data=lean_ob)
        # Risk ratio
exp(coef(adj_model_1))
vcovHC(adj_model_1, type="HC0") # Robust variance-covariance matrix (https://stat.ethz.ch/pipermail/r-help/2008-May/161610.html)

      # Confounder-adjusted + LBP model
E.out2 = multinom(butyrate_tertile ~ age + city + calories + physical_activity + fiber + LBP, data=lean_ob) # the exposure variable has more than 2 levels
ps2 = predict(E.out2, type = "probs")

        # Create stabilized weights, using a null model with E as the dependent variable
but_bool_ps2 = cbind(lean_ob$butyrate_tertile,ps2)
multinom_sptw <- function(x) {
  ifelse(x[,1]==1,mean(x[,1]==1,na.rm=T)/x[,2],ifelse(x[,1]==2,mean(x[,1]==2,na.rm=T)/x[,3],mean(x[,1]==3,na.rm=T)/x[,4]))
}
sptw2 = multinom_sptw(but_bool_ps2) # multinomial exposure variable

        # Fit a log binomial model (confounder-adjusted model)
adj_model_2 = glm(HT_bool~as.factor(butyrate_tertile), family=binomial(link="log"), weight=sptw2, data=lean_ob)
        # Risk ratio
exp(coef(adj_model_2))
vcovHC(adj_model_2, type="HC0") # Robust variance-covariance matrix (https://stat.ethz.ch/pipermail/r-help/2008-May/161610.html)


      # Per OTU richness tertiles ----
      # Count occurrences per levels or butyrate:
lean_ob %>% 
  group_by(richness_tertile) %>%
  summarise(no_rows = length(richness_tertile))

      # Fit a log binomial model (unadjusted model)
unadj_model = glm(HT_bool~as.factor(richness_tertile), family=binomial(link="log"), data=lean_ob)
      # Risk ratio
exp(coef(unadj_model))
vcovHC(unadj_model, type="HC0") # Robust variance-covariance matrix (https://stat.ethz.ch/pipermail/r-help/2008-May/161610.html)

      # Calculate denominators used in inverse probability weights (PA: if the exposure variable is binary, use glm with family=binomial(link=logit); if it has more than two levels, use multinom(E~Z1+Z2) or ordered logistic regression [polr])
        # Confounder-adjusted model
E.out1 = multinom(richness_tertile ~ age + city + calories + physical_activity + fiber, data=lean_ob) # the exposure variable has more than 2 levels
ps1 = predict(E.out1, type = "probs")

          # Create stabilized weights, using a null model with E as the dependent variable
but_bool_ps1 = cbind(lean_ob$richness_tertile,ps1)
multinom_sptw <- function(x) {
  ifelse(x[,1]==1,mean(x[,1]==1,na.rm=T)/x[,2],ifelse(x[,1]==2,mean(x[,1]==2,na.rm=T)/x[,3],mean(x[,1]==3,na.rm=T)/x[,4]))
}
sptw1 = multinom_sptw(but_bool_ps1) # multinomial exposure variable

          # Fit a log binomial model (confounder-adjusted model)
adj_model_1 = glm(HT_bool~as.factor(richness_tertile), family=binomial(link="log"), weight=sptw1, data=lean_ob)
          # Risk ratio
exp(coef(adj_model_1))
vcovHC(adj_model_1, type="HC0") # Robust variance-covariance matrix (https://stat.ethz.ch/pipermail/r-help/2008-May/161610.html)

        # Confounder-adjusted + LBP model
E.out2 = multinom(richness_tertile ~ age + city + calories + physical_activity + fiber + LBP, data=lean_ob) # the exposure variable has more than 2 levels
ps2 = predict(E.out2, type = "probs")

          # Create stabilized weights, using a null model with E as the dependent variable
but_bool_ps2 = cbind(lean_ob$richness_tertile,ps2)
multinom_sptw <- function(x) {
  ifelse(x[,1]==1,mean(x[,1]==1,na.rm=T)/x[,2],ifelse(x[,1]==2,mean(x[,1]==2,na.rm=T)/x[,3],mean(x[,1]==3,na.rm=T)/x[,4]))
}
sptw2 = multinom_sptw(but_bool_ps2) # multinomial exposure variable

          # Fit a log binomial model (confounder-adjusted model)
adj_model_2 = glm(HT_bool~as.factor(richness_tertile), family=binomial(link="log"), weight=sptw2, data=lean_ob)
          # Risk ratio
exp(coef(adj_model_2))
vcovHC(adj_model_2, type="HC0") # Robust variance-covariance matrix (https://stat.ethz.ch/pipermail/r-help/2008-May/161610.html)


  # Supplementary Table S3: SCFAs and stool consistency ----
# This table was constructed with multivariable-adjusted SCFA levels

# Calculating residuals from adjusted models
adjusted_SCFA = residuals(lm(log(0.1+SCFA)~age+city+calories+physical_activity+fiber, data = all.meta))
adjusted_acetate = residuals(lm(log(0.1+acetate)~age+city+calories+physical_activity+fiber, data = all.meta))
adjusted_propionate = residuals(lm(log(0.05+propionate)~age+city+calories+physical_activity+fiber, data = all.meta))
adjusted_butyrate = residuals(lm(log(0.01+butyrate)~age+city+calories+physical_activity+fiber, data = all.meta))
adjusted_isobutyrate = residuals(lm(log(0.001+isobutyrate)~age+city+calories+physical_activity+fiber, data = all.meta))

# Per total SCFA levels
aggregate(adjusted_SCFA ~ stool_consistency, FUN = mean, data = all.meta)
aggregate(adjusted_SCFA ~ stool_consistency, FUN = sem, data = all.meta)
summary(aov(adjusted_SCFA ~ stool_consistency, data = all.meta))

# Per acetate levels
aggregate(adjusted_acetate ~ stool_consistency, FUN = mean, data = all.meta)
aggregate(adjusted_acetate ~ stool_consistency, FUN = sem, data = all.meta)
summary(aov(adjusted_acetate ~ stool_consistency, data = all.meta))

# Per propionate levels
aggregate(adjusted_propionate ~ stool_consistency, FUN = mean, data = all.meta)
aggregate(adjusted_propionate ~ stool_consistency, FUN = sem, data = all.meta)
summary(aov(adjusted_propionate ~ stool_consistency, data = all.meta))

# Per butyrate levels
aggregate(adjusted_butyrate ~ stool_consistency, FUN = mean, data = all.meta)
aggregate(adjusted_butyrate ~ stool_consistency, FUN = sem, data = all.meta)
summary(aov(adjusted_butyrate ~ stool_consistency, data = all.meta))

# Per isobutyrate levels
aggregate(adjusted_isobutyrate ~ stool_consistency, FUN = mean, data = all.meta)
aggregate(adjusted_isobutyrate ~ stool_consistency, FUN = sem, data = all.meta)
summary(aov(adjusted_isobutyrate ~ stool_consistency, data = all.meta))


  # Supplementary Figure S4: PCoA with SCFA levels ----
# Calculate UniFrac distances
unifracs <- GUniFrac(microbio.rare, microbio.tree, alpha=c(0, 0.5, 1))$unifracs
dw <- unifracs[, , "d_1"]   # Weighted UniFrac
dw.dist = as.dist(dw)

du <- unifracs[, , "d_UW"]  # Unweighted UniFrac    
du.dist = as.dist(du)

# PCoA plots
e.pcoa.dw = cmdscale(dw.dist, k=5, eig = T)
e.PC1.dw = round(e.pcoa.dw$eig[1]/sum(e.pcoa.dw$eig), 4)* 100
e.PC2.dw = round(e.pcoa.dw$eig[2]/sum(e.pcoa.dw$eig), 4)* 100

e.pcoa.du = cmdscale(du.dist, k=5, eig = T)
e.PC1.du = round(e.pcoa.du$eig[1]/sum(e.pcoa.du$eig), 4)* 100
e.PC2.du = round(e.pcoa.du$eig[2]/sum(e.pcoa.du$eig), 4)* 100

# Adjusted SCFA levels
pcoa_table = data.frame(PC1.du = e.pcoa.du$points[, 1], PC2.du = e.pcoa.du$points[, 2], PC1.dw = e.pcoa.dw$points[, 1], PC2.dw = e.pcoa.dw$points[, 2], adjusted_SCFA = adjusted_SCFA, adjusted_acetate = adjusted_acetate, adjusted_propionate = adjusted_propionate, adjusted_butyrate = adjusted_butyrate, adjusted_isobutyrate = adjusted_isobutyrate, adjusted_bmi = adjusted_bmi, adjusted_waist = adjusted_waist, adjusted_diastolic_bp = adjusted_diastolic_bp, adj_SCFA_tertile = as.factor(all.meta$adj_SCFA_tertile), adj_acetate_tertile = as.factor(all.meta$adj_acetate_tertile), adj_propionate_tertile = as.factor(all.meta$adj_propionate_tertile), adj_butyrate_tertile = as.factor(all.meta$adj_butyrate_tertile), adj_isobutyrate_tertile = as.factor(all.meta$adj_isobutyrate_tertile))

scfa_tertiles_plot = ggplot(pcoa_table) +
  geom_point(aes(x=PC1.dw, y=PC2.dw, color=adj_SCFA_tertile)) +
  scale_size_area() +
  scale_colour_manual(values=cbPalette) +
  labs(colour="Adjusted total \nSCFA tertiles", x=paste("PCoA1",e.PC1.dw,"%"), y=paste("PCoA2",e.PC2.dw,"%"))

acetate_tertiles_plot = ggplot(pcoa_table) +
  geom_point(aes(x=PC1.dw, y=PC2.dw, color=adj_acetate_tertile)) +
  scale_size_area() +
  scale_colour_manual(values=cbPalette) +
  labs(colour="Adjusted \nacetate tertiles", x=paste("PCoA1",e.PC1.dw,"%"), y=paste("PCoA2",e.PC2.dw,"%"))

propionate_tertiles_plot = ggplot(pcoa_table) +
  geom_point(aes(x=PC1.dw, y=PC2.dw, color=adj_propionate_tertile)) +
  scale_size_area() +
  scale_colour_manual(values=cbPalette) +
  labs(colour="Adjusted \npropionate tertiles", x=paste("PCoA1",e.PC1.dw,"%"), y=paste("PCoA2",e.PC2.dw,"%"))

butyrate_tertiles_plot = ggplot(pcoa_table) +
  geom_point(aes(x=PC1.dw, y=PC2.dw, color=adj_butyrate_tertile)) +
  scale_size_area() +
  scale_colour_manual(values=cbPalette) +
  labs(colour="Adjusted \nbutyrate tertiles", x=paste("PCoA1",e.PC1.dw,"%"), y=paste("PCoA2",e.PC2.dw,"%"))

plot_grid(scfa_tertiles_plot, acetate_tertiles_plot, propionate_tertiles_plot, butyrate_tertiles_plot, nrow = 2, ncol = 2, labels="AUTO")


  # Supplementary Table S4: procrustes tests ----
# Procrustes tests with 10,000 permutations

# Unweighted unifrac
protest(du, pcoa_table$adjusted_SCFA, scale = TRUE, permutations = how(nperm = 10000))
protest(du, pcoa_table$adjusted_acetate, scale = TRUE, permutations = how(nperm = 10000))
protest(du, pcoa_table$adjusted_propionate, scale = TRUE, permutations = how(nperm = 10000))
protest(du, pcoa_table$adjusted_butyrate, scale = TRUE, permutations = how(nperm = 10000))
protest(du, pcoa_table$adjusted_isobutyrate, scale = TRUE, permutations = how(nperm = 10000))
protest(du, pcoa_table$adjusted_bmi, scale = TRUE, permutations = how(nperm = 10000))
protest(du, pcoa_table$adjusted_waist, scale = TRUE, permutations = how(nperm = 10000))
protest(du, pcoa_table$adjusted_diastolic_bp, scale = TRUE, permutations = how(nperm = 10000))

# Weighted unifrac
protest(dw, pcoa_table$adjusted_SCFA, scale = TRUE, permutations = how(nperm = 10000))
protest(dw, pcoa_table$adjusted_acetate, scale = TRUE, permutations = how(nperm = 10000))
protest(dw, pcoa_table$adjusted_propionate, scale = TRUE, permutations = how(nperm = 10000))
protest(dw, pcoa_table$adjusted_butyrate, scale = TRUE, permutations = how(nperm = 10000))
protest(dw, pcoa_table$adjusted_isobutyrate, scale = TRUE, permutations = how(nperm = 10000))
protest(dw, pcoa_table$adjusted_bmi, scale = TRUE, permutations = how(nperm = 10000))
protest(dw, pcoa_table$adjusted_waist, scale = TRUE, permutations = how(nperm = 10000))
protest(dw, pcoa_table$adjusted_diastolic_bp, scale = TRUE, permutations = how(nperm = 10000))
