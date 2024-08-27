library(pbkrtest)
library(psych)
library(sjPlot)
library(readxl)
library(reghelper)
library(merTools)
library(lmerTest)
library(lme4)
library(ggeffects)
library(emmeans)
library(effects)
library(dbplyr)
library(readxl)
library(parameters)
# library(ggplot2)
library(multcomp)
library(effectsize)
library(lattice)
library(plotly)
library(tidyverse)
library(papaja)
library(gt)
library(vtable)
library(plotrix)
library(BayesFactor)



# Function to detect outliers
check_outliers <- function(data, column_name) {
  column <- data[[column_name]]
  
  # Calculate the quartiles and IQR
  Q1 <- quantile(column, 0.25, na.rm = TRUE)
  Q3 <- quantile(column, 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  
  # Calculate the lower and upper bounds for outliers
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  
  # Find the outliers
  outliers <- which(column < lower_bound | column > upper_bound)
  
  # Plot the boxplot
  boxplot(column, main = "Boxplot with Outliers", xlab = column_name)
  
  # Return the indices of outliers
  return(outliers)
}


# Manipulation Checks ---- 
HFSdata <- read_excel("Path to the data folder/HFS_Data.xlsx")
HFSdata$group[HFSdata$group == 0] = "Nocebo"
HFSdata$group[HFSdata$group == 1] = "Control"
HFSdata$manipulationcheck_expectation[HFSdata$manipulationcheck_expectation == 0] = "Less"
HFSdata$manipulationcheck_expectation[HFSdata$manipulationcheck_expectation == 1] = "Same"
HFSdata$manipulationcheck_expectation[HFSdata$manipulationcheck_expectation == 2] = "Higher"
HFSdata$manipulationcheck_sensitivitymatch[HFSdata$manipulationcheck_sensitivitymatch==0] = "Less"
HFSdata$manipulationcheck_sensitivitymatch[HFSdata$manipulationcheck_sensitivitymatch==1] = "Match"
HFSdata$manipulationcheck_sensitivitymatch[HFSdata$manipulationcheck_sensitivitymatch==2] = "Higher"
HFSdata = HFSdata[1:64,]

# remove exclusions
HFSdata <- HFSdata[-c(4, 12, 15, 49), ]

Expecthigh = HFSdata$participant[HFSdata$manipulationcheck_expectation=="Higher"]
Expectlow = HFSdata$participant[HFSdata$manipulationcheck_expectation=="Less" | HFSdata$manipulationcheck_expectation=="Same"]

# Mean and SD Age
HFSdata$Age <- as.numeric(HFSdata$Age)
mean(HFSdata$Age, na.rm = TRUE)
sd(HFSdata$Age, na.rm = TRUE)
t.test(HFSdata$Age[HFSdata$group == "Control"], HFSdata$Age[HFSdata$group == "Nocebo"],
       var.equal = TRUE)
sd(HFSdata$Age[HFSdata$group == "Nocebo"])

# For detection threshold
HFSdata$numeric_threshold <- sub("^(\\d+\\.\\d+).*", "\\1", HFSdata$detection_threshold)
HFSdata$numeric_threshold <- as.numeric(HFSdata$numeric_threshold)
t.test(HFSdata$numeric_threshold[HFSdata$group == "Control"], HFSdata$numeric_threshold[HFSdata$group == "Nocebo"],
       var.equal = TRUE)
sd(HFSdata$numeric_threshold[HFSdata$group == "Nocebo"])


# For stimulation intensity
HFSdata$StimInt <- sub(".* (\\d+\\.\\d+)", "\\1", HFSdata$detection_threshold)
HFSdata$StimInt <- as.numeric(HFSdata$StimInt)

t.test(HFSdata$StimInt[HFSdata$group == "Control"], HFSdata$StimInt[HFSdata$group == "Nocebo"],
       var.equal = TRUE)
sd(HFSdata$StimInt[HFSdata$group == "Nocebo"])



# Conduct the Chi-Squared test and Fisher's exact test for the Expectation check
manip <- data.frame(
  manipulationcheck_expectation = HFSdata$manipulationcheck_expectation,
  group = HFSdata$group
)
manip_clean <- subset(manip, manipulationcheck_expectation != "-" & group != "-")
print(manip_clean)

chisq_result <- chisq.test(table(manip_clean$manipulationcheck_expectation, manip_clean$group))
print(chisq_result)

# Get the observed counts
contingency_table <- table(manip_clean$manipulationcheck_expectation, manip_clean$group)
observed_counts <- as.data.frame(contingency_table)
colnames(observed_counts) <- c("Manipulation_Check_Expectation", "Group", "Count")
print(observed_counts)

# Fisher's exact test instead
fisher_result <- fisher.test(table(manip_clean$manipulationcheck_expectation, manip_clean$group))
print(fisher_result)

axsize = 20
mytheme = theme(
  axis.title.x = element_text(size = axsize),
  axis.text.x = element_text(size = axsize-6),
  axis.title.y = element_text(size = axsize),
  axis.text.y = element_text(size = axsize-2),
  plot.title = element_text(size=30,hjust = 0.5),
  legend.text=element_text(size=15),
  legend.title=element_text(size=20))
# Plot the grouped bar plot
ggplot(observed_counts, aes(x = Manipulation_Check_Expectation, y = Count, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "A",
       x = "Expectations of Pinprick Intensity",
       y = "Count") + theme_blank() + theme(legend.position="top") +
  scale_fill_manual(values = c("blue", "red")) + mytheme






# Conduct the Chi-Squared test and Fisher's exact test for the sensitivity match
manip <- data.frame(
  sensitivity_match = HFSdata$manipulationcheck_sensitivitymatch,
  group = HFSdata$group
)
manip_clean <- subset(manip, sensitivity_match != "-" & group != "-")
print(manip_clean)

chisq_result <- chisq.test(table(manip_clean$sensitivity_match, manip_clean$group))
print(chisq_result)

# Get the observed counts
contingency_table <- table(manip_clean$sensitivity_match, manip_clean$group)
observed_counts <- as.data.frame(contingency_table)
colnames(observed_counts) <- c("sensitivity_match", "Group", "Count")
print(observed_counts)

# Fisher's exact test instead
fisher_result <- fisher.test(table(manip_clean$sensitivity_match, manip_clean$group))
print(fisher_result)


axsize = 20
mytheme = theme(
  axis.title.x = element_text(size = axsize),
  axis.text.x = element_text(size = axsize-6),
  axis.title.y = element_text(size = axsize),
  axis.text.y = element_text(size = axsize-2),
  plot.title = element_text(size=30,hjust = 0.5),
  legend.text=element_text(size=15),
  legend.title=element_text(size=20))
# Plot the grouped bar plot
ggplot(observed_counts, aes(x = sensitivity_match, y = Count, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "B",
       x = "Sensitivity Match",
       y = "Count") + theme_blank() + theme(legend.position="top") +
  ylim(0,30) + 
  scale_fill_manual(values = c("blue", "red")) + mytheme





# HFS Ratings analysis ---- 
HFSdata <- read_excel("Path to data folder /HFS_Data_Stacked.xlsx")
HFSdata$group[HFSdata$group == 0] = "Nocebo"
HFSdata$group[HFSdata$group == 1] = "Control"

HFSdata$participant[is.na(HFSdata$Ratings)] # No participants with missing ratings

## Condition
HFS <- lmer(Ratings ~ Train + group + (1+Train|participant), HFSdata, REML = TRUE, contrasts = list(group = "contr.sum"))
summary(HFS)
anova(HFS)
plot(HFS)
qqmath(HFS, id=0.05)
HFSplot <- ggpredict(HFS,terms = c('Train','group')) 
plot(HFSplot)
t_to_eta2(1.945,40.0690)


byTrainmean <- aggregate(Ratings ~ Train + group, data = HFSdata, 
                         FUN = function(x) c(mean = mean(x)))
byTrainse <- aggregate(Ratings ~ Train + group, data = HFSdata, 
                       FUN = function(x) c(se = std(x)))
byTrainmean$se <- byTrainse$Ratings
MFSconditionplot <- ggplot(data = byTrainmean, aes(x = Train, y = Ratings,color = group)) + labs(title="Perceived Intensity of MFS",y="Intensity (NRS)",x = "Train")+
  geom_line() + geom_point() + theme_blank() + geom_errorbar(data = byTrainmean,aes(ymin = byTrainmean$Ratings - byTrainmean$se,ymax = byTrainmean$Ratings + byTrainmean$se,width = 0.2)) + 
  ylim(c(30, 110)) + scale_x_continuous(breaks = c(1:12)) + mytheme + theme(legend.position = c(0.8, 0.8))
MFSconditionplot + scale_color_manual(values=c("red","black")) +  theme(legend.position = 'top')



# Pinprick data analysis ---- 
Pinprickdata <- read_excel("Path to data folder /PinprickData.xlsx")

# Some cleanup
Pinprickdata$Group[Pinprickdata$Group == 1] = "Control"
Pinprickdata$Group[Pinprickdata$Group == 0] = "Nocebo"
Pinprickdata <- Pinprickdata %>%
  pivot_longer(cols = c("T0Ratings", "T1Ratings"), names_to = "Time", values_to = "Ratings")
Pinprickdata$Time[Pinprickdata$Time == "T0Ratings"] = "T0"
Pinprickdata$Time[Pinprickdata$Time == "T1Ratings"] = "T1"

Pinprickdata <- Pinprickdata[Pinprickdata$ParticipantID != 27,] # Exclude participant 27


Pinprickdata$Expect <- ""
Pinprickdata$Expect[Pinprickdata$ParticipantID %in% Expecthigh] <- "High"
Pinprickdata$Expect[Pinprickdata$ParticipantID %in% Expectlow] <- "Low"

# Plot by expectations
participant_means <- Pinprickdata %>%
  group_by(ParticipantID, Time, Expect) %>%
  summarise(MeanRating = mean(Ratings, na.rm = TRUE), .groups = 'drop')
overall_means_sem <- participant_means %>%
  group_by(Time, Expect) %>%
  summarise(OverallMean = mean(MeanRating, na.rm = TRUE),
            SEM = sd(MeanRating, na.rm = TRUE) / sqrt(n()), .groups = 'drop')
ggplot() +
  geom_line(data = participant_means, aes(x = Time, y = MeanRating, group = interaction(ParticipantID), color = Expect), alpha = 0.5) +
  geom_line(data = overall_means_sem, aes(x = Time, y = OverallMean, color = Expect, group = Expect), size = 1) +
  geom_point(data = overall_means_sem, aes(x = Time, y = OverallMean, color = Expect), size = 2) +
  geom_errorbar(data = overall_means_sem, aes(x = Time, y = OverallMean, ymin = OverallMean - SEM, ymax = OverallMean + SEM), width = 0.2) +
  facet_wrap(~ Expect, scales = "free_y") +
  labs(title = "Mean Intensity Ratings by Time and Group",
       x = "Time", y = "Intensity Rating", color = "Arm") +
  theme_minimal() +
  ylim(c(0, 100)) +  
  theme(plot.title = element_text(hjust = 0.5))


# Look at individual data
participantdat = Pinprickdata[Pinprickdata$ParticipantID == 5,]
ggplot() +
  geom_line(data = participantdat, aes(x = Time, y = Ratings), alpha = 0.5) +
  geom_point(data = participantdat, aes(x = Time, y = Ratings), size = 2) +
  theme_minimal() +
  ylim(c(0, 100)) +  
  theme(plot.title = element_text(hjust = 0.5))


# Problem participants - 27 (Exclude)

mean(Pinprickdata$Ratings[Pinprickdata$Time == 'T1'])
sd(Pinprickdata$Ratings[Pinprickdata$Time == 'T1'])

# Fit the model
IntMod.basic <- lmer(Ratings ~ Time*Group + (1|ParticipantID),Pinprickdata,REML = TRUE, contrasts = list(Time = "contr.sum"))
summary(IntMod.basic)
anova(IntMod.basic)
plot(IntMod.basic)
resid = residuals(IntMod.basic)
hist(resid)
shapiro.test(resid)
qqmath(IntMod.basic, id=0.05)
emmeans(IntMod.basic, pairwise ~ Time,lmer.df = "satterthwaite")
t_to_eta2(1.342,209.20541)
F_to_eta2(0.0308,1, 2338)

t_to_eta2(1.342,209.20541)

F_to_eta2(.259,1,52627)

# Calculate the Bayes factor
Pinprickdata2 = Pinprickdata
missing = is.na(Pinprickdata2$Ratings)
Pinprickdata2 = Pinprickdata2[!missing,]
Pinprickdata2$ParticipantID = as.factor(Pinprickdata2$ParticipantID)
full_BF = lmBF(Ratings ~ Time*Group + ParticipantID ,data=Pinprickdata2, whichRandom = 'ParticipantID')
null_BF = lmBF(Ratings ~  Time + Group + ParticipantID ,data=Pinprickdata2, whichRandom = 'ParticipantID')
null_BF/full_BF

Pinprickdata$ParticipantID[Pinprickdata$Ratings == 0]
Pinprickdata2 <- Pinprickdata[Pinprickdata$Ratings != 0, ]
hist(log(Pinprickdata$Ratings+.01))
hist(Pinprickdata$Ratings)


hist(Pinprickdata$Ratings[Pinprickdata$Time == "T1"])



# Remove outliers
# Check for outliers for each Time point separately
T0 = Pinprickdata[Pinprickdata$Time == 'T0',]
T1 = Pinprickdata[Pinprickdata$Time == 'T1',]

T0outliers <- check_outliers(T0, "Ratings")
print(T0outliers)
T0$ParticipantID[T0outliers]
T0$Ratings[T0outliers] = NA

T1outliers <- check_outliers(T1, "Ratings")
print(T1outliers)
T1$ParticipantID[T1outliers]
T1$Ratings[T1outliers] = NA

PinprickdataOR = rbind(T0,T1)


# Check for outliers for each group separately
Control = Pinprickdata[Pinprickdata$Group == 'Control',]
Nocebo = Pinprickdata[Pinprickdata$Group == 'Nocebo',]

Coutliers <- check_outliers(Control, "Ratings")
print(Coutliers)
Control$ParticipantID[Coutliers]
Control$Ratings[Coutliers] = NA

Noutliers <- check_outliers(Nocebo, "Ratings")
print(Noutliers)
Nocebo$ParticipantID[Noutliers]
Nocebo$Ratings[Noutliers] = NA

PinprickdataOR = rbind(Control,Nocebo)


IntMod.basic <- lmer(Ratings ~ Time*Group + (1|ParticipantID),PinprickdataOR,REML = TRUE, contrasts = list(Time = "contr.sum"))
summary(IntMod.basic)
anova(IntMod.basic)
plot(IntMod.basic)
emms1 <- emmeans(IntMod.basic, ~ Group | Time)
con1 <- contrast(emms1, interaction = "pairwise")
pairs(con1, by = 'Group')
hyperalgesia <- ggpredict(IntMod.basic, terms = c('Time','Group'))
plot(hyperalgesia,add.data = FALSE)


df_summary <- Pinprickdata %>%
  group_by(Time, Group) %>%
  summarize(Mean = mean(Ratings),
            SE = sd(Ratings) / sqrt(n()),
            .groups = 'drop')

# Calculate upper and lower bounds for error bars
df_summary <- df_summary %>%
  mutate(upper = Mean + SE, lower = Mean - SE)

# Plot
ggplot(df_summary, aes(x = Time, y = Mean, color = Group, group = Group)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  labs(x = "Time", y = "Intensity (NRS)", title = "Ratings by Group Over Time") +
  theme_blank()+
  ylim(c(0, 75)) +
  labs(title = "Ratings by Group Over Time", y = "Intensity (NRS)", x = "Time") +
  mytheme + theme(text = element_text(size = axsize)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("black", "red"))



Intplot2 <- ggplot(data = Pinprickdata, aes(x = Time, y = Ratings, colour = Group)) +
  geom_boxplot(aes(middle = mean(Ratings)), width = .15, position = position_dodge(1.1), fatten = NULL, outlier.shape = NA) +
  theme_blank() +
  ylim(c(0, 110)) + scale_y_continuous(limits = c(0, 100)) +
  geom_violin(aes(x = Time, y = Ratings, fill = Group), alpha = .01, width = 1.1, show.legend = FALSE,trim = TRUE) +
  labs(title = "", y = "Intensity (NRS)", x = "Time") +
  mytheme +
  # theme(legend.position="top") +
  theme(text = element_text(size = axsize)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_summary(aes(x = Time, y = Ratings, group = Group), fun = mean, geom = "point", width = .5, size = 1.5, color = "red", fill = "red", position = position_dodge(1.1)) +
  scale_color_manual(values = c("black", "red")) +
  geom_point(data = Pinprickdata,  aes(x = Time, y = Ratings, colour = Group), shape=16, position=position_jitterdodge(dodge.width=1.1,jitter.width = .2), alpha = .1) + 
  scale_fill_manual(values = c("black", "red"))
Intplot2



# Area analysis ---- 
HFSdata <- read_excel("Path to data folder /HFS_Data.xlsx")
HFSdata$group[HFSdata$group == 0] = "Nocebo"
HFSdata$group[HFSdata$group == 1] = "Control"
HFSdata = HFSdata[1:64,]
HFSdata <- HFSdata[-c(4, 12, 15, 49), ]
HFSdata$numeric_threshold <- sub("^(\\d+\\.\\d+).*", "\\1", HFSdata$detection_threshold)
HFSdata$numeric_threshold <- as.numeric(HFSdata$numeric_threshold)


HFSdata$area_sensitivity_length <- as.numeric(HFSdata$area_sensitivity_length)

# Run ttest for the area
t.test(HFSdata$area_sensitivity_length ~ HFSdata$group, var.equal = TRUE)
t_to_eta2(-0.86279,58)

result <- ttestBF(formula = area_sensitivity_length ~ group, data = HFSdata)
print(result)

mean(HFSdata$area_sensitivity_length[HFSdata$group == "Nocebo"], na.rm = TRUE)
sd(HFSdata$area_sensitivity_length[HFSdata$group == "Nocebo"], na.rm = TRUE)
mean(HFSdata$area_sensitivity_length[HFSdata$group == "Control"], na.rm = TRUE)
sd(HFSdata$area_sensitivity_length[HFSdata$group == "Control"], na.rm = TRUE)


# Check for outliers for each group separately
control = HFSdata[HFSdata$group == 'Control',]
nocebo = HFSdata[HFSdata$group == 'Nocebo',]

Coutliers <- check_outliers(control, "area_sensitivity_length")
print(Coutliers)
control$participant[Coutliers]
control$area_sensitivity_length[Coutliers] = NA

Noutliers <- check_outliers(nocebo, "area_sensitivity_length")
print(Noutliers)
nocebo$participant[Noutliers]
nocebo$area_sensitivity_length[Noutliers] = NA

HFSdataOR = rbind(control,nocebo)

t.test(HFSdataOR$area_sensitivity_length ~ HFSdataOR$group, var.equal = TRUE)



# Plotting
axsize = 20
mytheme = theme(
  axis.title.x = element_text(size = axsize),
  axis.text.x = element_text(size = axsize-6),
  axis.title.y = element_text(size = axsize),
  axis.text.y = element_text(size = axsize-2),
  plot.title = element_text(size=25,hjust = 0.5),
  legend.text=element_text(size=15),
  legend.title=element_text(size=20))
summary_stats <- HFSdata %>%
  group_by(group) %>%
  summarise(mean_area_sensitivity_length = mean(area_sensitivity_length),
            se_area_sensitivity_length = sd(area_sensitivity_length) / sqrt(n()))

# Plot the means and standard errors
ggplot(summary_stats, aes(x = group, y = mean_area_sensitivity_length)) +
  geom_bar(stat = "identity", fill = "gray", color = "black") +  # Plot means as bars
  geom_errorbar(aes(ymin = mean_area_sensitivity_length - se_area_sensitivity_length, 
                    ymax = mean_area_sensitivity_length + se_area_sensitivity_length),  # Plot error bars
                width = 0.2, color = "black", size = 0.7) + theme_blank() +
  labs(title = "Spread of Hypersensitivity",
       x = "Group", y = "Mean Spread (cm)") + ylim(0,15) +
  mytheme



Areaplot <- ggplot(data = HFSdata, aes(x = group, y = area_sensitivity_length)) +
  geom_boxplot(aes(middle = mean(area_sensitivity_length)), width = .15, position = position_dodge(1.1), fatten = NULL, outlier.shape = NA) +
  theme_blank() +
  ylim(c(0, 20)) + scale_y_continuous(limits = c(0, 20)) +
  geom_violin(aes(x = group, y = area_sensitivity_length), alpha = .01, width = 1.1, show.legend = FALSE,trim = TRUE) +
  labs(title = "", y = "Mean Spread (cm)", x = "Group") +
  mytheme +
  # theme(legend.position="top") +
  theme(text = element_text(size = axsize)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_summary(aes(x = group, y = area_sensitivity_length), fun = mean, geom = "point", width = .5, size = 1.5, color = "red", fill = "red", position = position_dodge(1.1)) +
  scale_color_manual(values = c("black", "red")) +
  scale_fill_manual(values = c("black", "red"))
Areaplot


HFSdata$group <- as.factor(HFSdata$group)
# With points
Areaplot <- ggplot(data = HFSdata, aes(x = group, y = area_sensitivity_length)) +
  geom_boxplot(aes(middle = mean(area_sensitivity_length)), width = .15, position = position_dodge(1.1), fatten = NULL, outlier.shape = NA) +
  theme_blank() +
  ylim(c(0, 20)) + scale_y_continuous(limits = c(0, 20)) +
  geom_violin(aes(x = group, y = area_sensitivity_length), alpha = .01, width = 1.1, show.legend = FALSE,trim = TRUE) +
  labs(title = "", y = "Mean Spread (cm)", x = "Group") +
  mytheme +
  theme(text = element_text(size = axsize)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_summary(aes(x = group, y = area_sensitivity_length), fun = mean, geom = "point", width = .5, size = 1.5, color = "red", fill = "red", position = position_dodge(1.1)) +
  # geom_jitter(alpha = 0.3, color = "black", position = position_dodge(1.1), size = 1, width = 0.2) + # Adding individual points
  # geom_point(data = HFSdata,aes(x = group,y = area_sensitivity_length),shape=16, position=position_jitterdodge(dodge.width=1), alpha = .3) + 
  geom_point(data = HFSdata,aes(x = group,y = area_sensitivity_length),shape=16, position=position_jitter(width = .1), alpha = .3) + 
  scale_color_manual(values = c("black", "red")) +
  scale_fill_manual(values = c("black", "red"))
Areaplot



# Analysis of the resting state average pupil diameter - from Pre-registered processing
PupilRS <- read_excel("Path to data folder /RS_vals.xlsx") # - from Pre-registered processing
PupilRS <- read_excel("Path to data folder /RS_vals_Ketan.xlsx") # - from Emanuel's processing
PupilRS$Group[PupilRS$Group == 0] = "Nocebo"
PupilRS$Group[PupilRS$Group == 1] = "Control"
PupilRS <- PupilRS %>%
  pivot_longer(cols = c("RS1", "RS2"), names_to = "Time", values_to = "Values")


RSmod <- lmer(Values ~ Time*Group + (1|Participant),PupilRS,REML = TRUE)
summary(RSmod)
anova(RSmod)
plot(RSmod)
hyperalgesia <- ggpredict(RSmod, terms = c('Time','Group'))
plot(hyperalgesia,add.data = FALSE)

F_to_eta2(1.5490,1,47)

# Calculate the Bayes factor
PupilRS2 = PupilRS
missing = is.na(PupilRS2$Values)
PupilRS2 = PupilRS2[!missing,]
PupilRS2$Participant = as.factor(PupilRS2$Participant)
full_BF = lmBF(Values ~ Time*Group + Participant ,data=PupilRS2, whichRandom = 'Participant')
null_BF = lmBF(Values ~  Time + Group + Participant ,data=PupilRS2, whichRandom = 'Participant')
null_BF/full_BF


mean(PupilRS$Values[PupilRS$Group == "Nocebo" & PupilRS$Time == "RS1"])
sd(PupilRS$Values[PupilRS$Group == "Nocebo" & PupilRS$Time == "RS1"])

# Get estimated marginal means and standard errors
emmeans <- emmeans(RSmod, specs = ~ Time*Group)

# Convert to data frame
emmeans_df <- as.data.frame(emmeans)

# Plot
axsize = 20
mytheme = theme(
  axis.title.x = element_text(size = axsize),
  axis.text.x = element_text(size = axsize-6),
  axis.title.y = element_text(size = axsize),
  axis.text.y = element_text(size = axsize-2),
  plot.title = element_text(size=25,hjust = 0.5),
  legend.text=element_text(size=15),
  legend.title=element_text(size=20))
ggplot(emmeans_df, aes(x = Time, y = emmean, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black", alpha = 0.85) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2, position = position_dodge(width = 0.9)) +
  labs(title = "Resting Pupil Diameter",
       x = "Time",
       y = "Mean Pupil Diameter (a.u.)") +
  ylim(0, 1500) +
  theme_blank() + mytheme +  # Assuming mytheme contains your theme settings
  scale_fill_manual(values = c("black", "red"))  # Set fill colors for groups



Pupilplot <- ggplot(data = PupilRS, aes(x = Time, y = Values, colour = Group)) +
  geom_boxplot(aes(middle = mean(Values)), width = .2, position = position_dodge(1.0), fatten = NULL, outlier.shape = NA) +
  theme_blank() +
  ylim(c(0, 1700)) + scale_y_continuous(limits = c(0, 1700)) +
  geom_violin(aes(x = Time, y = Values, fill = Group), alpha = .01, width = 1.0, show.legend = FALSE,trim = TRUE) +
  labs(title = "Resting Pupil Diameter", y = "Mean Pupil Diameter (a.u.)", x = "Time") +
  mytheme +
  theme(text = element_text(size = axsize)) +
  theme(plot.title = element_text(hjust = 0.5),legend.position="top") +
  stat_summary(aes(x = Time, y = Values, group = Group), fun = mean, geom = "point", width = .5, size = 1.5, color = "red", fill = "red", position = position_dodge(1)) +
  geom_point(data = PupilRS,  aes(x = Time, y = Values, colour = Group), shape=16, position=position_jitterdodge(dodge.width=1.1,jitter.width = .2), alpha = .3) + 
  scale_color_manual(values = c("black", "red")) +
  scale_fill_manual(values = c("black", "red"))
Pupilplot




