library(tidyverse)

# Import and clean ####
df <- read_csv("~/Downloads/2021-03-06T1135_Grades-BIOL_1620_002_ _Spring_2021.csv")
df <- df[2:nrow(df),]


# Your-Choice totals ####
YC_Names <- names(df) %>% 
  grep(pattern="YOUR-CHOICE|First",value = TRUE) %>% 
  grep(pattern = "Midterm",invert = TRUE,value = TRUE) %>% 
  grep(pattern = "Final",invert = TRUE,value = TRUE)


yc <- df %>% 
  filter(Student != "Student, Test") %>% 
  select(contains(YC_Names)) 

yc[is.na(yc)] <- 0

yc$Sum <- rowSums(yc)

yc <- yc %>% 
  mutate(BarColor = case_when(Sum >= 50 ~ "Green",
                              Sum < 50 ~ "Red"))


# Midterm totals ####
yc %>% 
  arrange(Sum) %>% 
  ggplot(aes(x=1:nrow(yc),y=Sum,fill=BarColor)) +
  geom_col() +
  geom_hline(yintercept = 50,linetype=2) +
  annotate("text",label="Midterm requirement",x = (nrow(yc)/2)-10,y=55) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x="Student",y="Current Your-Choice Total") +
  scale_fill_manual(values = c("DarkGreen","Red"))

# Assign midterm totals to column for Canvas upload
df$`YOUR-CHOICE -- Midterm (5290555)`[df$Student != "Student, Test"] <- yc$Sum


# Exam 2 distribution
df$`Exam 2 (5290748)`
df %>% 
  ggplot(aes(x=`Exam 2 (5290748)`)) +
  geom_density(fill="DarkGreen") +
  lims(x=c(0,102)) +
  theme_bw() +
  # theme(axis.text.x = element_blank(),
  #       axis.title.x = element_blank()) +
  labs(x="Exam 2 scores")

# plot yc vs midterm score
ggplot(df,aes(x=`YOUR-CHOICE -- Midterm (5290555)`,y=`Exam 2 (5290748)`)) +
  geom_point(alpha=.5) +
  geom_smooth(method = "lm",formula = y~poly(x,3)) +
  theme_minimal() +
  labs(y="Exam 2 score",x="Total YC Points at midterm")

mod <- glm(data=df,
    formula = `Exam 2 (5290748)` ~ `YOUR-CHOICE -- Midterm (5290555)`)
summary(mod)

quizzes <- df %>% 
  select(contains("Quiz "))
quizzes[is.na(quizzes)] <- 0
rowSums(quizzes)

mod2 <- glm(df$`Exam 2 (5290748)` ~ rowSums(quizzes) + df$`YOUR-CHOICE -- Midterm (5290555)`)
summary(mod2)
plot(df$`Exam 2 (5290748)` ~ rowSums(quizzes))
ggplot(df, aes(x=rowSums(quizzes),y=`Exam 2 (5290748)`)) +
  geom_point() + geom_smooth(method=lm) + theme_minimal() +
  labs(x="Current quiz total",y="Exam 2 score")

summary(mod2)
# export
write.csv(df, "~/Desktop/1620_midterm_YC.csv",row.names = FALSE,col.names = TRUE)

# Final grades ####

# Find totals

# Add bonus YC points

# Add bonus extra points










# Define point ranges for letter grades
a.cutoff = 700 
bplus.cutoff = c(690,699)
b.cutoff = c(640,689)
cplus.cutoff = c(630,639)
c.cutoff = c(560,629)
dplus.cutoff = c(550,559)
d.cutoff = c(480,549)


final.points <- as.numeric(df$`Final Points`)

df %>% 
  mutate(Letter.Grade = case_when(final.points >= a.cutoff ~ "A",
                                  final.points >= bplus.cutoff[1] & final.points <= bplus.cutoff[2] ~ "B+",
                                  final.points >= b.cutoff[1] & final.points <= b.cutoff[2] ~ "B",
                                  final.points >= cplus.cutoff[1] & final.points <= cplus.cutoff[2] ~ "C+",
                                  final.points >= c.cutoff[1] & final.points <= c.cutoff[2] ~ "C",
                                  final.points >= dplus.cutoff[1] & final.points <= dplus.cutoff[2] ~ "D+",
                                  final.points >= d.cutoff[1] & final.points <= d.cutoff[2] ~ "D",
                                  final.points < d.cutoff[1] ~ "E"))

