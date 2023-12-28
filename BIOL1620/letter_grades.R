library(tidyverse)

gradebook <- "~/Downloads/2023-12-11T0922_Grades-BIOL_1620_002_ _2023_Fall_-_Full_Term.csv"

df <- read_csv(gradebook)
names(df)

# add up point categories
bonus <- 
  df %>% 
  select(starts_with("BONUS - "),
         starts_with("Pre-Course Survey "),
         starts_with("Post-Course Survey ")) %>% 
  rowSums(na.rm = TRUE)

exams <- 
  df %>% 
  select(starts_with("Exam "),
         contains("Final Exam ")) %>% 
  rowSums(na.rm = TRUE)

assignments <- 
  df %>% 
  select(starts_with("YOUR-CHOICE -- Final ")) %>% 
  rowSums(na.rm = TRUE)

quizzes <- 
  df %>% 
  select(contains(": Quiz "),
         contains("Pre-Test "),
         contains("Post-Test ")) %>% 
  rowSums(na.rm=TRUE)


# add together
total <- assignments + bonus + exams + quizzes


# add to gradebook
df$`Final Points` <- df$`Final Points` %>% as.numeric %>% floor

df$Student
# calculate grades
final_grades <- 
df %>% 
  mutate(lettergrade = case_when(`Final Points` >= 600 ~ "A",
                   `Final Points` >= 590 & `Final Points` < 600 ~ "B+",
                   `Final Points` >= 540 & `Final Points` < 590 ~ "B",
                   `Final Points` >= 530 & `Final Points` < 540 ~ "C+",
                   `Final Points` >= 460 & `Final Points` < 530 ~ "C",
                   `Final Points` >= 450 & `Final Points` < 460 ~ "D+",
                   `Final Points` >= 380 & `Final Points` < 450 ~ "D",
                   `Final Points` < 380 ~ "E")) %>% 
  pluck('lettergrade')

# plot final points for students
df %>% 
  janitor::clean_names() %>% 
  mutate(letter = factor(final_grades,
                         levels = c("A","B+","B","C+","C","D+","D","E"))) %>% 
  select(final_points,student,letter) %>% 
  filter(student != "Points Possible") %>%
  filter(student != "Student, Test") %>% 
  arrange(desc(final_points)) %>% 
  mutate(id=seq_along(student)) %>% 
  ggplot(aes(x=id,y=final_points,fill=letter)) +
  geom_col() +
  geom_errorbar(aes(ymin=final_points,ymax=final_points+100),alpha=.5) +
  geom_hline(yintercept = 600,linetype=2) +
  theme_bw() +
  labs(x="student id",y="current total points\n before final exam",
       fill='letter grade\n(current)') +
  scale_fill_viridis_d(direction = -1,end=.85)

