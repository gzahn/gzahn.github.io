library(tidyverse)
library(janitor)
library(splines)

df <- read_delim("./binf-data-skills/Case_Studies/01_extract_ITS1_region/key/summary_output_per_file.tsv",
           delim = "\t")
names(df) <- make_clean_names(names(df))

df <- df %>% 
  mutate(elapsed_seconds = as.numeric(elapsed_time))
  
dfl <- pivot_longer(df,c(seqs_in,seqs_out),names_to = "in_out",values_to = "sequence_count",names_prefix = "seqs_")


ggplot(dfl,aes(x=sequence_count,y=elapsed_seconds,color=in_out)) +
  geom_point() +
  geom_smooth(method = "lm",formula = y ~ ns(x,2)) +
  theme_bw() +
  labs(x="Sequence count",
       y="Elapsed time (seconds)",
       color="Sequence source\n(in or out)")

ggsave("./binf-data-skills/Media/ITSxpress_speed_stats.png")

glm(data=dfl,
    formula = elapsed_seconds ~ sequence_count * in_out) %>% summary()

df2 <- read_delim("./binf-data-skills/Assignments/keys/Assignment_6/Full_Summary_Output.tsv",
                  delim="\t")

df2

ggplot(df2,aes(x=Mean_Quality,y=`Elapsed time`,color=`Seqs in`)) +
  geom_point() +
  geom_smooth()
