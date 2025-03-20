# SETUP ####
library(googledrive)
library(googlesheets4)
library(tidyverse)
library(ggtext)
library(janitor)

# log in w/ zahn.ecology@gmail.com authentication
googledrive::drive_auth() 
3

theme_set(theme_bw() +
            theme(axis.text.y = element_text(face='bold',size=10),
                  legend.position = 'none',
                  strip.text = element_text(face='bold',size=14,color='white'),
                  strip.background = element_rect(fill = "gray30"),
                  axis.title = element_blank(),
                  axis.text.x = element_text(face='bold',size=12)))



# get data from google drive
lab_member_info <- "https://docs.google.com/spreadsheets/d/1y7mXtOH6TRtMFd87aa7lJRoqWm6fggFb-uhM81LCWrg/edit?usp=sharing"
dat <- googledrive::drive_get("Lab_Member_Info") %>% 
  googlesheets4::read_sheet() %>% 
  clean_names() %>% 
  mutate(across(ends_with("date"),as.Date))
3


# if no end date, set to current date
dat$end_date[is.na(dat$end_date)] <- Sys.Date() %>% as.Date()


dat %>% 
ggplot(aes(x=start_date,xend=end_date,
           y=name,yend=name,color=role)) +
  geom_segment(linewidth=10) +
  geom_vline(xintercept = as.Date(Sys.Date()), linetype=2,color='gray10') +
  scale_color_viridis_d(option = 'turbo',begin = .2) +
  scale_x_date(breaks = 'month') +
  theme(axis.text.x = element_text(angle = 90,hjust=0,vjust=.5))
ggsave("./assets/images/lab_member_gantt_chart.png",dpi=300,height = 3, width = 6)
