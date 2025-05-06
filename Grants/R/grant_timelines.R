# Build plot/table of grant timelines and alerts #

# SETUP ####
library(googledrive)
library(googlesheets4)
library(tidyverse)
library(ggtext)
options(scipen = 999)
# googledrive::drive_auth()

theme_set(theme_bw() +
            theme(axis.text.y = element_text(face='bold',size=10),
                  legend.position = 'right',
                  strip.text = element_text(face='bold',size=14,color='white'),
                  strip.background = element_rect(fill = "gray30"),
                  axis.title = element_blank(),
                  axis.text.x = element_text(face='bold',size=10,angle=90,hjust=1,vjust=.5)))

# LOAD DATA FROM GOOGLE DRIVE ####
# https://docs.google.com/spreadsheets/d/12I-JyaU1ROdT9-COqLNlGedOeWiHkWtr6iPSoTEdaTY/edit?gid=0#gid=0
dat <- googledrive::drive_get("Grant timelines") %>% 
  googlesheets4::read_sheet() %>% 
  mutate(due_date=as.Date(due_date))

# refactor 
newlevels <- 
  dat %>% 
  arrange(desc(due_date)) %>% 
  pluck("program")

dat <- 
  dat %>% 
  mutate(program = factor(program,levels=newlevels))

# weeks remaining
dat$weeks <- difftime(dat$due_date,Sys.Date(),units = 'weeks') %>% round(0) %>% as.character()


dat %>% 
  ggplot(aes(x=due_date,y=program,color=est_amount)) +
  geom_point(size=5) +
  geom_segment(aes(x=due_date,y=program,xend=due_date,yend=0),
               linetype=22) +
  scale_x_date(
    limits = as.Date(c(Sys.Date(), max(dat$due_date)+100)),
    date_breaks = "1 month",
    date_labels = "%b %Y"
  ) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
  scale_color_viridis_c(option = 'turbo',end=.8) +
  labs(color="Estimated amount") +
  geom_vline(xintercept = Sys.Date(),linetype=111) +
  geom_label(aes(x=due_date+100,y=program,label=paste(weeks,"weeks left")),color='black')
ggsave("./grant_due_dates_plot.png",dpi=300,height = nrow(dat),width = nrow(dat)*2)







