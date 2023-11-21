# SETUP ####
library(googledrive)
library(googlesheets4)
library(tidyverse)
library(ggtext)

theme_set(theme_bw() +
            theme(axis.text.y = element_text(face='bold',size=10),
                  legend.position = 'none',
                  strip.text = element_text(face='bold',size=14,color='white'),
                  strip.background = element_rect(fill = "gray30"),
                  axis.title = element_blank(),
                  axis.text.x = element_text(face='bold',size=12)))

plot_gantt <- function(x){
  x %>% 
    ggplot(aes(x=start_date,xend=end_date,
               y=task,yend=task,color=highlight)) +
    geom_segment(linewidth=5) +
    geom_vline(xintercept = Sys.Date(), linetype=2,color='gray10') +
    facet_wrap(~project,ncol = 1) +
    scale_y_discrete(position = 'right') +
    scale_color_manual(values=c("gray30","green4"))
}

# LOAD DATA FROM GOOGLE DRIVE ####
dat <- googledrive::drive_get("ResearchProjects") %>% 
  googlesheets4::read_sheet() %>% 
  mutate(start_date = as.Date(start_date),
         end_date = as.Date(end_date),
         task = factor(task, levels=rev(c("Planning","Securing funding","Data collection","Analysis","Writing",
                                      "Submission","Revisions","Completed"))),
         highlight = case_when(current_stage == task ~ TRUE, TRUE ~ FALSE))

dat$project <- factor(dat$project,levels = dat %>% 
                        select(project,next_step_target) %>% 
                        unique.data.frame() %>% arrange(next_step_target) %>% 
                        pluck('project'))


# GANTT CHARTS ####
full_gantt <- plot_gantt(dat)
current_gantt <- plot_gantt(dat %>% dplyr::filter(current_stage != "Completed" & !is.na(current_stage)))

current_gantt$data

# export objects
saveRDS(full_gantt,"./output/full_gantt.RDS")
saveRDS(current_gantt,"./output/current_gantt.RDS")


# TO-DO TASKS ####
current_projects <- 
dat %>% 
  dplyr::filter(current_stage != "Completed") %>% 
  dplyr::select(project,current_stage,lead_author,to_do) %>% 
  unique.data.frame()
row.names(current_projects) <- current_projects$project

current_projects[levels(forcats::fct_drop(current_projects$project)),] %>% 
  saveRDS("./output/current_projects.RDS")

