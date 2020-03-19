############## Gantt Chart for MRes Proposal ###############

# clear workspace
rm(list=ls())

# load packages
library(ggplot2)
require(dplyr)
library(forcats)

# create data frame of tasks and times
data <-data.frame(task = c('Probe analysis in lab',
                           'Literature review',
                           'Statistics and Modelling',
                           'Data Analysis',
                           'Final Submission Deadline',
                           'Data Exploration',
                           'Results and Final Write-Up',
                           'Introduction Write-Up',
                           'Methods Write-Up',
                           'HPC pipeline'),
                  start.date = as.Date(c("2020-03-02", 
                                         "2019-12-16", 
                                         "2020-02-03",
                                         "2020-03-02",
                                         "2020-08-26",
                                         "2019-12-30",
                                         "2020-06-22",
                                         "2019-12-30",
                                         "2020-03-16",
                                         "2019-12-16")),
                  end.date = as.Date(c("2020-03-29", 
                                       "2020-07-12", 
                                       "2020-03-02",
                                       "2020-06-29",
                                       "2020-08-27",
                                       "2020-02-03",
                                       "2020-07-30",
                                       "2020-03-02",
                                       "2020-05-15",
                                       "2019-12-20")),
                  type = c('Lab', 
                           'Reading',
                           'Taught module',
                           'Computational',
                           'Deadline',
                           'Computational',
                           'Writing',
                           'Writing',
                           'Writing',
                           'Computational'))

# plot Gantt chart

brks <- c(seq(from = as.Date("2019-12-09"), to = as.Date("2020-09-06"),
         by = "month"))
         
pdf("GanttChart.pdf", 9, 3)
data %>% mutate(task = fct_reorder(task, desc(start.date))) %>% # order tasks by start date in decending order
         ggplot(aes(x = start.date, 
                 xend = end.date, 
                 y = task,
                 yend = task,
                 group = start.date,
                 colour = type)) +
  geom_segment(size = 6) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        legend.text=element_text(size=11),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11)) +
  scale_x_date(breaks = brks, 
               date_labels = "%b",
               position = "top") +
  labs(x = NULL,
       y = NULL)
dev.off()
                  
                  