# Analyse results of survey on marine microfossil data guidelines

library(tidyverse)
library(brglm2)
library(broom)
library(waffle)
library(patchwork)
library(ggtext)

# prep data, remove prior to publication ####
library(googlesheets4)
library(readxl)

GoogleReps <- read_sheet('https://docs.google.com/spreadsheets/d/1Lka0AOU6q-7mZt2g7MhyBH8zPsip4wuuqladPrJLQ1I/edit?resourcekey#gid=463587945')
short_questions <- read_sheet('https://docs.google.com/spreadsheets/d/1Lka0AOU6q-7mZt2g7MhyBH8zPsip4wuuqladPrJLQ1I/edit?resourcekey#gid=463587945', sheet = 'Sheet1')

continents <- read_sheet('https://docs.google.com/spreadsheets/d/1FWfUm7wNjbAYR0vuQzd1seOMvpYc96bXxh8w2OYE43g/edit#gid=0')

names(GoogleReps) <- short_questions$short

# merge with response collected using MS forms
MSFormsReps <- read_xlsx('Marine micropalaeontological data requirements.xlsx') %>%
  select(all_of(short_questions$excelname[!is.na(short_questions$excelname)]))

names(MSFormsReps) <- short_questions$short[!is.na(short_questions$excelname)]

# MSFormsReps$Country %in% continents$Country

# merge 
reps <- bind_rows(GoogleReps, MSFormsReps)

reps <- reps %>% 
  mutate(
    Career = case_when(Career == 'Yes' ~ 'ECR',
                       Career == 'No' ~ 'ER'),
    Continent = continents$Continent[match(Country, continents$Country)]
    )

reps <- reps %>%
  select(-c("Timestamp", "Country", "First_name", "Last_name", "email"))

# View(reps)

write_csv(reps, "data/MiPaSurveyResponses.csv")

# save grouping
short_questions %>%
  select(short, group) %>%
  filter(!short %in% c("Country", "First_name", "Last_name", "email")) %>%
  drop_na() %>%
  write_csv(., "data/QuestionGrouping.csv")

# read pangaea assessment
PANGAEA_raw <- read_sheet("https://docs.google.com/spreadsheets/d/1VrPC7Y0OORy5ocENicCMWvS2GQYscy3cKuEZXJbA5Ps/edit#gid=0", skip = 6, n_max = 60)

PANGAEA_raw %>%
  drop_na(question) %>%
  mutate(across(everything(), ~gsub("n/a", NA_character_, .))) %>%
  select(-Section, -Description, -CC6) %>%
  write_csv(., "data/StatusLegacyData.csv")
  
rm(list = ls())

# publish from here on
#####

# read data #####
# read survey responses
responses <- read_csv("data/MiPaSurveyResponses.csv")

# read question grouping
qgrouping <- read_csv("data/QuestionGrouping.csv")

# read assessment of legacy datasets
legacy <- read_csv("data/StatusLegacyData.csv")
#####

# ANALYSIS #####
responses_long <- responses %>%
  pivot_longer(
    -c(Career, Fossil_group, Continent),
    names_to = "question",
    values_to = "answer"
    ) %>%
  mutate(
    qgroup = factor(qgrouping$group[match(question, qgrouping$short)]
    )
  ) 


# Demographics

# median and minimum number of answers per question
responses_long %>%
  drop_na(qgroup) %>%
  group_by(question) %>%
  reframe(
    n = sum(!is.na(answer))
    ) %>%
  reframe(
    median(n),
    min(n)
    )

# divsion by career stage
responses %>%
  group_by(Career) %>%
  summarise(n = n()) %>%
  mutate(
    per = round(n/sum(n)*100)
    )

# by continent
responses %>%
  group_by(Continent) %>%
  summarise(n = n()) %>%
  mutate(
    per = round(n/sum(n)*100)
    )

# by fossil group
responses %>%
  group_by(Fossil_group) %>%
  summarise(n = n()) %>%
  mutate(
    per = round(n/sum(n)*100)
    )

# Responses

# proportion of each category per question
# put "no" to "Desired" and "yes" to "Essential"
props <- responses_long %>%
  drop_na(qgroup) %>%
  group_by(question, qgroup) %>%
  reframe(
    Desired = sum(answer == "Desired" | answer == "No", na.rm =  TRUE),
    Recommended = sum(answer == "Recommended", na.rm =  TRUE),
    Essential = sum(answer == "Essential" | answer == "Yes", na.rm =  TRUE)
  ) %>%
  pivot_longer(
    -c(question, qgroup),
    names_to = "cat",
    values_to = "n"
    ) %>%
  group_by(question, qgroup) %>%
  mutate(
    p = n/sum(n),
    cat = factor(cat, levels = c("Essential", "Recommended", "Desired"))
  ) %>%
  ungroup()

# calculate ranking
rankings <- props %>%
  select(-n) %>%
  pivot_wider(
    names_from = cat, 
    values_from = p
    ) %>%
  mutate(
    ranking = case_when(
      Recommended == 0 ~ "yes",
      Desired > 0.5 ~ "desired",
      Recommended > 0.5 ~ "recommended",
      Essential > 0.5 ~ "essential",
      Desired + Recommended > Recommended + Essential & Desired + Recommended > 0.7 ~ "desired-recommended",
      Recommended + Essential > Desired + Recommended & Recommended + Essential > 0.7 ~ "recommended-essential",
      TRUE ~ "recommended"
    )
  )

# distribution of rankings
rankings %>%
  group_by(ranking) %>%
  reframe(n = n())

# Adjacent category logistic regression
aclr_in <- responses_long %>%
  filter(answer != "Yes", answer != "No") %>%
  drop_na() %>%
  mutate(
    Fossil_group = case_when(!grepl("foraminifera", Fossil_group) ~ "other",
                             TRUE ~ Fossil_group),
    across(c(Career, Fossil_group, Continent), as.factor),
    answer = factor(answer, labels = c("desired", "recommended", "essential"), ordered = TRUE)
    ) %>%
  group_by(answer, Career, Fossil_group) %>%
  reframe(n = n())

# set the contrast
contrasts(aclr_in$Fossil_group) <- contr.treatment(
  levels(aclr_in$Fossil_group),
  base = which(levels(aclr_in$Fossil_group) == "Benthic foraminifera")
  )

# make the model
acm <- bracl(answer ~ Career + Fossil_group, weights = n, data = aclr_in, parallel = FALSE, type = "ML")

# display results
tidy(acm) %>%
  mutate(
    odds_ratios = case_when(p.value <= 0.05 ~ round(exp(estimate), 1),
                                 TRUE ~ NA_real_)
    )

# Assess legacy data
# remove some unnecessary information
legacy_clean <- legacy %>%
  filter(!is.na(qgroup)) %>%
  pivot_longer(
    -c(question, qgroup),
    names_to = "fossil_group",
    values_to = "included"
    ) %>%
  mutate(
    fossil_group = gsub('[[:digit:]]+', '', fossil_group)
    )

# calculate proportion of attributes included and merge with rankings
legacy_assess <- legacy_clean %>%
  drop_na() %>%
  group_by(question, qgroup, included) %>%
  reframe(
    n = n()
    ) %>%
  group_by(question, qgroup) %>%
  mutate(t = sum(n),
         p = n/t
  ) %>%
  ungroup() %>%
  left_join(., 
            rankings %>%
              select(question, ranking),
            by = "question"
            ) 

# proportion of essential attributes included on average
legacy_assess %>%
  filter(ranking == "essential" & included == "yes") %>%
  summarise(mean(p))

# FIGURES #####
# set theme for figures
theme_set(theme_bw() + 
            theme(line = element_line(colour = "#464646"),
                  rect = element_rect(fill = NA, colour = NA),
                  text = element_text(colour = "#000000", size = 9),
                  axis.ticks = element_line(colour = "#464646"), 
                  legend.key = element_rect(colour = NA, fill = NA), 
                  panel.border = element_rect(colour = "#464646", fill = NA),
                  panel.grid = element_blank(),
                  plot.background = element_blank(),
                  panel.background = element_blank(),
                  strip.background = element_blank(),
                  plot.title = element_text(hjust = 0.5, size = 7),
                  plot.subtitle = element_text(hjust = 0.5))
)

# Demographics
careerWaffle <- responses %>%
  group_by(Career) %>%
  summarise(n = n()) %>%
  mutate(Career = case_when(Career == "ECR" ~ "ECR: <5 years since PhD",
                            TRUE ~ "ER: >5 years since PhD")) %>%
  ggplot(aes(fill = Career, values = n)) +
  geom_waffle(colour = 'white', n_rows = ceiling(sqrt(nrow(responses))), size = 0.25, na.rm = FALSE) +
  coord_equal() +
  scale_fill_brewer(palette = "Blues") +
  labs(title = 'Career stage',
       fill = NULL) +
  theme_void() +
  guides(fill = guide_legend(nrow = 2)) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)
  )


geoWaffle <- responses %>%
  group_by(Continent) %>%
  summarise(n = n()) %>%
  ggplot(aes(fill = Continent, values = n)) +
  geom_waffle(colour = 'white', n_rows = ceiling(sqrt(nrow(responses))), size = 0.25, na.rm = FALSE) +
  coord_equal() +
  scale_fill_brewer(palette = "Blues") +
  labs(title = 'Geographic origin',
       fill = NULL) +
  guides(fill = guide_legend(nrow = 3)) +
  theme_void() +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)
  )


groupWaffleInput <- responses %>%
  group_by(Fossil_group) %>%
  summarise(n = n()) %>%
  mutate(group = case_when(Fossil_group == 'Other' ~ 1,
                           TRUE ~ 0)) %>%
  group_by(group) %>%
  arrange(desc(n), .by_group = TRUE) %>%
  ungroup()

groupWaffle <- groupWaffleInput %>%
  mutate(Fossil_group = factor(Fossil_group, levels = groupWaffleInput$Fossil_group)) %>%
  ggplot(aes(fill = Fossil_group, values = n)) +
  geom_waffle(colour = 'white', n_rows = ceiling(sqrt(nrow(responses))), size = 0.25, na.rm = FALSE) +
  coord_equal() +
  scale_fill_brewer(palette = "Blues") +
  labs(title = 'Microfossil group',
       fill = NULL) +
  guides(fill = guide_legend(ncol = 2)) +
  theme_void() +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)
  )

demoPLot <- careerWaffle + geoWaffle + groupWaffle + plot_annotation(tag_levels = "a")

ggsave(filename = "~/Library/CloudStorage/GoogleDrive-texelwildlife@gmail.com/My Drive/manuscripts/mamipacs/figures/demo.png",
       plot = demoPLot,
       dpi = "retina", 
       units = "mm", 
       width = 360, 
       height = 180)

# Likert plots
# make input
likert_input <- props %>%
  group_by(question) %>%
  mutate(
    xmin = case_when(
      cat == "Desired" ~ -p[cat == "Recommended"]/2 - p,
      cat == "Recommended" ~ -p/2,
      cat == "Essential" ~ p[cat == "Recommended"]/2
    ),
    xmax = case_when(
      cat == "Desired" ~ -p[cat == "Recommended"]/2,
      cat == "Recommended" ~ p/2,
      cat == "Essential" ~ p[cat == "Recommended"]/2 + p
    )
  ) %>%
  ungroup() 

# define plotting order
ord_likert <- likert_input %>%
  group_by(qgroup, question) %>%
  reframe(Max = max(xmax)) %>%
  arrange(Max, .by_group = TRUE) %>%
  ungroup()

likert_input <- likert_input %>%
  mutate(question = factor(question, levels = ord_likert$question))

n_responses <- likert_input %>%
  group_by(question, qgroup) %>%
  reframe(n = sum(n)) %>%
  mutate(
    question = factor(question, levels = ord_likert$question)
    )

# attribute categories/plots
cats <- likert_input %>%
  group_by(qgroup) %>%
  reframe(n = n()/3)

rankings <- rankings %>%
  mutate(
    xpos = Recommended/2 + Essential,
    question = factor(question, levels = ord_likert$question)
  )

# make plots with height set by number of bars
map(cats$qgroup, function(x){
  plt <- likert_input %>%
    filter(qgroup == x) %>%
    ggplot(aes(y = question, colour = cat)) +
    geom_segment(aes(x = xmin, xend = xmax, yend = question), linewidth = 4) +
    scale_colour_brewer(palette = "Blues", guide = "none") +
    geom_text(data = rankings %>%
                filter(qgroup == x),
              aes(x = xpos, question, label = ranking),
              inherit.aes = FALSE,
              colour = "black",
              size = 2,
              hjust = 1,
              nudge_x = -0.02) +
    geom_text(data = n_responses %>%
                filter(qgroup == x),
              aes(x = 1, question, label = n), inherit.aes = FALSE, colour = "grey20", size = 1.5) +
    scale_x_continuous(limits = c(-0.7, 1), 
                       expand = expansion(mult = c(0.05, 0.05),
                                          add = c(0, 0)),
                       labels = function(x) paste0(x*100, "%")
    ) +
    scale_y_discrete(expand = expansion(mult = c(0, 0),
                                        add = c(0.6, 0.6))
    ) +
    labs(title = x,
         colour = NULL,
         y = NULL,
         x = NULL)
  
  ggsave(paste0("~/Library/CloudStorage/GoogleDrive-texelwildlife@gmail.com/My Drive/manuscripts/mamipacs/figures/", x, ".png"),
         plt, 
         dpi = "retina", 
         height = (cats$n[cats$qgroup == x] + 0.5)* 5 + 10, 
         width = 100, 
         units = "mm")
})

# Legacy data plot

# default attributes
pangaeaDefault <- c("Source", "Link to external ontology", "Contributor", "Links to ancilary data")

legacy_plot <- legacy_assess %>%
  mutate(question = ifelse(question %in% pangaeaDefault, paste0(question, "*"), question),
         y.label = ifelse(ranking == "essential" | ranking == "yes" | question == "Split", paste0("***", question, "***"), question),
         qgroup = factor(qgroup,
                         levels = cats$qgroup[c(6, 7, 5, 2, 8, 3, 4, 1)]
         )
  ) %>%
  ggplot(aes(p, y.label, fill = included)) +
  geom_bar(stat = "identity") +
  geom_text(aes(x = 1, y.label, label = t), inherit.aes = FALSE, colour = "grey20", size = 2, hjust = 0, nudge_x = 0.01) +
  facet_wrap(~qgroup, scale = "free_y", ncol = 2) +
  theme(axis.text.y = element_markdown()) +
  scale_fill_brewer(palette = "Blues") +
  scale_x_continuous(labels = scales::percent,
                     expand = expansion(mult = c(0.05, 0.1))) +
  labs(y = NULL,
       x = NULL,
       fill = "Included") +
  theme(legend.position = "bottom",
        legend.key.size = unit(4, "mm")
        )

ggsave(filename = "~/Library/CloudStorage/GoogleDrive-texelwildlife@gmail.com/My Drive/manuscripts/mamipacs/figures/PANGAEA_state.png",
       plot = legacy_plot,
       dpi = "retina", 
       units = "mm", 
       width = 160, 
       height = 200)
