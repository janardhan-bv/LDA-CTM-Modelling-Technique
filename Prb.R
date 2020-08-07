library(magrittr) # for %$% operator
library(tidyverse) # R is better when it is tidy
library(tidytext)  # for easy handling of text
library(ldatuning) # for mathematical hint of number of topics
library(topicmodels) # nomen est omen
library(knitr)
#install.packages("kableExtra")
library(kableExtra)

if (!exists("haiku_clean") || !exists("lemma_unique")){
  if (!file.exists("text_lemma.RData")){
    res <- tryCatch(download.file("http://bit.ly/haiku_lemma",
                                  "text_lemma.RData", mode = "wb"),
                    error=function(e) 1)
  }
  load("text_lemma.RData")
}

haiku_clean1 <- head(haiku_clean, 100)
lemma_unique1 <- head(lemma_unique, 1300)

haiku_clean1 <- haiku_clean1 %>%
  mutate (h_number = row_number())

haiku_lemma <- haiku_clean1 %>%
  unnest_tokens(word, text) %>%
  select (h_number, word) %>%
  left_join (lemma_unique1 %>%
               select(word:lemma, synonym))

#haiku_lemma <- haiku_lemma()

dtm_lemma <- haiku_lemma %>%
  count(h_number, lemma) %>%
  filter(!is.na(lemma)) %>%
  cast_dtm(h_number, lemma, n)

as.matrix(dtm_lemma)

control_list_gibbs <- list(
  burnin = 500,
  iter = 1000,
  seed = 0:4,
  nstart = 5,
  best = TRUE
)

system.time(
  topic_number_lemma <- FindTopicsNumber(
    dtm_lemma,
    topics = c(seq(from = 2, to = 9, by = 1), seq(10, 20, 2), seq(25, 50, 5)),
    metrics = c( "Griffiths2004", "CaoJuan2009", "Arun2010", "Deveaud2014"),
    method = "Gibbs",
    control = control_list_gibbs,
    mc.cores = 4L,
    verbose = TRUE
  )
)

FindTopicsNumber_plot(topic_number_lemma)

para <- tibble(k = c(2,3,9,12,50))

  lemma_tm <- para %>%
    mutate(lda = map(k,
                     function(k) LDA(
                       k=k,
                       x=dtm_lemma,
                       method="Gibbs",
                       control=control_list_gibbs
                     )
    )
    )


lemma_tm <- lemma_tm %>%
  mutate(lda_gamma = map(.x=lda,
                         .f=tidytext::tidy,
                         matrix="gamma"))

lemma_tm %>%
  unnest(lda_gamma) %>%
  group_by(k, document) %>%
  arrange(desc(gamma)) %>%
  slice(1) %>%
  #top_n(1, gamma) %>%
  ungroup() %>%
  ggplot(aes(x=gamma, fill=factor(k))) +
  geom_histogram(bins = 20) +
  scale_fill_discrete(name = "Number of\nTopics") +
  xlab("maximum gamma per document") +
  facet_wrap(~k) +
  geom_vline(aes(xintercept = 1/k),
             tibble(k=lemma_tm %$% unique(k)),
             color="darkred")

control_list_ctm <- list(
  seed = 5:9,
  nstart = 5,
  best = TRUE
)

system.time(
  lemma_tm <- lemma_tm %>%
    mutate(ctm = map(k,
                     function(k) CTM(
                       k=k,
                       x=dtm_lemma,
                       control=control_list_ctm
                     )
    )
    )
)

tidy_ctm_gamma  <- function(CTM_object){
  CTM_object %>%
    slot("gamma")  %>%
    as_tibble()  %>%    #changed value
    mutate (document = row_number()) %>%
    gather(topic, gamma, -document) %>%
    mutate(topic = strtoi(stringr::str_sub(topic,2)))
}

lemma_tm <- lemma_tm %>%
  mutate(ctm_gamma = map(.x=ctm, .f=tidy_ctm_gamma))

lemma_tm %>%
  unnest(ctm_gamma) %>%
  group_by(k, document) %>%
  arrange(desc(gamma)) %>%
  slice(1) %>%
  #top_n(1, ctm_gamma) %>%
  ungroup() %>%
  ggplot(aes(x=gamma, fill=factor(k))) +
  geom_histogram(bins = 15) +
  scale_fill_discrete(name = "Number of\nTopics") +
  facet_wrap(~k) +
  geom_vline(aes(xintercept = 1/k),
             tibble(k=lemma_tm %$% unique(k)),
             color="darkred")
lemma_tm <- lemma_tm %>%
  mutate(lda_beta = map(.x=lda, .f=tidytext::tidy, matrix = "beta"))

tidy_ctm_beta  <- function(CTM_object){
  Terms  <- CTM_object %>%
    slot("terms")
  
  CTM_object %>%
    slot("beta")  %>%
    as_tibble() %>%
    setNames(Terms) %>%
    mutate (topic = row_number()) %>%
    gather(term, beta, -topic) %>%
    mutate(beta = exp(beta))
}

lemma_tm <- lemma_tm %>%
  mutate(ctm_beta = map(.x=ctm, .f=tidy_ctm_beta))

lemma_tm  %>%
  unnest(ctm_beta) %>%
  filter(k==9) %>%
  group_by(topic) %>%
  top_n(5, beta) %>%
  arrange(topic, -beta)  %>%
  mutate(rank = paste0("rank ",row_number())) %>%
  ungroup() %>%
  select(-c(k, beta)) %>%
  spread(rank, term) %>%
  knitr::kable(format="html") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))

k9_ctm_gamma <- lemma_tm %>%
  unnest(ctm_gamma) %>%
  filter(k==9) %>%
  left_join(lemma_tm  %>%
              unnest(ctm_beta) %>%
              filter(k==9) %>%
              group_by(topic) %>%
              top_n(1, beta) %>%
              select(topic, term)) %>%
  rename(`top term`=term)

k9_ctm_gamma %>%
  select (-k)  %>%
  group_by(document) %>%
  arrange(desc(gamma)) %>%
  slice(1) %>%
  ungroup() %>%
  group_by(topic, `top term`) %>%
  summarize(n=n()) %>%
  ungroup() %>%
  mutate(percent=round(n/sum(n)*100, 1)) %>%
  knitr::kable(format="html") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))

## helper function to sort
## correlation matrix
reorder_cor <- function(x){
  ord <- corrplot::corrMatOrder(x)
  x[ord,ord]
}

## helper function to extract
## lower triangle of matrix
get_lower_tri<-function(x){
  x[upper.tri(x)] <- NA
  return(x)
}



## prepare data
cor_data <- k9_ctm_gamma %>%
  select(-topic) %>%
  spread('top term', gamma) %>%
  select(-c(k, document, lda, lda_gamma, ctm, lda_beta, ctm_beta)) %>%
  cor() %>%
  reorder_cor() %>%
  get_lower_tri() %>%
  as_tibble() %>%
  mutate(topic1 = forcats::as_factor(paste(names(.)))) %>%
  gather(topic2, correlation, - topic1) %>%
  mutate(topic2 = factor(topic2, levels=levels(topic1)))

#sapply(cor_data, class)

## create plot
cor_data %>%
  ggplot(aes(as.numeric(topic1), as.numeric(topic2), fill=correlation)) +
  geom_tile(color="white") +
  scale_fill_gradient2(low = "firebrick4", high = "dodgerblue4",
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       na.value = "white",
                       name="Pearson\nCorrelation")  +
  scale_x_continuous(
    breaks=1:15, labels = levels(cor_data$topic1), name="") +
  scale_y_continuous(
    breaks=1:15, labels = levels(cor_data$topic2), name="")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  coord_fixed()

beta_day_light <- lemma_tm  %>%
  unnest(ctm_beta) %>%
  filter(k==9)%>%
  left_join(lemma_tm  %>%
              unnest(ctm_beta) %>%
              filter(k==9) %>%
              group_by(topic) %>%
              top_n(1, beta) %>%
              rename(`top term`=term) %>%
              select(topic, `top term`)) %>%
  filter(`top term` %in% c("day", "claim")) %>% #light to claim
  select(-c(k, topic)) %>%
  spread(`top term`, beta) %>%
  filter(day > .001 | claim > .001) %>% #light to claim
  mutate(log_ratio = log2(day / claim)) #light to claim



beta_day_light %>%
  top_n(10, abs(log_ratio)) %>%
  mutate(term = forcats::fct_reorder(term,log_ratio)) %>%
  ggplot(aes(x=term, y=log_ratio, fill=log_ratio<0)) +
  scale_fill_manual(values=c("gold2", "deepskyblue3"),
                    name="topic",
                    labels=c("claim", "day")) +  #light to claim
  geom_col() +
  coord_flip()

dom_top <- k9_ctm_gamma %>%
  group_by(document) %>%
  top_n(1, gamma) %>%
  rename(dom_top = `top term`) %>%
  select(document, dom_top)

k9_ctm_gamma_wide <- k9_ctm_gamma %>%
  select(-c(k, topic)) %>%
  spread(`top term`, gamma)

haiku_gamma <- haiku_clean1 %>%
  left_join(dom_top, by = c("h_number" = "document")) %>%
  left_join(k9_ctm_gamma_wide, by = c("h_number" = "document"))

haiku_gamma %>%
  filter (dom_top == "claim") %>% #light to claim
  top_n (3, claim) %>% #light to claim
  mutate(claim = round(claim,2))  %>% #light to claim
  select (text, author, claim) %>%  #light to claim
  knitr::kable(format="html") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))

haiku_gamma %>%
  filter (dom_top == "day") %>%
  top_n (3, day) %>%
  mutate(day = round(day,2))  %>%
  select (text, author, day) %>%
  knitr::kable(format="html") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))


