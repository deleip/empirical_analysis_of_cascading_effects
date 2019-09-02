# Supplementary code to "Empirical Analysis of Cascading Effects using CCM, Internship Report"
# Once the section "Settings" is run, all other sections can be run without first runing the above 
# sections (assuming that the corresponding files were downloaded as well). Especially the sections
# "Run CCM on datasets" and "Include Climate [...]" can take quite a while and should be skipped if no 
# changes have been done.

# Everything can be reproduced from scratch using just the files:
# data_analysis.R
# utilities.R
# salmon_data.csv
# comtrade_banana_prepared.csv
# CL_FI_COUNTRY_GROUPS.csv
# CL_FI_SPECIES_GROUPS.csv
# TS_FI_CAPTURE.csv
# updated_climatedata.RData



######## Settings ########

library(tidyverse)
library(lubridate)
library(rEDM)
library(lattice)
library(reshape2)
library(future)
library(future.apply)
library(sna)
library(rlist)

source("utilities.R")
#

######## Datapreparation ########

dat_csv <- read_csv("salmon_data.csv")   # already preprocessed 

dat_csv <- dat_csv %>%
  ungroup() %>%
  mutate(Period = parse_date_time(Period, orders = "Ym"))

dat_csv$Period %>% unique() %>% length() # 120 time points!
dat_csv$link %>% unique() %>% length() # 4483 links

## Delete Bunkers, that's not a country. For more info: http://www.fao.org/faostat/en/#data/EM/metadata
dat_csv <- filter(dat_csv, importer_name != "Bunkers")
dat_csv <- filter(dat_csv, importer_name != "Ic")
dat_csv <- filter(dat_csv, exporter_name != "Ic")



# choosing all fresh and frozen commodities
dat <- dat_csv %>%
  filter(Commodity_Code %in% c('030212','030213', '030214', # fresh salmon
                               '030441',                      # fillets
                               '030311', '030312', '030313', '030319', '030322',    # frozen salmon
                               '030481'                           # fillets frozen
  ))   %>%
  select(-starts_with("diff"), -starts_with("Commodity")) %>%
  group_by(link, Period) %>%
  select(-n_obs) %>%
  mutate(
    total_netweight_tons = sum(max_netweight)/1000,
    total_value_us = sum(max_value_us)
  ) %>%
  select(-starts_with("max")) %>% unique() %>% ungroup()

dat <- dat %>%
  filter(Year > 2009 & Year < 2018)     # choosing right timeframe

dat <- dat %>% filter(exporter != importer)     # excluding trade within a country

save(dat, file = "./intermediate_results/prepared_data_fresh_frozen.RData")  

# get mapping between country names and country code
code_match <- dat %>% select(importer, importer_name) %>% distinct()
colnames(code_match) <- c("Code", "Country")
code_match2 <- dat %>% select(exporter, exporter_name) %>% distinct()
colnames(code_match2) <- c("Code", "Country")
code_match <- rbind(code_match, code_match2) %>% distinct()
code_match <- code_match[order(code_match$Code),]
rm(code_match2)

save(code_match, file = "./intermediate_results/country_codes.RData")



# choosing only fresh and frozen commodities
dat <- dat_csv %>%
  filter(Commodity_Code %in% c("030212", "030213", "030214", "030441") ) %>%
  select(-starts_with("diff"), -starts_with("Commodity")) %>%
  group_by(link, Period) %>%
  select(-n_obs) %>%
  mutate(
    total_netweight_tons = sum(max_netweight)/1000,
    total_value_us = sum(max_value_us)
  ) %>%
  select(-starts_with("max")) %>% unique() %>% ungroup()

dat <- dat %>%
  filter(Year > 2009 & Year < 2018)     # choosing right timeframe

dat <- dat %>% filter(exporter != importer)    # exculde export within a country

save(dat, file = "./intermediate_results/prepared_data_only_fresh.RData") 



######## Run CCM on datasets ########

# filters out timeseries with less than min_obs obsevations, and all that are "too linear" (at the moment
# checked by rho_{theta=0} < \rho_max - nonlin) and detrends the timeseries. Then ccm is run for all possible
# pairs of trade connections, those where rho is not converging are filtered out and than t-test are performed
# to decide which are strong/weak causally related. 
# Meta data is collected to know at which point how many timeseries/links are dropped.

load(file = "./intermediate_results/prepared_data_fresh_frozen.RData")
res <- get_links(dat, min_obs = 60, detrend = TRUE, nonlin = 0.05)
save(res, file = "./intermediate_results/ccm_res_fresh_frozen.RData")

load(file = "./intermediate_results/prepared_data_only_fresh.RData")
res <- get_links(dat, min_obs = 60, detrend = TRUE, nonlin = 0.05)
save(res, file = "./intermediate_results/ccm_res_only_fresh.RData")

######## Nonlinearity Visualization  #############

load(file = "./intermediate_results/ccm_res_fresh_frozen.RData")
load(file = "./intermediate_results/country_codes.RData")

df <- res$df
nonlinearity <- res$nonlinearity

interesting_shapes <- c(1, 19, 20, 288, 410) 

plot(nonlinearity[[1]]$theta, nonlinearity[[1]]$rho, ylim = c(0,1), col = 1, type = "l", xlab = "Degree of nonlinearity (Theta)", ylab = "Prediction Skill (Rho)")
c <- 1
for(i in interesting_shapes[2:5]){
  c <- c +1
  lines(nonlinearity[[i]]$theta, nonlinearity[[i]]$rho, ylim = c(0,1), col = c, type = "l")
}

l <- colnames(df) %>%  .[interesting_shapes+1] %>% as.data.frame() %>% 
  separate(1,  into = c("A", "B"), sep = "_") %>% codes_to_names(cols = 1:2, match = code_match)

legend("topright", legend = paste(l[,1], l[,2], sep = " -> "), col = 1:10, lty = 1)


######## Analysis Links for Main Dataset #######

load(file = "./intermediate_results/ccm_res_fresh_frozen.RData")

res$meta_information

links <- res$links
links_s <- links %>% filter(strong_detection == TRUE)
links_w <- links %>% filter(detection == TRUE)

# Total strong links: 3883
links_s %>% dim() %>% .[1]
# motif A -> B <- C:  61
links_s %>% filter(B == D) %>% dim() %>% .[1]
# motif A <- B -> C:  162
links_s %>% filter(A == C) %>% dim() %>% .[1]
# motif A -> B -> C: 60
links_s %>% filter(B == C) %>% dim() %>% .[1]
# motif A -> B <- A:  1
links_s %>% filter(D == A & C == B) %>% dim() %>% .[1]  
# motif A <- B <- C: 83
links_s %>% filter(A == D) %>% dim() %>% .[1]
# links w/o motif: 3518
links_s %>% filter(D != A & D != B & C != A & C != B) %>% dim() %>% .[1]  

# Total weak links: 8701
links_w %>% dim() %>% .[1]
# motif A -> B <- C: 171
links_w %>% filter(B == D) %>% dim() %>% .[1]
# motif A <- B -> C: 459
links_w %>% filter(A == C) %>% dim() %>% .[1] 
# motif A -> B -> C: 142
links_w %>% filter(B == C) %>% dim() %>% .[1]
# motif A -> B -> A: 3
links_w %>% filter(C == B & D == A) %>% dim() %>% .[1] 
# motif A <- B <- C: 173
links_w %>% filter(A == D) %>% dim() %>% .[1]
# links w/o motif: 7759
links_w %>% filter(D != A & D != B & C != A & C != B) %>% dim() %>% .[1]  



######## Geography as common factor? (For Main Dataset) #########

load(file = "./intermediate_results/prepared_data_fresh_frozen.RData")
load(file = "./intermediate_results/ccm_res_fresh_frozen.RData")
load(file = "./intermediate_results/country_codes.RData")

links <- res$links
links_s <- links %>% filter(strong_detection == TRUE)
links_w<- links %>% filter(detection == TRUE)

# strong links without explanation through underlying motive in trade network
links_s_wo_motif <- links_s %>% filter(B != D  & A != C & B != C & A != D) %>% select(lib_column, target_column, rho) 

links_s_wo_motif <- links_s_wo_motif %>%
  separate(lib_column, into = c("A", "B"), sep = "_", remove = FALSE) %>%
  separate(target_column, into = c("C", "D"), sep = "_", remove = FALSE)

links_s_wo_motif <- links_s_wo_motif %>% codes_to_names( c(2,3,5,6), match = code_match) %>% arrange(desc(rho))

links_s_wo_motif



######## Including Climate Dataset (For Main Dataset, Strong links) ########

load(file = "./intermediate_results/updated_climatedata.RData")
load(file = "./intermediate_results/ccm_res_fresh_frozen.RData")

links <- res$links
links_s <- links %>% filter(strong_detection == TRUE)
dim(links_s)    ## 3883  strong links



## find all feedback loops 

count <- 0
feedbacks <- list()
links_check <- links_s
for(i in 1:(dim(links_check)[1]-1)){
  links_loop <- links_check %>% filter(A == .[1,"C"]$C & B == .[1,"D"]$D & C == .[1,"A"]$A & D == .[1,"B"]$B) 
  if(dim(links_loop)[1]>0){
    if(dim(links_loop)[1]>1){print("???")}
    count <- count + 1
    #feedbacks[[count]] <- rbind(ind2[1,], ind_loop)
    feedbacks[[count]] <- links_loop %>% select(lib_column, target_column)
  }
  links_check <- links_check[-1,]
}

length(feedbacks)   ## 351 feedback loops


## get timeseries for those links

load(file = "./intermediate_results/prepared_data_fresh_frozen.RData")

### Prepare dataset for ccm
df <- dat %>% ungroup() %>% #skimr::skim()
  select(Period, link, total_netweight_tons) %>%
  group_by(Period, link) %>%
  spread(link, total_netweight_tons) %>%
  as.data.frame()

res <- lapply(feedbacks, function(x) {df %>% select(Period,x[1,1][[1]], x[1,2][[1]])})

res <- lapply(res, function(x){right_join(x, climdat)})

# calculate ccm for each listelement, and all combinations clim/trade, trade/clim and for trade/trade
res <- lapply(res, links_for_climdat, detrend = TRUE)

links <- lapply(res, function(x){return(x[[1]])})      

detection <- lapply(links, function(x){x <- x %>% filter(detection == TRUE); return(x)})
# strong_detection <- lapply(links, function(x){x <- x %>% filter(strong_detection == TRUE); return(x)})

# we are interested in cases where the trade connections predict a climate index, as that means that climate is causing the trade connection
detection <- lapply(detection, function(x){x <- x %>% filter(is.na(D)); return(x)})

# strong_detection <- lapply(strong_detection, function(x){x <- x %>% filter(is.na(D)); return(x)})

# and then in those that can be predicted by both of the trade connections, so that the climate index is a common driver
shared_driver <- lapply(detection, function(x){sel1 <- duplicated(x$C) %>% which()
sel2 <- duplicated(x$C, fromLast = T) %>% which()
sel <- c(sel1, sel2)
return(x[sel,])})

# shared_driver <- lapply(strong_detection, function(x){sel1 <- duplicated(x$C) %>% which()
#                                                sel2 <- duplicated(x$C, fromLast = T) %>% which()
#                                                sel <- c(sel1, sel2)
#                                                return(x[sel,])})

shared_driver_clean <- list.clean(shared_driver, function(x) dim(x)[1] == 0L)

length(shared_driver_clean)   ## 125 for strong

motifs <- lapply(shared_driver_clean, function(x){c(x$A, x$B) %>% unique() %>% length()}) %>% unlist()
sum(motifs != 4) ## 19 with motif for strong, leaving 106 between four different countries

number_climate_indices <- lapply(shared_driver_clean, function(x){return(dim(x)[1]/2)}) %>% unlist()
sum(number_climate_indices == 1) ## 68 share only one climate index as driver, 40 shaeres 2, 13 share 3, 3 share 4, and 1 shares 5 (strong)

load(file = "./intermediate_results/country_codes.RData")
shared_driver_clean <- lapply(shared_driver_clean, function(x){x <- codes_to_names(x, c(2,3), match = code_match)})

## remove these links from the link list created by CCM
links_rm <- lapply(shared_driver_clean, function(x){return(x$lib_column) %>% unique()}) %>% unlist()

load(file = "./intermediate_results/ccm_res_fresh_frozen.RData")
links <- res$links
links_s <- links %>% filter(strong_detection == TRUE)  ## 227 links

for(i in 1:(length(links_rm)/2)){
  links_s <- links_s[!(links_s$lib_column == links_rm[2*(i-1)+1] & links_s$target_column == links_rm[2*i]),]
  links_s <- links_s[!(links_s$lib_column == links_rm[2*i] & links_s$target_column == links_rm[2*(i-1)+1]),]
}


save(links_s, file = "./intermediate_results/climate_fresh_frozen_strong.RData" )

######## Including Climate Dataset (For Main Dataset, Weak links) ########

load(file = "updated_climatedata.RData")
load(file = "./intermediate_results/ccm_res_fresh_frozen.RData")

links <- res$links
links_w<- links %>% filter(detection == TRUE)

dim(links_w)    ## 8701  weak links

## find all feedback loops 

count <- 0
feedbacks <- list()
links_check <- links_w            ### this segment was run twice, once with links_check <- links_w, once with links_check <- links_s
for(i in 1:(dim(links_check)[1]-1)){
  links_loop <- links_check %>% filter(A == .[1,"C"]$C & B == .[1,"D"]$D & C == .[1,"A"]$A & D == .[1,"B"]$B) 
  if(dim(links_loop)[1]>0){
    if(dim(links_loop)[1]>1){print("???")}
    count <- count + 1
    #feedbacks[[count]] <- rbind(ind2[1,], ind_loop)
    feedbacks[[count]] <- links_loop %>% select(lib_column, target_column)
  }
  links_check <- links_check[-1,]
}

length(feedbacks)   ## weak: 1373 feedback lops


## get timeseries for those links

load(file = "./intermediate_results/prepared_data_fresh_frozen.RData")

### Prepare dataset for ccm
df <- dat %>% ungroup() %>% #skimr::skim()
  select(Period, link, total_netweight_tons) %>%
  group_by(Period, link) %>%
  spread(link, total_netweight_tons) %>%
  as.data.frame()


res <- lapply(feedbacks, function(x) {df %>% select(Period,x[1,1][[1]], x[1,2][[1]])})

res <- lapply(res, function(x){right_join(x, climdat)})

# calculate ccm for each listelement, and all combinations clim/trade, trade/clim and for trade/trade

res <- lapply(res, links_for_climdat, detrend = TRUE)

links <- lapply(res, function(x){return(x[[1]])})      

detection <- lapply(links, function(x){x <- x %>% filter(detection == TRUE); return(x)})

# we are interested in cases where the trade connections predict a climate index, as that means that climate is causing the trade connection
detection <- lapply(detection, function(x){x <- x %>% filter(is.na(D)); return(x)})

# and then in those that can be predicted by both of the trade connections, so that the climate index is a common driver
shared_driver <- lapply(detection, function(x){sel1 <- duplicated(x$C) %>% which()
                                                sel2 <- duplicated(x$C, fromLast = T) %>% which()
                                                sel <- c(sel1, sel2) 
                                                return(x[sel,])})

shared_driver_clean <- list.clean(shared_driver, function(x) dim(x)[1] == 0L)

length(shared_driver_clean)   ## 498 causal links have at least one climate index as shared driver and hence might 
                              ## not be causally related after all 

motifs <- lapply(shared_driver_clean, function(x){c(x$A, x$B) %>% unique() %>% length()}) %>% unlist()
sum(motifs != 4) ## 84 show some sort of motif (only 3 involved countries), while 414 are between four different countries (weak)

number_climate_indices <- lapply(shared_driver_clean, function(x){return(dim(x)[1]/2)}) %>% unlist()
sum(number_climate_indices == 1) ## 305 share only one climate index as driver, 135 share 2, 48 share 3, 9 share 4  and 1 shares 5 (weak)


load(file = "./intermediate_results/country_codes.RData")
shared_driver_clean <- lapply(shared_driver_clean, function(x){x <- codes_to_names(x, c(2,3), match = code_match)})


## remove these links from the link list created by CCM
links_rm <- lapply(shared_driver_clean, function(x){return(x$lib_column) %>% unique()}) %>% unlist()

load(file = "./intermediate_results/ccm_res_fresh_frozen.RData")
links <- res$links
links_w <- links %>% filter(detection == TRUE)  ## 2481 links

for(i in 1:(length(links_rm)/2)){
  links_w <- links_w[!(links_w$lib_column == links_rm[2*(i-1)+1] & links_w$target_column == links_rm[2*i]),]
  links_w <- links_w[!(links_w$lib_column == links_rm[2*i] & links_w$target_column == links_rm[2*(i-1)+1]),]
}


save(links_w, file = "climate_fresh_frozen_weak.RData" )



######## Links needing explanation for later comparison (Main Dataset) #######

load(file = "./intermediate_results/climate_fresh_frozen_strong.RData")
load(file = "./intermediate_results/climate_fresh_frozen_weak.RData")

links_s_wo_motif <- links_s %>% filter(B != D  & A != C & B != C & A != D) %>% select(lib_column, target_column, rho) 
timeseries_s <- c(links_s$lib_column, links_s$target_column) %>% unique
l <- length(timeseries_s)

wanted_links_s <- matrix(0, nrow = l, ncol = l)
colnames(wanted_links_s) <- timeseries_s
rownames(wanted_links_s) <- timeseries_s
for(i in 1:(dim(links_s_wo_motif)[1])){
  wanted_links_s[links_s_wo_motif$lib_column[i], links_s_wo_motif$target_column[i]] <- 1
}


links_w_wo_motif <- links_w %>% filter(B != D  & A != C & B != C & A != D) %>% select(lib_column, target_column, rho) 
timeseries_w <- c(links_w$lib_column, links_w$target_column) %>% unique
l <- length(timeseries_w)

wanted_links_w <- matrix(0, nrow = l, ncol = l)
colnames(wanted_links_w) <- timeseries_w
rownames(wanted_links_w) <- timeseries_w
for(i in 1:(dim(links_w_wo_motif)[1])){
  wanted_links_w[links_w_wo_motif$lib_column[i], links_w_wo_motif$target_column[i]] <- 1
}

save(wanted_links_s, wanted_links_w, timeseries_w, timeseries_s, file = "./intermediate_results/need_explanation_fresh_frozen.RData" )

######## Trasitivety Analysis (Main Dataset, Strong links) ########

load(file = "./intermediate_results/climate_fresh_frozen_strong.RData")
load(file = "./intermediate_results/climate_fresh_frozen_weak.RData")
load(file = "./intermediate_results/need_explanation_fresh_frozen.RData")   


dim(links_s)    ## 3633 strong links left

sum(wanted_links_w)  # 6931 need explanation

links_s_motif <- links_s %>% filter(B == D  | A == C | B == C | A == D) %>% select(lib_column, target_column, rho)
dim(links_s_motif)  ## 327 strong links with motif

m1 <- matrix(0,nrow = length(timeseries_s), ncol = length(timeseries_s))
colnames(m1) <- timeseries_s
rownames(m1) <- timeseries_s
for(i in 1:(dim(links_s_motif)[1])){
  m1[links_s_motif$lib_column[i], links_s_motif$target_column[i]] <- links_s_motif$rho[i]
}

m_all <- list()
null_all <- list()
m_all[[1]] <- m1
null_all[[1]] <- m1
null_all[[1]][m1 != 0] <- 1
m <- list()
n <- list()
n[[1]] <- null_all[[1]]
number_links <- list()
added_links <- list()
m[[1]] <- m1
number_links[[1]] <- sum(m1!=0) - dim(m1)[1] 
for(i in 1:20){
  n[[i+1]] <- n[[i]] %*% n[[1]]
  m[[i+1]] <- m[[i]] %*% m1
  m_all[[i+1]] <- m_all[[i]]
  m_all[[i+1]][m_all[[i]] == 0] <- m[[i+1]][m_all[[i]] == 0]
  null_all[[i+1]] <- null_all[[i]]
  null_all[[i+1]][null_all[[i]] == 0] <- n[[i+1]][null_all[[i]] == 0]
  # average rho for new links that can be reached by multiple paths of length i+1
  m_all[[i+1]][null_all[[i]] == 0 & null_all[[i+1]] != 0] <- m_all[[i+1]][null_all[[i]] == 0 & null_all[[i+1]] != 0]/
    null_all[[i+1]][null_all[[i]] == 0 & null_all[[i+1]] != 0]
  number_links[[i+1]] <- sum(m_all[[i+1]]!=0)
  added_links[[i]] <- number_links[[i+1]] - number_links[[i]]
}


#### Are the new causal links part of the ones we couldn't explain?

links_expl <- list()
for(c in 1:20){
  tmp <- matrix(0,length(timeseries_s),length(timeseries_s))
  tmp[m_all[[c+1]]!=0 & wanted_links_s!=0] <- 1
  links_expl[[c]] <- tmp
}

number_explained_links <- lapply(links_expl, function(x){sum(x) }) %>% unlist()
new_explained_links <- diff(number_explained_links)
new_explained_links <- c(number_explained_links[1], new_explained_links)
number_unexplained_links <- sum(wanted_links_s) - number_explained_links
number_extra_links_through_paths <- unlist(number_links)[2:21] - number_explained_links
new_extra_links <- diff(number_extra_links_through_paths)
new_extra_links <- c(number_extra_links_through_paths[1], new_extra_links)

visualize1 <- rbind(number_explained_links, rep(sum(wanted_links_s), 20)-number_explained_links)
colnames(visualize1) <- 2:21
barplot(visualize1[,1:14], col = c("slategray4","lightcoral"), xlab = "Maximal Pathlength", ylab = "Number of explained and unexplained links")

visualize2 <- rbind(number_explained_links, number_extra_links_through_paths)
colnames(visualize2) <- 2:21
barplot(visualize2[,1:14], col = c("slategray4","lightcoral"), xlab = "Maximal Pathlength", ylab = "Number of wanted and additional links")


##### ROC Curve:

# false negative rate = sum false negative / sum condition positive
# false positive rate = sum false positive / sum condition negative
# true negative rate = sum true negative / sum condition negative
# true positive rate = sum true positive / sum condition positive

# condition positiv is links without underlying motif
# resulted positive is links we get through paths
# condition negative is pairs of trade connection that CCM didn't link causally
# false negative is a unexplained link of CCM that we don't explain by paths


load(file = "./intermediate_results/ccm_res_fresh_frozen.RData")

condition_pos <- sum(wanted_links_s)
condition_neg <- res$meta_information$`Number of possible causal links` - length(links_s)

fnr <- number_unexplained_links / condition_pos
fpr <- number_extra_links_through_paths / condition_neg
tnr <- (condition_neg - number_extra_links_through_paths) / condition_neg
tpr <- number_explained_links / condition_pos

  
plot(fpr, tpr, type = "b", xlim = c(0,0.3), ylim = c(0,0.3), col = 3, xlab = "False Positive Rate", ylab = "True Positive Rate")
lines(c(0,0.3), c(0,0.3), col = "gray")

max(tpr)


### Do paths that explain causal links have a higher rho then spurious ones?

load(file = "./intermediate_results/ccm_res_fresh_frozen.RData")

m_all_expl <- lapply(m_all, function(x){x[wanted_links_s == 0] <- 0; 
                                            x[m_all[[1]] != 0] <- 0;
                                            return(x)})
m_all_spurious <- lapply(m_all, function(x){x[wanted_links_s != 0] <-0;  x[m_all[[1]] != 0] <- 0; return(x)})

# mean rho for paths of length 2 that explain a causal link: 
sum(m_all_expl[[2]]) / sum(m_all_expl[[2]]!=0)   ## 0.09272209

# mean rho for paths of length 2 that explain a causal link: 
sum(m_all_spurious[[2]]) / sum(m_all_spurious[[2]]!=0)   # 0.06590845

m_all_expl[[2]] %>% as.vector() %>% sort(decreasing = TRUE) %>% head()
### [1] 0.2760557 0.2403009 0.2181777 0.2123437 0.2079930 0.2079028


m_all_spurious[[2]] %>% as.vector() %>% sort(decreasing = TRUE) %>% head()
### [1] 0.2823791 0.2592526 0.2590621 0.2410574 0.2356279 0.2318002


######## Analysis Links for only fresh #######

load(file = "./intermediate_results/ccm_res_only_fresh.RData")

res$meta_information

links <- res$links
links_s <- links %>% filter(strong_detection == TRUE)
links_w <- links %>% filter(detection == TRUE)

# Total strong links: 1508
links_s %>% dim() %>% .[1]
# motif A -> B <- C: 20
links_s %>% filter(B == D) %>% dim() %>% .[1]
# motif A <- B -> C:  90
links_s %>% filter(A == C) %>% dim() %>% .[1]
# motif A -> B <- A: 0
links_s %>% filter(D == A & C == B) %>% dim() %>% .[1]
# motif A -> B -> C:  18
links_s %>% filter(B == C) %>% dim() %>% .[1]
# motif A <- B <- C: 27
links_s %>% filter(A == D) %>% dim() %>% .[1]
# links w/o motif: 1353
links_s %>% filter(D != A & D != B & C != A & C != B) %>% dim() %>% .[1]


# Total weak links: 3189
links_w %>% dim() %>% .[1]
# motif A -> B <- C: 69
links_w %>% filter(B == D) %>% dim() %>% .[1]
# motif A <- B -> C: 294
links_w %>% filter(A == C) %>% dim() %>% .[1]
# motif A -> B -> C:  51
links_w %>% filter(B == C) %>% dim() %>% .[1]
# motif A -> B -> A: 0
links_w %>% filter(C == B & D == A) %>% dim() %>% .[1]
# motif A <- B <- C: 60
links_w %>% filter(A == D) %>% dim() %>% .[1]
# links w/o motif: 2715
links_w %>% filter(D != A & D != B & C != A & C != B) %>% dim() %>% .[1]


######## Including Climate Dataset (Only fresh, Strong links) ########

load(file = "updated_climatedata.RData")

load(file = "./intermediate_results/ccm_res_only_fresh.RData")


links <- res$links
links_s <- links %>% filter(strong_detection == TRUE)
dim(links_s)    ## 1508  strong links


## find all feedback loops 

count <- 0
feedbacks <- list()
links_check <- links_s
for(i in 1:(dim(links_check)[1]-1)){
  links_loop <- links_check %>% filter(A == .[1,"C"]$C & B == .[1,"D"]$D & C == .[1,"A"]$A & D == .[1,"B"]$B) 
  if(dim(links_loop)[1]>0){
    if(dim(links_loop)[1]>1){print("???")}
    count <- count + 1
    feedbacks[[count]] <- links_loop %>% select(lib_column, target_column)
  }
  links_check <- links_check[-1,]
}

length(feedbacks)   ## 152 feedback loops


## get timeseries for those links

load(file = "./intermediate_results/prepared_data_fresh_frozen.RData")

### Prepare dataset for ccm
df <- dat %>% ungroup() %>% #skimr::skim()
  select(Period, link, total_netweight_tons) %>%
  group_by(Period, link) %>%
  spread(link, total_netweight_tons) %>%
  as.data.frame()

res <- lapply(feedbacks, function(x) {df %>% select(Period,x[1,1][[1]], x[1,2][[1]])})

res <- lapply(res, function(x){right_join(x, climdat)})

# calculate ccm for each listelement, and all combinations clim/trade, trade/clim and for trade/trade
res <- lapply(res, links_for_climdat, detrend = TRUE)

links <- lapply(res, function(x){return(x[[1]])})      

detection <- lapply(links, function(x){x <- x %>% filter(detection == TRUE); return(x)})

# we are interested in cases where the trade connections predict a climate index, as that means that climate is causing the trade connection
detection <- lapply(detection, function(x){x <- x %>% filter(is.na(D)); return(x)})

# and then in those that can be predicted by both of the trade connections, so that the climate index is a common driver
shared_driver <- lapply(detection, function(x){sel1 <- duplicated(x$C) %>% which()
                                                  sel2 <- duplicated(x$C, fromLast = T) %>% which()
                                                  sel <- c(sel1, sel2)
                                                  return(x[sel,])})

shared_driver_clean <- list.clean(shared_driver, function(x) dim(x)[1] == 0L)

length(shared_driver_clean)  ## 50 for strong

motifs <- lapply(shared_driver_clean, function(x){c(x$A, x$B) %>% unique() %>% length()}) %>% unlist()
sum(motifs != 4) ## 11 with motif for strong, leaving 39 between four different countries

number_climate_indices <- lapply(shared_driver_clean, function(x){return(dim(x)[1]/2)}) %>% unlist()
sum(number_climate_indices == 4) ## 38 share only one climate index as driver, 6 shaere 2, 3 share 3, 3 share 4 (strong)

load(file = "./intermediate_results/country_codes.RData")
shared_driver_clean <- lapply(shared_driver_clean, function(x){x <- codes_to_names(x, c(2,3), match = code_match)})

## remove these links from the link list created by CCM
links_rm <- lapply(shared_driver_clean, function(x){return(x$lib_column) %>% unique()}) %>% unlist()

load(file = "./intermediate_results/ccm_res_only_fresh.RData")
links <- res$links
links_s <- links %>% filter(strong_detection == TRUE)  ## 1508 links

for(i in 1:(length(links_rm)/2)){
  links_s <- links_s[!(links_s$lib_column == links_rm[2*(i-1)+1] & links_s$target_column == links_rm[2*i]),]
  links_s <- links_s[!(links_s$lib_column == links_rm[2*i] & links_s$target_column == links_rm[2*(i-1)+1]),]
}

save(links_s, file = "./intermediate_results/climate_only_fresh_strong.RData" )



######## Including Climate Dataset (Only fresh, Weak links) ########

load(file = "updated_climatedata.RData")
load(file = "./intermediate_results/ccm_res_only_fresh.RData")

links <- res$links
links_w<- links %>% filter(detection == TRUE)

dim(links_w)    ## 3189  weak links

## find all feedback loops 

count <- 0
feedbacks <- list()
links_check <- links_w            ### this segment was run twice, once with links_check <- links_w, once with links_check <- links_s
for(i in 1:(dim(links_check)[1]-1)){
  links_loop <- links_check %>% filter(A == .[1,"C"]$C & B == .[1,"D"]$D & C == .[1,"A"]$A & D == .[1,"B"]$B) 
  if(dim(links_loop)[1]>0){
    if(dim(links_loop)[1]>1){print("???")}
    count <- count + 1
    #feedbacks[[count]] <- rbind(ind2[1,], ind_loop)
    feedbacks[[count]] <- links_loop %>% select(lib_column, target_column)
  }
  links_check <- links_check[-1,]
}

length(feedbacks)   ## weak: 636 feedback lops


## get timeseries for those links

load(file = "./intermediate_results/prepared_data_fresh_frozen.RData")

### Prepare dataset for ccm
df <- dat %>% ungroup() %>% #skimr::skim()
  select(Period, link, total_netweight_tons) %>%
  group_by(Period, link) %>%
  spread(link, total_netweight_tons) %>%
  as.data.frame()


res <- lapply(feedbacks, function(x) {df %>% select(Period,x[1,1][[1]], x[1,2][[1]])})

res <- lapply(res, function(x){right_join(x, climdat)})

# calculate ccm for each listelement, and all combinations clim/trade, trade/clim and for trade/trade

res <- lapply(res, links_for_climdat, detrend = TRUE)

links <- lapply(res, function(x){return(x[[1]])})      

detection <- lapply(links, function(x){x <- x %>% filter(detection == TRUE); return(x)})

# we are interested in cases where the trade connections predict a climate index, as that means that climate is causing the trade connection
detection <- lapply(detection, function(x){x <- x %>% filter(is.na(D)); return(x)})

# and then in those that can be predicted by both of the trade connections, so that the climate index is a common driver
shared_driver <- lapply(detection, function(x){sel1 <- duplicated(x$C) %>% which()
sel2 <- duplicated(x$C, fromLast = T) %>% which()
sel <- c(sel1, sel2) 
return(x[sel,])})

shared_driver_clean <- list.clean(shared_driver, function(x) dim(x)[1] == 0L)

length(shared_driver_clean)   ## 222 causal links have at least one climate index as shared driver and hence might 
                              ## not be causally related after all 

motifs <- lapply(shared_driver_clean, function(x){c(x$A, x$B) %>% unique() %>% length()}) %>% unlist()
sum(motifs == 4) ## 57 show some sort of motif (only 3 involved countries), while 165 are between four different countries (weak)

number_climate_indices <- lapply(shared_driver_clean, function(x){return(dim(x)[1]/2)}) %>% unlist()
sum(number_climate_indices == 5) ## 162 share only one climate index as driver, 39 share 2, 16 share 3, 5 share 4  and 1 shares 5 (weak)


load(file = "./intermediate_results/country_codes.RData")
shared_driver_clean <- lapply(shared_driver_clean, function(x){x <- codes_to_names(x, c(2,3), match = code_match)})


## remove these links from the link list created by CCM
links_rm <- lapply(shared_driver_clean, function(x){return(x$lib_column) %>% unique()}) %>% unlist()

load(file = "./intermediate_results/ccm_res_fresh_frozen.RData")
links <- res$links
links_w <- links %>% filter(detection == TRUE)  ## 2481 links

for(i in 1:(length(links_rm)/2)){
  links_w <- links_w[!(links_w$lib_column == links_rm[2*(i-1)+1] & links_w$target_column == links_rm[2*i]),]
  links_w <- links_w[!(links_w$lib_column == links_rm[2*i] & links_w$target_column == links_rm[2*(i-1)+1]),]
}


save(links_w, file = "./intermediate_results/climate_only_fresh_weak.RData" )





######## Links needing explanation for later comparison (Only fresh) #######

load(file = "./intermediate_results/climate_only_fresh_strong.RData")
load(file = "./intermediate_results/climate_only_fresh_weak.RData")


links_s_wo_motif <- links_s %>% filter(B != D  & A != C & B != C & A != D) %>% select(lib_column, target_column, rho) 
timeseries_s <- c(links_s$lib_column, links_s$target_column) %>% unique
l <- length(timeseries_s)

wanted_links_s <- matrix(0, nrow = l, ncol = l)
colnames(wanted_links_s) <- timeseries_s
rownames(wanted_links_s) <- timeseries_s
for(i in 1:(dim(links_s_wo_motif)[1])){
  wanted_links_s[links_s_wo_motif$lib_column[i], links_s_wo_motif$target_column[i]] <- 1
}


links_w_wo_motif <- links_w %>% filter(B != D  & A != C & B != C & A != D) %>% select(lib_column, target_column, rho) 
timeseries_w <- c(links_w$lib_column, links_w$target_column) %>% unique
l <- length(timeseries_w)

wanted_links_w <- matrix(0, nrow = l, ncol = l)
colnames(wanted_links_w) <- timeseries_w
rownames(wanted_links_w) <- timeseries_w
for(i in 1:(dim(links_w_wo_motif)[1])){
  wanted_links_w[links_w_wo_motif$lib_column[i], links_w_wo_motif$target_column[i]] <- 1
}

save(wanted_links_s, wanted_links_w, timeseries_w, timeseries_s, file = "./intermediate_results/need_explanation_only_fresh.RData" )


######## Trasitivety Analysis (Only fresh salmon, Strong links) ########

load(file = "./intermediate_results/climate_only_fresh_strong.RData")
load(file = "./intermediate_results/need_explanation_only_fresh.RData")   

dim(links_s)    ## 1408 strong links left

sum(wanted_links_s)  # 1275
links_s_motif <- links_s %>% filter(B == D  | A == C | B == C | A == D) %>% select(lib_column, target_column, rho)
dim(links_s_motif)  ## 133 strong links with motif

m1 <- matrix(0,nrow = length(timeseries_s), ncol = length(timeseries_s))
colnames(m1) <- timeseries_s
rownames(m1) <- timeseries_s
for(i in 1:(dim(links_s_motif)[1])){
  m1[links_s_motif$lib_column[i], links_s_motif$target_column[i]] <- links_s_motif$rho[i]
}

m_all <- list()
null_all <- list()
m_all[[1]] <- m1
null_all[[1]] <- m1
null_all[[1]][m1 != 0] <- 1
m <- list()
n <- list()
n[[1]] <- null_all[[1]]
number_links <- list()
added_links <- list()
m[[1]] <- m1
number_links[[1]] <- sum(m1!=0) - dim(m1)[1] 
for(i in 1:20){
  n[[i+1]] <- n[[i]] %*% n[[1]]
  m[[i+1]] <- m[[i]] %*% m1
  m_all[[i+1]] <- m_all[[i]]
  m_all[[i+1]][m_all[[i]] == 0] <- m[[i+1]][m_all[[i]] == 0]
  null_all[[i+1]] <- null_all[[i]]
  null_all[[i+1]][null_all[[i]] == 0] <- n[[i+1]][null_all[[i]] == 0]
  m_all[[i+1]][null_all[[i]] == 0 & null_all[[i+1]] != 0] <- m_all[[i+1]][null_all[[i]] == 0 & null_all[[i+1]] != 0]/
    null_all[[i+1]][null_all[[i]] == 0 & null_all[[i+1]] != 0]
  number_links[[i+1]] <- sum(m_all[[i+1]]!=0)
  added_links[[i]] <- number_links[[i+1]] - number_links[[i]]
}


#### Are the new causal links part of the ones we couldn't explain?

links_expl <- list()
for(c in 1:20){
  tmp <- matrix(0,length(timeseries_s),length(timeseries_s))
  tmp[m_all[[c+1]]!=0 & wanted_links_s!=0] <- 1
  links_expl[[c]] <- tmp
}

number_explained_links <- lapply(links_expl, function(x){sum(x) }) %>% unlist()
new_explained_links <- diff(number_explained_links)
new_explained_links <- c(number_explained_links[1], new_explained_links)
number_unexplained_links <- sum(wanted_links_w) - number_explained_links
number_extra_links_through_paths <- unlist(number_links)[2:21] - number_explained_links
new_extra_links <- diff(number_extra_links_through_paths)
new_extra_links <- c(number_extra_links_through_paths[1], new_extra_links)

visualize1 <- rbind(number_explained_links, rep(2062, 20)-number_explained_links)
colnames(visualize1) <- 2:21
barplot(visualize1[,1:14], col = c("slategray4","lightcoral"), xlab = "Maximal Pathlength", ylab = "Number of explained and unexplained links")

visualize2 <- rbind(number_explained_links, number_extra_links_through_paths)
colnames(visualize2) <- 2:21
barplot(visualize2[,1:14], col = c("slategray4","lightcoral"), xlab = "Maximal Pathlength", ylab = "Number of wanted and additional links")


### Do paths that explain causal links have a higher rho then spurious ones?

m_all_expl <- lapply(m_all, function(x){x[wanted_links_s == 0] <- 0; 
x[m_all[[1]] != 0] <- 0;
return(x)})
m_all_spurious <- lapply(m_all, function(x){x[wanted_links_s != 0] <-0;  x[m_all[[1]] != 0] <- 0; return(x)})


# mean rho for paths of length 2 that explain a causal link: 
sum(m_all_expl[[2]]) / sum(m_all_expl[[2]]!=0)   # 0.1059094   

# mean rho for paths of length 2 that explain a causal link: 
sum(m_all_spurious[[2]]) / sum(m_all_spurious[[2]]!=0)   # 0.07200075

m_all_expl[[2]] %>% as.vector() %>% sort(decreasing = TRUE) %>% head()
# [1] 0.30694391 0.16840620 0.11674479 0.09952947 0.09740515 0.07657036

m_all_spurious[[2]] %>% as.vector() %>% sort(decreasing = TRUE) %>% head()
# [1] 0.2163682 0.1866366 0.1721670 0.1707119 0.1654798 0.1640326




######## Getting country dependencies from A->B<-C motif ########


load(file = "./intermediate_results/ccm_res_fresh_frozen.RData")
load(file = "./intermediate_results/country_codes.RData")

links <- res$links

links_w <- links %>% filter(detection == T)
links_s <- links_w %>% filter(strong_detection == T)

motif_s <- links_s %>% filter(B == D) %>% codes_to_names(match = code_match)
motif_w <- links_w %>% filter(B == D) %>% codes_to_names(match = code_match)

# Visualize strong links
involved_countries <- c(motif_s$A, motif_s$C) %>% unique()
m_s <- matrix(NA, length(involved_countries), length(involved_countries))
colnames(m_s) <- involved_countries
rownames(m_s) <- involved_countries
for(i in 1:(dim(motif_s)[1])){
  m_s[motif_s$A[i], motif_s$C[i]] <- motif_s$rho[i]
}

colnames(m_s)[colnames(m_s) == "Russian Federation"] <- "Russia"
colnames(m_s)[colnames(m_s) == "United States of America"] <- "USA"
colnames(m_s)[colnames(m_s) == "United Kingdom"] <- "UK"
rownames(m_s) <- colnames(m_s)
levelplot(m_s, scales=list(x=list(rot=40), alternating = 2), xlab ="", ylab = "", 
          #     main = "Strong links between salmon capture data of countries", 
          col.regions = heat.colors(100)[80:1],at=seq(0, 1, length.out= 80))


# Visualize weak links
involved_countries <- c(motif_w$A, motif_w$C) %>% unique()
m_w <- matrix(NA, length(involved_countries), length(involved_countries))
colnames(m_w) <- involved_countries
rownames(m_w) <- involved_countries
for(i in 1:(dim(motif_w)[1])){
  m_w[motif_w$A[i], motif_w$C[i]] <- motif_w$rho[i]
}

colnames(m_w)[colnames(m_w) == "Russian Federation"] <- "Russia"
colnames(m_w)[colnames(m_w) == "United States of America"] <- "USA"
colnames(m_w)[colnames(m_w) == "United Kingdom"] <- "UK"
rownames(m_w) <- colnames(m_w)
levelplot(m_w, scales=list(x=list(rot=40), alternating = 2), xlab ="", ylab = "", 
          #     main = "Weak links between salmon capture data of countries", 
          col.regions = heat.colors(100)[80:1],at=seq(0, 1, length.out= 80))










######## Preparation and CCM banana data ########

data <- read.csv(file = "comtrade_banana_prepared.csv")    ## preprocessed data

data <- data %>%
  ungroup() %>%
  mutate(Period = parse_date_time(Period, orders = "Ym"))

data$Period %>% unique() %>% length() # 113 time points!
data$link %>% unique() %>% length() # 2697 links

## Delete Bunkers, that's not a country. For more info: http://www.fao.org/faostat/en/#data/EM/metadata
data <- filter(data, importer_name != "Bunkers")
data <- filter(data, importer_name != "Ic")
data <- filter(data, exporter_name != "Ic")

data <- data %>% 
  select(-starts_with("diff"), -starts_with("Commodity")) %>%
  group_by(link, Period) %>%
  select(-n_obs) %>%
  mutate(
    total_netweight_tons = sum(max_netweight)/1000,
    total_value_us = sum(max_value_us)
  ) %>%
  select(-starts_with("max")) %>% unique() %>% ungroup()

dat <- data %>%
  filter(Year > 2009 & Year < 2018)

save(dat, file = "./intermediate_results/prepared_data_bananas.RData")

res <- get_links(dat, min_obs = 50, detrend = TRUE, nonlin = 0.05)

save(res, file = "./intermediate_results/ccm_bananas.RData")

######## Analysis of links in banana dataset ########

load(file = "./intermediate_results/ccm_bananas.RData")

res$meta_information

links <- res$links
links_s <- links %>% filter(strong_detection == TRUE)
links_w <- links %>% filter(detection == TRUE)

# Total strong links: 92
links_s %>% dim() %>% .[1]
# motif A -> B <- C: 5
links_s %>% filter(B == D) %>% dim() %>% .[1]
# motif A <- B -> C: 5
links_s %>% filter(A == C) %>% dim() %>% .[1]
# motif A -> B -> C: 2
links_s %>% filter(B == C) %>% dim() %>% .[1]
# motif A -> B <- A: 0
links_s %>% filter(D == A & C == B) %>% dim() %>% .[1]  
# motif A <- B <- C: 4
links_s %>% filter(A == D) %>% dim() %>% .[1]
# links w/o motif: 76
links_s %>% filter(D != A & D != B & C != A & C != B) %>% dim() %>% .[1]  

# Total weak links: 318
links_w %>% dim() %>% .[1]
# motif A -> B <- C: 12
links_w %>% filter(B == D) %>% dim() %>% .[1]
# motif A <- B -> C: 27
links_w %>% filter(A == C) %>% dim() %>% .[1] 
# motif A -> B -> C: 7
links_w %>% filter(B == C) %>% dim() %>% .[1]
# motif A -> B -> A: 1
links_w %>% filter(C == B & D == A) %>% dim() %>% .[1] 
# motif A <- B <- C: 12
links_w %>% filter(A == D) %>% dim() %>% .[1]
# links w/o motif: 261
links_w %>% filter(D != A & D != B & C != A & C != B) %>% dim() %>% .[1]  




######## Including Climate Dataset for banana (Strong links) ########

load(file = "updated_climatedata.RData")

load(file = "./intermediate_results/ccm_bananas.RData")

links <- res$links
links_s <- links %>% filter(strong_detection == TRUE)
dim(links_s)    ## 92  strong links


## find all feedback loops 

count <- 0
feedbacks <- list()
links_check <- links_s
for(i in 1:(dim(links_check)[1]-1)){
  links_loop <- links_check %>% filter(A == .[1,"C"]$C & B == .[1,"D"]$D & C == .[1,"A"]$A & D == .[1,"B"]$B) 
  if(dim(links_loop)[1]>0){
    if(dim(links_loop)[1]>1){print("???")}
    count <- count + 1
    feedbacks[[count]] <- links_loop %>% select(lib_column, target_column)
  }
  links_check <- links_check[-1,]
}

length(feedbacks)   ## 1 feedback loops


## get timeseries for those links

load(file = "./intermediate_results/prepared_data_bananas.RData")


### Prepare dataset for ccm
df <- dat %>% ungroup() %>% #skimr::skim()
  select(Period, link, total_netweight_tons) %>%
  group_by(Period, link) %>%
  spread(link, total_netweight_tons) %>%
  as.data.frame()

res <- lapply(feedbacks, function(x) {df %>% select(Period,x[1,1][[1]], x[1,2][[1]])})

res <- lapply(res, function(x){right_join(x, climdat)})

# calculate ccm for each listelement, and all combinations clim/trade, trade/clim and for trade/trade
res <- lapply(res, links_for_climdat, detrend = TRUE)

links <- lapply(res, function(x){return(x[[1]])})      

detection <- lapply(links, function(x){x <- x %>% filter(detection == TRUE); return(x)})

# we are interested in cases where the trade connections predict a climate index, as that means that climate is causing the trade connection
detection <- lapply(detection, function(x){x <- x %>% filter(is.na(D)); return(x)})

# and then in those that can be predicted by both of the trade connections, so that the climate index is a common driver
shared_driver <- lapply(detection, function(x){sel1 <- duplicated(x$C) %>% which()
                                                sel2 <- duplicated(x$C, fromLast = T) %>% which()
                                                sel <- c(sel1, sel2)
                                                return(x[sel,])})

shared_driver_clean <- list.clean(shared_driver, function(x) dim(x)[1] == 0L)

length(shared_driver_clean)  ## 1 for strong

load(file = "./intermediate_results/country_codes.RData")
shared_driver_clean <- lapply(shared_driver_clean, function(x){x <- codes_to_names(x, c(2,3), match = code_match)})

## remove these links from the link list created by CCM
links_rm <- lapply(shared_driver_clean, function(x){return(x$lib_column) %>% unique()}) %>% unlist()

load(file = "./intermediate_results/ccm_bananas.RData")
links <- res$links
links_s <- links %>% filter(strong_detection == TRUE)  ## 1508 links

for(i in 1:(length(links_rm)/2)){
  links_s <- links_s[!(links_s$lib_column == links_rm[2*(i-1)+1] & links_s$target_column == links_rm[2*i]),]
  links_s <- links_s[!(links_s$lib_column == links_rm[2*i] & links_s$target_column == links_rm[2*(i-1)+1]),]
}

save(links_s, file = "./intermediate_results/climate_bananas_strong.RData" )






######## Including Climate Dataset for banana (Weak links) ########

load(file = "updated_climatedata.RData")
load(file = "./intermediate_results/ccm_bananas.RData")

links <- res$links
links_w<- links %>% filter(detection == TRUE)

dim(links_w)    ## 318  weak links

## find all feedback loops 

count <- 0
feedbacks <- list()
links_check <- links_w            ### this segment was run twice, once with links_check <- links_w, once with links_check <- links_s
for(i in 1:(dim(links_check)[1]-1)){
  links_loop <- links_check %>% filter(A == .[1,"C"]$C & B == .[1,"D"]$D & C == .[1,"A"]$A & D == .[1,"B"]$B) 
  if(dim(links_loop)[1]>0){
    if(dim(links_loop)[1]>1){print("???")}
    count <- count + 1
    #feedbacks[[count]] <- rbind(ind2[1,], ind_loop)
    feedbacks[[count]] <- links_loop %>% select(lib_column, target_column)
  }
  links_check <- links_check[-1,]
}

length(feedbacks)   ## weak: 24 feedback lops


## get timeseries for those links

load(file = "./intermediate_results/prepared_data_bananas.RData")

dat <- dat %>% 
  group_by(link) %>%
  mutate(log_weight = log1p(total_netweight_tons)) %>%
  mutate(norm_weight = (total_netweight_tons - mean(total_netweight_tons, na.rm = TRUE)) / sd(total_netweight_tons, na.rm = TRUE),
         norm_log_weight = (log_weight - mean(log_weight, na.rm = TRUE)) / sd(log_weight, na.rm = TRUE))

### Prepare dataset for ccm
df <- dat %>% ungroup() %>% #skimr::skim()
  select(Period, link, norm_weight) %>%
  group_by(Period, link) %>%
  spread(link, norm_weight) %>%
  as.data.frame()

res <- lapply(feedbacks, function(x) {df %>% select(Period,x[1,1][[1]], x[1,2][[1]])})

res <- lapply(res, function(x){right_join(x, climdat)})

# calculate ccm for each listelement, and all combinations clim/trade, trade/clim and for trade/trade

res <- lapply(res, links_for_climdat, detrend = TRUE)

links <- lapply(res, function(x){return(x[[1]])})      

detection <- lapply(links, function(x){x <- x %>% filter(detection == TRUE); return(x)})

# we are interested in cases where the trade connections predict a climate index, as that means that climate is causing the trade connection
detection <- lapply(detection, function(x){x <- x %>% filter(is.na(D)); return(x)})

# and then in those that can be predicted by both of the trade connections, so that the climate index is a common driver
shared_driver <- lapply(detection, function(x){sel1 <- duplicated(x$C) %>% which()
                                                sel2 <- duplicated(x$C, fromLast = T) %>% which()
                                                sel <- c(sel1, sel2) 
                                                return(x[sel,])})

shared_driver_clean <- list.clean(shared_driver, function(x) dim(x)[1] == 0L)

length(shared_driver_clean)   ## 6 causal links have at least one climate index as shared driver

load(file = "./intermediate_results/country_codes.RData")
shared_driver_clean <- lapply(shared_driver_clean, function(x){x <- codes_to_names(x, c(2,3), match = code_match)})


## remove these links from the link list created by CCM
links_rm <- lapply(shared_driver_clean, function(x){return(x$lib_column) %>% unique()}) %>% unlist()

load(file = "./intermediate_results/ccm_bananas.RData")
links <- res$links
links_w <- links %>% filter(detection == TRUE)  ## 2481 links

for(i in 1:(length(links_rm)/2)){
  links_w <- links_w[!(links_w$lib_column == links_rm[2*(i-1)+1] & links_w$target_column == links_rm[2*i]),]
  links_w <- links_w[!(links_w$lib_column == links_rm[2*i] & links_w$target_column == links_rm[2*(i-1)+1]),]
}


save(links_w, file = "./intermediate_results/climate_bananas_weak.RData" )





######## Links needing explanation for later comparison (bananas) ########

load(file = "./intermediate_results/climate_bananas_strong.RData")
load(file = "./intermediate_results/climate_bananas_weak.RData")

links_s_wo_motif <- links_s %>% filter(B != D  & A != C & B != C & A != D) %>% select(lib_column, target_column, rho) 
timeseries_s <- c(links_s$lib_column, links_s$target_column) %>% unique
l <- length(timeseries_s)

wanted_links_s <- matrix(0, nrow = l, ncol = l)
colnames(wanted_links_s) <- timeseries_s
rownames(wanted_links_s) <- timeseries_s
for(i in 1:(dim(links_s_wo_motif)[1])){
  wanted_links_s[links_s_wo_motif$lib_column[i], links_s_wo_motif$target_column[i]] <- 1
}


links_w_wo_motif <- links_w %>% filter(B != D  & A != C & B != C & A != D) %>% select(lib_column, target_column, rho) 
timeseries_w <- c(links_w$lib_column, links_w$target_column) %>% unique
l <- length(timeseries_w)

wanted_links_w <- matrix(0, nrow = l, ncol = l)
colnames(wanted_links_w) <- timeseries_w
rownames(wanted_links_w) <- timeseries_w
for(i in 1:(dim(links_w_wo_motif)[1])){
  wanted_links_w[links_w_wo_motif$lib_column[i], links_w_wo_motif$target_column[i]] <- 1
}


save(wanted_links_s, wanted_links_w, timeseries_w, timeseries_s, file = "./intermediate_results/need_explanation_bananas.RData" )



######## Trasitivety Analysis (Banana Dataset, Weak links) ########

# for strong links no new connections are explained

load(file = "./intermediate_results/climate_bananas_weak.RData")
load(file = "./intermediate_results/need_explanation_bananas.RData")   

dim(links_w)    ## 306 weak links left

sum(wanted_links_w)  # 249

links_w_motif <- links_w %>% filter(B == D  | A == C | B == C | A == D) %>% select(lib_column, target_column, rho)
dim(links_w_motif)  ## 57 strong links with motif

m1 <- matrix(0,nrow = length(timeseries_w), ncol = length(timeseries_w))
colnames(m1) <- timeseries_w
rownames(m1) <- timeseries_w
for(i in 1:(dim(links_w_motif)[1])){
  m1[links_w_motif$lib_column[i], links_w_motif$target_column[i]] <- links_w_motif$rho[i]
}

m_all <- list()
null_all <- list()
m_all[[1]] <- m1
null_all[[1]] <- m1
null_all[[1]][m1 != 0] <- 1
m <- list()
n <- list()
n[[1]] <- null_all[[1]]
number_links <- list()
added_links <- list()
m[[1]] <- m1
number_links[[1]] <- sum(m1!=0) - dim(m1)[1] 
for(i in 1:20){
  n[[i+1]] <- n[[i]] %*% n[[1]]
  m[[i+1]] <- m[[i]] %*% m1
  m_all[[i+1]] <- m_all[[i]]
  m_all[[i+1]][m_all[[i]] == 0] <- m[[i+1]][m_all[[i]] == 0]
  null_all[[i+1]] <- null_all[[i]]
  null_all[[i+1]][null_all[[i]] == 0] <- n[[i+1]][null_all[[i]] == 0]
  # average rho for new links that can be reached by multiple paths of length i+1
  m_all[[i+1]][null_all[[i]] == 0 & null_all[[i+1]] != 0] <- m_all[[i+1]][null_all[[i]] == 0 & null_all[[i+1]] != 0]/
    null_all[[i+1]][null_all[[i]] == 0 & null_all[[i+1]] != 0]
  number_links[[i+1]] <- sum(m_all[[i+1]]!=0)
  added_links[[i]] <- number_links[[i+1]] - number_links[[i]]
}


#### Are the new causal links part of the ones we couldn't explain?

links_expl <- list()
for(c in 1:20){
  tmp <- matrix(0,length(timeseries_w),length(timeseries_w))
  tmp[m_all[[c+1]]!=0 & wanted_links_w!=0] <- 1
  links_expl[[c]] <- tmp
}

number_explained_links <- lapply(links_expl, function(x){sum(x) }) %>% unlist()
new_explained_links <- diff(number_explained_links)
new_explained_links <- c(number_explained_links[1], new_explained_links)
number_unexplained_links <- sum(wanted_links_w) - number_explained_links
number_extra_links_through_paths <- unlist(number_links)[2:21] - number_explained_links
new_extra_links <- diff(number_extra_links_through_paths)
new_extra_links <- c(number_extra_links_through_paths[1], new_extra_links)

visualize1 <- rbind(number_explained_links, rep(sum(wanted_links_w), 20)-number_explained_links)
colnames(visualize1) <- 2:21
barplot(visualize1[,1:14], col = c("slategray4","lightcoral"), xlab = "Maximal Pathlength", ylab = "Number of explained and unexplained links")

visualize2 <- rbind(number_explained_links, number_extra_links_through_paths)
colnames(visualize2) <- 2:21
barplot(visualize2[,1:14], col = c("slategray4","lightcoral"), xlab = "Maximal Pathlength", ylab = "Number of wanted and additional links")




######## Preparation FAO dataset ########

fao <- read.csv(file = "TS_FI_CAPTURE.csv")
code_match_countries <- read.csv(file = "CL_FI_COUNTRY_GROUPS.csv") %>% select(Code = UN_Code, Country = Name_En)
code_match_species <- read.csv(file = "CL_FI_SPECIES_GROUPS.csv") %>% select(Code = X3Alpha_Code, Species = Name_En) %>%
                        filter(Code %in% c("SAL", "PIN", "CHU", "ONC", "CHE", "SOC", "CHI", "COH", "ORC", "SWI"))

fao <- fao %>%  select(-SYMBOL) %>% filter(SPECIES %in% c("SAL", "PIN", "CHU", "ONC", "CHE", "SOC", "CHI", "COH", "ORC", "SWI")) 
## only 28 countries cpture salmon according to FAO
fao_countries <- fao %>% select(COUNTRY) %>% unique() %>% .[,1]

save(fao_countries, file = "./intermediate_results/fao_countries.RData" )

code_match_countries[,2] <- as.character(code_match_countries[,2])
code_match_species[,2] <- as.character(code_match_species[,2])
code_match_species[,1] <- as.character(code_match_species[,1])

fao_prep <- fao %>% select(-FISHING_AREA, -UNIT, -SPECIES)
fao_prep <- fao_prep %>% group_by(COUNTRY, YEAR) %>% mutate(total_catch = sum(QUANTITY)) %>% select(-QUANTITY) %>% unique() %>% ungroup()

# 833: Isle of Man reports same amount each year and messes up CC
fao_prep <- fao_prep %>% filter(COUNTRY != 833) 

save(fao_prep, file = "./intermediate_results/fao_salmon_data_prepared.RData")
save(code_match_countries,code_match_species, file = "./intermediate_results/fao_code_match.RData")

######## CCM on FAO Dataset #########

load(file = "./intermediate_results/fao_salmon_data_prepared.RData")

res <- get_links_fao(fao_prep, min_obs = 40, detrend = detrend, nonlin = 0.01)

save(res, file = "./intermediate_results/ccm_fao.RData")

######## Visualization of FAO causal links ########

load(file = "./intermediate_results/ccm_fao.RData")
load(file = "./intermediate_results/fao_code_match.RData")

links <- res$links

links_s <- links %>% filter(strong_detection == T)   ## 10 strong links
links_w <- links %>% filter(detection == T)    ## 12 weak links

links_s <- codes_to_names(links_s, cols = c(1,2), match = code_match_countries)
links_w <- codes_to_names(links_w, cols = c(1,2), match = code_match_countries)

# Visualize strong links
involved_countries <- c(links_s$lib_column, links_s$target_column) %>% unique()
m_s <- matrix(NA, length(involved_countries), length(involved_countries))
colnames(m_s) <- involved_countries
rownames(m_s) <- involved_countries
for(i in 1:(dim(links_s)[1])){
  m_s[links_s$lib_column[i], links_s$target_column[i]] <- links_s$rho[i]
}

colnames(m_s)[colnames(m_s) == "United Kingdom"] <- "UK"
rownames(m_s) <- colnames(m_s)
levelplot(m_s, scales=list(x=list(rot=40), alternating = 2), xlab ="", ylab = "", 
          #     main = "Strong links between salmon capture data of countries", 
          col.regions = heat.colors(100)[80:1],at=seq(0, 1, length.out= 80))

# Visualize weak links
involved_countries <- c(links_w$lib_column, links_w$target_column) %>% unique()
m_w <- matrix(NA, length(involved_countries), length(involved_countries))
colnames(m_w) <- involved_countries
rownames(m_w) <- involved_countries
for(i in 1:(dim(links_w)[1])){
  m_w[links_w$lib_column[i], links_w$target_column[i]] <- links_w$rho[i]
}

colnames(m_w)[colnames(m_w) == "United Kingdom"] <- "UK"
colnames(m_w)[colnames(m_w) == "United States of America"] <- "USA"
rownames(m_w) <- colnames(m_w)
levelplot(m_w, scales=list(x=list(rot=40),alternating=2), xlab ="", ylab = "", 
          #   main = "Weak links between salmon capture data of countries", 
          col.regions = heat.colors(100)[80:1],at=seq(0, 1, length.out= 80))






######## Reducing trade network to 28 FAO countries as exporter ########

load(file = "./intermediate_results/prepared_data_fresh_frozen.RData")
load(file = "./intermediate_results/fao_countries.RData")

dat <- dat %>% filter(exporter %in% fao_countries)
# reduces dataset from 69813 to 36847 links

save(dat, file = "./intermediate_results/prepared_data_fresh_frozen_fao_countries.RData")

# res <- get_links(dat, min_obs = 60)
# save(res, file = "ccm_fresh_frozen_fao_countries.RData")
load(file = "./intermediate_results/ccm_fresh_frozen_fao_countries.RData")

links <- res$links

links_s <- links %>% filter(strong_detection == T)
links_w <- links %>% filter(detection == T)

# Total strong links: 490
links_s %>% dim() %>% .[1]
# motif A -> B <- C: 4
links_s %>% filter(B == D) %>% dim() %>% .[1]
# motif A <- B -> C: 24
links_s %>% filter(A == C) %>% dim() %>% .[1]
# motif A -> B -> C: 11
links_s %>% filter(B == C) %>% dim() %>% .[1]
# motif A -> B <- A: 0
links_s %>% filter(D == A & C == B) %>% dim() %>% .[1]  
# motif A <- B <- C: 9
links_s %>% filter(A == D) %>% dim() %>% .[1]
# links w/o motif: 442
links_s %>% filter(D != A & D != B & C != A & C != B) %>% dim() %>% .[1]  

# Total weak links: 1134
links_w %>% dim() %>% .[1]
# motif A -> B <- C: 23
links_w %>% filter(B == D) %>% dim() %>% .[1]
# motif A <- B -> C: 68
links_w %>% filter(A == C) %>% dim() %>% .[1] 
# motif A -> B -> C: 22
links_w %>% filter(B == C) %>% dim() %>% .[1]
# motif A -> B -> A: 0
links_w %>% filter(C == B & D == A) %>% dim() %>% .[1] 
# motif A <- B <- C: 16
links_w %>% filter(A == D) %>% dim() %>% .[1]
# links w/o motif: 1005
links_w %>% filter(D != A & D != B & C != A & C != B) %>% dim() %>% .[1]  


load(file = "./intermediate_results/country_codes.RData")

motif_s <- links_s %>% filter(B == D) %>% codes_to_names(match = code_match)
motif_w <- links_w %>% filter(B == D) %>% codes_to_names(match = code_match)

# Visualize strong links
involved_countries <- c(motif_s$A, motif_s$C) %>% unique()
m_s <- matrix(NA, length(involved_countries), length(involved_countries))
colnames(m_s) <- involved_countries
rownames(m_s) <- involved_countries
for(i in 1:(dim(motif_s)[1])){
  m_s[motif_s$A[i], motif_s$C[i]] <- motif_s$rho[i]
}

colnames(m_s)[colnames(m_s) == "Russian Federation"] <- "Russia"
colnames(m_s)[colnames(m_s) == "United States of America"] <- "USA"
colnames(m_s)[colnames(m_s) == "United Kingdom"] <- "UK"
rownames(m_s) <- colnames(m_s)
levelplot(m_s, scales=list(x=list(rot=40), alternating = 2), xlab ="", ylab = "", 
          #     main = "Strong links between salmon capture data of countries", 
          col.regions = heat.colors(100)[80:1],at=seq(0, 1, length.out= 80))


# Visualize weak links
involved_countries <- c(motif_w$A, motif_w$C) %>% unique()
m_w <- matrix(NA, length(involved_countries), length(involved_countries))
colnames(m_w) <- involved_countries
rownames(m_w) <- involved_countries
for(i in 1:(dim(motif_w)[1])){
  m_w[motif_w$A[i], motif_w$C[i]] <- motif_w$rho[i]
}

colnames(m_w)[colnames(m_w) == "Russian Federation"] <- "Russia"
colnames(m_w)[colnames(m_w) == "United States of America"] <- "USA"
colnames(m_w)[colnames(m_w) == "United Kingdom"] <- "UK"
rownames(m_w) <- colnames(m_w)
levelplot(m_w, scales=list(x=list(rot=40), alternating = 2), xlab ="", ylab = "", 
          #     main = "Weak links between salmon capture data of countries", 
          col.regions = heat.colors(100)[80:1],at=seq(0, 1, length.out= 80))






