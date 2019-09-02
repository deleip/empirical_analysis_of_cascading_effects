
######## Add timelags ########

add_timelags <- function(data, columns , number_of_lags, cut = TRUE){
  if(number_of_lags < 1){
    return(data)
  }
  names <- colnames(data) 
  data <- as.data.frame(data)
  d <- dim(data)[[1]]
  if(dim(data)[[2]] == 1){
    columns <- 1
  }
  else if(!is.numeric(columns)){
    columns <- match(columns, colnames(data))
  }
  columns <- sort(columns, decreasing = TRUE)
  for(i in columns){
    lags <- matrix(NA, d, number_of_lags)
    colnames(lags) <- 1:number_of_lags
    for(j in 1:number_of_lags){
      lags[-(1:j),j] <- data[1:(d-j), i]
      colnames(lags)[j] <- paste(names[i], j, sep= "")
    }
    if(i== dim(data)[[2]]){
      data <- cbind(data[,1:i],lags)
    }
    else{
      data <- cbind(data[,1:i],lags,data[,(i+1):(dim(data)[[2]])])
    }
  }
  if(as.logical(cut)){
    data <- data[(number_of_lags+1):d,]
  }
  return(data)
}


######## Replace country codes by country names ########

codes_to_names <- function(data, cols = c(2,3,5,6), match = code_match){
  for(j in cols){
    tmp <- colnames(data)[j]
    colnames(data)[j] <- "Country"
    for(i in 1:(dim(data)[1])){
      data[i,j] <- match[match$Code == data[i,j][[1]], 2]
    }
    colnames(data)[j] <- tmp
  }
  return(data)
}
######## Function to apply CCM on trade data ########
get_links <- function(data, meta = TRUE, min_obs = 60, detrend = TRUE, nonlin = 0.1){
  
  meta_information <- list()
  metanames <- c()
  
  if(meta == TRUE){
    meta_information[[length(meta_information)+1]] <- data %>% select(link) %>% distinct() %>% dim() %>% .[[1]]
    metanames <- c(metanames, "Total trade connections in the begining")
  }
  
  
  ## There is too little observations for some links, filter and drop:
  if(meta == TRUE){
    drop <- data %>% group_by(link) %>% mutate(n_obs = n()) %>% ungroup () %>% select(exporter_name, importer_name, n_obs) %>% distinct()
    drop <- sum(drop$n_obs<=min_obs)
    
    meta_information[[length(meta_information)+1]] <- drop
    metanames <- c(metanames, "Trade connections dropped due to too few observations")
  }
  
  data <- data %>%
    group_by(link) %>%
    mutate(n_obs = n()) %>%
    filter(n_obs > min_obs)
  
  
  ### Prepare dataset for ccm
  df <- data %>% ungroup() %>% #skimr::skim()
    select(Period, link, total_netweight_tons) %>%
    group_by(Period, link) %>%
    spread(link, total_netweight_tons) %>%
    as.data.frame()
  
  
  
  period <- as.data.frame(df[,1])
  df <- df[,-1]
  
  # detrending
  if(detrend == TRUE){
    df <- cbind(apply(df, 2, function(x){l <- 1:length(x);
    lm_res <- lm(x~l);
    pred <- lm_res$coefficients[[1]] + l *lm_res$coefficients[[2]];
    x <- x - pred;
    return(x)}))
  }
  
  
  # Normalizing
  df <- cbind(apply(df, 2, function(x) {x <- (x-mean(x, na.rm = TRUE))/sd(x, na.rm = T); return(x)}))
  
  
  
  
  df <- cbind(period, df)
  
  df_return <- df
  colnames(df)[1] <- "Period"
  
  
  ## get best embedding dimensions
  emb <- list()
  for (i in 2:dim(df)[2]){
    emb[[i-1]] <- df[c(1,i)] %>%
      simplex()
    emb[[i-1]]$link <- colnames(df)[i]
  }
  
  # best embedding dimension:
  bestE <- emb %>% map_dbl(~ .$E[which.max(.$rho)])
  bestE_df <- data_frame(
    bestE = bestE,
    link = colnames(df)[-1]
  )
  
  
  ## Nonlinearity test
  nonlinearity <- list()
  for (i in 2:dim(df)[2]){
    nonlinearity[[i-1]] <- df[c(1,i)] %>% s_map(., E = bestE[i-1])
  }
  
  nonlinearity_return <- nonlinearity
  
  
  # check nonlinearity
  linear <- map(.x = nonlinearity, .f = ~(max(.x$rho)<(.x$rho[1]+nonlin))) %>% unlist() %>% which()
  neg <- map(.x = nonlinearity, .f = ~(max(.x$rho)<0)) %>% unlist() %>% which()
  
  ## delete timeseries that seem to be linear or have negative rho from dataset
  del <- c(linear, neg) %>% unique()
  df_lin <- df[,del+1]
  df <- df[,-(del+1)]
  
  if(meta == TRUE){
    meta_information[[length(meta_information)+1]] <- length(linear)
    meta_information[[length(meta_information)+1]] <- length(neg)
    meta_information[[length(meta_information)+1]] <- c(linear, neg) %>% unique() %>% length()
    meta_information[[length(meta_information)+1]] <- (length(colnames(df))-1) 
    meta_information[[length(meta_information)+1]] <- (length(colnames(df))-1) * (length(colnames(df))-2)
    metanames <- c(metanames, "Trade connections dropped due to linearity")
    metanames <- c(metanames, "Trade connections dropped due to negative rho")
    metanames <- c(metanames, "Trade connections dropped due to linearity or negative rho (might be that for some connections both applied)")
    metanames <- c(metanames, "Trade connections left")
    metanames <- c(metanames, "Number of possible causal links")
  }  
  
  
  
  
  ######## CCM for all combinations
  
  ind <- crossing(lib_column = colnames(df)[-1], target_column = colnames(df)[-1])
  ind <- ind %>% filter(lib_column != target_column)
  
  ind <- ind %>% left_join(., bestE_df, by = c('lib_column' = 'link'))
  ind <- ind %>% rename(E = bestE)
  
  
  param <- rlang::as_list(ind)
  
  # parallelized version
  rho_list <- future_apply(ind, 1, function(x) {ccm(lib_column = x[1], target_column = x[2], E = as.numeric(x[3]), block = df,
                                                    lib_sizes = seq(5, dim(df)[1], by = 10), replace = FALSE, silent = TRUE,
                                                    random_libs = FALSE)})
  
  rho_list_means <- lapply(rho_list, ccm_means)
  
  
  # In some cases ccm can't make predictions and hence the t-test would fail. Therefore these are excluded
  not_working <- lapply(rho_list_means, function(x) {return(x$rho)}) 
  not_working <- lapply(not_working, function(x) {return(mean(x, na.rm = FALSE))})
  not_working <- which(is.na(not_working))
  
  
  if(length(not_working)>0){
    ind <- ind[-not_working, ]
    rho_list <- rho_list[-not_working]
    rho_list_means <- rho_list_means[-not_working]  
  }
  
  
  
  
  # checking convergance by demanding at least 0.1 increase in rho between shortest and longest library size
  
  increasing <- map(.x = rho_list_means, .f = ~ .x$rho[length(.x$rho)]-.x$rho[1]) %>% unlist()
  sign_incr <- which(increasing < 0.1)
  if(length(sign_incr)>0){
    rho_list <- rho_list[-sign_incr]
    rho_list_means <- rho_list_means[-sign_incr]
    ind <- ind[-sign_incr,]   
  }
  
  
  if(meta == TRUE){
    meta_information[[length(meta_information)+1]] <- length(sign_incr)
    metanames <- c(metanames, "Links dropped due to non-convergence")
  }
  
  
  
  
  # get correlations between the trade timeseries
  corr_list <-  map2(
    .x = ind$lib_column, .y = ind$target_column ,
    .f = ~ ccf(x = df[.x], y = df[.y],
               type = "correlation", plot = FALSE,
               na.action = na.pass, demean = FALSE)$acf
  )  
  
  
  ## t-test for each relationship: testing if we can say with significance that in mean the prediction skill is positive
  list_for_ttest <-  lapply(rho_list, function(x) {
    m <- max(x$lib_size)
    return(x[x$lib_size >= m-10, ])
  }
  )
  
  t_tests <-  map(.x = list_for_ttest, safely(
    .f = ~ t.test(.x$rho, alternative = "greater", mu = 0, na.action = na.exclude))
  )
  
  t_tests <- transpose(t_tests)
  
  fail <- t_tests$error %>% map_lgl(is_null) %>% unlist()
  if(sum(!fail)!=0){ ## TODO all shoul be true. I've had a few occasions where there was one false, what then?
    warning("Error in t_test", immediate. = TRUE)
  } 
  
  
  corr <-  map_dbl(.x = corr_list, function(x) max(abs(x)))
  
  for(i in 1:length(corr)){
    list_for_ttest[[i]]$corr <- corr[i]
  }
  
  t_tests2 <-  map(.x = list_for_ttest, safely(
    .f = ~ t.test(.x$rho, alternative = "greater", mu = .x$corr[1], na.action = na.exclude))
  )
  
  t_tests2 <- transpose(t_tests2)
  
  fail2 <- t_tests2$error %>% map_lgl(is_null) %>% unlist()
  if(sum(!fail2)!=0){ ## TODO all shoul be true. I've had a few occasions where there was one false, what then?
    warning("Error in t_test2", immediate. = TRUE)
  }
  
  
  # put information together
  
  ind <- ind %>%
    mutate(
      rho = map_dbl(.x = rho_list, .f = ~ mean(.x$rho, na.rm = TRUE)),
      rho_t = map_dbl(.x = t_tests$result, function(x) x$estimate ),
      p_value = map_dbl(.x = t_tests$result, function(x) x$p.value ),
      p_value2 = map_dbl(.x = t_tests2$result, function(x) x$p.value ),
      corr = map_dbl(.x = corr_list, function(x) max(abs(x))),
      detection = ifelse(p_value <0.05 & rho > 0.1, TRUE, FALSE),
      strong_detection = ifelse(p_value2 < 0.05 & p_value <0.05 & rho > 0.1, TRUE, FALSE)
    )
  
  if(meta == TRUE){
    meta_information[[length(meta_information)+1]] <- sum(!ind$detection)
    metanames <- c(metanames, "Links dropped due to insignificance (weak)")
    meta_information[[length(meta_information)+1]] <- sum(!ind$strong_detection)
    metanames <- c(metanames, "Links dropped due to insignificance (strong)")
  } 
  
  
  significant_weak <- which(ind$detection)
  significant_strong <- which(ind$strong_detection)    
  
  
  if(meta == TRUE){
    meta_information[[length(meta_information)+1]] <- length(significant_weak)
    meta_information[[length(meta_information)+1]] <- length(significant_strong)
    metanames <- c(metanames, "Total weak links in the end")
    metanames <- c(metanames, "Total strong links in the end")
  } 
  
  
  ind <- ind %>%
    separate(lib_column, into = c("A", "B"), sep = "_", remove = FALSE) %>%
    separate(target_column, into = c("C", "D"), sep = "_", remove = FALSE)
  
  
  names(meta_information) <- metanames
  
  return(list(links = ind, meta_information = meta_information, nonlinearity = nonlinearity_return, df = df_return))
}



######## Function to apply CCM on trade data including climate ########
links_for_climdat <- function(df, meta = TRUE, detrend = TRUE){
  
  meta_information <- list()
  metanames <- c()
  
  period <- as.data.frame(df[,1])
  df <- df[,-1]
  
  # detrending trade dataset
  if(detrend == TRUE){
    df_detrend <- cbind(apply(df[,1:2], 2, function(x){l <- 1:length(x);
    lm_res <- lm(x~l);
    pred <- lm_res$coefficients[[1]] + l *lm_res$coefficients[[2]];
    x <- x - pred;
    return(x)}))
  }
  df[,1:2] <- df_detrend
  
  
  # Normalizing
  df <- cbind(apply(df, 2, function(x) {x <- (x-mean(x, na.rm = TRUE))/sd(x, na.rm = T); return(x)}))
  
  df <- cbind(period, df)
  
  df_return <- df
  colnames(df)[1] <- "Period"
  
  
  ## get best embedding dimensions
  emb <- list()
  for (i in 2:dim(df)[2]){
    emb[[i-1]] <- df[c(1,i)] %>%
      simplex()
    emb[[i-1]]$link <- colnames(df)[i]
  }
  
  # best embedding dimension:
  bestE <- emb %>% map_dbl(~ .$E[which.max(.$rho)])
  bestE_df <- data_frame(
    bestE = bestE,
    link = colnames(df)[-1]
  )
  
  
  
  ######## CCM for all combinations
  
  ind <- crossing(lib_column = colnames(df)[-1], target_column = colnames(df)[-1])
  ind <- ind %>% filter(lib_column != target_column) %>% filter((lib_column %in% colnames(df)[2:3] | target_column %in% colnames(df)[2:3]))
  
  
  ind <- ind %>% left_join(., bestE_df, by = c('lib_column' = 'link'))
  ind <- ind %>% rename(E = bestE)
  
  
  # parallelized version
  rho_list <- future_apply(ind, 1, function(x) {ccm(lib_column = x[1], target_column = x[2], E = as.numeric(x[3]), block = df,
                                                    lib_sizes = seq(5, dim(df)[1], by = 10), replace = FALSE, silent = TRUE,
                                                    random_libs = FALSE)})
  
  rho_list_means <- lapply(rho_list, ccm_means)
  
  
  # In some cases ccm can't make predictions and hence the t-test would fail. Therefore these are excluded 
  not_working <- lapply(rho_list_means, function(x) {return(x$rho)}) 
  not_working <- lapply(not_working, function(x) {return(mean(x, na.rm = FALSE))})
  not_working <- which(is.na(not_working))
  
  if(length(not_working)>0){
    ind <- ind[-not_working, ]
    rho_list <- rho_list[-not_working]
    rho_list_means <- rho_list_means[-not_working]  
  }
  
  
  increasing <- map(.x = rho_list_means, .f = ~ .x$rho[length(.x$rho)]-.x$rho[1]) %>% unlist()
  sign_incr <- which(increasing < 0.1)
  if(length(sign_incr)>0){
    rho_list <- rho_list[-sign_incr]
    rho_list_means <- rho_list_means[-sign_incr]
    ind <- ind[-sign_incr,]   
  }
  
  
  if(meta == TRUE){
    meta_information[[length(meta_information)+1]] <- length(sign_incr)
    metanames <- c(metanames, "Links dropped due to non-convergence")
  }
  
  
  
  corr_list <-  map2(
    .x = ind$lib_column, .y = ind$target_column ,
    .f = ~ ccf(x = df[.x], y = df[.y],
               type = "correlation", plot = FALSE,
               na.action = na.pass, demean = FALSE)$acf
  )  
  
  rho_list <- lapply(rho_list, function(x) {x[x == -Inf] <- NA; return(x)})
  
  
  ## t-test for each relationship:
  list_for_ttest <-  lapply(rho_list, function(x) {
    m <- max(x$lib_size)
    return(x[x$lib_size >= m-10, ])
  }
  )
  
  t_tests <-  map(.x = rho_list, safely(
    .f = ~ t.test(.x$rho, alternative = "greater", mu = 0, na.action = na.exclude))
  )
  
  
  t_tests <- transpose(t_tests)
  fail <- t_tests$error %>% map_lgl(is_null) %>% unlist()
  if(sum(!fail)!=0){ ## all should be true
    warning("Error in t_test", immediate. = TRUE)
  } 
  
  corr <-  map_dbl(.x = corr_list, function(x) max(abs(x)))
  
  for(i in 1:length(corr)){
    list_for_ttest[[i]]$corr <- corr[i]
  }
  
  t_tests2 <-  map(.x = list_for_ttest, safely(
    .f = ~ t.test(.x$rho, alternative = "greater", mu = .x$corr[1], na.action = na.exclude))
  )
  
  t_tests2 <- transpose(t_tests2)
  
  fail2 <- t_tests2$error %>% map_lgl(is_null) %>% unlist()
  if(sum(!fail2)!=0){ ## TODO all shoul be true. I've had a few occasions where there was one false, what then?
    warning("Error in t_test2", immediate. = TRUE)
  }
  
  
  ind <- ind %>%
    mutate(
      rho = map_dbl(.x = rho_list, .f = ~ mean(.x$rho, na.rm = TRUE)),
      rho_t =(map(.x = t_tests$result, function(x) x$estimate) %>% unlist()),
      p_value = map_dbl(.x = t_tests$result, function(x) x$p.value ),
      p_value2 = map_dbl(.x = t_tests2$result, function(x) x$p.value ),
      corr = map_dbl(.x = corr_list, function(x) max(abs(x))),
      detection = ifelse(p_value <0.05 & rho > 0.1, TRUE, FALSE),
      strong_detection = ifelse(p_value2 < 0.05 & p_value <0.05 & rho > 0.1, TRUE, FALSE)
    )
  
  
  if(meta == TRUE){
    meta_information[[length(meta_information)+1]] <- sum(!ind$detection)
    metanames <- c(metanames, "Links dropped due to insignificance (weak)")
    meta_information[[length(meta_information)+1]] <- sum(!ind$strong_detection)
    metanames <- c(metanames, "Links dropped due to insignificance (strong)")
  } 
  
  
  significant_weak <- which(ind$detection)
  significant_strong <- which(ind$strong_detection)    
  
  
  if(meta == TRUE){
    meta_information[[length(meta_information)+1]] <- length(significant_weak)
    meta_information[[length(meta_information)+1]] <- length(significant_strong)
    metanames <- c(metanames, "Total weak links in the end")
    metanames <- c(metanames, "Total strong links in the end")
  } 
  
  
  ind <- ind %>%
    separate(lib_column, into = c("A", "B"), sep = "_", remove = FALSE) %>%
    separate(target_column, into = c("C", "D"), sep = "_", remove = FALSE)
  
  names(meta_information) <- metanames
  
  return(list(links = ind, meta_information = meta_information))
}
