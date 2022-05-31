# function to get variables from a cmdstanr summary table
extract_from_summary <- function(summary_table, param, indices = NULL, mappings = NULL, troubleshoot = FALSE) {
  
  d <- summary_table
  
  if (!is.null(mappings)) {
    indices <- names(mappings)
    
    if (length(indices) > 1) {
      
      d <- d %>%
        filter(str_detect(variable, str_c("^", param, "\\["))) %>%
        separate(variable, into = indices, sep = ",", remove = troubleshoot) 
      
      for (i in 1:length(indices)) {
        d <- d %>%
          mutate(!!indices[i] := str_extract(eval(parse(text = indices[i])), "[0-9]+")) %>%
          mutate(!!indices[i] := plyr::mapvalues(
            eval(parse(text = indices[i])),
            mappings[[i]]$number %>% as.character(),
            mappings[[i]]$value
          ))
      }
      
    } else if (length(indices) == 1) {
      
      d <- d %>%
        filter(str_detect(variable, param)) %>%
        mutate(!!indices := str_extract(variable, "\\[[0-9]+") %>% str_extract("[0-9]+")) %>%
        mutate(!!indices := plyr::mapvalues(
          eval(parse(text = indices)),
          mappings[[1]]$number %>% as.character(),
          mappings[[1]]$value
        ))
    }
  
    return(d)
    
  } else if (!is.null(indices)) {
    
    d <- d %>%
      filter(str_detect(variable, str_c("^", param, "\\["))) %>%
      separate(variable, into = indices, sep = ",", remove = troubleshoot) %>%
      mutate(across(all_of(indices), ~str_extract(.x, "[0-9]+")))
    
    return(d)
    
  } else {
    stop("Provide either indices (as a vector) or mappings (named list of dataframes with columns 'name' and 'value')")
  }
  
}


