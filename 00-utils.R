library(tidyverse)
library(broom)
library(kableExtra)
library(knitr)
library(lme4)
library(lmerTest)
library(MuMIn)
library(rlang)

# function that transforms factors into factors and change contrasts
# to sum contrasts
prepare_factor_groups <- function(my_data) {
  for (variable in colnames(my_data)) {
    if (class(my_data[[variable]]) == "character") {
      my_data[[variable]] <- as.factor(my_data[[variable]])
      contrasts(my_data[[variable]]) <- "contr.sum"
    }
  }
  return(my_data)
}


# function that takes in log odds and convert it to probabilities
logit2prob <- function(logit) {
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

# function to get all level estimates from a glmer model
get_all_estimates <- function(my_variable, my_model, my_data) {
  selected_variable <- my_model %>%
    filter(grepl(my_variable, term)) %>%
    select(estimate)
  
  how_many_levels <- length(selected_variable$estimate) + 1
  terms <- paste0(my_variable, c(1:how_many_levels))
  
  estimates <- as.data.frame(rowSums((t(t(
    contrasts(my_data[[my_variable]])) * selected_variable[[1]]))))
  colnames(estimates) <- "estimate"
  estimates$term <- terms
  estimates$levels <- rownames(estimates)
  estimates$factor_group <- my_variable
  estimates$weight <- logit2prob(estimates$estimate)
  estimates$range <- diff(range(estimates$weight))
  
  return(estimates)
}

# function to get count (i.e., n) and proportion (i.e., % dep_variable)
get_count_proportion <- function(my_factor_group, my_data, my_dep_variable) {
  my_dep_variable <- sym(my_dep_variable)
  my_factor_group <- sym(my_factor_group)
  count_df <- my_data %>%
    group_by(!!my_factor_group) %>%
    summarise(proportion = mean(!!my_dep_variable) * 100,
              n = n()) %>%
    ungroup()
  colnames(count_df)[1] <- "levels"
  count_df$factor_group <- as_name(my_factor_group)
  return(count_df)
}

# function to get the whole output table
get_results_table <- function(my_data, my_model_df, 
                              my_factor_group, my_dep_variable) {
  count_and_proportion <- get_count_proportion(my_factor_group,
                                               my_data,
                                               my_dep_variable)
  all_estimates <- get_all_estimates(my_factor_group, my_model_df, my_data)
  
  count_and_proportion$levels <- as.character(count_and_proportion$levels)
  all_estimates$levels <- as.character(all_estimates$levels)
  
  all_results <- left_join(count_and_proportion,
                           all_estimates,
                           by = c("levels", "factor_group"))
  return(all_results)
}

# function to get logistic regression results
logistic_regression <- function(my_formula, my_data) {
  # transform factors into factors and change contrasts to sum contrasts
  my_data <- prepare_factor_groups(my_data)
  
  
  # check if mixed effects
  formula_string <- paste(as.character(my_formula), collapse = " ")
  if (grepl("\\|", formula_string)) {
    mixed_effects <- TRUE
  } else {
    mixed_effects <- FALSE
  }
  
  if (mixed_effects) {
    model <- glmer(my_formula,
                   data = my_data,
                   family = "binomial")
  } else {
    model <- glm(my_formula,
                   data = my_data,
                   family = "binomial")
  }

  
  model_df <- tidy(model)
  model_df$weight <- logit2prob(model_df$estimate)
  
  input <- model_df %>%
    filter(term == "(Intercept)") %>%
    select(weight)
  
  input <- input$weight
  variance_explained <- as.data.frame(suppressWarnings(r.squaredGLMM(model)))
  n <- nobs(model)
  
  complete <- all_factor_groups(my_formula, my_data, model_df)
  
  nice_table <- get_output_table(complete)
  
  model_results <- list(model = model,
                        model_df = model_df,
                        complete = complete,
                        nice_table = nice_table,
                        input = input,
                        variance_explained = variance_explained,
                        n = n)
  return(model_results)
}

all_factor_groups <- function(my_formula, my_data, my_model_df) {
  # get estimates for hidden terms
  my_terms <- terms(my_formula)
  my_factor_groups <- colnames(attr(my_terms, "factors"))
  
  all_new_estimates <- data.frame()
  for (factor_group in my_factor_groups) {
    if (!grepl("\\|", factor_group)) {
      these_new_estimates <- get_results_table(my_data, my_model_df, 
                                               factor_group, 
                                               all.vars(my_formula)[1])
      all_new_estimates <- bind_rows(all_new_estimates, these_new_estimates)
    }
  }
  return(all_new_estimates)
}

get_output_table <- function(model_df) {
  table2print <- arrange(model_df, -range, -weight)
  table2print$range <- format(table2print$range, digits = 2)
  
  how_many_rows <- nrow(table2print)
  copy_of_factor_group <- table2print$factor_group
  
  for (i in c(1:how_many_rows)) {
    if (i > 1) {
      if (copy_of_factor_group[i] == copy_of_factor_group[i - 1]) {
        table2print$factor_group[i] <- ""
      }
    }
  }
  
  for (i in c(how_many_rows:1)) {
    if (i < how_many_rows) {
      if (copy_of_factor_group[i] == copy_of_factor_group[i + 1]) {
        table2print$range[i] <- ""
      }
    }
  }
  
  table2print <- table2print %>%
    select(factor_group, levels, n, proportion, 
           "logodds" = estimate, weight, range)
  
  return(kable(table2print, digits = 2, format = "rst"))
}

######################### SET UP AND DOWN #################################

get_best_model <- function(step_results, this_step){
  # filter data per step (-1 indicates all steps)
  if (this_step > 0) {
    step_results <- step_results %>%
      filter(step == this_step)
  } 
  
  # filter data to keep only significant results
  step_clean <- step_results[step_results[, ncol(step_results)] <= .05, ]
  
  # get top result by r squared total
  best_step_model <- step_clean %>%
    arrange(-r_squared_total) %>%
    head(n = 1)
  
  # get formula for top model
  return(best_step_model$model_formula)
}

update_fixed_terms <- function(fixed_terms, best_formula, step_type) {
  new_fixed_terms <- c()
  if (step_type == "up") {
    for (term in fixed_terms) {
      if (!grepl(term, best_formula)) {
        new_fixed_terms <- c(new_fixed_terms, term)
      }
    }
  } else {
    for (term in fixed_terms) {
      if (grepl(term, best_formula)) {
        new_fixed_terms <- c(new_fixed_terms, term)
      }
    }
  }
  
  return(new_fixed_terms)
}

step_up_formulas <- function(fixed_terms, base_formula) {
  all_formulas <- c()
  factor_count <- length(fixed_terms)
  for (i in c(1:factor_count)) {
    this_formula <- paste(base_formula, "+", fixed_terms[i])
    this_formula <- gsub("~ 1 \\+", "~", this_formula)
    all_formulas <- c(all_formulas, this_formula)
  }
  return(all_formulas)
}

reorganize_formula <- function(my_formula_string) {
  if (grepl("\\|", my_formula_string)) {
    my_formula_string_clean <- gsub("\\s+\\|\\s+", "\\|", my_formula_string)
    formula_parts <- str_split(my_formula_string_clean, pattern = " ")[[1]]
    where_random_effects <- grepl("\\|", formula_parts)
    random_effects <- formula_parts[where_random_effects]
    not_random_effects <- formula_parts[!where_random_effects]
    new_formula_part1 <- paste(not_random_effects, collapse = " ")
    new_formula_part2 <- paste(random_effects, collapse = " + ")
    new_formula <- paste(new_formula_part1, "+", new_formula_part2)
    new_formula <- gsub("\\+ \\+", "\\+", new_formula)
    return(new_formula)
  } else {
    return(my_formula_string)
  }
  
}

step_down_formulas <- function(fixed_terms, base_formula) {
  all_formulas <- c()
  formula_parts <- str_split(base_formula, pattern = " ")[[1]]
  factor_count <- length(fixed_terms)
  for (term in fixed_terms) {
    this_formula <- paste(formula_parts[!formula_parts == term], collapse = " ")
    this_formula <- gsub("\\+ \\+", "\\+", this_formula)
    this_formula <- gsub("~ \\+", "~", this_formula)
    this_formula <- gsub("\\s\\+$", "", this_formula)
    all_formulas <- c(all_formulas, this_formula)
  }
  return(all_formulas)
}


step_one_way <- function(my_formula, my_data, step_type = "up") {
  dep_var <- all.vars(my_formula)[1]
  my_terms <- terms(my_formula)
  my_factor_groups <- colnames(attr(my_terms, "factors"))
  
  fixed_terms <- c()
  random_effects <- c()
  
  for (term in my_factor_groups) {
    mixed_model <- FALSE
    if (grepl("\\|", term)) {
      this_random_effect <- paste0("(", term, ")")
      random_effects <- c(random_effects, this_random_effect)
      mixed_model <- TRUE
    } else {
      fixed_terms <- c(fixed_terms, term)
    }
  }
  
  random_effects_string <- paste(random_effects, collapse = " + ")
  
  # get results for first model (no predictors)
  if (step_type == "up") {
    base_formula <- paste(dep_var, "~ 1")
  } else if (step_type == "down") {
    all_fixed_terms <- paste(fixed_terms, collapse = " + ")
    base_formula <- paste(dep_var, "~", all_fixed_terms)
  } else {
    stop(paste(step_type, 'is not a valid step_type. Enter a valid step_type: "up" or "down"'))
    return()
  }
  
  if (mixed_model) {
    if (grepl("\\|", base_formula)) {
      first_formula <- formula(base_formula)
    } else {
      first_formula <- formula(paste(base_formula, "+", random_effects_string))
    }
    first_model <- glmer(first_formula,
                         data = my_data,
                         family = "binomial")
  } else {
    first_formula <- formula(base_formula)
    first_model <- glm(first_formula,
                       data = my_data,
                       family = "binomial")
  }

  all_models_for_ref <- list(first_model)
  
  all_models <- data.frame()
  best_model_formula_string <- base_formula

  
  how_many_steps <- length(fixed_terms)
  if (step_type == "down") {
    how_many_steps <- how_many_steps - 1
  }
  
  for (step in c(1:how_many_steps)) {
    if (step_type == "up") {
      all_formulas <- step_up_formulas(fixed_terms, best_model_formula_string)
    } else {
      all_formulas <- step_down_formulas(fixed_terms, best_model_formula_string)
    }
    
    formula_count <- length(all_formulas)
    for (i in c(1:formula_count)) {
      if (mixed_model) {
        if (grepl("\\|", all_formulas[i])) {
          a_formula <- paste(all_formulas[i])
        } else {
          a_formula <- paste(all_formulas[i], "+", random_effects_string)
        }
        
        a_formula <- reorganize_formula(a_formula)
  
        finished_formula <- formula(a_formula)
        current_model <- glmer(finished_formula,
                               data = my_data,
                               family = "binomial")
      } else {
        a_formula <- paste(all_formulas[i])
        finished_formula <- formula(a_formula)
        current_model <- glm(finished_formula,
                             data = my_data,
                             family = "binomial")
      }
      
      variance_explained <- as.data.frame(suppressWarnings(r.squaredGLMM(current_model)))
      
      previous_model <- all_models_for_ref[[step]]
      models_anova <- anova(previous_model, current_model, test="Chisq")
      model_comparison <- as.data.frame(models_anova)[2, ]
      
      comparison_string <- gsub("\\n", " -- ", paste(attr(models_anova, "heading"), collapse = " "))
      
      current_model_df <- data.frame(step = step,
                                     model_formula = a_formula,
                                     r_squared_fixed = variance_explained$R2m[1],
                                     r_squared_total = variance_explained$R2c[1],
                                     comparison = comparison_string)
      
      current_model_df$model_formula <- as.character(current_model_df$model_formula)
      current_model_df <- bind_cols(current_model_df, model_comparison)
      all_models <- bind_rows(all_models, current_model_df)
    }
    
    this_best_model_string <- get_best_model(all_models, step)
    if (length(this_best_model_string) != 0) {
      best_model_formula_string <- this_best_model_string
    }

    fixed_terms <- update_fixed_terms(fixed_terms, best_model_formula_string, step_type)
    
    best_model_formula <- formula(best_model_formula_string)
    
    # get results for best model
    if (mixed_model) {
      best_step_model <- glmer(best_model_formula,
                               data = my_data,
                               family = "binomial")
    } else {
      best_step_model <- glm(best_model_formula,
                             data = my_data,
                             family = "binomial")
    }
    
    all_models_for_ref[[step+1]] <- best_step_model
    
  }

  best_of_all <- get_best_model(all_models, -1)
  
  all_models <- all_models %>%
    mutate(best = if_else(model_formula == best_of_all, "best", ""))
  
  return(all_models)
  
}

step_up_and_down <- function(my_formula, my_data) {
  df_up <- step_one_way(my_formula, my_data, "up")
  df_up$step_type <- "up"
  
  df_down <- step_one_way(my_formula, my_data, "down")
  df_down$step_type <- "down"
  
  best_up <- df_up %>%
    filter(best == "best")
  
  best_down <- df_down %>%
    filter(best == "best")
  
  if (best_up$model_formula == best_down$model_formula) {
    match <- "Step up and step down best models match."
  } else {
    match <- "Step up and step down best models DO NOT match."
  }
  
  warning(paste0("STEP UP best model: ",  best_up$model_formula, 
        "\nSTEP DOWN best model: ",  best_down$model_formula, 
        "\n", match))
  
  return(bind_rows(df_up, df_down))
  
}
