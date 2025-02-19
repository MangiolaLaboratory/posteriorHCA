#' posterior_test: Perform posterior predictive analysis
#'
#' This function performs posterior predictive sampling and returns
#' a table and a density plot.
#'
#' @param proportions A data frame of observed cell type proportions.
#' @param sex A string indicating sex ("male", "female", or "unknown").
#' @param age_bin A string for age bin category.
#' @param disease_groups A string indicating disease group.
#' @param ethnicity_groups A string indicating ethnicity group.
#' @param assay_groups A string indicating assay group.
#' @param tissue_groups A string indicating tissue group.
#' @param path_to_model Path to the `.rds` posterior model file.
#' @return A list with a result table and a density plot.
#' @export
#'
#' @import dplyr tidyr purrr ggplot2 sccomp dittoSeq
#' @importFrom scales trans_new
posterior_test <-
  function(
    proportions = NULL,
    sex = NULL,
    age_bin = NULL,
    disease_groups = NULL,
    ethnicity_groups = NULL,
    assay_groups = NULL,
    tissue_groups = NULL,
    load_model_to_global_env = T
    ){

    # Get the path to the model file (cached or downloaded)
    path_to_model <- get_model_ready()

    # Check if model is loaded in the wanted environment
    if (load_model_to_global_env){
      if (exists("estimates_age_bins___L2", envir = .GlobalEnv)){
        message('Model file loaded: Global Env.')
      }else{
        if (exists("estimates_age_bins___L2", envir = environment())){
          message('Model file loaded in local Env. Moving it to Global Env.')
          estimates_age_bins___L2 <<- estimates_age_bins___L2
          message('Model file loaded: Global Env.')
        }else{
          message('Load model file... (to Global Env)')
          estimates_age_bins___L2 <<- readRDS(path_to_model)
          message('Model file loaded: Global Env.')
        }
      }
    }else{
      if (exists("estimates_age_bins___L2")){
        message('Model file loaded.')
      }else{
        message('Load model file... (to local Env)')
        estimates_age_bins___L2 <- readRDS(path_to_model)
        message('Model file loaded: local Env.')
      }
    }

    # if (!exists("estimates_age_bins___L2", envir = if (load_model_to_global_env) .GlobalEnv else environment())) {
    #   message('Load model...')
    #
    #   if (load_model_to_global_env) {
    #     estimates_age_bins___L2 <<- readRDS(path_to_model)
    #   } else {
    #     estimates_age_bins___L2 <- readRDS(path_to_model)
    #   }
    # }
    #
    # message('Model loaded!')

    # Define valid values
    valid_cell_types <- estimates_age_bins___L2 %>% attr('truncation_df2') %>% pull(L2) %>% unique() # refer to L2 cell type annotation
    valid_sex <- estimates_age_bins___L2 %>% attr('truncation_df2') %>% pull(sex) %>% unique()
    valid_age <- estimates_age_bins___L2 %>% attr('truncation_df2') %>% pull(age_bin) %>% unique()
    valid_disease_groups <- estimates_age_bins___L2 %>% attr('truncation_df2') %>% pull(disease_groups) %>% unique()
    valid_ethnicity_groups <- estimates_age_bins___L2 %>% attr('truncation_df2') %>% pull(ethnicity_groups) %>% unique()
    valid_tissue_groups <- estimates_age_bins___L2 %>% attr('truncation_df2') %>% pull(tissue_groups) %>% unique()

    # Validate input arguments

    # Check proportions
    if (!is.null(proportions)) {
      if (!is.data.frame(proportions) || !ncol(proportions) %in% c(2,3)) {
        stop("Error: 'proportions' must be a data frame with either two (single sample) or three (multiple samples) columns.")
      }
      if (ncol(proportions) == 2){
        colnames(proportions) <- c("cell_type", "proportion")
      }else{
        colnames(proportions) <- c("sample_id", "cell_type", "proportion")
      }

      if (!all(proportions$cell_type %in% valid_cell_types)) {
        stop("Error: The 'cell_type' column must contain valid cell types.")
      }
      if (proportions %>% group_by(sample_id) %>% reframe(d = duplicated(cell_type)) %>% pull(d) %>% any) {
        stop("Error: The 'cell_type' column must not contain duplicate values.")
      }
      if (!is.numeric(proportions$proportion) || any(proportions$proportion < 0) || any(proportions$proportion >= 1)) {
        stop("Error: The 'proportion' column must be numeric and within the range [0,1).")
      }
    }

    # Check sex
    if (!is.null(sex)) {
      if (length(sex) != 1 || !sex %in% valid_sex) {
        stop("Error: 'sex' must be a single value and one of 'female', 'male', or 'unknown'.")
      }
    }

    # Check age
    if (!is.null(age_bin)) {
      if (length(age_bin) != 1 || !age_bin %in% valid_age) {
        stop("Error: 'age' must be a single value and one of the predefined age categories.")
      }
    }

    # Check disease_groups
    if (!is.null(disease_groups)) {
      if (length(disease_groups) != 1 || !disease_groups %in% valid_disease_groups) {
        stop("Error: 'disease_groups' must be a single valid value.")
      }
    }

    # Check ethnicity_groups
    if (!is.null(ethnicity_groups)) {
      if (length(ethnicity_groups) != 1 || !ethnicity_groups %in% valid_ethnicity_groups) {
        stop("Error: 'ethnicity_groups' must be a single valid value.")
      }
    }

    # Check tissue_groups
    if (!is.null(tissue_groups)) {
      if (length(tissue_groups) != 1 || !tissue_groups %in% valid_tissue_groups) {
        stop("Error: 'tissue_groups' must be a single valid value.")
      }
    }

    message("All inputs of metadata are valid.")

    # Create sample metadata tibble
    sample_metadata <- tidyr::tibble(
      sample_id = "dummy_sample_id_for_query", # this column name must be sample_id to match HCA data trained by model
      age_bin = if (!is.null(age_bin)) age_bin else NA,
      sex = if (!is.null(sex)) sex else NA,
      disease_groups = if (!is.null(disease_groups)) disease_groups else NA,
      ethnicity_groups = if (!is.null(ethnicity_groups)) ethnicity_groups else NA,
      tissue_groups = if (!is.null(tissue_groups)) tissue_groups else NA
    ) %>% dplyr::select(where(~!any(is.na(.))))


    # Generate formula dynamically based on user input as sub-formula of full formula used in posterior model
    full_formula  <- estimates_age_bins___L2 %>% attr('formula_composition')
    environment(full_formula) <- environment() # reset environment


    sub_formula <- get_sub_formula(
      full_formula, factor_names = sample_metadata %>% select(!sample_id) %>% colnames
    )

    # browser()

    message("Start to generative sampling form posterior...")

    # Perform predictive sampling
    predict_res = estimates_age_bins___L2 %>% sccomp::sccomp_predict(
      formula_composition = sub_formula %>% as.formula,
      new_data = sample_metadata,
      summary_instead_of_draws = F
    ) %>% rename(cell_type = L2)

    # Create summary statistics table for posterior distribution
    dist_by_cell_type = predict_res %>%
      group_by(cell_type) %>%
      reframe(
        mean = mean(proportion),
        lower = quantile(proportion, probs = 0.025),
        upper = quantile(proportion, probs = 0.975)
      )

    # Define arcsine square root transformation
    arcsine_sqrt_trans <- scales::trans_new(
      "arcsine_sqrt",
      transform = function(x) asin(sqrt(x)),
      inverse = function(x) (sin(x))^2
    )

    # If user provided observed proportions, filter out unwanted cell types
    if (!is.null(proportions)){
      predict_res <-
        predict_res %>% filter(cell_type %in% (proportions$cell_type %>% unique))
    }

    # Generate density plot for predicted proportions
    dist_plot <-
      ggplot(predict_res, aes(x = proportion)) +
      geom_density(alpha = 0.3) +
      facet_wrap(~cell_type) +  # Facet by cell_group
      scale_x_continuous(
        trans = arcsine_sqrt_trans,
        name = "Proportion (Arcsine-Sqrt Scaled)"
      ) +
      labs(
        # title = "Density of Proportions by Cell Group",
        x = "Proportion",
        y = "Density"
      )+
      theme_minimal() +
      theme(strip.text = element_text(size = 12, face = "bold"))

    # If observed proportions exist, compare against predictions
    if (!is.null(proportions)){

      dist_by_cell_type <-
        proportions %>%
        right_join(
          x = predict_res,
          by = 'cell_type',
          suffix = c('_sampled', '_observed'),
          relationship = "many-to-many"
        ) %>%
        group_by(sample_id_observed, cell_type) %>%
        reframe(
          proportion_observed = proportion_observed %>% unique(),
          p_value = 2 * pmin(
            mean(proportion_observed > proportion_sampled), # this calculate the quantile of observed proportion against sampling
            1- mean(proportion_observed > proportion_sampled)
          )
        ) %>% left_join(
          y = dist_by_cell_type,
          by = 'cell_type',
          relationship = "many-to-many"
        )

      # add dash line for observed proportion against density plot
      dist_plot <- dist_plot +
        geom_vline(
          data = dist_by_cell_type,
          aes(xintercept = proportion_observed, color = sample_id_observed),
          # color = "red",
          linetype = "dashed",
          linewidth = 0.5
        )  +
        scale_color_manual(values = dittoSeq::dittoColors()) +
        labs(color = "Sample ID") +  # Legend label
        theme(
          legend.position = "bottom",  # Move legend below the plot
          legend.title = element_text(size = 12, face = "bold"),  # Style legend title
          legend.text = element_text(size = 10)  # Style legend text
        )
    }

    # Print results for immediate feedback
    print(dist_by_cell_type)
    print(dist_plot)

    # Return results as a list
    return(
      list(
        result_table = dist_by_cell_type,
        plot = dist_plot
      )
    )

}
