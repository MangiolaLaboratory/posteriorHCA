#' composition_posterior_test: Perform posterior predictive analysis
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
#' @param load_model_to_global_env Flag whether load model file to global env.
#' @return A list with a result table and a density plot.
#' @export
#'
#' @import dplyr tidyr purrr ggplot2 sccomp dittoSeq ggrepel
#' @importFrom scales trans_new
composition_posterior_test <-
  function(
    proportions = NULL,
    sex = NULL,
    age_decade = NULL,
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
      if (exists("sccomp_est", envir = .GlobalEnv)){
        message('Model file loaded: Global Env.')
      }else{
        if (exists("sccomp_est", envir = environment())){
          message('Model file loaded in local Env. Moving it to Global Env.')
          sccomp_est <<- sccomp_est
          message('Model file loaded: Global Env.')
        }else{
          message('Load model file... (to Global Env)')
          sccomp_est <<- readRDS(path_to_model)
          message('Model file loaded: Global Env.')
        }
      }
    }else{
      if (exists("sccomp_est")){
        message('Model file loaded.')
      }else{
        message('Load model file... (to local Env)')
        sccomp_est <- readRDS(path_to_model)
        message('Model file loaded: local Env.')
      }
    }

    # if (!exists("sccomp_est", envir = if (load_model_to_global_env) .GlobalEnv else environment())) {
    #   message('Load model...')
    #
    #   if (load_model_to_global_env) {
    #     sccomp_est <<- readRDS(path_to_model)
    #   } else {
    #     sccomp_est <- readRDS(path_to_model)
    #   }
    # }
    #
    # message('Model loaded!')

    # Define valid values
    valid_cell_types <- sccomp_est %>% attr('count_data') %>% pull(L3) %>% unique() # refer to L2 cell type annotation
    valid_sex <- sccomp_est %>% attr('count_data') %>% pull(sex) %>% unique()
    valid_age <- sccomp_est %>% attr('count_data') %>% pull(age_decade) %>% unique()
    valid_disease_groups <- sccomp_est %>% attr('count_data') %>% pull(disease_groups___altered) %>% unique()
    valid_ethnicity_groups <- sccomp_est %>% attr('count_data') %>% pull(ethnicity_groups_imputed) %>% unique()
    valid_tissue_groups <- sccomp_est %>% attr('count_data') %>% pull(tissue_groups) %>% unique()

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
    if (!is.null(age_decade)) {
      if (length(age_decade) != 1 || !age_decade %in% valid_age) {
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
      age_decade = if (!is.null(age_decade)) as.character(age_decade) else NA,
      sex = if (!is.null(sex)) sex else NA,
      disease_groups___altered = if (!is.null(disease_groups)) disease_groups else NA,
      ethnicity_groups_imputed = if (!is.null(ethnicity_groups)) ethnicity_groups else NA,
      tissue_groups = if (!is.null(tissue_groups)) tissue_groups else NA
    ) %>% dplyr::select(where(~!any(is.na(.))))


    # Generate formula dynamically based on user input as sub-formula of full formula used in posterior model
    full_formula  <- sccomp_est %>% attr('formula_composition')
    environment(full_formula) <- environment() # reset environment


    sub_formula <- get_sub_formula(
      full_formula, factor_names = sample_metadata %>% select(!sample_id) %>% colnames
    )

    # browser()

    message("Start to generative sampling form posterior...")

    # Perform predictive sampling
    predict_res = sccomp_est %>% sccomp::sccomp_predict(
      formula_composition = sub_formula %>% as.formula,
      new_data = sample_metadata,
      summary_instead_of_draws = F
    ) %>% rename(cell_type = L3)

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

    # browser()

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
          Empirical_Confidence = 2 * pmin(
            mean(proportion_observed > proportion_sampled), # this calculate the quantile of observed proportion against sampling
            1- mean(proportion_observed > proportion_sampled)
          )
        ) %>% left_join(
          y = dist_by_cell_type,
          by = 'cell_type',
          relationship = "many-to-many"
        )

      # add dash line for observed proportion against density plot
      dist_plot <-
        dist_plot +
        geom_vline(
          data = dist_by_cell_type,
          aes(xintercept = proportion_observed, color = sample_id_observed),
          linetype = "dashed",
          linewidth = 0.5,
          show.legend = TRUE  # Ensure only vlines contribute to the legend
        )  +
        scale_color_manual(values = dittoSeq::dittoColors()) +
        labs(color = "EC: Empirical Confidence, Sample ID:") +  # Legend label
        theme(
          legend.position = "bottom",  # Move legend below the plot
          legend.title = element_text(size = 12, face = "bold"),  # Style legend title
          legend.text = element_text(size = 10)  # Style legend text
        ) +
        geom_text_repel(
          data = dist_by_cell_type,
          aes(
            x = proportion_observed,  #
            y = Inf,  # keep text to right top corner
            label = paste0("EC:", signif(Empirical_Confidence, 2)),
            color = sample_id_observed  # Match text color to vertical line color
          ),
          size = 4,
          direction = "y",  # Spread text vertically
          segment.color = NA,  # Remove repelling line
          inherit.aes = FALSE,
          show.legend = FALSE  # Prevents text from modifying the legend
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

#' Query a stored gene-level brms model
#'
#' Loads the fit, builds a covariate grid with [build_newdata_grid()],
#' then draws with [expr_draws()].
#'
#' Metadata rules (same as [build_newdata_grid()]):
#' * `NA` — marginalise over every level the model knows
#' * one value — fix that level
#' * several values — marginalise over only those levels
#'
#' @param cell_type A string indicating the cell type.
#' @param gene A gene symbol or Ensembl id (e.g. `"ADRB2"` or `"ENSG00000169252"`).
#' @param age_decade,sex,disease_groups,ethnicity_groups,assay_groups,tissue_groups
#'   Metadata choices. See the rules above.
#' @param version `"latest"` or a pinned container/version (e.g. `"V1"`).
#' @param quantity `"linpred"` for log(μ), `"predict"` for posterior
#'   predicted counts, `"epred"` for expected counts.
#' @param collapse `"mean"` averages grid profiles within each draw.
#'   `"pool"` stacks every draw x profile. `"sample"` picks one profile
#'   per draw.
#' @param ndraws Number of posterior draws, or `NULL` for all.
#' @param seed Optional RNG seed.
#' @return A list with `pred` (data frame of draws), `summary`, `plot`
#'   (from [plot_hca_draws()]), plus `grid`, `quantity`, and `collapse`.
#'   For `quantity = "predict"` or `"epred"`, the density plot uses a
#'   log1p-scaled x-axis.
#' @export
#'
#' @import ggplot2
expr_predict <- function(
  cell_type,
  gene,
  age_decade = NA,
  sex = NA,
  disease_groups = NA,
  ethnicity_groups = NA,
  assay_groups = NA,
  tissue_groups = NA,
  version = "latest",
  quantity = c("linpred", "predict", "epred"),
  collapse = c("mean", "pool", "sample"),
  ndraws = NULL,
  seed = NULL
) {
  quantity <- match.arg(quantity)
  collapse <- match.arg(collapse)

  fit <- load_expr_fit(
    cell_type = cell_type,
    gene = gene,
    version = version
  )

  grid <- build_newdata_grid(
    fit,
    age_decade = age_decade,
    sex = sex,
    disease_groups = disease_groups,
    ethnicity_groups = ethnicity_groups,
    assay_groups = assay_groups,
    tissue_groups = tissue_groups
  )

  out <- expr_draws(
    fit,
    newdata = grid,
    quantity = quantity,
    collapse = collapse,
    ndraws = ndraws,
    seed = seed
  )

  pred_df <- data.frame(value = out$draws)

  peak_location <- tryCatch(
    {
      dens <- stats::density(out$draws)
      dens$x[which.max(dens$y)]
    },
    error = function(e) NA_real_
  )

  density_plot <- plot_hca_draws(out)

  list(
    pred = pred_df,
    summary = list(
      mean = mean(out$draws),
      median = stats::median(out$draws),
      peak_location = peak_location
    ),
    plot = density_plot,
    draws = out$draws,
    grid = out$grid,
    quantity = out$quantity,
    collapse = out$collapse,
    cell_type = out$cell_type,
    gene_ensg = out$gene_ensg,
    gene_symbol = out$gene_symbol
  )
}
