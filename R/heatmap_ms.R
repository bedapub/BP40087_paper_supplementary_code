tidy_column_to_rownames <- function(df, column) {
  df %>%
    tibble::column_to_rownames(names(dplyr::select(df, {{ column }})))
}


mk_significant_cell_fun <- function (m_significant, text = "X", fontsize = 10)
{
  stopifnot(is.matrix(m_significant))
  cell_fun <- function(j, i, x, y, width, height, fill) {
    if (m_significant[i, j]) {
      grid::grid.text(text, x, y, gp = grid::gpar(fontsize = fontsize))
    }
  }
}
pd_heatmap_plot <- function(dat,
                            top_table,
                            unit = c("ABS", "%"),
                            top_flag = TOPFLAG,
                            timepoint = TIMEPOINT,
                            baseline_timepoint = "C1_D1",
                            cfb = CFB,
                            effect = Effect,
                            subject = USUBJID,
                            group = NULL,
                            test = testName,
                            row_label = {{ test }},
                            row_split_var = NULL,
                            test_annot = dplyr::vars(),
                            result = AVAL,
                            show_baseline = TRUE,
                            show_n_subjects = TRUE,
                            row_annotation_col = NULL,
                            show_row_names = TRUE,
                            column_annotation_col = NULL,
                            use_cfb_effect = FALSE,
                            show_cell_values = FALSE,
                            column_prefix = "",
                            column_title = NULL)
{
  unit <- match.arg(unit)

  top_table <- top_table %>%
    dplyr::select({{ test }}, {{ timepoint }}, {{ top_flag }}, {{ effect }})

  test_names <- dat %>%
    dplyr::distinct({{ test }}, {{ row_label }})

  suppressMessages(
    assertthat::assert_that(nrow(janitor::get_dupes(test_names, everything())) == 0)
  )

  # subset data for heatmap
  hm_dat <- dat %>%
    dplyr::select(
      {{ test }}, !!!test_annot, {{ timepoint }}, {{ group }},
      {{ subject }}, {{ result }}, {{ cfb }}
    ) %>%
    dplyr::filter(!is.na({{ cfb }}), !is.na({{ result }})) %>%
    dplyr::group_by({{ test }}, !!!test_annot, {{ timepoint }}, {{ group }}) %>%
    dplyr::summarize(
      {{ cfb }} := median({{ cfb }}),
      {{ result }} := median({{ result }}),
      .n_subjects = dplyr::n_distinct({{ subject }}),
      .groups = "drop"
    ) %>%
    dplyr::mutate({{ timepoint }} := droplevels({{ timepoint }}))

  test_anno_df <- dplyr::distinct(hm_dat, {{ test }}, !!!test_annot)

  # make sure we do not use any columns due to missing value
  if (rlang::quo_is_null(rlang::enquo(group))) {
    hm_dat_complete <- tidyr::crossing(
      {{ test }} := unique(dplyr::pull(hm_dat, {{ test }})),
      {{ timepoint }} := unique(dplyr::pull(hm_dat, {{ timepoint }}))
    )
  } else {
    hm_dat_complete <- tidyr::crossing(
      {{ test }} := unique(dplyr::pull(hm_dat, {{ test }})),
      {{ timepoint }} := unique(dplyr::pull(hm_dat, {{ timepoint }})),
      {{ group }} := unique(dplyr::pull(hm_dat, {{ group }}))
    )
  }
  hm_dat_complete <- hm_dat_complete %>%
    dplyr::left_join(test_anno_df, by = names(dplyr::select(hm_dat, {{ test }})))

  hm_dat <- dplyr::left_join(hm_dat_complete, hm_dat,
                             by = names(hm_dat_complete)
  ) %>%
    dplyr::mutate(column = {{ timepoint }}) %>%
    dplyr::arrange({{ test }}, {{ timepoint}})

  # number of column groups
  n_groups <- dat %>%
    dplyr::distinct({{ group }}) %>%
    nrow()

  if (n_groups > 1) {
    hm_dat <- hm_dat %>%
      dplyr::mutate(column = forcats::fct_cross({{ group }}, column))
  }

  n_sub_cols <- nrow(dplyr::distinct(hm_dat, {{ group }}))

  # significance "matrix"
  m_significant <- hm_dat %>%
    dplyr::filter({{ timepoint }} != baseline_timepoint) %>%
    dplyr::left_join(dplyr::filter(top_table, {{ top_flag }}),
                     by = names(dplyr::select(hm_dat, {{ test }}, {{ timepoint }}))
    ) %>%
    dplyr::mutate({{ top_flag }} := replace_na({{ top_flag }}, FALSE)) %>%
    dplyr::select({{ test }}, column, {{ top_flag }}) %>%
    tidyr::pivot_wider(names_from = "column", values_from = {{ top_flag }}) %>%
    tidy_column_to_rownames({{ test }}) %>%
    as.matrix()

  # effect matrix
  m <- hm_dat %>%
    dplyr::filter({{ timepoint }} != baseline_timepoint) %>%
    dplyr::select({{ test }}, {{ timepoint }}, column, {{ cfb }})


  if (use_cfb_effect || !rlang::quo_is_null(rlang::enquo(group))) {
    message("Effects are estimated based on median change from baseline")
    m <- m %>%
      dplyr::rename(.effect := {{ cfb }})
  } else {
    # in case when no group is provided, use effects from the top table
    stopifnot(nrow(top_table) == nrow(top_table))
    top_eff <- dplyr::select(top_table, -{{ top_flag }})
    top_eff_cols <- names(dplyr::select(top_eff, {{ test }}, {{ timepoint }}))
    m <- m %>%
      dplyr::select(-{{ cfb }}) %>%
      dplyr::left_join(top_eff, top_eff_cols) %>%
      dplyr::rename(.effect := {{ effect }})

    # backtransform
    if (unit == "ABS") {
      m$.effect <- log2(m$.effect)
    } else {
      # percentages
      m$.effect <- log(m$.effect)
    }
  }

  m <- m %>%
    dplyr::select(-{{ timepoint }}) %>%
    tidyr::pivot_wider(names_from = "column", values_from = .effect) %>%
    tidy_column_to_rownames({{ test }}) %>%
    as.matrix()

  if (unit == "%") {
    # transform log(odds ratio) to log2(odds ratio) for visualization
    m <- m / log(2)
  }

  row_labels <- tibble::tibble({{ test }} := rownames(m)) %>%
    dplyr::left_join(test_names, by = names(dplyr::select(test_names, {{ test }}))) %>%
    dplyr::pull({{ row_label }}) %>%
    str_replace(" ABS$", "")

  max_effect <- max(abs(m), na.rm = TRUE)
  max_effect_int <- ceiling(max_effect)
  breaks <- seq.int(-max_effect_int, max_effect_int, length.out = min(
    max_effect_int * 2 + 1,
    7
  ))
  m_pal <- circlize::colorRamp2(c(min(breaks), 0, max(breaks)), c("blue", "white", "red"))

  # number of observations
  if (show_n_subjects) {
    n_subjects <- hm_dat %>%
      dplyr::filter({{ timepoint }} != baseline_timepoint) %>%
      dplyr::group_by(column) %>%
      dplyr::summarize(n_subjects = max(.n_subjects), .groups = "drop") %>%
      dplyr::arrange(column)

    n_subjects_anno <- ComplexHeatmap::anno_text(n_subjects$n_subjects)
  } else {
    n_subjects_anno <- NULL
  }

  col_anno <- tibble::tibble(column = colnames(m)) %>%
    dplyr::left_join(dplyr::distinct(hm_dat, column, {{ timepoint }}, {{ group }}), by = "column")

  row_anno <- tibble::tibble({{ test }} := rownames(m)) %>%
    dplyr::left_join(dplyr::distinct(hm_dat, {{ test }}, !!!test_annot), by = names(dplyr::select(hm_dat, {{ test }})))

  if (!rlang::quo_is_null(rlang::enquo(row_split_var))) {
    set_row_split <- row_anno %>%
      dplyr::left_join(dat %>% dplyr::distinct({{ test }}, {{ row_split_var }}),
                       by = names(dplyr::select(row_anno, {{ test }}))
      ) %>%
      dplyr::pull({{ row_split_var }})
  } else {
    set_row_split <- NULL
  }

  # make sure group color are consistent across columns & baseline
  if (!rlang::quo_is_null(rlang::enquo(group))) {
    if (is.null(column_annotation_col)) {
      column_annotation_col <- list()
    }
    group_column <- names(dplyr::select(hm_dat, {{ group }}))
    if (is.null(column_annotation_col[[group_column]])) {
      group_levels <- levels(as.factor(hm_dat[[group_column]]))
      column_annotation_col[[group_column]] <- viridis::turbo(
        length(group_levels)
      ) %>%
        set_names(group_levels)
    }
  }

  column_ha <- ComplexHeatmap::HeatmapAnnotation(
    df = tidy_column_to_rownames(dplyr::select(col_anno, -{{ timepoint }}), column),
    n_subjects = n_subjects_anno,
    col = column_annotation_col
  )

  if (ncol(row_anno) > 1) {
    row_ha <- ComplexHeatmap::rowAnnotation(
      df = tidy_column_to_rownames(row_anno, {{ test }}),
      col = row_annotation_col
    )
  } else {
    row_ha <- NULL
  }

  if (show_baseline) {
    m_bl <- hm_dat %>%
      dplyr::filter({{ timepoint }} == baseline_timepoint) %>%
      dplyr::select({{ test }}, column, {{ result }}) %>%
      tidyr::pivot_wider(names_from = "column", values_from = {{ result }}) %>%
      dplyr::select(-{{ test }})

    if (unit == "ABS") {
      m_bl_t <- as.matrix(log2(m_bl + 1))
      mn <- min(m_bl_t)
      mx <- max(m_bl_t)
      blcol <- circlize::colorRamp2(seq(mn, mx, length.out = 10), viridis::viridis(10))
      at <- seq(mn, mx, length.out = 4)
      labels <- sprintf("%.0f", 2**at - 1)
    } else { # percent
      m_bl_t <- as.matrix(logit_mod(m_bl / 100))
      mn <- min(m_bl_t)
      mx <- max(m_bl_t)
      blcol <- circlize::colorRamp2(logit_mod(c(0, 0.1, 0.5, 0.9, 1)), viridis::viridis(5))
      at <- logit_mod(c(0.01, 0.1, 0.5, 0.9, 0.99))
      labels <- c("1%", "10%", "50%", "90%", "99%")
    }

    if (!rlang::quo_is_null(rlang::enquo(group))) {
      bl_col_anno <- tibble::tibble(column = colnames(m_bl)) %>%
        dplyr::left_join(dplyr::distinct(hm_dat, column, {{ timepoint }}, {{ group }}), by = "column")
      bl_column_ha <- ComplexHeatmap::HeatmapAnnotation(
        df = tidy_column_to_rownames(dplyr::select(bl_col_anno, -{{ timepoint }}), column),
        show_annotation_name = FALSE,
        col = column_annotation_col
      )
    } else {
      bl_column_ha <- NULL
    }

    colnames(m_bl_t) <- str_c(column_prefix, colnames(m_bl_t))
    hm_bl <- ComplexHeatmap::rowAnnotation(baseline = m_bl_t,
                                           annotation_legend_param = list(
                                             baseline=list(
                                       title = "baseline (cells/uL)",
                                       at = at,
                                       labels = labels,
                                       direction = "horizontal"
                                     )
    ), col = list(baseline=blcol))
  } else {
    hm_bl <- NULL
  }

  if (show_cell_values) {
    vals <- as.matrix(m)
    if (unit == "ABS") {
      vals <- 2**vals
    }
    cell_fn <- mk_values_cell_fun(vals)
  } else if (n_sub_cols == 1) {
    column_ha <- NULL
    cell_fn <- mk_significant_cell_fun(m_significant, text="∗", fontsize = 6)
  } else {
    stop("not available")
  }
  if (unit == "ABS") {
    name <- "Fold-change*"
    labels <- dplyr::if_else(breaks >= 0, as.character(2**breaks), stringr::str_c("1/", 2**(-breaks)))
  } else { # percentages
    name <- "Log2 odds ratio*"
    labels <- sprintf("%.2g", breaks)
  }

  hm <- ComplexHeatmap::Heatmap(m,
                                col = m_pal,
                                name = name,
                                row_labels = row_labels,
                                column_labels = str_c(column_prefix, colnames(m)),
                                cluster_rows = FALSE,
                                cluster_columns = FALSE,
                                show_column_names = TRUE,
                                show_row_names = show_row_names,
                                column_split = dplyr::pull(col_anno, {{ timepoint }}),
                                row_split = set_row_split,
                                #row_title_gp = grid::gpar(fontsize = 20),
                                cell_fun = cell_fn,
                                heatmap_legend_param = list(
                                  at = breaks,
                                  labels = labels,
                                  direction = "horizontal"
                                ),
                                row_names_side = "left",
                                left_annotation = hm_bl,
                                top_annotation = column_ha,
                                column_title = column_title
  )
  #hm_bl +
  hm
}
