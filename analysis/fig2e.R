# Author Iakov Davydov and Eveline Nueesch
library(tidyverse)

# reading internal parameters (such as bucket UUID)
# RS and EX are ADaM (analysis data model) domains
source("R/fig2e_internal.R")

respcolors <- c(
  "CR" = "#009E73",
  "PR" = "#0072B2",
  "SD" = "#F0E442",
  "PD" = "#D55E00",
  "Missing" = "grey",
  "CR (Complete response)" = "#009E73",
  "PR (Partial response)" = "#0072B2",
  "PD (Progressive disease)" = "#D55E00",
  "SD (Stable disease)" = "#F0E442"
)
group_shape <- c(
  "RO7122290 single agent" = 23,
  "RO7122290 + atezolizumab" = 21
)

# response table
resp <- rs.df %>%
    filter(PARAMCD == "BOR" & ((PART01 == "PART A" & PARTTP == "PAQWC") | (PART01 == "PART B" & PARTTP == "PBQ3WC") | (PART01 == "IMAG" & PARTTP == "PBQ3WC"))) %>%
    mutate(RESP = factor(AVALC, levels = c("PD", "SD", "PR", "CR"))) %>%
    select(USUBJID, RESP, PT) %>%
  distinct() %>%
  mutate(PatientID=str_c("Pt", PT))

# read rna-seq
e <- as.environment(list())
aws.s3::s3load(
  params$rnaseq_data$object,
  params$rnaseq_data$collection,
  envir = e
)
ae <- e$AEeset
rm(e)

ae <- ae[Biobase::fData(ae)$GeneOrSignature == "Gene",]
ae <- ae[,replace_na(ae$Selected == 1, FALSE)]

Biobase::pData(ae) <- Biobase::pData(ae) %>%
  left_join(resp, by=c("Category:PatientID"="PatientID"))

arm_info <- ex.df %>%
  select(USUBJID, ARMCD) %>%
  distinct() %>%
  transmute(PatientID=str_replace(USUBJID, "^BP40087-", "Pt"), ARMCD)

Biobase::pData(ae) <- Biobase::pData(ae) %>%
  left_join(arm_info, by=c("Category:PatientID"="PatientID"))
ae_bl <- ae[,Biobase::pData(ae)$VisitCode == 801]

ae_bl %>%
  Biobase::pData() %>%
  group_by(USUBJID) %>%
  filter(n()>1) %>%
  nrow() %>%
  {stopifnot(. == 0)}

# gene of interest
goi <- c("TNFRSF9")

gids <- which(!is.na(Biobase::fData(ae_bl)$Gene) & Biobase::fData(ae_bl)$Gene %in% goi)

ae_goi <- ae_bl[gids,]

ae_goi_genes <- ae_goi %>%
  Biobase::exprs() %>%
  t() %>%
  as_tibble() %>%
  set_names(Biobase::fData(ae_goi)$Gene)

ae_goi_df <- bind_cols(ae_goi_genes,
          Biobase::pData(ae_goi))

ae_goi_df_long <- ae_goi_df %>%
  pivot_longer(1:nrow(ae_goi), names_to="gene", values_to="expr") %>%
  mutate(GROUP = case_when(
    ARMCD == "PARTAQW" ~ 1,
    ARMCD == "PARTBQ3W" ~ 2,
    TRUE ~ NA_real_
  )) %>%
  mutate(GROUP = factor(GROUP, levels = c("1", "2"), labels = c("RO7122290 single agent", "RO7122290 + atezolizumab")))

geneBORplt <- . %>%
  # transform from log2(tpm + 5e-4) to log2(tpm+1e-3)
  mutate(expr = log2(2**expr - 5e-4 + 5e-3)) %>%
  filter(!is.na(RESP)) %>%
  mutate(RESP = factor(RESP, levels = names(respcolors[1:4]))) %>%
  group_by(gene) %>%
  nest() %>%
  mutate(
    plt = map(
      data,
      ~ ggplot(.x) +
        ggtitle(gene) +
        aes(RESP, expr) +
        geom_boxplot(aes(fill = RESP), outlier.shape = NA, alpha = 0.8) +
        xlab("Best overall response (BOR)") +
        ylab("Log2(TPM+5e-3)") +
        geom_hline(
          yintercept = log2(5e-3 + 1), col = "blue", linetype = "dashed"
        ) +
        geom_hline(yintercept = log2(5e-3), col = "red", linetype = "dashed") +
        guides(fill = guide_legend(title = "BOR")) +
        geom_jitter(
          aes(shape = GROUP),
          height = 0,
          width = 0.2,
          size = 3,
          alpha = 0.8,
          color = "black",
          fill = "#525252"
        ) +
        scale_fill_manual(values = respcolors[1:4], name = "") +
        scale_shape_manual(values = group_shape, name = "") +
        scale_x_discrete(limits = rev) +
        theme_bw()
    )
  )
set.seed(42)
ae_goi_df_long %>%
  geneBORplt() %>%
  pull(plt)
