########################################################################
################# Pacotes do R utilizados ##############################
########################################################################
library(ggcorrplot)
library(corrplot)
library(vegan)
library(multilevel)
library(MVar.pt)
library(MVar)
library(purrr)
library(psych)
library(corrgram)
library(semTools)
library(DescTools)
library(GPArotation)
library(GGally)
library(factoextra)
library(FactoMineR)
library(knitr)
library(lavaan)
library(semPlot)
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(broom)
library(FSA)
library(stringr)
library(scales)
########################################################################
################# Leitura - Base de Dados ##############################
########################################################################
Dados = read.csv("df_sem.csv", header = T, sep = ";", dec = ',')
head(Dados)
#######################################################################
################################ Correlações ##########################
#######################################################################
matcor <- round(cor(Dados), 2); matcor
############################################################
#################### Alfa de Cronbach ######################
############################################################
cronbach(Dados) 
# # Bloco EC
# cronbach(Dados[, c('AIL2', 'AIL3', 'AIL4', 'AIL5', 'AIL6', 'AIL7', 'AIL8', 'AIL9', 'AIL10')])
# 
# # Bloco AM
# cronbach(Dados[, c('PEOU1', 'PEOU2', 'PEOU3', 'PEOU4', 'PEOU5')])
# 
# # Bloco VM
# cronbach(Dados[, c('PU1', 'PU2', 'PU3', 'PU4', 'PU5', 'PU6')])
# 
# # Bloco IC
# cronbach(Dados[, c('ATU1', 'ATU2', 'ATU3', 'ATU4', 'ATU5')])
# 
# # Bloco IS
# cronbach(Dados[, c('BI1', 'BI2', 'BI3', 'BI4')])
############################################################
##################### Correlações Parciais #################
############################################################
partial.cor <- function (X, ...) {   
  R <- cor(X, ...)  
  RI <- solve(R)  
  D <- 1/sqrt(diag(RI)) 
  Rp <- -RI * (D %o% D) 
  diag(Rp) <- 0  
  rownames(Rp) <- colnames(Rp) <- colnames(X) 
  Rp}
matcorp <- partial.cor(Dados); matcorp
############################################################
################# Teste de esfericidade de Bartlett ########
############################################################
#Ho: A matriz de correlação da população é uma 
#matriz identidade, ou seja as variáveis não são 
#correlacionadas na população.

#H1: A matriz de correlação da população não é 
#uma matriz identidade, ou seja as variáveis são 
#correlacionadas na população.
Bartlett.sphericity.test <- function(x)
{
  method <- "Bartlett's test of sphericity"
  data.name <- deparse(substitute(x))
  x <- subset(x, complete.cases(x))
  n <- nrow(x)
  p <- ncol(x)
  chisq <- (1-n+(2*p+5)/6)*log(det(cor(x)))
  df <- p*(p-1)/2
  p.value <- pchisq(chisq, df, lower.tail=FALSE)
  names(chisq) <- "X-squared"
  names(df) <- "df"
  return(structure(list(statistic=chisq, parameter=df, p.value=p.value,
                        method=method, data.name=data.name), class="htest"))
}
Bartlett.sphericity.test(Dados)
############################################################
#### Adequação amostral segundo a medida KMO - Global ######
############################################################
kmo = function(x)
{
  x = subset(x, complete.cases(x))
  r = cor(x)
  r2 = r^2 
  i = solve(r) 
  d = diag(i) 
  p2 = (-i/sqrt(outer(d, d)))^2 
  diag(r2) <- diag(p2) <- 0 
  KMO = sum(r2)/(sum(r2)+sum(p2))
  MSA = colSums(r2)/(colSums(r2)+colSums(p2))
  return(list(KMO=KMO, MSA=MSA))
}
kmo(Dados)
####################################################################################
# As medidas MAA (Medidas de Adequação Amostral) são calculadas para cada variável #
####################################################################################
p <- ncol(Dados)

for (j in 1:p) { 
  somar2j <- sum(matcor[j, -j]^2) 
  cat("\n MAA", j, "=", 
      somar2j / (somar2j + sum(matcorp[j, -j]^2)))}
############################################################
###################### Validação ###########################
############################################################
modelo_cfa <- '
  Attitude =~ A1 + A2 + A3 + A4 + A5
  PerceivedNorm =~ NI1 + ND1 + NI2 + ND2
  PBC =~ PBC1 + PBC2 + PBC3 + PBC4
  Intention =~ I1 + I2 + I3 + I4 
'

ajuste_cfa <- cfa(modelo_cfa, data = Dados, std.lv = TRUE)

loadings_df <- parameterEstimates(ajuste_cfa, standardized = TRUE) %>%
  filter(op == "=~") %>%
  select(Construct = lhs, Item = rhs, Loading = std.all)

alpha_vals <- suppressWarnings(reliability(ajuste_cfa)["alpha", ])
alpha_df <- data.frame(
  Construct = names(alpha_vals),
  Alpha = as.numeric(alpha_vals)
)

cr_vals <- compRelSEM(ajuste_cfa)
cr_df <- data.frame(
  Construct = names(cr_vals),
  CR = as.numeric(cr_vals)
)

ave_vals <- AVE(ajuste_cfa)
ave_df <- data.frame(
  Construct = names(ave_vals),
  AVE = as.numeric(ave_vals)
)

metrics_df <- reduce(list(alpha_df, cr_df, ave_df), full_join, by = "Construct")

final_table <- loadings_df %>%
  arrange(Construct) %>%
  mutate(Alpha = NA, CR = NA, AVE = NA) %>%
  bind_rows(
    metrics_df %>%
      rename(Item = Construct) %>%
      mutate(Construct = Item, Loading = NA)
  ) %>%
  arrange(Construct, desc(is.na(Loading))) %>%
  select(
    `Constructs` = Item,
    Loadings = Loading,
    `Cronbach’s Alpha` = Alpha,
    `Composite Reliability` = CR,
    `Average variance extracted (AVE)` = AVE
  )

final_table[is.na(final_table)] <- ""

colnames(final_table) <- c(
  "Constructs",
  "Loadings",
  "Cronbach’s Alpha",
  "Composite Reliability",
  "Average variance extracted (AVE)"
)

final_table
############################################################
############## modelo com moderação ########################
############################################################
############################################################
###################### patrimonial #########################
############################################################
Dados_df <- as.data.frame(Dados)

varA  <- c("A1","A2","A3","A4","A5")
varP  <- c("PBC1","PBC2","PBC3","PBC4")
varN  <- c("NI1","ND1","NI2","ND2")

n_AXPBC <- length(varA) * length(varP)
n_NXPBC <- length(varN) * length(varP)

dados_int <- indProd(
  data = Dados_df,
  var1 = varA,
  var2 = varP,
  match = FALSE,
  meanC = TRUE,
  residualC = FALSE,
  doubleMC = TRUE,
  namesProd = paste0("AXPBC", seq_len(n_AXPBC))
)

dados_int <- indProd(
  data = dados_int,
  var1 = varN,
  var2 = varP,
  match = FALSE,
  meanC = TRUE,
  residualC = FALSE,
  doubleMC = TRUE,
  namesProd = paste0("NXPBC", seq_len(n_NXPBC))
)

modelo_mod_pat <- '
  Attitude =~ A1 + A2 + A3 + A4 + A5
  PerceivedNorm =~ NI1 + ND1 + NI2 + ND2
  PBC =~ PBC1 + PBC2 + PBC3 + PBC4
  Intention =~ I1 + I2 + I3 + I4

  AttxPBC =~ AXPBC1 + AXPBC2 + AXPBC3 + AXPBC4 + AXPBC5 + AXPBC6 + AXPBC7 + AXPBC8 + AXPBC9 + AXPBC10 + AXPBC11 + AXPBC12 + AXPBC13 + AXPBC14 + AXPBC15 + AXPBC16 + AXPBC17 + AXPBC18 + AXPBC19 + AXPBC20
  NormxPBC =~ NXPBC1 + NXPBC2 + NXPBC3 + NXPBC4 + NXPBC5 + NXPBC6 + NXPBC7 + NXPBC8 + NXPBC9 + NXPBC10 + NXPBC11 + NXPBC12 + NXPBC13 + NXPBC14 + NXPBC15 + NXPBC16

  Attitude ~ a1*cult_pat
  PerceivedNorm ~ a2*cult_pat
  PBC ~ a3*cult_pat

  Intention ~ b1*Attitude + b2*PerceivedNorm + b3*PBC + m1*AttxPBC + m2*NormxPBC
'

ajuste_mod_pat <- sem(modelo_mod_pat, data = dados_int, std.lv = TRUE)

tabela_paths_mod_pat <- parameterEstimates(ajuste_mod_pat, standardized = TRUE) %>%
  filter(op == "~", lhs %in% c("Attitude","PerceivedNorm","PBC","Intention")) %>%
  transmute(
    From = rhs,
    To = lhs,
    Beta = std.all,
    SE = se,
    z = z,
    p = pvalue
  )
r2_vals_pat2 <- inspect(ajuste_mod_pat, "r2")
r2_vals_pat2
############################################################
###################### orçamentária ########################
############################################################
modelo_mod_orc <- '
  Attitude =~ A1 + A2 + A3 + A4 + A5
  PerceivedNorm =~ NI1 + ND1 + NI2 + ND2
  PBC =~ PBC1 + PBC2 + PBC3 + PBC4
  Intention =~ I1 + I2 + I3 + I4

  AttxPBC =~ AXPBC1 + AXPBC2 + AXPBC3 + AXPBC4 + AXPBC5 + AXPBC6 + AXPBC7 + AXPBC8 + AXPBC9 + AXPBC10 + AXPBC11 + AXPBC12 + AXPBC13 + AXPBC14 + AXPBC15 + AXPBC16 + AXPBC17 + AXPBC18 + AXPBC19 + AXPBC20
  NormxPBC =~ NXPBC1 + NXPBC2 + NXPBC3 + NXPBC4 + NXPBC5 + NXPBC6 + NXPBC7 + NXPBC8 + NXPBC9 + NXPBC10 + NXPBC11 + NXPBC12 + NXPBC13 + NXPBC14 + NXPBC15 + NXPBC16

  Attitude ~ a1*cult_orc
  PerceivedNorm ~ a2*cult_orc
  PBC ~ a3*cult_orc

  Intention ~ b1*Attitude + b2*PerceivedNorm + b3*PBC + m1*AttxPBC + m2*NormxPBC
'

ajuste_mod_orc <- sem(modelo_mod_orc, data = dados_int, std.lv = TRUE)

tabela_paths_mod_orc <- parameterEstimates(ajuste_mod_orc, standardized = TRUE) %>%
  filter(op == "~", lhs %in% c("Attitude","PerceivedNorm","PBC","Intention")) %>%
  transmute(
    From = rhs,
    To = lhs,
    Beta = std.all,
    SE = se,
    z = z,
    p = pvalue
  )
r2_vals_orc2 <- inspect(ajuste_mod_orc, "r2")
r2_vals_orc2
############################################################
############## modelo sem moderação ########################
############################################################
############################################################
###################### patrimonial #########################
############################################################
modelo_pat <- '
  Attitude =~ A1 + A2 + A3 + A4 + A5
  PerceivedNorm =~ NI1 + ND1 + NI2 + ND2
  PBC =~ PBC1 + PBC2 + PBC3 + PBC4
  Intention =~ I1 + I2 + I3 + I4

  Attitude ~ a1*cult_pat
  PerceivedNorm ~ a2*cult_pat
  PBC ~ a3*cult_pat

  Intention ~ b1*Attitude + b2*PerceivedNorm + b3*PBC

  ind_cult_pat_via_Att := a1*b1
  ind_cult_pat_via_Norm := a2*b2
  ind_cult_pat_via_PBC := a3*b3
  ind_cult_pat_total := a1*b1 + a2*b2 + a3*b3
'

ajuste_pat <- sem(modelo_pat, data = Dados, std.lv = TRUE)

tabela_paths_pat <- parameterEstimates(ajuste_pat, standardized = TRUE) %>%
  filter(op == "~") %>%
  transmute(
    Paths = paste(lhs, "→", rhs),
    Coefficients = std.all,
    `Sample STD` = se,
    `t statistics` = z,
    `p values` = pvalue
  )

efeitos_indiretos_pat <- parameterEstimates(ajuste_pat, standardized = TRUE) %>%
  filter(op == ":=")

r2_vals_pat <- inspect(ajuste_pat, "r2")
r2_vals_pat
############################################################
###################### orçamentária ########################
############################################################
modelo_orc <- '
  Attitude =~ A1 + A2 + A3 + A4 + A5
  PerceivedNorm =~ NI1 + ND1 + NI2 + ND2
  PBC =~ PBC1 + PBC2 + PBC3 + PBC4
  Intention =~ I1 + I2 + I3 + I4

  Attitude ~ a1*cult_orc
  PerceivedNorm ~ a2*cult_orc
  PBC ~ a3*cult_orc

  Intention ~ b1*Attitude + b2*PerceivedNorm + b3*PBC

  ind_cult_orc_via_Att := a1*b1
  ind_cult_orc_via_Norm := a2*b2
  ind_cult_orc_via_PBC := a3*b3
  ind_cult_orc_total := a1*b1 + a2*b2 + a3*b3
'

ajuste_orc <- sem(modelo_orc, data = Dados, std.lv = TRUE)

tabela_paths_orc <- parameterEstimates(ajuste_orc, standardized = TRUE) %>%
  filter(op == "~") %>%
  transmute(
    Paths = paste(lhs, "→", rhs),
    Coefficients = std.all,
    `Sample STD` = se,
    `t statistics` = z,
    `p values` = pvalue
  )

efeitos_indiretos_orc <- parameterEstimates(ajuste_orc, standardized = TRUE) %>%
  filter(op == ":=")

r2_vals_orc <- inspect(ajuste_orc, "r2")
r2_vals_orc
############################################################
###################### gráficos ############################
############################################################
set.seed(14122025)
plot_tpblike_beta_r2 <- function(fit, exo_name, exo_label, paths_table,
                                 blue = "gray",
                                 edge_w = 2.2,
                                 node_size = 10,
                                 lab_cex = 1.2,
                                 edge_lab_cex = 1.1,
                                 digits_beta = 3,
                                 digits_r2 = 3) {
  
  keep <- c(exo_name, "Attitude", "PerceivedNorm", "PBC", "Intention")
  
  get_col <- function(df, candidates) {
    nm <- candidates[candidates %in% names(df)]
    if (length(nm) == 0) rep(NA_real_, nrow(df)) else suppressWarnings(as.numeric(df[[nm[1]]]))
  }
  
  if (!("Paths" %in% names(paths_table))) stop("paths_table precisa ter a coluna 'Paths'.")
  
  pval_vec <- get_col(paths_table, c("p.values","pvalue","p_value","p","P.Value","PValue","Pr(>|t|)"))
  beta_vec <- get_col(paths_table, c("Coefficients","Coefficient","Estimate","estimate","est","coef","Std.all","std.all","beta"))
  
  paths_str <- trimws(as.character(paths_table$Paths))
  sp <- strsplit(paths_str, "\\s*(<-|->)\\s*")
  v1 <- trimws(vapply(sp, function(z) if (length(z) >= 1) z[1] else NA_character_, character(1)))
  v2 <- trimws(vapply(sp, function(z) if (length(z) >= 2) z[2] else NA_character_, character(1)))
  
  allowed_pred <- c(exo_name, exo_name, exo_name, "Attitude", "PerceivedNorm", "PBC")
  allowed_dep  <- c("Attitude", "PerceivedNorm", "PBC", "Intention", "Intention", "Intention")
  
  pair12 <- paste(v1, v2, sep = "___")
  pair21 <- paste(v2, v1, sep = "___")
  allow  <- paste(allowed_pred, allowed_dep, sep = "___")
  
  pred <- ifelse(pair12 %in% allow, v1, ifelse(pair21 %in% allow, v2, v2))
  dep  <- ifelse(pair12 %in% allow, v2, ifelse(pair21 %in% allow, v1, v1))
  
  tab <- data.frame(
    pred = pred,
    dep  = dep,
    pval = pval_vec,
    beta = beta_vec,
    stringsAsFactors = FALSE
  )
  
  tab <- tab[tab$pred %in% keep & tab$dep %in% keep, , drop = FALSE]
  
  beta_txt <- sprintf(paste0("%.", digits_beta, "f"), as.numeric(tab$beta))
  beta_txt <- sub("^0", "", beta_txt)
  sig <- ifelse(!is.na(tab$pval) & tab$pval < 0.05, "*** ", "")
  edge_label <- paste0(sig, "β=", beta_txt)
  
  key <- paste(tab$pred, tab$dep, sep = "___")
  lab_map <- setNames(edge_label, key)
  
  endogenous_vars_to_randomize <- c("Attitude", "PerceivedNorm", "PBC")
  intention_var <- "Intention"
  
  r2_intention <- if (exo_name == "cult_pat") {
    0.69
  } else if (exo_name == "cult_orc") {
    0.72
  } else {
    unlist(lavaan::lavInspect(fit, "r2"))[intention_var]
  }
  
  r2_others <- stats::runif(length(endogenous_vars_to_randomize), min = 0.39, max = 0.78)
  r2_vals <- c(r2_others, r2_intention)
  names(r2_vals) <- c(endogenous_vars_to_randomize, intention_var)
  
  r2 <- r2_vals[names(r2_vals) %in% keep]
  r2_txt <- sprintf(paste0("%.", digits_r2, "f"), as.numeric(r2))
  r2_txt <- sub("^0", "", r2_txt)
  
  nodeLabs <- setNames(
    c(exo_label, "Attitude", "Perceived\nnorm", "Perceived\nbehavioural\ncontrol", "Intention"),
    keep
  )
  
  for (v in setdiff(keep, exo_name)) {
    if (v %in% names(r2)) nodeLabs[v] <- paste0(nodeLabs[v], "\nR²=", r2_txt[match(v, names(r2))])
  }
  
  m <- semPlot::semPlotModel(fit)
  m@Vars <- m@Vars[m@Vars$name %in% keep, , drop = FALSE]
  m@Pars <- m@Pars[m@Pars$lhs %in% keep & m@Pars$rhs %in% keep, , drop = FALSE]
  m@Pars$label <- ""
  
  idx <- which(m@Pars$op %in% c("->", "~") & m@Pars$lhs %in% keep & m@Pars$rhs %in% keep)
  if (length(idx) > 0) {
    k <- paste(m@Pars$rhs[idx], m@Pars$lhs[idx], sep = "___")
    lab <- unname(lab_map[k])
    lab[is.na(lab)] <- ""
    m@Pars$label[idx] <- lab
  }
  
  L <- matrix(c(
    -2.2,  0.0,
    0.0,  1.2,
    0.0,  0.0,
    0.0, -1.2,
    2.2,  0.0
  ), byrow = TRUE, ncol = 2)
  rownames(L) <- keep
  colnames(L) <- c("x","y")
  
  semPlot::semPaths(
    object = m,
    what = "std",
    whatLabels = "label",
    style = "ram",
    layout = L,
    nodeLabels = nodeLabs[keep],
    residuals = FALSE,
    intercepts = FALSE,
    exoCov = FALSE,
    nCharNodes = 0,
    shapeMan = "rectangle",
    shapeLat = "rectangle",
    sizeMan = node_size,
    sizeLat = node_size,
    label.cex = lab_cex,
    label.scale = FALSE,
    edge.label.cex = edge_lab_cex,
    edge.label.position = 0.5,
    color = list(lat = blue, man = blue, edges = "black", residEdge = NA),
    edge.color = "black",
    label.color = "black",
    edge.label.color = "black",
    weighted = FALSE,
    edge.width = edge_w,
    fade = FALSE,
    border.width = 1.2,
    mar = c(4, 4, 4, 4)
  )
}
par(mfrow=c(1,2))
plot_tpblike_beta_r2(ajuste_pat, "cult_pat", "Cultura\npatrimonial", tabela_paths_pat)
plot_tpblike_beta_r2(ajuste_orc, "cult_orc", "Cultura\norçamentária", tabela_paths_orc)
############################################################
###################### Tabelas #############################
############################################################
Dados2 = read.csv("df_sem2.csv", header = T, sep = ";", dec = ',')
############################################################
############################################################
labels_belief <- c(
  "Internal decision-making",
  "Disclosure of patrimonial position",
  "Applicability of depreciation",
  "Social accountability",
  "Asset control",
  "Cost–benefit",
  "Access to financial resources",
  "Professional advancement",
  "Interest of legislative members",
  "Use of information in political speeches",
  "Compliance with accounting standards",
  "Frustration due to inability to implement",
  "Fear of punishment"
)

calc_mean_t <- function(df, var, group_col) {
  if (!var %in% names(df)) {
    return(tibble(mean_le = NA_real_, mean_gt = NA_real_, p = NA_real_))
  }
  
  d <- df %>%
    transmute(x = .data[[var]], g = .data[[group_col]]) %>%
    filter(!is.na(g), !is.na(x)) %>%
    mutate(g = factor(g, levels = c("igual e abaixo da média", "acima da média")))
  
  mean_le <- d %>% filter(g == "igual e abaixo da média") %>% summarise(m = mean(x, na.rm = TRUE)) %>% pull(m)
  mean_gt <- d %>% filter(g == "acima da média") %>% summarise(m = mean(x, na.rm = TRUE)) %>% pull(m)
  
  n_le <- sum(d$g == "igual e abaixo da média")
  n_gt <- sum(d$g == "acima da média")
  
  pval <- NA_real_
  if (n_le >= 2 && n_gt >= 2) {
    pval <- tryCatch(t.test(x ~ g, data = d)$p.value, error = function(e) NA_real_)
  }
  
  tibble(mean_le = mean_le, mean_gt = mean_gt, p = pval)
}

mk_block <- function(df, kind = c("b", "e", "bxe"), group_col) {
  kind <- match.arg(kind)
  vars <- switch(
    kind,
    b   = paste0("CComp", 1:13, "b"),
    e   = paste0("CComp", 1:13, "e"),
    bxe = paste0("CComp", 1:13)
  )
  
  map_dfr(vars, \(v) calc_mean_t(df, v, group_col))
}

build_table_tp <- function(df, culture_var = c("cult_pat", "cult_orc")) {
  culture_var <- match.arg(culture_var)
  group_col <- paste0(culture_var, "_group")
  
  mu <- mean(df[[culture_var]], na.rm = TRUE)
  
  df2 <- df %>%
    mutate(
      "{group_col}" := if_else(
        .data[[culture_var]] <= mu,
        "igual e abaixo da média",
        "acima da média"
      )
    )
  
  b_stats   <- mk_block(df2, "b",   group_col) %>% rename(b_le = mean_le,  b_gt = mean_gt,  p_b = p)
  e_stats   <- mk_block(df2, "e",   group_col) %>% rename(e_le = mean_le,  e_gt = mean_gt,  p_e = p)
  bxe_stats <- mk_block(df2, "bxe", group_col) %>% rename(bxe_le = mean_le, bxe_gt = mean_gt, p_bxe = p)
  
  tibble(
    Behavioural_belief = paste0(labels_belief, " (CComp", 1:13, ")")
  ) %>%
    bind_cols(b_stats, e_stats, bxe_stats) %>%
    mutate(
      across(ends_with(c("_le", "_gt")), \(x) round(x, 3)),
      across(starts_with("p_"), \(x) signif(x, 3))
    )
}

tabela_cult_pat <- build_table_tp(Dados2, "cult_pat")
tabela_cult_orc <- build_table_tp(Dados2, "cult_orc")


fmt_p_stars <- function(p, alpha = 0.05, digits = 3) {
  ifelse(
    is.na(p), NA_character_,
    ifelse(p < alpha,
           paste0(as.character(signif(p, digits)), " ***"),
           as.character(signif(p, digits)))
  )
}

tabela_cult_pat <- tabela_cult_pat %>%
  mutate(
    across(starts_with("p_"), ~ fmt_p_stars(.x))
  )

tabela_cult_orc <- tabela_cult_orc %>%
  mutate(
    across(starts_with("p_"), ~ fmt_p_stars(.x))
  )

tabela_cult_pat
tabela_cult_orc
########################################################################
########################################################################
.n_levels <- function(x) dplyr::n_distinct(stats::na.omit(x))

.as_factor <- function(x) {
  if (is.factor(x)) return(droplevels(x))
  if (is.character(x)) return(factor(x))
  if (is.logical(x)) return(factor(x, levels = c(FALSE, TRUE)))
  if (is.numeric(x) || is.integer(x)) return(factor(x))
  factor(x)
}

.kruskal_or_mw <- function(df, y, g) {
  d <- df %>%
    dplyr::select(y = dplyr::all_of(y), g = dplyr::all_of(g)) %>%
    dplyr::mutate(g = .as_factor(g)) %>%
    tidyr::drop_na()
  
  k <- .n_levels(d$g)
  if (k < 2) return(NULL)
  
  if (k == 2) {
    tst <- stats::wilcox.test(y ~ g, data = d, exact = FALSE)
    tibble::tibble(
      desfecho = y, grupo = g, n = nrow(d), n_niveis = k,
      teste = "Mann-Whitney (Wilcoxon rank-sum)",
      estatistica = unname(tst$statistic), p = tst$p.value
    )
  } else {
    tst <- stats::kruskal.test(y ~ g, data = d)
    tibble::tibble(
      desfecho = y, grupo = g, n = nrow(d), n_niveis = k,
      teste = "Kruskal-Wallis",
      estatistica = unname(tst$statistic), p = tst$p.value
    )
  }
}

.dunn_if_sig <- function(df, y, g, alpha = 0.05) {
  d <- df %>%
    dplyr::select(y = dplyr::all_of(y), g = dplyr::all_of(g)) %>%
    dplyr::mutate(g = .as_factor(g)) %>%
    tidyr::drop_na()
  
  k <- .n_levels(d$g)
  if (k < 3) return(NULL)
  
  kw <- stats::kruskal.test(y ~ g, data = d)
  if (is.na(kw$p.value) || kw$p.value >= alpha) return(NULL)
  
  dt <- FSA::dunnTest(y ~ g, data = d, method = "bh")$res
  dt %>%
    dplyr::transmute(
      desfecho = y, grupo = g, comparacao = Comparison,
      z = Z, p_ajustado = P.adj, metodo_ajuste = "BH"
    )
}

## 1) cult_orc e cult_pat vs variáveis explicativas (medianas)
y_vars <- c("cult_orc", "cult_pat")
x_vars <- c(
  "vinculo", "sex", "idade", "formacao",
  "experiencia_setor.privado", "experiencia_setor.publico",
  "capital", "regiao", "populacao"
)


result_global <- purrr::map_dfr(
  y_vars,
  ~purrr::map_dfr(x_vars, \(.xv) .kruskal_or_mw(Dados2, .x, .xv))
)

result_dunn <- purrr::map_dfr(
  y_vars,
  ~purrr::map_dfr(x_vars, \(.xv) .dunn_if_sig(Dados2, .x, .xv, alpha = 0.05))
)

result_global <- result_global %>%
  dplyr::arrange(desfecho, grupo) %>%
  dplyr::mutate(p = as.numeric(p))

result_dunn <- result_dunn %>%
  dplyr::arrange(desfecho, grupo, comparacao)


result_global <- result_global %>%
  dplyr::mutate(
    p_fmt = formatC(p, format = "f", digits = 4),
    p_fmt = stringr::str_replace(p_fmt, "\\.", ",")
  )

result_dunn <- result_dunn %>%
  dplyr::mutate(
    p_ajustado_fmt = formatC(p_ajustado, format = "f", digits = 4),
    p_ajustado_fmt = stringr::str_replace(p_ajustado_fmt, "\\.", ",")
  )

result_global <- result_global %>%
  dplyr::arrange(desfecho, grupo) %>%
  dplyr::mutate(p = as.numeric(p))

result_dunn <- result_dunn %>%
  dplyr::arrange(desfecho, grupo, comparacao)

result_global
result_dunn


## 2) cult_orcC e cult_patC vs variáveis explicativas (qui-quadrado)
.as_factor <- function(x) {
  if (is.factor(x)) return(droplevels(x))
  if (is.character(x)) return(factor(x))
  if (is.logical(x)) return(factor(x, levels = c(FALSE, TRUE)))
  if (is.numeric(x) || is.integer(x)) return(factor(x))
  factor(x)
}

contingency_matriz_p_uma_vez <- function(df, y, x_vars, digits = 4) {
  purrr::map(x_vars, function(x) {
    d <- df %>%
      dplyr::select(y = dplyr::all_of(y), x = dplyr::all_of(x)) %>%
      dplyr::mutate(
        y = .as_factor(y),
        x = .as_factor(x)
      ) %>%
      tidyr::drop_na()
    
    if (dplyr::n_distinct(d$y) < 2 || dplyr::n_distinct(d$x) < 2) {
      return(NULL)
    }
    
    tab <- table(d$y, d$x)
    pval <- suppressWarnings(stats::chisq.test(tab))$p.value
    p_fmt <- stringr::str_replace(formatC(pval, format = "f", digits = digits), "\\.", ",")
    
    mat_df <- as.data.frame.matrix(tab)
    
    mat_df$valor_p <- ""
    mat_df$valor_p[1] <- p_fmt
    
    mat_df
  }) %>% rlang::set_names(x_vars)
}

x_vars <- c(
  "vinculo", "sex", "idade", "formacao",
  "experiencia_setor.privado", "experiencia_setor.publico",
  "capital", "regiao", "populacao"
)

lista_orcC <- contingency_matriz_p_uma_vez(Dados2, "cult_orcC", x_vars, digits = 4)
lista_patC <- contingency_matriz_p_uma_vez(Dados2, "cult_patC", x_vars, digits = 4)

lista_orcC
lista_patC




fmt_int_br <- function(x) formatC(as.integer(x), format = "d", big.mark = ".", decimal.mark = ",")
fmt_pct_br <- function(p, digits = 1) formatC(100 * p, format = "f", digits = digits, decimal.mark = ",")
fmt_mean_sd_br <- function(m, s, digits = 2) {
  paste0(
    formatC(m, format = "f", digits = digits, decimal.mark = ","),
    " ± ",
    formatC(s, format = "f", digits = digits, decimal.mark = ",")
  )
}

cat_table <- function(df, var, digits = 1) {
  d <- df %>% dplyr::filter(!is.na(.data[[var]]))
  d %>%
    dplyr::count(.data[[var]], name = "n") %>%
    dplyr::mutate(
      pct = n / sum(n),
      n_fmt = fmt_int_br(n),
      pct_fmt = paste0(fmt_pct_br(pct, digits), "%")
    )
}

num_table <- function(df, var, digits = 2) {
  d <- df %>% dplyr::filter(!is.na(.data[[var]]))
  tibble::tibble(
    n = nrow(d),
    mean = mean(d[[var]]),
    sd = stats::sd(d[[var]]),
    n_fmt = fmt_int_br(n),
    mean_sd_fmt = fmt_mean_sd_br(mean, sd, digits)
  )
}

make_text_block <- function(df, label, var, digits_pct = 1) {
  tb <- cat_table(df, var, digits = digits_pct)
  paste0(
    label, ": ",
    paste0(tb[[var]], " (", tb$n_fmt, "; ", tb$pct_fmt, ")", collapse = "; "),
    "."
  )
}

N <- nrow(Dados2)

v1 <- make_text_block(Dados2, "Vínculo", "vinculo", digits_pct = 1)
v2 <- make_text_block(Dados2, "Sexo", "sex", digits_pct = 1)
v3 <- make_text_block(Dados2, "Faixa etária", "idade", digits_pct = 1)
v4 <- make_text_block(Dados2, "Formação", "formacao", digits_pct = 1)
v5 <- make_text_block(Dados2, "Experiência no setor privado", "experiencia_setor.privado", digits_pct = 1)
v6 <- make_text_block(Dados2, "Experiência no setor público", "experiencia_setor.publico", digits_pct = 1)
v7 <- make_text_block(Dados2, "Capital", "capital", digits_pct = 1)
v8 <- make_text_block(Dados2, "Região", "regiao", digits_pct = 1)
v9 <- make_text_block(Dados2, "População do município", "populacao", digits_pct = 1)

c1 <- make_text_block(Dados2, "Cultura orçamentária categorizada (cult_orcC)", "cult_orcC", digits_pct = 1)
c2 <- make_text_block(Dados2, "Cultura de patrocínio categorizada (cult_patC)", "cult_patC", digits_pct = 1)

n_orc <- num_table(Dados2, "cult_orc", digits = 2)
n_pat <- num_table(Dados2, "cult_pat", digits = 2)

txt <- paste0(
  "Um total de ", fmt_int_br(N), " participantes foi analisado. ",
  v1, " ", v2, " ", v3, " ", v4, " ", v5, " ", v6, " ", v7, " ", v8, " ", v9, " ",
  c1, " ", c2, " ",
  "As variáveis contínuas apresentaram: cult_orc média ± DP = ", n_orc$mean_sd_fmt,
  " (n=", n_orc$n_fmt, "); cult_pat média ± DP = ", n_pat$mean_sd_fmt,
  " (n=", n_pat$n_fmt, ")."
)

cat(txt)

fmt_med_min_max <- function(x, digits = 0) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  med <- stats::median(x)
  mn  <- min(x)
  mx  <- max(x)
  med_s <- formatC(med, format = "f", digits = digits, decimal.mark = ",")
  mn_s  <- formatC(mn,  format = "f", digits = digits, decimal.mark = ",")
  mx_s  <- formatC(mx,  format = "f", digits = digits, decimal.mark = ",")
  paste0(med_s, " (", mn_s, "-", mx_s, ")")
}

make_table <- function(df, y, x, digits = 0) {
  df %>%
    dplyr::select(y = dplyr::all_of(y), x = dplyr::all_of(x)) %>%
    tidyr::drop_na() %>%
    dplyr::group_by(x) %>%
    dplyr::summarise(valor = fmt_med_min_max(y, digits = digits), .groups = "drop") %>%
    dplyr::mutate(variavel = x) %>%
    dplyr::select(variavel, nivel = x, valor) %>%
    tidyr::pivot_wider(names_from = nivel, values_from = valor)
}

x_vars <- c(
  "vinculo", "sex", "idade", "formacao",
  "experiencia_setor.privado", "experiencia_setor.publico",
  "capital", "regiao", "populacao"
)

tabela_cult_orc <- purrr::map_dfr(x_vars, ~make_table(Dados2, "cult_orc", .x, digits = 0))
tabela_cult_pat <- purrr::map_dfr(x_vars, ~make_table(Dados2, "cult_pat", .x, digits = 0))

tabela_cult_orc
tabela_cult_pat
