# Gráficos do painel NGS (cohort.maf) — coorte de 93 pacientes, 10 genes do painel
# Tumor_Sample_Barcode do MAF = ID_REDCAP do banco conferido no REDCap
#   (Banco_Atualizado_REDCAP_19112026.xlsx)
# Filtro: só os 10 genes do painel (Tabela Suplementar 1). Variantes codificantes +
#         sítio de splicing entram nas análises; Intron, Splice_Region e 5'Flank
#         ficam como "silenciosas" e não entram.
# Figuras e tabelas em português: as classes de variante são renomeadas antes do
#   read.maf e os textos fixos do maftools são traduzidos por traducao_maftools.R.
# Uso: Rscript graficos_maf.R   (rodar a partir desta pasta)

invisible(Sys.setlocale("LC_ALL", "en_US.UTF-8"))   # acentos nos rótulos e nos CSV
suppressPackageStartupMessages({
  library(maftools); library(readxl); library(dplyr); library(ggplot2); library(tidyr)
})
source("traducao_maftools.R"); traduz_maftools()
options(OutDec = ",")   # vírgula decimal nos eixos das figuras (CSV continuam com ponto)

# Análise de sensibilidade: Rscript graficos_maf.R 74  -> exclui o paciente (ID_REDCAP) 74
excluir <- commandArgs(TRUE)
out <- if (length(excluir)) paste0("Figuras_MAF_sem_", paste(excluir, collapse = "_")) else "Figuras_MAF"
dir.create(out, showWarnings = FALSE)
dir.create(file.path(out, "lollipop"), showWarnings = FALSE)

genes_painel <- c("ASXL1","DNMT3A","FLT3","IDH1","IDH2","KMT2A","SF3B1","SRSF2","TET2","TP53")

# Executa um bloco de gráfico; se falhar, avisa e segue para o próximo
fig <- function(nome, expr, w = 3000, h = 2000) {
  png(file.path(out, nome), width = w, height = h, res = 300)
  ok <- tryCatch({ force(expr); TRUE },
                 error = function(e) { message("  [falhou] ", nome, ": ", conditionMessage(e)); FALSE })
  dev.off()
  if (!ok) unlink(file.path(out, nome)) else message("  ok: ", nome)
}

# Salva CSV com cabeçalhos em português (UTF-8)
nomes_pt <- c(Hugo_Symbol = "Gene", Tumor_Sample_Barcode = "ID_REDCAP", MutatedSamples = "Pacientes_mutados",
  AlteredSamples = "Pacientes_alterados", total = "Total_variantes", clusters = "Agrupamentos",
  muts_in_clusters = "Variantes_em_agrupamentos", clusterScores = "Escore_agrupamento",
  protLen = "Tamanho_proteina_aa", zscore = "Escore_z", pval = "p", fdr = "FDR", adjPval = "FDR",
  fract_muts_in_clusters = "Fracao_em_agrupamentos", or = "OR", ci.up = "IC95_superior", ci.low = "IC95_inferior",
  pct = "Percentual", hotspot = "Alteracao_proteica", n_mut = "Variantes", n_genes = "Genes_mutados",
  Diagnostico = "Diagnóstico", Blastos = "Blastos_pct", Cariotipo = "Cariótipo",
  Kt_complexo = "Cariótipo_complexo", Obito = "Óbito")
salva_csv <- function(df, nome) {
  df <- as.data.frame(df); nm <- names(df)
  nm <- ifelse(nm %in% names(nomes_pt), nomes_pt[nm], nm)
  nm <- sub("^n_mut_", "Variantes_", sub("^n_genes_", "Genes_mutados_", nm))
  names(df) <- nm
  write.csv(df, file.path(out, nome), row.names = FALSE, fileEncoding = "UTF-8")
}

# ---- Dados clínicos -------------------------------------------------------
num <- function(x) suppressWarnings(as.numeric(gsub(",", ".", as.character(x))))
na_txt <- function(x) { x <- as.character(x); x[grepl("informado|seguimento", x, ignore.case = TRUE)] <- NA; x }

# IPSS-M em dois grupos (Bernard et al., 2022): VL/L/ML = menor risco; MH/H/VH = maior risco
grupo_risco <- function(ipssm) {
  g <- dplyr::case_when(
    ipssm == "LMA" ~ "LMA-s",
    ipssm %in% c("Very Low", "Low", "Moderate Low") ~ "SMD baixo risco",
    ipssm %in% c("Moderate High", "High", "Very High") ~ "SMD alto risco")
  factor(g, levels = c("SMD baixo risco", "SMD alto risco", "LMA-s"))
}

# readxl falha em caminhos do OneDrive com acento: lê de uma cópia temporária
xlsx_tmp <- tempfile(fileext = ".xlsx")
file.copy("Banco_Atualizado_REDCAP_19112026.xlsx", xlsx_tmp)
ipssm_pt <- c("Very Low" = "Muito baixo", "Low" = "Baixo", "Moderate Low" = "Moderado baixo",
              "Moderate High" = "Moderado alto", "High" = "Alto", "Very High" = "Muito alto", "LMA" = "LMA-s")
aus_pres <- function(x) dplyr::case_when(grepl("^Aus", x) ~ "Ausente", grepl("^Pres", x) ~ "Presente",
                                         grepl("aplica", x) ~ "Não se aplica", TRUE ~ NA_character_)
clin <- read_excel(xlsx_tmp) |>
  filter(NGS == "Sim") |>
  transmute(Tumor_Sample_Barcode = as.character(ID_REDCAP),
            Grupo   = grupo_risco(Grupo_de_risco),
            IPSS_M  = factor(ipssm_pt[Grupo_de_risco], levels = ipssm_pt),
            Diagnostico = ifelse(Grupo_de_risco == "LMA", "LMA-s", "SMD"),
            Sexo, Idade = num(Idade), Blastos = num(Blastos),
            Cariotipo = na_txt(Resultado), Kt_complexo = aus_pres(Kt_complexo),
            Obito = na_txt(Obito))
stopifnot(!anyNA(clin$IPSS_M))
clin <- clin[!clin$Tumor_Sample_Barcode %in% excluir, ]
if (length(excluir)) cat("Excluídos:", excluir, "\n")
stopifnot(!anyNA(clin$Grupo))
n_grupo <- table(clin$Grupo)
cat("Pacientes por grupo:\n"); print(n_grupo)

cores_grupo <- setNames(c("#4C72B0", "#DD8452", "#C44E52"), levels(clin$Grupo))
cores_ipssm <- setNames(c("#1A9850","#91CF60","#D9EF8B","#FEE08B","#FC8D59","#D73027","#7B3294"), levels(clin$IPSS_M))
cores_dx    <- c(SMD = "#4C72B0", `LMA-s` = "#C44E52")
# maftools 2.28 falha ao gerar cores automáticas: definir uma paleta para cada anotação
pal <- function(x, cols) { u <- sort(unique(na.omit(as.character(x)))); setNames(cols[seq_along(u)], u) }
cores_cat <- list(Cariotipo   = c(Normal = "#BDBDBD", Alterado = "#E7298A"),
                  Kt_complexo = c(Ausente = "#BDBDBD", "Não se aplica" = "#F0F0F0", Presente = "#E7298A"),
                  Obito       = c("Não" = "#BDBDBD", Sim = "#252525"))
stopifnot(all(na.omit(clin$Cariotipo) %in% names(cores_cat$Cariotipo)),
          all(na.omit(clin$Kt_complexo) %in% names(cores_cat$Kt_complexo)),
          all(na.omit(clin$Obito) %in% names(cores_cat$Obito)))

# ---- MAF (só os 10 genes) -------------------------------------------------
maf_df <- read.delim("cohort.maf", comment.char = "#", check.names = FALSE,
                     colClasses = c(Tumor_Sample_Barcode = "character"))
maf_df <- maf_df[!maf_df$Tumor_Sample_Barcode %in% excluir, ]
stopifnot(all(maf_df$Tumor_Sample_Barcode %in% clin$Tumor_Sample_Barcode))
cat("Variantes no MAF:", nrow(maf_df), "| nos 10 genes:", sum(maf_df$Hugo_Symbol %in% genes_painel), "\n")
maf_df <- maf_df[maf_df$Hugo_Symbol %in% genes_painel, ]
# classes de variante em português (as não codificantes ficam como estão e são descartadas)
vc <- maf_df$Variant_Classification
maf_df$Variant_Classification <- ifelse(vc %in% names(classes_pt), classes_pt[vc], vc)

maf <- read.maf(maf_df, clinicalData = as.data.frame(clin), vc_nonSyn = unname(classes_pt), verbose = FALSE)
mut <- subsetMaf(maf, genes = genes_painel, fields = c("HGVSp_Short", "t_vaf"), includeSyn = FALSE, mafObj = FALSE)
cat("Variantes retidas (não silenciosas):", nrow(mut), "| pacientes com variante:",
    length(unique(mut$Tumor_Sample_Barcode)), "de", nrow(clin), "\n\n")

sink("/dev/null"); write.mafSummary(maf, basename = file.path(out, "T_maf")); sink()   # resumos por gene/amostra

# ---- Panorama -------------------------------------------------------------
fig("F01_resumo_maf.png",
    plotmafSummary(maf, color = cores_classes, rmOutlier = FALSE, addStat = "median", dashboard = TRUE, top = 10))

fig("F02_oncoplot_grupo.png", w = 3200, h = 2000,
    oncoplot(maf, genes = genes_painel, colors = cores_classes, clinicalFeatures = "Grupo",
             annotationColor = list(Grupo = cores_grupo), sortByAnnotation = TRUE,
             cohortSize = nrow(clin), showTitle = FALSE))

# Faixas de anotação mais altas (evita sobreposição dos nomes).
# O maftools só aceita nomes sintáticos do R (sem acento, espaço ou hífen)
# (acentos são aceitos; espaço e hífen não -> "IPSS_M" e "Cariótipo_complexo")
anno <- data.frame(Tumor_Sample_Barcode = clin$Tumor_Sample_Barcode,
                   Grupo = clin$Grupo, IPSS_M = clin$IPSS_M, Cariótipo = clin$Cariotipo,
                   Cariótipo_complexo = clin$Kt_complexo, Óbito = clin$Obito)
cores_anno <- list(Grupo = cores_grupo, IPSS_M = cores_ipssm, Cariótipo = cores_cat$Cariotipo,
                   Cariótipo_complexo = cores_cat$Kt_complexo, Óbito = cores_cat$Obito)
fig("F03_oncoplot_clinico.png", w = 3600, h = 3000,
    oncoplot(maf, genes = genes_painel, colors = cores_classes, annotationDat = anno,
             clinicalFeatures = names(cores_anno), annotationColor = cores_anno,
             anno_height = 3, annotationFontSize = 1.1, legendFontSize = 1.1, legend_height = 5,
             annoBorderCol = "white", sortByAnnotation = TRUE,
             cohortSize = nrow(clin), showTitle = FALSE))

fig("F04_oncoplot_vaf.png", w = 3200, h = 2000,
    oncoplot(maf, genes = genes_painel, colors = cores_classes, clinicalFeatures = "Grupo",
             annotationColor = list(Grupo = cores_grupo), sortByAnnotation = TRUE,
             cohortSize = nrow(clin), showTitle = FALSE,
             leftBarData = aggregate(t_vaf ~ Hugo_Symbol, data = data.frame(
               Hugo_Symbol = as.character(mut$Hugo_Symbol), t_vaf = as.numeric(mut$t_vaf)), FUN = median)))

fig("F05_oncostrip.png", w = 3000, h = 1600,
    oncostrip(maf, genes = genes_painel, clinicalFeatures = "Grupo",
              annotationColor = list(Grupo = cores_grupo), sortByAnnotation = TRUE))

fig("F06_titv.png", plotTiTv(titv(maf, plot = FALSE, useSyn = FALSE)))

fig("F07_vaf_maftools.png", plotVaf(maf, vafCol = "t_vaf", top = 10))

# ---- Lollipops (um por gene) ----------------------------------------------
for (g in genes_painel) {
  fig(file.path("lollipop", paste0("L_", g, ".png")), w = 3000, h = 1400,
      lollipopPlot(maf, gene = g, colors = cores_classes, AACol = "HGVSp_Short", showMutationRate = TRUE,
                   labelPos = "all", labPosSize = 0.7, cBioPortal = FALSE))
}

# ---- Interações, vias e drogas --------------------------------------------
fig("F08_interacoes_somaticas.png", w = 2400, h = 2400,
    somaticInteractions(maf, genes = genes_painel, pvalue = c(0.05, 0.1), showSum = TRUE))

fig("F09a_processos_biologicos.png", w = 2400, h = 1800, pathways(maf, pathdb = "smgbp", plotType = "bar"))
fig("F09b_vias_sinalizacao.png",   w = 2400, h = 1800, pathways(maf, pathdb = "sigpw", plotType = "bar"))
fig("F10_interacoes_drogas.png", w = 2400, h = 1800, drugInteractions(maf, fontSize = 0.75))

od <- tryCatch(oncodrive(maf, AACol = "HGVSp_Short", minMut = 3, pvalMethod = "zscore"),
               error = function(e) { message("  [oncodrive] ", conditionMessage(e)); NULL })
if (!is.null(od)) {
  salva_csv(od, "T_oncodrive.csv")
  fig("F11_oncodrive.png", plotOncodrive(res = od, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 0.6))
}

# ---- Comparação SMD vs. LMA-s ---------------------------------------------
smd <- subsetMaf(maf, clinQuery = "Diagnostico == 'SMD'")
lma <- subsetMaf(maf, clinQuery = "Diagnostico == 'LMA-s'")
n_smd <- sum(clin$Diagnostico == "SMD"); n_lma <- sum(clin$Diagnostico == "LMA-s")

cmp <- mafCompare(m1 = smd, m2 = lma, m1Name = "SMD", m2Name = "LMA-s", minMut = 3)
salva_csv(cmp$results, "T_mafCompare_SMD_vs_LMAs.csv")
print(cmp$results[, c("Hugo_Symbol","SMD","LMA-s","pval","or","adjPval")])

fig("F12_forest_SMD_vs_LMAs.png", w = 2600, h = 1800, forestPlot(mafCompareRes = cmp, pVal = 1, geneFontSize = 0.9))
fig("F13_cooncoplot_SMD_vs_LMAs.png", w = 3600, h = 2000,
    coOncoplot(m1 = smd, m2 = lma, colors = cores_classes, m1Name = "SMD",
               m2Name = "LMA-s", genes = genes_painel, removeNonMutated = FALSE,
               geneNamefont = 0.7, gene_mar = 1.6))
# atenção: no coOncoplot os percentuais são calculados só entre os pacientes COM variante
fig("F14_cobarplot_SMD_vs_LMAs.png", w = 2600, h = 1800,
    coBarplot(m1 = smd, m2 = lma, m1Name = "SMD", m2Name = "LMA-s", genes = genes_painel))
for (g in intersect(c("TET2", "DNMT3A", "KMT2A", "SF3B1"), genes_painel)) {
  fig(file.path("lollipop", paste0("L2_", g, "_SMD_vs_LMAs.png")), w = 3000, h = 2000,
      lollipopPlot2(m1 = smd, m2 = lma, gene = g, colors = cores_classes, AACol1 = "HGVSp_Short", AACol2 = "HGVSp_Short",
                    m1_name = "SMD", m2_name = "LMA-s"))
}

# ---- Enriquecimento clínico (genes x grupo de risco) ----------------------
ce <- tryCatch(clinicalEnrichment(maf, clinicalFeature = "Grupo", minMut = 3),
               error = function(e) { message("  [clinicalEnrichment] ", conditionMessage(e)); NULL })
if (!is.null(ce)) {
  salva_csv(ce$groupwise_comparision, "T_enriquecimento_grupo.csv")
  fig("F15_enriquecimento_grupo.png", w = 2800, h = 1800, plotEnrichmentResults(ce, pVal = 0.1))
}

# ---- Gráficos próprios (ggplot): denominadores com os 93 pacientes --------
freq <- mut |>
  distinct(Tumor_Sample_Barcode, Hugo_Symbol) |>
  left_join(clin, by = "Tumor_Sample_Barcode") |>
  count(Hugo_Symbol, Grupo) |>
  complete(Hugo_Symbol = genes_painel, Grupo = levels(clin$Grupo), fill = list(n = 0)) |>
  mutate(Grupo = factor(Grupo, levels(clin$Grupo)),
         pct = 100 * n / as.numeric(n_grupo[as.character(Grupo)]))
salva_csv(freq, "T_freq_gene_grupo.csv")
ordem <- freq |> group_by(Hugo_Symbol) |> summarise(t = sum(n)) |> arrange(t) |> pull(Hugo_Symbol)
ggsave(file.path(out, "F16_freq_gene_grupo.png"), width = 7, height = 5, dpi = 300, device = grDevices::png,
  ggplot(freq, aes(factor(Hugo_Symbol, ordem), pct, fill = Grupo)) +
    geom_col(position = position_dodge(0.8), width = 0.75) + coord_flip() +
    scale_fill_manual(values = cores_grupo) +
    labs(x = NULL, y = "Pacientes com variante (%)", fill = "Grupo") + theme_minimal(base_size = 12))

hot <- mut |>
  filter(!is.na(HGVSp_Short), HGVSp_Short != "") |>
  mutate(hotspot = paste0(Hugo_Symbol, ":", HGVSp_Short)) |>
  distinct(Tumor_Sample_Barcode, hotspot) |>
  count(hotspot, sort = TRUE) |> filter(n >= 2)
salva_csv(hot, "T_hotspots.csv")
ggsave(file.path(out, "F17_hotspots.png"), width = 7, height = max(3, 0.3 * nrow(hot) + 1), dpi = 300, device = grDevices::png,
  ggplot(hot, aes(reorder(hotspot, n), n)) + geom_col(fill = "#4C72B0") + coord_flip() +
    labs(x = NULL, y = "Número de pacientes") + theme_minimal(base_size = 12))

vaf <- mut |> left_join(clin, by = "Tumor_Sample_Barcode") |> mutate(vaf = 100 * as.numeric(t_vaf))
ggsave(file.path(out, "F18_vaf_gene_grupo.png"), width = 8, height = 5, dpi = 300, device = grDevices::png,
  ggplot(vaf, aes(reorder(Hugo_Symbol, vaf, median), vaf)) +
    geom_boxplot(outlier.shape = NA, fill = "grey92") +
    geom_jitter(aes(colour = Grupo), width = 0.15, size = 1.7, alpha = 0.8) +
    scale_colour_manual(values = cores_grupo) + coord_flip() +
    labs(x = NULL, y = "VAF (%)", colour = "Grupo") + theme_minimal(base_size = 12))

por_pac <- clin |>
  left_join(mut |> group_by(Tumor_Sample_Barcode) |>
              summarise(n_mut = n(), n_genes = n_distinct(Hugo_Symbol)),
            by = "Tumor_Sample_Barcode") |>
  mutate(across(c(n_mut, n_genes), ~ coalesce(.x, 0L)))
salva_csv(por_pac, "T_mutacoes_por_paciente_individual.csv")
resumo <- por_pac |> group_by(Grupo) |>
  summarise(n = n(), across(c(n_mut, n_genes),
            list(media = mean, dp = sd, mediana = median,
                 q1 = ~ quantile(.x, .25), q3 = ~ quantile(.x, .75), min = min, max = max)))
salva_csv(resumo, "T_mutacoes_por_paciente.csv")
kw_mut <- kruskal.test(n_mut ~ Grupo, data = por_pac)$p.value
kw_gen <- kruskal.test(n_genes ~ Grupo, data = por_pac)$p.value
cat("\nKruskal-Wallis: n_mut p =", signif(kw_mut, 3), "| n_genes p =", signif(kw_gen, 3), "\n")
ggsave(file.path(out, "F19_mutacoes_por_paciente.png"), width = 6, height = 4.5, dpi = 300, device = grDevices::png,
  ggplot(por_pac, aes(Grupo, n_mut, fill = Grupo)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.15, height = 0.1, size = 1.5) +
    scale_fill_manual(values = cores_grupo, guide = "none") +
    labs(x = NULL, y = "Variantes por paciente",
         subtitle = paste0("Kruskal-Wallis p = ", format(signif(kw_mut, 2), decimal.mark = ","))) +
    theme_minimal(base_size = 12))

# ---- Fisher com denominador real (93): maftools ignora pacientes sem variante ----
mutados <- mut |> distinct(Tumor_Sample_Barcode, Hugo_Symbol)
fisher_n93 <- function(col) {
  bind_rows(lapply(genes_painel, function(g) {
    tem <- clin$Tumor_Sample_Barcode %in% mutados$Tumor_Sample_Barcode[mutados$Hugo_Symbol == g]
    tb  <- table(factor(tem, c(TRUE, FALSE)), clin[[col]])
    data.frame(Gene = g, t(setNames(paste0(tb["TRUE", ], "/", colSums(tb)), colnames(tb))),
               p = fisher.test(tb)$p.value, check.names = FALSE)
  })) |> mutate(FDR = p.adjust(p, "BH")) |> arrange(p)
}
f_dx <- fisher_n93("Diagnostico"); f_gr <- fisher_n93("Grupo")
salva_csv(f_dx, "T_fisher_n93_SMD_vs_LMAs.csv")
salva_csv(f_gr, "T_fisher_n93_tres_grupos.csv")
cat("\nFisher (n = 93) SMD vs LMA-s:\n"); print(f_dx, row.names = FALSE)
cat("\nFisher (n = 93) tres grupos:\n"); print(f_gr, row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("\nPronto. Saídas em", normalizePath(out), "\n")
