# ============================================================
# Reanálise da PARTE CLÍNICA do metiloma com o banco conferido no REDCap (19/11/2026).
# Rótulos, tabelas e figuras em português.
#
# A matriz de metilação, o QC, os filtros e os agrupamentos são os mesmos de
# analysis_epicv2_pronon.qmd (mesmos parâmetros e mesma semente). Só a informação
# clínica muda. Cada teste é rodado com o pdata ANTERIOR e com o REDCap, lado a lado.
#
# Uso (de dentro desta pasta):  Rscript 01_reanalise_clinica_PT.R
# Entradas: sdfs_betas.rds, pdata_pronon_REDCap_PT.csv, ../pdata_pronon.csv
# Saída:    saida_clinica_REDCap/{tabelas,figuras}
# ============================================================
invisible(Sys.setlocale("LC_ALL", "en_US.UTF-8"))
suppressPackageStartupMessages({
  library(sesame); library(limma); library(DMRcate); library(ggplot2)
  library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
  library(pheatmap); library(EpiDISH); library(ConsensusClusterPlus)
})

out <- "saida_clinica_REDCap"
dir.create(file.path(out, "tabelas"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "figuras"), recursive = TRUE, showWarnings = FALSE)
tab <- function(df, nome) write.csv(df, file.path(out, "tabelas", paste0(nome, ".csv")),
                                    row.names = FALSE, fileEncoding = "UTF-8")
fig <- function(g, nome, w = 8, h = 6)
  ggsave(file.path(out, "figuras", paste0(nome, ".png")), g, width = w, height = h, dpi = 300,
         device = grDevices::png)

cores_grupo <- c("SMD baixo risco" = "#3E7CB1", "SMD alto risco" = "#2E9E4F", "LMA-s" = "#E4572E")

# ---------------- Dados clínicos ----------------
pd  <- read.csv("pdata_pronon_REDCap_PT.csv", encoding = "UTF-8")
pd$Grupo <- factor(pd$Grupo, levels = names(cores_grupo))
ant <- read.csv("../pdata_pronon.csv")                       # pdata anterior, para comparação
ant <- ant[match(pd$Barcode, ant$Barcode), ]

# ---------------- Metilação: mesmo QC e filtros do .qmd ----------------
obj   <- readRDS("sdfs_betas.rds")
betas <- obj$betas_norm
det   <- 1 - colMeans(is.na(betas))
boas  <- names(det)[det >= 0.80]
boas  <- boas[!pd$ID_amostra[match(boas, pd$Barcode)] %in% c(26, 12)]   # outliers
b     <- betas[, boas]
pd    <- pd[match(boas, pd$Barcode), ];  ant <- ant[match(boas, ant$Barcode), ]
cat("Amostras analisadas:", ncol(b), "\n"); print(table(pd$Grupo))

b_det  <- b[rowMeans(is.na(b)) <= 0.05, ]
M_snp  <- rmSNPandCH(minfi::logit2(b_det), rmcrosshyb = TRUE)
ann    <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
M_auto <- M_snp[!(rownames(M_snp) %in% rownames(ann)[ann$chr %in% c("chrX", "chrY")]), ]
b_col  <- betasCollapseToPfx(b_det[rownames(M_auto), ])
M_col  <- minfi::logit2(b_col); M_col <- M_col[is.finite(rowSums(M_col)), ]
b_col  <- b_col[rownames(M_col), ]
cat("CpGs finais:", nrow(M_col), "\n")

# ---------------- Agrupamentos (idênticos aos anteriores) ----------------
top  <- names(sort(apply(M_col, 1, var), decreasing = TRUE))[1:2000]
mat  <- M_col[top, ]
hc   <- hclust(dist(t(mat)), method = "ward.D2")
pd$Cluster_hierarquico <- factor(paste0("H", cutree(hc, 3)[pd$Barcode]))

mat_cc <- sweep(as.matrix(mat), 1, rowMeans(as.matrix(mat)))
set.seed(1234)
cc <- ConsensusClusterPlus(mat_cc, maxK = 6, reps = 1000, pItem = 0.8, pFeature = 1,
        clusterAlg = "hc", distance = "pearson", innerLinkage = "ward.D2",
        finalLinkage = "ward.D2", seed = 1234, plot = "png",
        title = file.path(out, "figuras", "consenso"))
pd$Cluster_consenso <- factor(paste0("C", cc[[3]]$consensusClass[pd$Barcode]))

tab(data.frame(Cluster = rownames(table(pd$Cluster_hierarquico, pd$Grupo)),
               as.data.frame.matrix(table(pd$Cluster_hierarquico, pd$Grupo)), check.names = FALSE),
    "T1_cluster_hierarquico_x_grupo")
tab(data.frame(Cluster = rownames(table(pd$Cluster_consenso, pd$Grupo)),
               as.data.frame.matrix(table(pd$Cluster_consenso, pd$Grupo)), check.names = FALSE),
    "T2_cluster_consenso_x_grupo")
tab(pd[, c("ID_amostra", "ID_REDCAP", "Grupo", "Cluster_hierarquico", "Cluster_consenso")],
    "T0_clusters_por_amostra")

# ---------------- Clusters × clínica: anterior vs REDCap ----------------
# rótulo em português = c(coluna REDCap, coluna no pdata anterior)
num_v <- list("Blastos (%)"        = c("Blastos_pct", "Blast_pct"),
              "Idade"              = c("Idade", "Age"),
              "Hemoglobina"        = c("Hemoglobina", "Hemoglobin"),
              "Plaquetas"          = c("Plaquetas", "Platelets"),
              "Neutrófilos segmentados" = c("Segmentados", "Segmented_neutrophils"))
cat_v <- list("Grupo diagnóstico"  = c("Grupo", "Sample_type"),
              "Sexo"               = c("Sexo", "Sex"),
              "Cariótipo"          = c("Cariotipo", "Karyotype_status"),
              "Cariótipo complexo" = c("Cariotipo_complexo", "Karyotype_Complex"),
              "Celularidade"       = c("Celularidade", "Cellularity"),
              "Dependência transfusional" = c("Dependencia_transfusional", "Transfusion_dependence"),
              "Óbito"              = c("Obito", "Death"),
              "Perda do Y"         = c("Perda_Y", "Loss_of_Y"))
limpa <- function(x) { x <- as.character(x)
  x[x %in% c("Perda de seguimento", "Not_informed")] <- NA        # desfecho/dado desconhecido
  if (any(x %in% c("MDS_VERY_HIGH"))) x[x == "MDS_VERY_HIGH"] <- "MDS_HIGH"
  x }
testa <- function(df, col, cl, numerica) {
  x <- if (numerica) as.numeric(df[[col]]) else limpa(df[[col]]); ok <- !is.na(x)
  if (numerica) return(c(n = sum(ok), p = kruskal.test(x[ok] ~ droplevels(cl[ok]))$p.value))
  tb <- table(x[ok], droplevels(cl[ok])); if (nrow(tb) < 2) return(c(n = sum(ok), p = NA))
  set.seed(1234); c(n = sum(ok), p = fisher.test(tb, simulate.p.value = TRUE, B = 1e5)$p.value) }

compara <- function(cl) {
  r <- rbind(
    do.call(rbind, lapply(names(num_v), function(v) data.frame(Variável = v, Teste = "Kruskal-Wallis",
      t(setNames(testa(ant, num_v[[v]][2], cl, TRUE), c("n_anterior", "p_anterior"))),
      t(setNames(testa(pd,  num_v[[v]][1], cl, TRUE), c("n_REDCap",   "p_REDCap")))))),
    do.call(rbind, lapply(names(cat_v), function(v) data.frame(Variável = v, Teste = "Fisher (Monte Carlo)",
      t(setNames(testa(ant, cat_v[[v]][2], cl, FALSE), c("n_anterior", "p_anterior"))),
      t(setNames(testa(pd,  cat_v[[v]][1], cl, FALSE), c("n_REDCap",   "p_REDCap")))))))
  r$p_anterior <- signif(r$p_anterior, 2); r$p_REDCap <- signif(r$p_REDCap, 2)
  r$Mudou <- ifelse((r$p_anterior < 0.05) != (r$p_REDCap < 0.05), "SIM", "")
  r[order(r$p_REDCap), ] }
res_h  <- compara(pd$Cluster_hierarquico); tab(res_h,  "T3_cluster_hierarquico_vs_clinica")
res_cc <- compara(pd$Cluster_consenso);    tab(res_cc, "T4_cluster_consenso_vs_clinica")
cat("\n== Hierárquico × clínica ==\n"); print(res_h, row.names = FALSE)
cat("\n== Consenso × clínica ==\n");    print(res_cc, row.names = FALSE)

# Distribuição das variáveis que interessam, por cluster (para a descrição na tese)
desc <- function(cl, nome) {
  d <- do.call(rbind, lapply(levels(cl), function(k) { s <- pd[cl == k, ]
    data.frame(Cluster = k, n = nrow(s),
      Blastos_mediana = median(s$Blastos_pct, na.rm = TRUE),
      Blastos_IIQ = paste(quantile(s$Blastos_pct, c(.25, .75), na.rm = TRUE), collapse = "–"),
      Transfusao_Sim = sprintf("%d/%d", sum(s$Dependencia_transfusional == "Sim", na.rm = TRUE),
                               sum(!is.na(s$Dependencia_transfusional))),
      Obito_Sim = sprintf("%d/%d", sum(s$Obito == "Sim", na.rm = TRUE),
                          sum(s$Obito %in% c("Sim", "Não")))) }))
  tab(d, nome); print(d, row.names = FALSE) }
cat("\n"); desc(pd$Cluster_hierarquico, "T5_descricao_cluster_hierarquico")
desc(pd$Cluster_consenso, "T6_descricao_cluster_consenso")

# ---------------- PCA: R² de cada componente ----------------
pcs <- prcomp(t(M_col))$x[, 1:5]
r2  <- function(x) { ok <- !is.na(x); sapply(1:5, function(k) summary(lm(pcs[ok, k] ~ x[ok]))$r.squared) }
vars_pc <- list("Idade" = pd$Idade, "Blastos (%)" = pd$Blastos_pct, "Hemoglobina" = pd$Hemoglobina,
                "Plaquetas" = pd$Plaquetas, "Grupo diagnóstico" = pd$Grupo, "Sexo" = factor(pd$Sexo),
                "Celularidade" = factor(pd$Celularidade), "BeadChip" = factor(pd$Chip))
r2t <- data.frame(Variável = names(vars_pc), round(t(sapply(vars_pc, r2)), 3), check.names = FALSE)
names(r2t)[-1] <- paste0("PC", 1:5); tab(r2t, "T7_R2_componentes_principais"); print(r2t, row.names = FALSE)

g <- ggplot(data.frame(PC1 = pcs[, 1], PC2 = pcs[, 2], Blastos = pd$Blastos_pct, Grupo = pd$Grupo),
            aes(PC1, PC2, color = Blastos, shape = Grupo)) + geom_point(size = 3) +
  scale_color_viridis_c(option = "C", na.value = "grey70") +
  labs(title = "PCA colorida pelo percentual de blastos", color = "Blastos (%)", shape = "Grupo") +
  theme_minimal()
fig(g, "F1_PCA_por_blastos")

# ---------------- Heatmap com anotação em português ----------------
anot <- data.frame(Grupo = pd$Grupo, `Blastos (%)` = pd$Blastos_pct, Cariótipo = pd$Cariotipo,
                   `Dependência transfusional` = pd$Dependencia_transfusional,
                   Óbito = ifelse(pd$Obito == "Perda de seguimento", NA, pd$Obito),
                   `Cluster consenso` = pd$Cluster_consenso, check.names = FALSE)
rownames(anot) <- pd$Barcode
cores_anot <- list(Grupo = cores_grupo,
  Cariótipo = c(Normal = "#BDD7E7", Alterado = "#08519C"),
  `Dependência transfusional` = c(Não = "#FEE0D2", Sim = "#CB181D"),
  Óbito = c(Não = "#D9D9D9", Sim = "#252525"),
  `Cluster consenso` = c(C1 = "#8DA0CB", C2 = "#FC8D62", C3 = "#66C2A5"))
pheatmap(mat[, pd$Barcode], scale = "row", show_rownames = FALSE, show_colnames = FALSE,
         annotation_col = anot, annotation_colors = cores_anot, clustering_method = "ward.D2",
         main = "2.000 sítios CpG mais variáveis (valor M, padronizado por sítio)",
         filename = file.path(out, "figuras", "F2_heatmap_top2000.png"), width = 10, height = 7)

png(file.path(out, "figuras", "F3_dendrograma.png"), width = 2400, height = 1300, res = 220)
plot(hc, labels = paste0(pd$ID_amostra, "-", c("SMD baixo risco" = "BR", "SMD alto risco" = "AR",
     "LMA-s" = "LMA")[as.character(pd$Grupo)]), main = "Agrupamento hierárquico (Ward)", xlab = "", sub = "",
     ylab = "Distância")
invisible(dev.off())

# ---------------- Sensibilidade: blastos como covariável ----------------
ok <- !is.na(pd$Blastos_pct)
dd <- droplevels(pd[ok, ]); dd$Grupo <- factor(dd$Grupo, labels = make.names(levels(dd$Grupo)))
Db <- model.matrix(~ 0 + Grupo + Sexo + factor(Chip) + Blastos_pct, data = dd)
colnames(Db) <- make.names(colnames(Db))
Cb <- makeContrasts(LMAs_vs_baixo = GrupoLMA.s - GrupoSMD.baixo.risco, levels = Db)
fb <- eBayes(contrasts.fit(lmFit(M_col[, ok], Db), Cb))
sens <- data.frame(Análise = "SMD baixo risco vs. LMA-s, ajustado por blastos",
                   n = sum(ok), n_por_grupo = paste(table(dd$Grupo), collapse = "/"),
                   DMPs_FDR_0.05 = sum(decideTests(fb, p.value = 0.05) != 0))
tab(sens, "T8_sensibilidade_blastos"); print(sens, row.names = FALSE)

# ---------------- EpiDISH: composição celular estimada ----------------
fr <- epidish(beta.m = b_col, ref.m = centDHSbloodDMC.m, method = "RPC")$estF
ep <- data.frame(ID_amostra = pd$ID_amostra, Grupo = pd$Grupo, Blastos = pd$Blastos_pct,
                 round(fr[pd$Barcode, ], 4), check.names = FALSE)
tab(ep, "T9_composicao_celular_EpiDISH")
kw <- kruskal.test(ep$Neutro ~ ep$Grupo)$p.value
sp <- cor.test(ep$Neutro, ep$Blastos, method = "spearman", exact = FALSE)
epi <- data.frame(Teste = c("Neutrófilos entre grupos (Kruskal-Wallis)", "Neutrófilos × blastos (Spearman)"),
                  n = c(nrow(ep), sum(!is.na(ep$Blastos))),
                  rho = c(NA, round(unname(sp$estimate), 2)), p = signif(c(kw, sp$p.value), 2))
tab(epi, "T10_EpiDISH_testes"); print(epi, row.names = FALSE)
g <- ggplot(ep, aes(Blastos, Neutro, color = Grupo)) + geom_point(size = 3) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "grey40", linetype = "dashed") +
  scale_color_manual(values = cores_grupo) +
  labs(x = "Blastos (%)", y = "Fração estimada de neutrófilos (EpiDISH)",
       title = sprintf("Neutrófilos estimados × blastos (Spearman ρ = %.2f; p = %.3f)",
                       sp$estimate, sp$p.value)) + theme_minimal()
fig(g, "F4_EpiDISH_neutrofilos_vs_blastos")

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("\nConcluído. Saída em", out, "\n")
