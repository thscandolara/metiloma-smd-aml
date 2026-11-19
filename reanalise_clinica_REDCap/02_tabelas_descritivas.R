# ============================================================
# Tabelas descritivas da coorte (banco conferido no REDCap, 19/11/2026)
#   Tabela 1: coorte do painel NGS (93 pacientes)
#   Tabela 2: coorte do metiloma (43 pacientes com array)
# Gera:
#   Tabelas_descritivas_coorte.xlsx   -> para a aluna (tabelas com coluna de p em branco
#                                        + dados prontos para o programa de estatística)
#   GABARITO_testes_tabelas.xlsx      -> p-valores calculados no R, para conferência
# Uso (de dentro desta pasta): Rscript 02_tabelas_descritivas.R
# ============================================================
invisible(Sys.setlocale("LC_ALL", "en_US.UTF-8"))
suppressMessages({ library(readxl); library(openxlsx) })

xlsx <- "../dados/Banco_Atualizado_REDCAP_19112026.xlsx"   # ajuste para o caminho local do banco clínico (não versionado)
tmp  <- tempfile(fileext = ".xlsx"); file.copy(xlsx, tmp)   # readxl falha com acento no caminho
red  <- as.data.frame(read_excel(tmp))

faltante <- function(x) { x <- trimws(as.character(x))
  x[grepl("informad|ignorad|realizad", x, ignore.case = TRUE) | x == ""] <- NA; x }
aus_pres <- function(x) { x <- faltante(x)
  ifelse(is.na(x), NA, ifelse(grepl("^Aus", x), "Ausente", ifelse(grepl("^Pres", x), "Presente", x))) }
grupos <- c("SMD baixo risco", "SMD alto risco", "LMA-s")
ipssm_pt <- c("Very Low" = "Muito baixo", "Low" = "Baixo", "Moderate Low" = "Moderado baixo",
              "Moderate High" = "Moderado alto", "High" = "Alto", "Very High" = "Muito alto")
subtipo_pt <- c(
  "MDS with low blasts (MDS-LB)" = "SMD-BB",
  "MDS with low blasts and SF3B1 mutationa (MDS-SF3B1)" = "SMD-SF3B1",
  "MDS with low blasts and isolated 5q deletion (MDS-5q)" = "SMD-5q",
  "MDS with increased blasts (MDS-IB) 1" = "SMD-AB1",
  "MDS with increased blasts (MDS-IB) 2" = "SMD-AB2",
  "MDS, hypoplasticb (MDS-h)" = "SMD-h")

d <- data.frame(
  ID_REDCAP = red$ID_REDCAP,
  Metiloma  = red$Methylation_Array,
  Grupo     = factor(ifelse(red$Grupo_de_risco == "LMA", "LMA-s",
                ifelse(red$Grupo_de_risco %in% c("Very Low", "Low", "Moderate Low"),
                       "SMD baixo risco", "SMD alto risco")), levels = grupos),
  `IPSS-M` = factor(ipssm_pt[red$Grupo_de_risco], levels = ipssm_pt),
  `Subtipo OMS 2022` = factor(subtipo_pt[trimws(red$Subtipos_SMD)], levels = subtipo_pt),
  Sexo      = red$Sexo,
  Idade     = as.numeric(red$Idade),
  `Faixa etária` = factor(ifelse(as.numeric(red$Idade) < 60, "< 60 anos", "≥ 60 anos")),
  Celularidade   = faltante(red$Celularidade),
  Fibrose        = aus_pres(red$Fibrose),
  `Sideroblastos em anel` = aus_pres(red$Sideroblastos),
  Diseritropoese = aus_pres(red$Diseritropoese),
  Disgranulopoese = aus_pres(red$Disgranulopoese),
  Dismegacariopoese = aus_pres(red$Dismegacariopoese),
  `Blastos (%)`  = as.numeric(red$Blastos),
  `Blastos, categorias` = factor(as.character(cut(as.numeric(red$Blastos), c(-Inf, 4.999, 9.999, 19.999, Inf),
                          c("< 5%", "5–9%", "10–19%", "≥ 20%"))), levels = c("< 5%", "5–9%", "10–19%", "≥ 20%")),
  `Hemoglobina (g/dL)` = as.numeric(red$Hemoglobina),
  `Neutrófilos segmentados (/µL)` = as.numeric(red$Segmentados),
  `Plaquetas (/µL)` = as.numeric(red$Plaquetas),
  Cariótipo = faltante(red$Resultado),
  `Cariótipo complexo` = ifelse(grepl("aplica", red$Kt_complexo), NA, aus_pres(red$Kt_complexo)),
  `Dependência transfusional` = faltante(red$Dependencia_transfusional),
  Óbito = ifelse(red$Obito == "Perda de seguimento", NA, faltante(red$Obito)),
  check.names = FALSE)

num_vars <- c("Idade", "Blastos (%)", "Hemoglobina (g/dL)", "Neutrófilos segmentados (/µL)", "Plaquetas (/µL)")
cat_vars <- c("Sexo", "Faixa etária", "IPSS-M", "Subtipo OMS 2022", "Celularidade", "Fibrose",
              "Sideroblastos em anel", "Diseritropoese", "Disgranulopoese", "Dismegacariopoese",
              "Blastos, categorias", "Cariótipo", "Cariótipo complexo", "Dependência transfusional", "Óbito")
so_smd <- c("IPSS-M", "Subtipo OMS 2022")          # não se aplicam à LMA-s
ordem  <- c("Sexo", "Idade", "Faixa etária", "IPSS-M", "Subtipo OMS 2022", "Blastos (%)", "Blastos, categorias",
            "Hemoglobina (g/dL)", "Neutrófilos segmentados (/µL)", "Plaquetas (/µL)", "Celularidade", "Fibrose",
            "Sideroblastos em anel", "Diseritropoese", "Disgranulopoese", "Dismegacariopoese",
            "Cariótipo", "Cariótipo complexo", "Dependência transfusional", "Óbito")

fmt_num <- function(x) if (all(is.na(x))) "—" else {
  q <- quantile(x, c(.5, .25, .75), na.rm = TRUE); dig <- if (max(x, na.rm = TRUE) > 1000) 0 else 1
  f <- function(v) formatC(v, format = "f", digits = dig, big.mark = ".", decimal.mark = ",")
  sprintf("%s (%s–%s)", f(q[1]), f(q[2]), f(q[3])) }
fmt_cat <- function(x, nivel) { ok <- !is.na(x); if (!any(ok)) return("—")
  k <- sum(x[ok] == nivel); sprintf("%d (%s%%)", k, formatC(100 * k / sum(ok), format = "f", digits = 1, decimal.mark = ",")) }

# Tabela descritiva + gabarito de testes
monta <- function(dd) {
  linhas <- list(); gab <- list()
  cab <- c("Variável", sprintf("Total (n = %d)", nrow(dd)),
           sprintf("%s (n = %d)", grupos, as.integer(table(dd$Grupo))), "p", "Teste", "Dados faltantes")
  for (v in ordem) {
    x <- dd[[v]]; nf <- sum(is.na(x))
    if (v %in% num_vars) {
      linhas[[length(linhas) + 1]] <- c(paste0(v, ", mediana (IIQ)"), fmt_num(x),
        sapply(grupos, function(g) fmt_num(x[dd$Grupo == g])), "", "Kruskal-Wallis", nf)
      ok <- !is.na(x)
      gab[[length(gab) + 1]] <- data.frame(Variável = v, Teste = "Kruskal-Wallis", n = sum(ok),
        p = kruskal.test(x[ok] ~ dd$Grupo[ok])$p.value, p_qui_quadrado = NA)
    } else {
      grp <- if (v %in% so_smd) grupos[1:2] else grupos
      niveis <- if (is.factor(x)) levels(droplevels(x)) else sort(unique(na.omit(x)))
      teste <- if (v %in% so_smd) "— (classificação define o grupo)" else "Qui-quadrado / Fisher"
      linhas[[length(linhas) + 1]] <- c(paste0(v, ", n (%)"), "", "", "", "", "", teste, nf)
      for (nv in niveis) linhas[[length(linhas) + 1]] <- c(paste0("   ", nv), fmt_cat(x, nv),
        sapply(grupos, function(g) if (g %in% grp) fmt_cat(x[dd$Grupo == g], nv) else "—"), "", "", "")
      if (!v %in% so_smd) {
        ok <- !is.na(x); tb <- table(x[ok], droplevels(dd$Grupo[ok]))
        set.seed(1234)
        gab[[length(gab) + 1]] <- data.frame(Variável = v, Teste = "Fisher (R, Monte Carlo)", n = sum(ok),
          p = fisher.test(tb, simulate.p.value = TRUE, B = 1e5)$p.value,
          p_qui_quadrado = suppressWarnings(chisq.test(tb)$p.value))
      }
    }
  }
  t <- as.data.frame(do.call(rbind, linhas), stringsAsFactors = FALSE); names(t) <- cab
  g <- do.call(rbind, gab); g$p <- signif(g$p, 2); g$p_qui_quadrado <- signif(g$p_qui_quadrado, 2)
  list(tabela = t, gabarito = g) }

c93 <- d
c43 <- d[d$Metiloma == "Sim", ]
t1 <- monta(c93); t2 <- monta(c43)

# ---------- Planilhas de dados para o programa de estatística ----------
# (a) formato longo: uma linha por paciente (JASP, jamovi, SPSS, R)
longo <- function(dd) dd[, c("ID_REDCAP", "Grupo", ordem)]
# (b) formato em colunas: uma coluna por grupo, para cada variável numérica (GraphPad Prism: "Column")
colunas <- function(dd, v) { l <- split(dd[[v]], dd$Grupo); m <- max(lengths(l))
  as.data.frame(lapply(l, function(z) c(z, rep(NA, m - length(z)))), check.names = FALSE) }
# (c) tabelas de contingência: contagens categoria × grupo (Prism: "Contingency")
conting <- function(dd, v) { x <- dd[[v]]; ok <- !is.na(x)
  tb <- as.data.frame.matrix(table(x[ok], droplevels(dd$Grupo[ok])))
  data.frame(Categoria = rownames(tb), tb, check.names = FALSE) }

wb <- createWorkbook()
est_tit <- createStyle(textDecoration = "bold", fontSize = 12)
est_cab <- createStyle(textDecoration = "bold", fgFill = "#DDEBF7", border = "bottom", wrapText = TRUE)
est_p   <- createStyle(fgFill = "#FFF2CC")
escreve_tab <- function(nome, titulo, t, nota) {
  addWorksheet(wb, nome); writeData(wb, nome, titulo, startRow = 1); addStyle(wb, nome, est_tit, 1, 1)
  writeData(wb, nome, t, startRow = 3, headerStyle = est_cab)
  addStyle(wb, nome, est_p, rows = 4:(nrow(t) + 3), cols = 6, gridExpand = TRUE)
  writeData(wb, nome, nota, startRow = nrow(t) + 5)
  setColWidths(wb, nome, 1, 38); setColWidths(wb, nome, 2:5, 20); setColWidths(wb, nome, 6:8, 18)
  freezePane(wb, nome, firstActiveRow = 4) }
nota <- c("IIQ = intervalo interquartil (1º–3º quartil). Percentuais calculados sobre os pacientes com o dado disponível.",
          "Óbito: pacientes com perda de seguimento contados como dado faltante.",
          "Cariótipo complexo: calculado apenas entre os pacientes com cariótipo avaliável.",
          "A coluna p (amarela) deve ser preenchida com o resultado do teste indicado na coluna Teste.")
escreve_tab("Tabela 1 - Painel (93)", "Tabela 1. Características clínicas e laboratoriais dos 93 pacientes com painel de sequenciamento",
            t1$tabela, nota)
escreve_tab("Tabela 2 - Metiloma (43)", "Tabela 2. Características clínicas e laboratoriais dos 43 pacientes com array de metilação",
            t2$tabela, c(nota, "Destes 43, 36 entraram na análise de metilação (5 excluídos por qualidade e 2 por perfil discrepante)."))

addWorksheet(wb, "Dados 93 (1 linha por paciente)"); writeData(wb, 3, longo(c93), headerStyle = est_cab)
addWorksheet(wb, "Dados 43 (1 linha por paciente)"); writeData(wb, 4, longo(c43), headerStyle = est_cab)

for (co in list(list("93", c93), list("43", c43))) {
  aba <- paste0("Prism colunas ", co[[1]]); addWorksheet(wb, aba); r <- 1
  for (v in num_vars) { writeData(wb, aba, v, startRow = r); addStyle(wb, aba, est_tit, r, 1)
    x <- colunas(co[[2]], v); writeData(wb, aba, x, startRow = r + 1, headerStyle = est_cab); r <- r + nrow(x) + 4 }
  aba <- paste0("Prism contingência ", co[[1]]); addWorksheet(wb, aba); r <- 1
  for (v in setdiff(cat_vars, so_smd)) { writeData(wb, aba, v, startRow = r); addStyle(wb, aba, est_tit, r, 1)
    x <- conting(co[[2]], v); writeData(wb, aba, x, startRow = r + 1, headerStyle = est_cab); r <- r + nrow(x) + 4 }
}

# ---------- Passo a passo no GraphPad Prism ----------
passos <- c(
  "COMO PREENCHER A COLUNA p NO GRAPHPAD PRISM",
  "",
  "Antes de começar: faça um arquivo do Prism para a Tabela 1 (93) e outro para a Tabela 2 (43). Cada variável vira uma tabela de dados dentro do arquivo.",
  "Use as abas que terminam em 93 para a Tabela 1 e as que terminam em 43 para a Tabela 2.",
  "",
  "PARTE A. VARIÁVEIS NUMÉRICAS (idade, blastos, hemoglobina, neutrófilos, plaquetas) -> teste de Kruskal-Wallis",
  "1. No Prism: File > New > New Data Table and Graph. Escolha o tipo 'Column'.",
  "2. Em Data table, marque 'Enter or import data into a new table'. Em Options, marque 'Enter replicate values, stacked into columns'. Clique em Create.",
  "3. Na aba 'Prism colunas 93' (ou 43), copie as três colunas da variável, junto com os títulos (SMD baixo risco, SMD alto risco, LMA-s).",
  "4. No Prism, clique na primeira célula do título do Grupo A e cole (Ctrl+V / Cmd+V). As três colunas devem aparecer como Grupo A, B e C. Células vazias são normais: os grupos têm tamanhos diferentes.",
  "5. Clique em Analyze > Column analyses > 'One-way ANOVA (and nonparametric or mixed)' > OK.",
  "6. Na aba Experimental Design: marque 'No matching or pairing' e, em 'Assume Gaussian distribution?', marque 'No. Use nonparametric test'.",
  "7. Na aba Nonparametric test (ou Options), confirme 'Kruskal-Wallis test'. Na aba Multiple Comparisons pode deixar 'No multiple comparisons'. Clique em OK.",
  "8. Nos resultados, procure a linha 'P value' do bloco 'Kruskal-Wallis test'. Esse é o p da tabela.",
  "9. (Conferência) Analyze > Column analyses > Descriptive statistics mostra mediana e quartis. Os valores devem bater com a tabela.",
  "10. Para a próxima variável: no painel da esquerda, clique em '+ New...' > New Data Table, e repita os passos 1 a 8.",
  "",
  "PARTE B. VARIÁVEIS CATEGÓRICAS (sexo, cariótipo, transfusão, óbito etc.) -> qui-quadrado",
  "1. File > New > New Data Table and Graph. Escolha o tipo 'Contingency'. Marque 'Enter or import data into a new table' e clique em Create.",
  "2. Na aba 'Prism contingência 93' (ou 43), copie o bloco da variável inteiro (coluna Categoria + as três colunas de contagem, com os títulos).",
  "3. No Prism, clique na célula de título das linhas (canto superior esquerdo) e cole. As categorias ficam nas linhas e os grupos nas colunas. Só aparecem números inteiros (contagens), nunca porcentagens.",
  "4. Analyze > Contingency table analyses > 'Chi-square (and Fisher's exact) test' > OK.",
  "5. O Prism só faz o teste exato de Fisher em tabelas 2 x 2. Como aqui sempre há três grupos, ele vai fazer o qui-quadrado. Isso é esperado; deixe assim.",
  "6. Nos resultados, copie o 'P value'. Se aparecer um aviso de que algumas caselas têm valor esperado menor que 5, anote quais variáveis tiveram esse aviso.",
  "",
  "PARTE C. COMO ESCREVER O p NA TABELA",
  "- Use vírgula decimal e três casas: p = 0,049. Se o Prism mostrar algo como < 0,0001, escreva p < 0,001.",
  "- IPSS-M e Subtipo OMS 2022 não têm teste: eles são a própria definição dos grupos.",
  "- Anote embaixo da tabela: 'Variáveis numéricas: teste de Kruskal-Wallis; variáveis categóricas: teste do qui-quadrado.'",
  "",
  "DÚVIDAS COMUNS",
  "- 'Meu p ficou um pouco diferente do da orientadora.' Para o qui-quadrado, diferenças pequenas podem acontecer se outra versão do teste foi usada. Mande a variável e o valor que encontrou.",
  "- 'O Prism não aceitou colar.' Confira se você escolheu o tipo certo de tabela (Column para números, Contingency para categorias).",
  "- 'Apareceu uma coluna a mais.' Verifique se não copiou a coluna em branco entre os blocos.")
addWorksheet(wb, "Como fazer no Prism", tabColour = "#FFC000")
writeData(wb, "Como fazer no Prism", data.frame(Passo = passos), colNames = FALSE)
addStyle(wb, "Como fazer no Prism", createStyle(wrapText = TRUE, valign = "top"), rows = seq_along(passos), cols = 1)
addStyle(wb, "Como fazer no Prism", createStyle(textDecoration = "bold", fontSize = 12),
         rows = grep("^[A-ZÇÃÉÚ ]{8,}|^PARTE|^COMO|^DÚVIDAS", passos), cols = 1, stack = TRUE)
setColWidths(wb, "Como fazer no Prism", 1, 120)
worksheetOrder(wb) <- c(length(names(wb)), seq_len(length(names(wb)) - 1))
activeSheet(wb) <- 1
saveWorkbook(wb, "Tabelas_descritivas_coorte.xlsx", overwrite = TRUE)

g <- createWorkbook()
for (x in list(list("Gabarito 93", t1$gabarito), list("Gabarito 43", t2$gabarito))) {
  addWorksheet(g, x[[1]]); writeData(g, x[[1]], x[[2]], headerStyle = est_cab); setColWidths(g, x[[1]], 1:5, c(32, 24, 6, 10, 16)) }
saveWorkbook(g, "GABARITO_testes_tabelas.xlsx", overwrite = TRUE)

cat("== Tabela 1 (93) ==\n"); print(t1$tabela, row.names = FALSE)
cat("\n== Gabarito 93 ==\n"); print(t1$gabarito, row.names = FALSE)
cat("\n== Gabarito 43 ==\n"); print(t2$gabarito, row.names = FALSE)
