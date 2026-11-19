# ============================================================
# Tradução para o português dos textos que o maftools desenha nas figuras.
# O maftools não tem opção de idioma: os rótulos (eixos, legendas, títulos) estão
# escritos dentro das funções. Este arquivo troca apenas esses textos de exibição,
# sem mexer em nenhum cálculo. Uso: source("traducao_maftools.R") depois de library(maftools).
# ============================================================

# Classes de variante: renomeadas no próprio MAF antes do read.maf (legenda do oncoplot/lollipop)
classes_pt <- c(
  Missense_Mutation      = "Missense",
  Nonsense_Mutation      = "Nonsense",
  Frame_Shift_Del        = "Deleção frameshift",
  Frame_Shift_Ins        = "Inserção frameshift",
  In_Frame_Del           = "Deleção in-frame",
  In_Frame_Ins           = "Inserção in-frame",
  Splice_Site            = "Sítio de splicing",
  Translation_Start_Site = "Perda do códon inicial",
  Nonstop_Mutation       = "Perda do códon de parada")
cores_classes <- c(
  "Missense" = "#33A02C", "Nonsense" = "#E31A1C", "Deleção frameshift" = "#1F78B4",
  "Inserção frameshift" = "#6A3D9A", "Deleção in-frame" = "#B2DF8A", "Inserção in-frame" = "#CAB2D6",
  "Sítio de splicing" = "#FF7F00", "Perda do códon inicial" = "#FB9A99",
  "Perda do códon de parada" = "#A6CEE3", "Múltiplas variantes" = "black")

# Textos fixos dentro das funções do maftools -> português
# (só textos de exibição; nomes de colunas internas não são tocados)
textos_pt <- c(
  "Multi_Hit"                              = "Múltiplas variantes",
  "TMB"                                    = "Nº de variantes",
  " (N = "                                 = " (com variante: ",
  "Domains"                                = "Domínios",
  "No. of samples"                         = "Nº de amostras",
  "Altered in "                            = "Alterado em ",
  "%) of "                                 = "%) de ",
  " samples."                              = " amostras.",
  " : [Somatic Mutation Rate: "            = " : [Taxa de mutação: ",
  "Fraction of variants within clusters"   = "Fração das variantes em agrupamentos",
  "-log10 (fdr)"                           = "-log10 (FDR)",
  "> 3 (Mutually exclusive)"               = "> 3 (Exclusividade mútua)",
  ">3 (Co-occurence)"                      = "> 3 (Coocorrência)",
  "-log10(P-value)"                        = "-log10 (valor de p)",
  "-log10(fdr)"                            = "-log10 (FDR)",
  " P < "                                  = " p < ",
  " fdr < "                                = " FDR < ",
  "Variant Classification"                 = "Classe da variante",
  "Variant Type"                           = "Tipo de variante",
  "SNV Class"                              = "Classe de SNV",
  "Variants per sample"                    = "Variantes por amostra",
  "Mean: "                                 = "Média: ",
  "Median: "                               = "Mediana: ",
  "Variant Classification \nsummary"       = "Resumo das classes\nde variante",
  "No. of alterations"                     = "Nº de alterações",
  "Percent of cases"                       = "Percentual de casos",
  "Number of cases"                        = "Número de casos",
  "P-value"                                = "Valor de p",
  "Fraction of pathway affected"           = "Fração da via alterada",
  "Druggable categories"                   = "Categorias com fármacos",
  "# Genes"                                = "Nº de genes")

traduz_maftools <- function(textos = textos_pt) {
  ns <- asNamespace("maftools")
  alvo <- paste0("\"", names(textos), "\"")
  for (f in ls(ns)) {
    fun <- get(f, envir = ns)
    if (!is.function(fun) || is.primitive(fun)) next
    cod <- deparse(fun, width.cutoff = 500L, control = c("keepNA", "keepInteger", "niceNames", "showAttributes"))
    ach <- vapply(alvo, function(a) any(grepl(a, cod, fixed = TRUE)), logical(1))
    if (!any(ach)) next
    for (i in which(ach)) cod <- gsub(alvo[i], paste0("\"", textos[[i]], "\""), cod, fixed = TRUE)
    novo <- tryCatch(eval(parse(text = cod, encoding = "UTF-8"), envir = ns), error = function(e) NULL)
    if (is.null(novo)) { message("  [tradução] não foi possível traduzir ", f); next }
    environment(novo) <- environment(fun)
    utils::assignInNamespace(f, novo, ns = "maftools")
    # funções exportadas também têm uma cópia no pacote anexado (library): trocar lá também
    pkg <- "package:maftools"
    if (pkg %in% search() && exists(f, envir = as.environment(pkg), inherits = FALSE)) {
      env <- as.environment(pkg); unlockBinding(f, env); assign(f, novo, envir = env); lockBinding(f, env)
    }
  }
  invisible(TRUE)
}
