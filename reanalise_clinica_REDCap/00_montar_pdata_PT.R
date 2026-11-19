# ============================================================
# Monta o pdata do metiloma em PORTUGUÊS a partir do banco conferido no REDCap
# (Banco_Atualizado_REDCAP_19112026.xlsx). Chip/posição/barcode vêm do pdata original.
# Uso: Rscript 00_montar_pdata_PT.R      (rodar de dentro desta pasta)
# Saída: pdata_pronon_REDCap_PT.csv
# ============================================================
suppressMessages(library(readxl))
invisible(Sys.setlocale("LC_ALL", "en_US.UTF-8"))   # acentos corretos no CSV

xlsx  <- "../dados/Banco_Atualizado_REDCAP_19112026.xlsx"   # ajuste para o caminho local do banco clínico (não versionado)
tmp   <- tempfile(fileext = ".xlsx"); file.copy(xlsx, tmp)   # readxl falha com acentos no caminho
red   <- as.data.frame(read_excel(tmp))
red   <- red[red$Methylation_Array == "Sim" & !is.na(red$ID_DILUICAO), ]
pd    <- read.csv("../pdata_pronon.csv")
stopifnot(nrow(red) == 43, setequal(red$ID_DILUICAO, pd$Sample_ID))
red   <- red[match(pd$Sample_ID, red$ID_DILUICAO), ]

# "Dado não informado", "Ignorado" etc. -> NA
faltante <- function(x) { x <- trimws(as.character(x))
  x[grepl("informad|ignorad|realizad", x, ignore.case = TRUE) | x == ""] <- NA; x }
aus_pres <- function(x) { x <- faltante(x)
  ifelse(is.na(x), NA, ifelse(grepl("^Aus", x), "Ausente", ifelse(grepl("^Pres", x), "Presente", x))) }

grupo <- c(Low = "SMD baixo risco", High = "SMD alto risco",
           "Very High" = "SMD alto risco", LMA = "LMA-s")[red$Grupo_de_risco]
risco <- c(Low = "Baixo", High = "Alto", "Very High" = "Muito alto", LMA = NA)[red$Grupo_de_risco]
subtipo_pt <- c(
  "MDS with low blasts (MDS-LB)"                        = "SMD com baixa contagem de blastos (SMD-BB)",
  "MDS with low blasts and SF3B1 mutationa (MDS-SF3B1)" = "SMD com baixa contagem de blastos e mutação de SF3B1 (SMD-SF3B1)",
  "MDS with increased blasts (MDS-IB) 1"                = "SMD com aumento de blastos 1 (SMD-AB1)",
  "MDS with increased blasts (MDS-IB) 2"                = "SMD com aumento de blastos 2 (SMD-AB2)",
  "MDS, hypoplasticb (MDS-h)"                           = "SMD hipoplásica (SMD-h)")
sub_raw <- trimws(red$Subtipos_SMD)
subtipo <- ifelse(grepl("Leucemia", sub_raw), "LMA secundária", subtipo_pt[sub_raw])

blastos <- as.numeric(red$Blastos)
cario   <- faltante(red$Cariotipo)
sexo_dna <- read.csv("CHECAGEM_sexo_metilacao_vs_planilha.csv")
sexo_dna <- sexo_dna[, c(grep("Sample_ID", names(sexo_dna))[1], grep("METILACAO", names(sexo_dna))[1])]
names(sexo_dna) <- c("ID", "Sexo_DNA")

obito <- faltante(red$Obito)

out <- data.frame(
  Barcode               = pd$Barcode,
  Chip                  = pd$Sentrix_ID,
  Posicao_chip          = pd$Sentrix_Position,
  ID_amostra            = pd$Sample_ID,
  ID_REDCAP             = red$ID_REDCAP,
  Grupo                 = factor(grupo, levels = c("SMD baixo risco", "SMD alto risco", "LMA-s")),
  Risco_IPSSM           = risco,
  Subtipo_OMS           = subtipo,
  Sexo                  = red$Sexo,
  Sexo_DNA              = c(Male = "Masculino", Female = "Feminino")[sexo_dna$Sexo_DNA[match(pd$Sample_ID, sexo_dna$ID)]],
  Idade                 = as.numeric(red$Idade),
  Faixa_etaria          = ifelse(as.numeric(red$Idade) < 60, "<60", ">=60"),
  Celularidade          = faltante(red$Celularidade),
  Fibrose               = aus_pres(red$Fibrose),
  Sideroblastos_em_anel = aus_pres(red$Sideroblastos),
  Diseritropoese        = aus_pres(red$Diseritropoese),
  Disgranulopoese       = aus_pres(red$Disgranulopoese),
  Dismegacariopoese     = aus_pres(red$Dismegacariopoese),
  Blastos_pct           = blastos,
  Categoria_blastos     = as.character(cut(blastos, c(-Inf, 4.999, 9.999, 19.999, Inf),
                                           c("<5", "5-9", "10-19", ">=20"))),
  Hemoglobina           = as.numeric(red$Hemoglobina),
  Segmentados           = as.numeric(red$Segmentados),
  Plaquetas             = as.numeric(red$Plaquetas),
  Cariotipo_ISCN        = cario,
  Cariotipo             = faltante(red$Resultado),
  Cariotipo_complexo    = ifelse(grepl("aplica", red$Kt_complexo), "Não se aplica", aus_pres(red$Kt_complexo)),
  Perda_Y               = ifelse(is.na(cario), NA, ifelse(grepl("-Y", cario), "Sim", "Não")),
  Dependencia_transfusional = faltante(red$Dependencia_transfusional),
  Obito                 = obito,
  FISH                  = faltante(red$FISH),
  Data_diagnostico      = red$Data_do_diagnostico,
  Data_obito            = red$Data_do_obito,
  check.names = FALSE)

write.csv(out, "pdata_pronon_REDCap_PT.csv", row.names = FALSE, fileEncoding = "UTF-8")
cat("pdata_pronon_REDCap_PT.csv:", nrow(out), "amostras,", ncol(out), "colunas\n")
for (v in c("Grupo","Risco_IPSSM","Subtipo_OMS","Sexo","Celularidade","Fibrose","Sideroblastos_em_anel",
            "Diseritropoese","Disgranulopoese","Dismegacariopoese","Categoria_blastos","Cariotipo",
            "Cariotipo_complexo","Perda_Y","Dependencia_transfusional","Obito","FISH")) {
  cat("\n--", v, "\n"); print(table(out[[v]], useNA = "ifany")) }
