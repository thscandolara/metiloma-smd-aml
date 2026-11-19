# metiloma-smd-aml

Pipeline reprodutível de análise de metilação de DNA (Illumina Infinium MethylationEPIC v2.0, EPICv2, hg38) comparando **leucemia mieloide aguda (AML)** e **neoplasia mielodisplásica (MDS)** de alto e baixo risco, em medula óssea.

## Objetivo
Explorar diferenças de metilação (DMPs e DMRs) entre AML e MDS e caracterizar as fontes de variação. Dado o tamanho amostral pequeno, o estudo é **exploratório e gerador de hipóteses**, não confirmatório.

## Pipeline (`analysis.qmd`)
1. **Pré-processamento (sesame):** `openSesame` (QCDPB — noob + correção de dye bias + detecção pOOBAH).
2. **QC de amostra:** fração de sondas detectadas (pOOBAH); exclusão por baixa detecção e **remoção de outliers** (densidade aberrante / PCA).
3. **Filtragem de sondas** (com contagem a cada passo): detecção ≥95%, remoção de SNP e cross-reactive (`rmSNPandCH`), remoção de chrX/chrY, e colapso de réplicas EPICv2 pela média.
4. **Fontes de variação:** PCA/MDS e associação com variáveis clínicas/técnicas; efeito de **lote (chip)** modelado como covariável.
5. **Análise diferencial (limma):** contrastes **AML vs MDS (agregado)**, **AML vs MDS-high**, **AML vs MDS-low** e **MDS-high vs MDS-low**. DMPs com `logFC`, **delta-beta**, p nominal e FDR; tabelas de **candidatos exploratórios** (P<0,001 e |Δβ|≥0,10); diagnósticos (histograma de p, QQ-plot, distribuição de |Δβ|, volcano).
6. **DMRs (DMRcate).**
7. **Metilação por região genômica** (promotor, 5'UTR, corpo do gene, 3'UTR, intergênico) e por **relação com ilha CpG** (island/shore/shelf/open sea); distribuição array vs estudo.
8. **Exploração não supervisionada:** heatmap anotado, clustering hierárquico, **consensus clustering** (`ConsensusClusterPlus`); associação dos clusters com variáveis clínicas; caracterização dos clusters (CpGs e listas de genes por cluster, exportadas para enriquecimento externo).
9. **Sensibilidade e diagnósticos:** com/sem outliers, modelos secundários com covariáveis, confundimento lote×grupo, estimativa de composição celular (EpiDISH).

## Reanálise clínica com o banco conferido no REDCap (nov/2026) — `reanalise_clinica_REDCap/`
A parte de metilação (QC, filtros, limma, DMRcate, agrupamentos) não depende dos dados clínicos e não muda. Estes scripts refazem só o que usa a clínica, com o banco conferido no REDCap e tudo em português:
1. `00_montar_pdata_PT.R` — monta o pdata em português (grupos SMD baixo risco / SMD alto risco / LMA-s; IPSS-M; achados medulares; cariótipo; dependência transfusional; óbito, com perda de seguimento como dado ausente).
2. `01_reanalise_clinica_PT.R` — repete QC/filtros/agrupamentos com os mesmos parâmetros e semente do `analysis.qmd`, testa clusters × clínica (Kruskal-Wallis; Fisher com Monte Carlo, B = 100.000) lado a lado com o pdata anterior, R² dos PCs, heatmap anotado, sensibilidade com blastos e EpiDISH.
3. `02_tabelas_descritivas.R` — tabelas descritivas da coorte (93 do painel; 43 do array), com planilhas prontas para o GraphPad Prism e gabarito dos testes.

## Painel NGS (arquivo MAF) — `painel_ngs/`
- `graficos_maf.R` — maftools nos 10 genes do painel (93 pacientes, variantes codificantes). Testes entre grupos com o denominador completo (pacientes sem variante incluídos). `Rscript graficos_maf.R 74` refaz tudo sem o paciente 74 (análise de sensibilidade).
- `traducao_maftools.R` — figuras em português: classes de variante renomeadas e textos fixos do maftools traduzidos (só rótulos; nenhum cálculo é alterado).

## Estrutura do repositório
- `analysis.qmd` — pipeline narrado (Quarto)
- `METHODS.md` — decisões metodológicas e justificativas com referências
- `reanalise_clinica_REDCap/` — reanálise clínica com o banco do REDCap (em português)
- `painel_ngs/` — análise do painel de sequenciamento (maftools)
- `output/tabelas/` e `output/figuras/` — tabelas (CSV) e figuras (PNG/PDF)
- `data/`, `.gitignore`, `LICENSE`, `CITATION.cff`

## Como reproduzir
R (>= 4.5, Bioconductor) e Quarto. Pacotes: `maftools`, `readxl`, `openxlsx`, `dplyr`, `tidyr`, `sesame`, `minfi`, `DMRcate`, `limma`, `ConsensusClusterPlus`, `pheatmap`, `EpiDISH`, `IlluminaHumanMethylationEPICv2manifest`, `IlluminaHumanMethylationEPICv2anno.20a1.hg38`. Ajuste os caminhos no bloco `setup` e rode `quarto render analysis.qmd`. O HTML gerado é autocontido (compartilhável).

## Disponibilidade de dados
Arquivos IDAT e dados clínicos identificáveis **não são versionados** (ver `.gitignore`), por razões éticas e de privacidade.

## Licença
Código sob licença MIT (ver `LICENSE`).
