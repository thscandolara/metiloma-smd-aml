# metiloma-smd-aml

Pipeline reprodutível de análise de metilação de DNA (Illumina Infinium MethylationEPIC v2.0, EPICv2, hg38) comparando **leucemia mieloide aguda (AML)** e **síndrome mielodisplásica (MDS)** de alto e baixo risco, em medula óssea.

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

## Estrutura do repositório
- `analysis.qmd` — pipeline narrado (Quarto)
- `METHODS.md` — decisões metodológicas e justificativas com referências
- `output/tabelas/` e `output/figuras/` — tabelas (CSV) e figuras (PNG/PDF)
- `data/`, `.gitignore`, `LICENSE`, `CITATION.cff`

## Como reproduzir
R (>= 4.5, Bioconductor) e Quarto. Pacotes: `sesame`, `minfi`, `DMRcate`, `limma`, `ConsensusClusterPlus`, `pheatmap`, `EpiDISH`, `IlluminaHumanMethylationEPICv2manifest`, `IlluminaHumanMethylationEPICv2anno.20a1.hg38`. Ajuste os caminhos no bloco `setup` e rode `quarto render analysis.qmd`. O HTML gerado é autocontido (compartilhável).

## Disponibilidade de dados
Arquivos IDAT e dados clínicos identificáveis **não são versionados** (ver `.gitignore`), por razões éticas e de privacidade.

## Licença
Código sob licença MIT (ver `LICENSE`).
