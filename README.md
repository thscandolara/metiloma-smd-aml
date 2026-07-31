# metiloma-smd-aml

Pipeline reprodutível de análise de metilação de DNA (Illumina Infinium MethylationEPIC v2.0, EPICv2, hg38) comparando **leucemia mieloide aguda (AML)** e **neoplasia mielodisplásica (MDS)** de alto e baixo risco.

## Objetivo
Explorar diferenças de metilação (DMPs e DMRs) entre AML e MDS. Dado o tamanho amostral pequeno (n≈37), o estudo é **exploratório**: o foco é caracterizar as fontes de variação e descrever padrões, não construir classificadores.

## Pipeline (visão geral)
`sesame::openSesame` (noob + correção de dye bias + detecção pOOBAH) → QC (SNPs de controle, densidade, MDS) → filtro de sondas (detecção, SNP, cross-reactive, chrX/Y) → colapso de réplicas EPICv2 → **checagem de efeito de lote e composição celular (PCA × variáveis)** → `limma` (DMPs, contrastes pareados) → `DMRcate` (DMRs).

## Estrutura do repositório
- `analysis.qmd` — pipeline narrado em Quarto (documento principal)
- `data/` — metadados clínicos de-identificados e dicionário (idats **não** versionados)
- `results/` — figuras e tabelas geradas
- `CITATION.cff`, `LICENSE`, `.gitignore`

## Como reproduzir
1. Instalar R (>= 4.5) e os pacotes Bioconductor listados no topo do `analysis.qmd`.
2. Colocar os idats em `data/idats/` (não versionados — ver Disponibilidade de dados).
3. Renderizar: `quarto render analysis.qmd`.

## Disponibilidade de dados
Os arquivos `.idat` e os dados clínicos identificáveis de pacientes **não são versionados** (ver `.gitignore`), por razões éticas e de privacidade. Uma versão de-identificada dos metadados acompanha o repositório; dados brutos ficam sob acesso controlado / depósito apropriado (ex.: GEO).

## Referências principais
- Peters TJ et al. (2024) Characterisation and reproducibility of the HumanMethylationEPIC v2.0 BeadChip. *BMC Genomics* 25:251.
- Cabezón M et al. (2021) DNA methylation in high-risk MDS/sAML. *Clinical Epigenetics*.
- Figueroa ME et al. (2010) DNA methylation subtypes of AML. *Cancer Cell*.

## Licença
Código sob licença MIT (ver `LICENSE`).
