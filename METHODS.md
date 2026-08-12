# Decisões metodológicas do pipeline (metiloma AML vs MDS, EPICv2)

Justificativa de cada escolha, com a base (padrão do pacote, literatura ou julgamento) e a compatibilidade com a literatura. Serve de "Métodos" para tese e relatório. Escolhas marcadas como *julgamento* são parâmetros ajustáveis, a reportar explicitamente.

## 1. Array e anotação
EPICv2, hg38 (`IlluminaHumanMethylationEPICv2manifest`, `...anno.20a1.hg38`). Base: array utilizado; Peters et al. 2024.

## 2. Pré-processamento (sesame, `openSesame` "QCDPB")
noob + correção de dye bias + detecção **pOOBAH** + qualityMask, por amostra. pOOBAH (sinal out-of-band) é mais robusto que o detection p-value baseado em controles negativos, que é anticonservador em EPIC (Zhou et al. 2018; Heiss et al. 2019). Avaliação sistemática aponta sesame como melhor e quantile como pior (Welsh et al. 2023). Trade-off: retém menos sondas que o minfi, porém mais confiáveis.

## 3. QC de amostra e remoção de outliers
Exclusão por detecção < 80% (>20% de sondas mascaradas). Adicionalmente, remoção de outliers identificados por densidade de beta aberrante e por posição extrema no PCA (duas amostras: uma MDS_LOW com distribuição atípica; um AML extremo). *Julgamento; limiares ajustáveis.*

## 4. Filtragem de sondas (com contagem por etapa)
Detecção ≥95% (não `complete.cases`), remoção de SNP-próximas e cross-reactive (`rmSNPandCH`), remoção de chrX/chrY (via anotação) e colapso de réplicas EPICv2 pela média. Base: Chen et al. 2013; Hop et al. 2020; Peters et al. 2024; vignette EPICv2 do DMRcate.

## 5. Escala estatística
Modelagem em M-value; beta para visualização (Du et al. 2010).

## 6. Modelo, sexo e efeito de lote
Design sem intercepto ajustado por **sexo registrado** (não predito — LOY é comum em mieloides e distorce a predição) e por **chip/lote como covariável**. O lote é modelado no design em vez de removido com ComBat, pois a remoção prévia não reduz corretamente os graus de liberdade, gerando p-valores anticonservadores (Nygaard et al. 2016). Verificou-se, por tabela cruzada, que grupo e chip não estão confundidos. O **% de blastos não é ajustado** (está na definição de AML, na via causal).

## 7. Contrastes
Dois modelos, ambos ajustados por sexo e lote: (a) **AML vs MDS** com fator de 2 níveis (MDS agregado, ponderado pelo n — doença vs doença, sem pesos); (b) modelo de 3 níveis do qual se extraem os pareados **AML vs MDS-high**, **AML vs MDS-low** e **MDS-high vs MDS-low**. Um modelo único por grupo aproveita a variância compartilhada via eBayes, com mais poder que subdividir os dados (Smyth 2004; Ritchie et al. 2015).

## 8. DMPs, tamanhos de efeito e diagnósticos
Tabelas completas com `logFC`, **delta-beta** (diferença de beta médio), p nominal e FDR; **candidatos exploratórios** por limiares pré-especificados (P<0,001 e |Δβ|≥0,10), explicitamente não confirmatórios; diagnósticos de poder/dispersão (histograma de p, QQ-plot, distribuição de |Δβ|, volcano).

## 9. DMRs
DMRcate (`cpg.annotate` com `epicv2Filter="mean"`), FDR de semeadura mantido em 0,05 (rigoroso). Como DMRs são construídas a partir de CpGs significativos, na ausência de DMPs não há DMRs; a exploração de regiões nominais exigiria relaxar o FDR (opção não adotada). Peters et al. 2015; 2024.

## 10. Região genômica e ilha CpG
% de sondas metiladas por grupo gênico e por relação com ilha CpG (rótulos reais do EPICv2: 5UTR/3UTR/exon/intron/TSS; Island/Shore/Shelf/OpenSea); distribuição das regiões no array vs no conjunto de maior sinal.

## 11. Exploração não supervisionada e subtipos
Sondas mais variáveis → clustering hierárquico (Ward) e heatmap anotado; **consensus clustering** (`ConsensusClusterPlus`) para número/estabilidade de subtipos (Monti et al. 2003; Jian et al. 2023). Associação dos clusters com variáveis clinicopatológicas (Kruskal-Wallis/Fisher; exploratório, sem correção múltipla). Caracterização dos clusters por limma cluster-vs-resto (genoma todo, para reduzir circularidade), com exportação de CpGs e listas de genes por cluster para enriquecimento externo. Abordagem análoga à definição de "epitipos" de AML (Figueroa et al. 2010; Giacopelli et al. 2021).

## 12. Sensibilidade e confundidores
Refit com e sem outliers; modelos secundários com covariáveis (ex.: blastos); quantificação do confundimento lote×grupo (Cramér's V; MDS por grupo e chip); estimativa de composição celular (EpiDISH, referência de sangue — aproximação, pois medula não tem referência padrão).

## 13. Enquadramento
Estudo exploratório. A ausência de DMP/DMR genome-wide reflete poder limitado (n pequeno, heterogeneidade de blastos), não ausência de biologia (Cabezón et al. 2021). O achado central é a estratificação metilômica por gravidade/desfecho (blastos, dependência transfusional, óbito), independente do diagnóstico.

## Referências
- Zhou W et al. NAR 2018;46(20):e123 (SeSAMe/pOOBAH).
- Welsh H et al. Clin Epigenetics 2023 (normalização EPIC).
- Heiss JA, Just AC. Clin Epigenetics 2019;11:15.
- Chen YA et al. Epigenetics 2013;8(2):203-209.
- Hop PJ et al. NAR Genom Bioinform 2020;2(4).
- Peters TJ et al. BMC Genomics 2024;25:251; e Epigenetics & Chromatin 2015;8:6 (DMRcate).
- Du P et al. BMC Bioinformatics 2010;11:587 (M vs beta).
- Ritchie ME et al. NAR 2015;43(7):e47 (limma); Smyth GK. Stat Appl Genet Mol Biol 2004;3.
- Nygaard V, Rødland EA, Hovig E. Biostatistics 2016;17(1):29-39 (batch no modelo).
- Monti S et al. Machine Learning 2003 (consensus clustering).
- Figueroa ME et al. Cancer Cell 2010; Giacopelli B et al. Genome Research 2021; Jian J et al. Clin Exp Med 2023 (subtipos/epitipos de AML).
- Cabezón M et al. Clin Epigenetics 2021 (MDS/sAML).
