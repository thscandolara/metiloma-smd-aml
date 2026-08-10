# Decisões metodológicas do pipeline (metiloma AML vs MDS, EPICv2)

Documento de justificativa de cada escolha do pipeline, com a base (padrão do pacote,
literatura, ou julgamento) e a **compatibilidade com a literatura**. Serve de "Métodos"
para a tese e para o relatório PRONON. Escolhas marcadas como *julgamento* são
parâmetros ajustáveis que devem ser reportados explicitamente, não fatos.

## 1. Array e anotação
**Decisão:** EPICv2, genoma hg38 (`IlluminaHumanMethylationEPICv2manifest`, `...anno.20a1.hg38`).
**Base:** array utilizado; build nativa do EPICv2. Peters et al. 2024 (BMC Genomics 25:251).

## 2. Pré-processamento — `openSesame(prep = "QCDPB")`
**Decisão:** qualityMask (Q) + inferência de canal (C) + dye bias não-linear (D) + pOOBAH (P) + noob (B), tudo por amostra.
**Por quê:** ecossistema único, nativo de EPICv2; pOOBAH mais robusto que o `detectionP` do minfi; noob corrige fundo.
**Compatibilidade com a literatura:** ALTA. Welsh et al. 2023 (Clinical Epigenetics) avaliaram sistematicamente métodos de normalização em EPIC e encontraram o SeSAMe (com rodada extra de pOOBAH, "SeSAMe 2") como o **melhor**, e os métodos **quantile como os piores** — validando tanto o uso do sesame quanto o abandono do quantile. pOOBAH: Zhou et al. 2018 (NAR). noob: Triche et al. 2013 (NAR).
**Trade-off consciente:** abre mão da functional normalization (Funnorm; Fortin et al. 2014). O sesame normaliza dentro da amostra (noob + dye bias), sem forçar distribuições iguais (ao contrário do quantile), o que é adequado a dado com diferença global.

## 3. QC de amostra — excluir detecção < 80% (>20% de sondas mascaradas)
**Decisão:** fração de sondas mascaradas por pOOBAH como métrica de qualidade da amostra; excluir amostras acima de 20% (equivalente a `sesameQC_calcStats` frac_dt < 0.80). Removidas: Pacientes 4, 9, 17, 20, 24.
**Por quê:** análogo ao mean detection p-value do minfi. Limiar de 0,20 escolhido pela quebra natural no dado (amostras ruins com 31–55% mascarado) preservando n e o menor grupo (MDS_LOW). *Julgamento — parâmetro ajustável.*
**Compatibilidade com a literatura:** os p-valores de detecção convencionais são insuficientes; abordagens baseadas em background/out-of-band são preferíveis (Heiss et al. 2019, Clinical Epigenetics; Zhou et al. 2018). O limiar exato é definido pelo usuário por convenção.

## 4. QC de sonda — manter sondas detectadas em ≥95% das amostras (≤5% NA)
**Decisão:** após excluir amostra ruim, remover sondas com >5% de NA (pOOBAH) entre as boas. NÃO usar `complete.cases` (descartaria sonda falha em uma única amostra — rígido demais com pOOBAH).
**Por quê:** meio-termo comum em EWAS; evita perda excessiva de sondas por causa de poucas amostras.
**Compatibilidade com a literatura:** filtro por detection p-value é passo essencial e recomendado (Heiss et al. 2019). Muitas sondas EPIC têm baixa reprodutibilidade (ICC<0.5), sobretudo com beta perto de 0/1 e baixa variância (Welsh et al. 2023) — reforça filtrar sonda não confiável. *O corte de 5% é convenção ajustável.*

## 5. Filtro de sondas — `rmSNPandCH` (SNP + cross-reactive + cromossomos sexuais)
**Decisão:** remover sondas próximas a SNP, cross-reactive (lista EPICv2) e de chrX/chrY.
**Por quê:** SNP → sinal de genótipo, não metilação; cross-reactive → hibridização em múltiplos locais (sinal ambíguo); XY → sinal de sexo (com desbalanço de sexo, gera falso-positivo); LOY fica preservado no cariótipo.
**Compatibilidade com a literatura:** ALTA. Chen et al. 2013 (Epigenetics): ~6% das sondas co-hibridizam e podem gerar falso sinal autossômico associado a sexo; sondas polimórficas refletem SNP. Hop et al. 2020 (NAR G&B): cautela com cross-reatividade em EWAS. EPICv2: Peters et al. 2024.

## 6. Colapso de réplicas — `epicv2Filter = "mean"`
**Decisão:** colapsar réplicas EPICv2 pela média (dentro do `cpg.annotate`).
**Por quê:** EPICv2 mede vários CpGs com múltiplas sondas; DMRcate precisa de 1 sonda/CpG senão o kernel enviesa.
**Compatibilidade com a literatura:** vignette EPICv2 do DMRcate (Peters TJ, 2025); Peters et al. 2024.

## 7. Escala estatística — M-value para testar, beta para visualizar
**Decisão:** modelagem em M-value; gráficos em beta.
**Por quê:** M-value é aproximadamente homocedástico/normal (válido para modelo linear); beta é intuitivo mas heterocedástico.
**Compatibilidade com a literatura:** Du et al. 2010 (BMC Bioinformatics).

## 8. Agrupamento — juntar MDS_VERY_HIGH em MDS_HIGH
**Decisão:** fundir very_high em high.
**Por quê:** very_high com n≈1, inviável como grupo; ambos são alto risco (IPSS-R). *Julgamento clínico.*

## 9. Modelo — `~ 0 + Group + Sex`, sexo REGISTRADO (não predito)
**Decisão:** ajustar por sexo registrado; `getSex`/`inferSex` só como checagem de troca de amostra.
**Por quê:** sexo afeta metilação autossômica e os grupos são desbalanceados por sexo; a **perda do cromossomo Y (LOY)**, comum em neoplasia mieloide, distorce a predição de sexo por metilação — então o registrado é a verdade confiável.
**Compatibilidade com a literatura:** ajuste por sexo é prática padrão em EWAS; cross-reatividade sexo-associada reforça o cuidado (Chen et al. 2013).

## 10. Contrastes — três pareados, de UM modelo único
**Decisão:** AML vs HIGH, AML vs LOW, HIGH vs LOW, todos do mesmo `fit`/`eBayes`.
**Por quê:** estimativa de variância compartilhada via eBayes → mais poder que três análises separadas (crucial com n pequeno); pareado é mais interpretável que "AML vs MDS agregado".
**Compatibilidade com a literatura:** limma — Ritchie et al. 2015 (NAR); eBayes — Smyth 2004.

## 11. DMPs (limma) e DMRs (DMRcate)
**Compatibilidade com a literatura:** Ritchie et al. 2015; Peters et al. 2015 (Epigenetics & Chromatin) + suporte EPICv2 (Peters et al. 2024).

## 12. Enquadramento exploratório + checagem de confundidores
**Decisão:** n pequeno → estudo exploratório; checar PCA × variáveis clínicas (lote/composição) antes de interpretar; SVA/RUVm como opção.
**Por quê:** a doença é eixo menor de variação; composição celular e lote são confundidores reconhecidos.
**Compatibilidade com a literatura:** em MDS/sAML o não-supervisionado costuma não separar por subtipo (Cabezón et al. 2021, Clinical Epigenetics); SVA — Leek & Storey 2007; RUVm — Maksimovic et al. 2015 (missMethyl).

## Discrepância minfi vs sesame na contagem de sondas (nota)
O menor número de sondas retidas no sesame decorre do uso do **pOOBAH** (detecção baseada
em sinal out-of-band) e do **qualityMask**, mais estritos e específicos que o detection
p-value do minfi, que é **anticonservador** em arrays EPIC (super-detecta; Zhou et al. 2018;
Heiss et al. 2019). É trade-off qualidade × quantidade: menos sondas, porém mais confiáveis.
Os dois pipelines não são comparáveis em contagem de sonda.

## Parâmetros ajustáveis (a reportar explicitamente)
- Limiar de QC de amostra (0,20 de falha / detecção 0,80).
- Limiar de QC de sonda (≤5% de NA).
- Cross-reactive: dropar vs manter+remapear (aqui: dropar).
- Ajuste por % de blastos / SVA / RUVm (a decidir conforme o PCA × clínico).

## Referências
- Zhou W, Triche TJ, Laird PW, Shen H. SeSAMe: reducing artifactual detection of DNA methylation by Infinium BeadChips in genomic deletions. Nucleic Acids Res. 2018;46(20):e123.
- Welsh H, et al. A systematic evaluation of normalization methods and probe replicability using Infinium EPIC methylation data. Clin Epigenetics. 2023.
- Heiss JA, Just AC. Improved filtering of DNA methylation microarray data by detection p values and its impact on downstream analyses. Clin Epigenetics. 2019;11:15.
- Chen YA, et al. Discovery of cross-reactive probes and polymorphic CpGs in the Illumina Infinium HumanMethylation450 microarray. Epigenetics. 2013;8(2):203-9.
- Hop PJ, et al. Cross-reactive probes on Illumina DNA methylation arrays: a large study on ALS. NAR Genom Bioinform. 2020;2(4).
- Peters TJ, et al. Characterisation and reproducibility of the HumanMethylationEPIC v2.0 BeadChip. BMC Genomics. 2024;25:251.
- Peters TJ, et al. De novo identification of differentially methylated regions in the human genome. Epigenetics Chromatin. 2015;8:6.
- Du P, et al. Comparison of Beta-value and M-value methods for quantifying methylation levels. BMC Bioinformatics. 2010;11:587.
- Triche TJ, et al. Low-level processing of Illumina Infinium DNA methylation BeadArrays. Nucleic Acids Res. 2013;41(7):e90.
- Fortin JP, et al. Functional normalization of 450k methylation array data. Genome Biol. 2014;15:503.
- Ritchie ME, et al. limma powers differential expression analyses. Nucleic Acids Res. 2015;43(7):e47.
- Cabezón M, et al. Different methylation signatures in high-risk MDS and secondary AML. Clin Epigenetics. 2021.
- Leek JT, Storey JD. Capturing heterogeneity in gene expression studies by surrogate variable analysis. PLoS Genet. 2007.
- Maksimovic J, et al. Removing unwanted variation in a differential methylation analysis (RUVm). Nucleic Acids Res. 2015 / missMethyl.
