# Projeto de Análise de RNA-Seq: C3KO vs WT (Pseudomonas aeruginosa)

Repositório destinado ao pipeline completo de transcriptómica (Tarefa 2 da disciplina de Bioinformática). O objetivo é avaliar o impacto da ausência do gene C3 na resposta imunitária de macrófagos alveolares murinos perante infeção bacteriana.

## 1. Metadados do Projeto
*   **Base de Dados:** [NCBI GEO - GSE281001](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE281001)
*   **Tipo de Dados:** Bulk RNA-Seq (Expression profiling by high throughput sequencing).
*   **Organismo:** *Mus musculus* (camundongo / genoma de referência GRCm39).
*   **Desenho Experimental:** 2 condições (Murganhos WT vs. Murganhos C3KO, ambos expostos a *Pseudomonas aeruginosa*), com 4 réplicas biológicas por condição, totalizando 8 amostras.
*   **Amostras:** Selecionadas e descarregadas via NCBI SRA Run Selector (tabela `SraRunTable.csv`).

---

## 2. Preparação do Ambiente e Aquisição de Dados
Ativação do ambiente Conda e extração dos identificadores para download das leituras brutas:

```bash
conda activate bioinfo
conda install -c conda-forge -c bioconda fastqc multiqc -y
cd "/home/mateo/Bioinformatica-Geral-Mateo/Tarefa 2"   

# Extrair lista de amostras e baixar leituras
awk -F',' 'NR>1 {print $1}' SraRunTable.csv | head -n 8 > SRR_list.txt
mkdir -p reads
cd reads                                  

fastq-dump --split-files --gzip SRR31221361
fastq-dump --split-files --gzip SRR31221363
fastq-dump --split-files --gzip SRR31221365
fastq-dump --split-files --gzip SRR31221367
fastq-dump --split-files --gzip SRR31221369
fastq-dump --split-files --gzip SRR31221371
fastq-dump --split-files --gzip SRR31221373
fastq-dump --split-files --gzip SRR31221375
cd ..
```

---

## 3. Controlo de Qualidade (QC)
Avaliação da qualidade bruta das leituras utilizando FastQC e agregação com MultiQC:

```bash
mkdir -p qc
fastqc reads/*.fastq.gz -o qc/ -t 4
multiqc qc/ -o qc/
```

---

## 4. Limpeza das Leituras (Trimming)
Remoção de adaptadores e bases de baixa qualidade utilizando o `fastp` (executado individualmente para cada amostra emparelhada):

```bash
conda install -c bioconda fastp -y
mkdir -p reads_trimmed

fastp -i reads/SRR31221361_1.fastq.gz -I reads/SRR31221361_2.fastq.gz -o reads_trimmed/SRR31221361_1_trimmed.fastq.gz -O reads_trimmed/SRR31221361_2_trimmed.fastq.gz -h reads_trimmed/SRR31221361_fastp.html -w 4
fastp -i reads/SRR31221363_1.fastq.gz -I reads/SRR31221363_2.fastq.gz -o reads_trimmed/SRR31221363_1_trimmed.fastq.gz -O reads_trimmed/SRR31221363_2_trimmed.fastq.gz -h reads_trimmed/SRR31221363_fastp.html -w 4
fastp -i reads/SRR31221365_1.fastq.gz -I reads/SRR31221365_2.fastq.gz -o reads_trimmed/SRR31221365_1_trimmed.fastq.gz -O reads_trimmed/SRR31221365_2_trimmed.fastq.gz -h reads_trimmed/SRR31221365_fastp.html -w 4
fastp -i reads/SRR31221367_1.fastq.gz -I reads/SRR31221367_2.fastq.gz -o reads_trimmed/SRR31221367_1_trimmed.fastq.gz -O reads_trimmed/SRR31221367_2_trimmed.fastq.gz -h reads_trimmed/SRR31221367_fastp.html -w 4
fastp -i reads/SRR31221369_1.fastq.gz -I reads/SRR31221369_2.fastq.gz -o reads_trimmed/SRR31221369_1_trimmed.fastq.gz -O reads_trimmed/SRR31221369_2_trimmed.fastq.gz -h reads_trimmed/SRR31221369_fastp.html -w 4
fastp -i reads/SRR31221371_1.fastq.gz -I reads/SRR31221371_2.fastq.gz -o reads_trimmed/SRR31221371_1_trimmed.fastq.gz -O reads_trimmed/SRR31221371_2_trimmed.fastq.gz -h reads_trimmed/SRR31221371_fastp.html -w 4
fastp -i reads/SRR31221373_1.fastq.gz -I reads/SRR31221373_2.fastq.gz -o reads_trimmed/SRR31221373_1_trimmed.fastq.gz -O reads_trimmed/SRR31221373_2_trimmed.fastq.gz -h reads_trimmed/SRR31221373_fastp.html -w 4
fastp -i reads/SRR31221375_1.fastq.gz -I reads/SRR31221375_2.fastq.gz -o reads_trimmed/SRR31221375_1_trimmed.fastq.gz -O reads_trimmed/SRR31221375_2_trimmed.fastq.gz -h reads_trimmed/SRR31221375_fastp.html -w 4
```

---

## 5. Download do Genoma e Indexação
Download do genoma de referência (FASTA) e anotação (GTF) do Ensembl (Release 112) e construção do índice do `HISAT2`:

```bash
conda install -c bioconda hisat2 samtools -y
mkdir -p genoma
cd genoma

wget -q --show-progress http://ftp.ensembl.org/pub/release-112/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
wget -q --show-progress http://ftp.ensembl.org/pub/release-112/gtf/mus_musculus/Mus_musculus.GRCm39.112.gtf.gz

gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
hisat2-build -p 4 Mus_musculus.GRCm39.dna.primary_assembly.fa index_rato
cd ..
```

---

## 6. Alinhamento ao Genoma
Alinhamento das leituras ao genoma e conversão para BAM ordenado (taxa de alinhamento obtida: ~93-96%):

```bash
mkdir -p alignments

cat << 'EOF' > rodar_hisat2.sh
hisat2 -p 4 -x genoma/index_rato -1 reads_trimmed/SRR31221361_1_trimmed.fastq.gz -2 reads_trimmed/SRR31221361_2_trimmed.fastq.gz | samtools view -bS - | samtools sort -@ 4 -o alignments/SRR31221361.bam
hisat2 -p 4 -x genoma/index_rato -1 reads_trimmed/SRR31221363_1_trimmed.fastq.gz -2 reads_trimmed/SRR31221363_2_trimmed.fastq.gz | samtools view -bS - | samtools sort -@ 4 -o alignments/SRR31221363.bam
hisat2 -p 4 -x genoma/index_rato -1 reads_trimmed/SRR31221365_1_trimmed.fastq.gz -2 reads_trimmed/SRR31221365_2_trimmed.fastq.gz | samtools view -bS - | samtools sort -@ 4 -o alignments/SRR31221365.bam
hisat2 -p 4 -x genoma/index_rato -1 reads_trimmed/SRR31221367_1_trimmed.fastq.gz -2 reads_trimmed/SRR31221367_2_trimmed.fastq.gz | samtools view -bS - | samtools sort -@ 4 -o alignments/SRR31221367.bam
hisat2 -p 4 -x genoma/index_rato -1 reads_trimmed/SRR31221369_1_trimmed.fastq.gz -2 reads_trimmed/SRR31221369_2_trimmed.fastq.gz | samtools view -bS - | samtools sort -@ 4 -o alignments/SRR31221369.bam
hisat2 -p 4 -x genoma/index_rato -1 reads_trimmed/SRR31221371_1_trimmed.fastq.gz -2 reads_trimmed/SRR31221371_2_trimmed.fastq.gz | samtools view -bS - | samtools sort -@ 4 -o alignments/SRR31221371.bam
hisat2 -p 4 -x genoma/index_rato -1 reads_trimmed/SRR31221373_1_trimmed.fastq.gz -2 reads_trimmed/SRR31221373_2_trimmed.fastq.gz | samtools view -bS - | samtools sort -@ 4 -o alignments/SRR31221373.bam
hisat2 -p 4 -x genoma/index_rato -1 reads_trimmed/SRR31221375_1_trimmed.fastq.gz -2 reads_trimmed/SRR31221375_2_trimmed.fastq.gz | samtools view -bS - | samtools sort -@ 4 -o alignments/SRR31221375.bam
EOF

bash rodar_hisat2.sh
```

---

## 7. Quantificação Genética
Contagem da expressão de genes utilizando o `featureCounts` (pacote `subread`):

```bash
conda install -c bioconda subread -y
gunzip -f genoma/Mus_musculus.GRCm39.112.gtf.gz
mkdir -p contagens

featureCounts -T 4 -p -a genoma/Mus_musculus.GRCm39.112.gtf -o contagens/contagens_genes.txt alignments/*.bam
```

---

## 8. Análise de Expressão Diferencial e Vias (GSEA)
Instalação das dependências para análise estatística no R (DESeq2, anotações, visualizações avançadas e ClusterProfiler) e configuração do GitHub para ignorar ficheiros pesados:

```bash
conda install -c bioconda -c conda-forge bioconductor-deseq2 bioconductor-enhancedvolcano r-ggplot2 r-pheatmap r-httpgd bioconductor-org.mm.eg.db bioconductor-annotationdbi bioconductor-clusterprofiler bioconductor-enrichplot r-gdtools -y

# Criação do escudo para evitar upload de ficheiros pesados
cat << 'EOF' > .gitignore
reads/
reads_trimmed/
alignments/
genoma/
*.bam
*.fastq
*.fastq.gz
*.gtf
*.fasta
*.fna
EOF
```

A modelagem matemática, anotação de genes e geração dos gráficos finais (PCA, Volcano, Heatmap e Enriquecimento de Vias) foram realizadas utilizando o script `analise_expressao.R`:

```R
# 1. Preparação da Matriz de Contagens
dados <- read.table("contagens/contagens_genes.txt", header = TRUE, row.names = 1, skip = 1)
matriz_contagens <- dados[, 6:ncol(dados)]
colnames(matriz_contagens) <- gsub("alignments\\.|\\.bam$", "", colnames(matriz_contagens))

# 2. Delineamento Experimental (DESeq2)
library(DESeq2)
grupos <- factor(c("C3KO", "C3KO", "C3KO", "C3KO", "WT", "WT", "WT", "WT"))
grupos <- relevel(grupos, ref = "WT")
design_exp <- data.frame(row.names = colnames(matriz_contagens), condicao = grupos)

dds <- DESeqDataSetFromMatrix(countData = matriz_contagens, colData = design_exp, design = ~ condicao)
dds <- DESeq(dds)
resultados <- results(dds)
vsd <- vst(dds, blind = FALSE)

# 3. Anotação dos Genes (Ensembl para Gene Symbol)
library(AnnotationDbi)
library(org.Mm.eg.db)
resultados_df <- as.data.frame(resultados)
resultados_df$Gene_Symbol <- mapIds(org.Mm.eg.db, keys = rownames(resultados_df), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
resultados_df$Gene_Name <- mapIds(org.Mm.eg.db, keys = rownames(resultados_df), column = "GENENAME", keytype = "ENSEMBL", multiVals = "first")
write.csv(resultados_df, file="contagens/resultados_anotados_finais.csv")

# 4. Gráficos Base: PCA e Volcano Plot Limpo
library(ggplot2)
library(EnhancedVolcano)
pca_plot <- plotPCA(vsd, intgroup = "condicao") + ggtitle("PCA - WT vs C3KO") + theme_minimal()

nomes_grafico <- ifelse(is.na(resultados_df[["Gene_Symbol"]]), rownames(resultados_df), resultados_df[["Gene_Symbol"]])
top_20_genes <- nomes_grafico[order(resultados_df[["padj"]])][1:20]

volcano_plot_limpo <- EnhancedVolcano(resultados_df, lab = nomes_grafico, x = 'log2FoldChange', y = 'pvalue',
    title = 'Expressão Diferencial: C3KO vs WT', subtitle = 'Infeção por P. aeruginosa',
    pCutoff = 0.05, FCcutoff = 1.0, pointSize = 1.5, labSize = 4.0, selectLab = top_20_genes,
    drawConnectors = TRUE, widthConnectors = 0.5, legendPosition = 'right')

pdf("contagens/volcano_plot_apresentacao.pdf", width=10, height=8)
print(volcano_plot_limpo)
dev.off()

# 5. Top 30 Genes: Heatmap
library(pheatmap)
genes_top30 <- head(rownames(resultados_df[order(resultados_df$padj), ]), 30)
matriz_top30 <- assay(vsd)[genes_top30, ]
rownames(matriz_top30) <- ifelse(is.na(resultados_df[genes_top30, "Gene_Symbol"]), rownames(matriz_top30), resultados_df[genes_top30, "Gene_Symbol"])

pdf("contagens/heatmap_top30.pdf", width=8, height=10)
pheatmap(matriz_top30, scale = "row", annotation_col = as.data.frame(colData(dds)[,"condicao", drop=FALSE]), main = "Top 30 Genes (C3KO vs WT)")
dev.off()

# 6. Gene Set Enrichment Analysis (GSEA)
library(clusterProfiler)
library(enrichplot)

res_gsea <- resultados_df[!is.na(resultados_df$log2FoldChange), ]
lista_ranqueada <- sort(setNames(res_gsea$log2FoldChange, rownames(res_gsea)), decreasing = TRUE)

gsea_resultado <- gseGO(geneList = lista_ranqueada, OrgDb = org.Mm.eg.db, keyType = "ENSEMBL", ont = "BP", pvalueCutoff = 0.05, eps = 0) 

pdf("contagens/gsea_dotplot.pdf", width=12, height=8)
dotplot(gsea_resultado, showCategory=10, split=".sign") + facet_grid(.~.sign) + ggtitle("GSEA: Vias Biológicas em C3KO")
dev.off()

pdf("contagens/gsea_top3_vias.pdf", width=10, height=8)
gseaplot2(gsea_resultado, geneSetID = 1:3, pvalue_table = TRUE)
dev.off()