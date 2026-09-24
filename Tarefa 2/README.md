 # Projeto de Análise de RNA-Seq: C3KO vs WT (Pseudomonas aeruginosa)

Repositório destinado ao pipeline completo de transcriptómica (Tarefa 2 da disciplina de Bioinformática). O objetivo é avaliar o impacto da ausência do gene C3 na resposta imunitária de macrófagos alveolares murinos perante infeção bacteriana.

## 1. Metadados do Projeto
*   **Base de Dados:** [NCBI GEO - GSE281001](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE281001)
*   **Tipo de Dados:** Bulk RNA-Seq (Expression profiling by high throughput sequencing).
*   **Organismo:** *Mus musculus* (camundongo / genoma de referência GRCm39).
*   **Desenho Experimental:** 2 condições (Murganhos WT vs. Murganhos C3KO, ambos expostos a *Pseudomonas aeruginosa*), com 4 réplicas biológicas por condição, totalizando 8 amostras.
*   **Amostras:** Selecionadas e descarregadas via NCBI SRA Run Selector (tabela `SraRunTable.csv`).                       

## 2. Preparação do Ambiente e Aquisição de Dados
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

## 3. Controle de Qualidade (QC)
mkdir -p qc
fastqc reads/*.fastq.gz -o qc/ -t 4
multiqc qc/ -o qc/

## 4. Limpeza das leituras (Trimming)
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

## 5. Download do genoma e indexação
conda install -c bioconda hisat2 samtools -y
mkdir -p genoma
cd genoma

wget -q --show-progress [http://ftp.ensembl.org/pub/release-112/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz](http://ftp.ensembl.org/pub/release-112/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz)
wget -q --show-progress [http://ftp.ensembl.org/pub/release-112/gtf/mus_musculus/Mus_musculus.GRCm39.112.gtf.gz](http://ftp.ensembl.org/pub/release-112/gtf/mus_musculus/Mus_musculus.GRCm39.112.gtf.gz)

gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
hisat2-build -p 4 Mus_musculus.GRCm39.dna.primary_assembly.fa index_rato
cd ..

## 6. Alinhamento do genoma
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

## 7. Quantificação genética
conda install -c bioconda subread -y
gunzip -f genoma/Mus_musculus.GRCm39.112.gtf.gz
mkdir -p contagens

featureCounts -T 4 -p -a genoma/Mus_musculus.GRCm39.112.gtf -o contagens/contagens_genes.txt alignments/*.bam

## 8. Análise de Expressão Diferencial e Versionamento
conda install -c bioconda -c conda-forge bioconductor-deseq2 bioconductor-enhancedvolcano r-ggplot2 r-pheatmap r-httpgd bioconductor-org.mm.eg.db bioconductor-annotationdbi -y

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

# 1. Ler o ficheiro gerado pelo featureCounts e isolar as contagens
dados <- read.table("contagens/contagens_genes.txt", header = TRUE, row.names = 1, skip = 1)
matriz_contagens <- dados[, 6:ncol(dados)]

# 2. Limpar os nomes das colunas para remover texto extra
colnames(matriz_contagens) <- gsub("alignments\\.", "", colnames(matriz_contagens))
colnames(matriz_contagens) <- gsub("\\.bam$", "", colnames(matriz_contagens))

# 3. Criar o delineamento experimental (Metadata)
library(DESeq2)
grupos <- factor(c("C3KO", "C3KO", "C3KO", "C3KO", "WT", "WT", "WT", "WT"))
grupos <- relevel(grupos, ref = "WT")
design_exp <- data.frame(row.names = colnames(matriz_contagens), condicao = grupos)

# 4. Executar a Expressão Diferencial
dds <- DESeqDataSetFromMatrix(countData = matriz_contagens, colData = design_exp, design = ~ condicao)
dds <- DESeq(dds)
resultados <- results(dds)

# 5. Traduzir identificadores Ensembl para Nomes de Genes Oficiais
library(AnnotationDbi)
library(org.Mm.eg.db)

resultados_df <- as.data.frame(resultados)
resultados_df$Gene_Symbol <- mapIds(org.Mm.eg.db, keys = rownames(resultados_df), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
resultados_df$Gene_Name <- mapIds(org.Mm.eg.db, keys = rownames(resultados_df), column = "GENENAME", keytype = "ENSEMBL", multiVals = "first")
resultados_df <- resultados_df[, c("Gene_Symbol", "Gene_Name", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]

# 6. Transformar dados e gerar PCA Plot
library(EnhancedVolcano)
library(ggplot2)

vsd <- vst(dds, blind = FALSE)
pca_plot <- plotPCA(vsd, intgroup = "condicao") +
  ggtitle("Análise de Componentes Principais (PCA) - WT vs C3KO") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(pca_plot)

# 7. Gerar Volcano Plot
volcano_plot <- EnhancedVolcano(resultados,
    lab = rownames(resultados),
    x = 'log2FoldChange',
    y = 'pvalue',
    title = 'Expressão Diferencial: C3KO vs WT',
    subtitle = 'Infeção por Pseudomonas aeruginosa',
    pCutoff = 0.05,
    FCcutoff = 1.0,
    pointSize = 2.0,
    labSize = 3.0,
    legendPosition = 'right')

print(volcano_plot)

# 8. Guardar ficheiros finais
write.csv(resultados_df, file="contagens/resultados_anotados_finais.csv")
pdf("contagens/volcano_plot.pdf", width=10, height=8)
print(volcano_plot)
dev.off()



