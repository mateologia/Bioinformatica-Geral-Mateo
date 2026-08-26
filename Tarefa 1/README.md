# Pipeline de Alinhamento - Bioinformática (Tarefa 1)

**Organismo de Referência:** *Pseudomonas fluorescens* SBW25 (AM181176.4)  
**Organismo das Leituras (Reads):** *Pseudomonas aeruginosa* (SRR40334418 - Illumina Paired-end)

## 1. Resolução de Problemas com o Banco de Dados (SRA)
Durante o download das leituras (reads), enfrentei desafios comuns de bioinformática:
* **Erro 404 (Not Found):** Códigos iniciados com `ERR` (do ENA) falharam ao sincronizar com o SRA Toolkit. A solução foi buscar um identificador `SRR` nativo do NCBI.
* **Erro 403 (Access Denied):** O identificador `SRR38877806` retornou erro de embargo, indicando que os dados ainda não haviam sido publicados publicamente pelos autores.
* **Solução e Mapeamento Cruzado (Cross-mapping):** Para contornar os erros e enriquecer a discussão biológica no visualizador IGV, utilizei dados públicos antigos e sem embargo de uma espécie do mesmo gênero (*Pseudomonas aeruginosa*). Isso permitiu visualizar SNPs e divergências evolutivas contra o genoma de referência.

## 2. Preparação do Ambiente e Download do Genoma
```bash
mkdir -p ref reads bam
cd ref
wget -O genome.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=AM181176.4&rettype=fasta&retmode=text"
```

## 3. Download das Leituras (reads)
```bash
cd ../reads
fastq-dump -X 100000 --split-files SRR40334418
```

## 4. Indexação e Alinhamento (BWA-MEM e Samtools)
```bash
cd ..
# Indexação da referência para BWA e IGV
bwa index ref/genome.fasta
samtools faidx ref/genome.fasta

# Alinhamento, conversão para BAM e ordenação
bwa mem -t 4 ref/genome.fasta reads/SRR40334418_1.fastq reads/SRR40334418_2.fastq | samtools view -bS - | samtools sort -o bam/aligned.sorted.bam -

# Indexação do arquivo BAM resultante
samtools index bam/aligned.sorted.bam
```
## 5. Validação do Alinhamento
```bash
samtools view -F 4 bam/aligned.sorted.bam | head -n 5
```
