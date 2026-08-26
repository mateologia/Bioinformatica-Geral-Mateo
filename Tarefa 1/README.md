# Pipeline de Alinhamento - Bioinformática (Tarefa 1)

**Organismo de Referência:** *Pseudomonas fluorescens* SBW25 (AM181176.4) 

**Organismo das Leituras (Reads):** *Pseudomonas aeruginosa* (SRR40334418 - Illumina Paired-end)

## 1. Configuração do Ambiente e Ferramentas
Criação do diretório de trabalho, inicialização do repositório local e instalação das ferramentas de bioinformática (BWA, Samtools e IGV) diretamente pelo terminal do Linux, garantindo que todas as dependências estejam prontas:
```bash
# Criação do diretório de trabalho e navegação
mkdir -p /home/mateo/tarefa1
cd /home/mateo/tarefa1


# Inicialização do repositório Git
git init
git branch -M main
mkdir -p ref reads bam

# Instalação dos softwares necessários (Ubuntu/Linux Mint)
sudo apt update
sudo apt install bwa samtools igv
```

## 2. Genoma de referência
```bash
cd /home/mateo/tarefa1/ref
wget -O genome.fasta "[https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=AM181176.4&rettype=fasta&retmode=text](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=AM181176.4&rettype=fasta&retmode=text)"
```

## 3. Download das leituras
```bash
cd /home/mateo/tarefa1/reads
fastq-dump -X 100000 --split-files SRR40334418
```

## 4. Indexação e Alinhamento Otimizado (BWA-MEM e Samtools)
Os comandos foram adaptados para utilizar pipes (|), processando o alinhamento e a conversão simultaneamente sem gerar arquivos .sam intermediários, economizando espaço em disco:

```bash
cd /home/mateo/tarefa1
# Indexação da referência para BWA e IGV
bwa index ref/genome.fasta
samtools faidx ref/genome.fasta

# Alinhamento, conversão direta para BAM e ordenação
bwa mem -t 4 ref/genome.fasta reads/SRR40334418_1.fastq reads/SRR40334418_2.fastq | samtools view -bS - | samtools sort -o bam/aligned.sorted.bam -

# Indexação do arquivo BAM resultante
samtools index bam/aligned.sorted.bam
```

## 5. Validação do alinhamento
```bash
samtools view -F 4 bam/aligned.sorted.bam | head -n 5
```

## 6. Comando de visualização do IGV
```bash
igv &
```

## 7. Erros
Pipeline Otimizado e Resolução de Problemas (SRA)
Novo Organismo de Referência: Pseudomonas fluorescens SBW25 (AM181176.4)

Novo Organismo das Leituras (Reads): Pseudomonas aeruginosa (SRR40334418 - Illumina Paired-end)

Durante a obtenção das novas leituras, enfrentei desafios:

Erro 404 (Not Found): Códigos iniciados com ERR (do ENA) falharam ao sincronizar com o SRA Toolkit. A solução foi buscar um identificador SRR nativo do NCBI.

Erro 403 (Access Denied): O identificador SRR38877806 retornou erro de embargo, indicando que os dados ainda não haviam sido publicados publicamente pelos autores.

Solução e Mapeamento Cruzado (Cross-mapping): Para contornar os erros e deixar o trabalho mais interessante (biologicamente) no visualizador IGV, utilizei dados públicos antigos e sem embargo de uma espécie do mesmo gênero (Pseudomonas aeruginosa). Isso permitiu visualizar SNPs e divergências evolutivas contra o genoma de referência.

