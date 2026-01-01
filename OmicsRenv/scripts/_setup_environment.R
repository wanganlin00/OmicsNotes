# scripts/setup_environment.R
# Bioinformatics 项目环境设置

cat("=== 设置 Bioinformatics 项目环境 ===\n")


# 加载 renv
if (!require("renv")) {
    install.packages("renv")
}

# 使用 pak 进行快速包安装
if (!require("pak")) {
    install.packages("pak")
}

# 加载元数据

# pak::meta_update()
# pak::meta_list()

cat("---------------------安装核心分析包---------------------\n")


# 安装数据分析包 some dependencies : broom, dbplyr, dtplyr ,readxl,tinytex,
pak::pak("tidyverse")
pak::pak("writexl")
pak::pak("arrow")
pak::pak("DT")
pak::pak("gtsummary")
pak::pak("collapse")      
pak::pak("furrr")


# 会话信息报告
pak::pak("devtools")
pak::pak("styler")

# 生存分析    some dependencies :ggpubr，rstatix，ggsci，ggsignif，corrplot
pak::pak("survminer")

# 注释数据库
pak::pak("msigdb")
pak::pak("AnnotationHub")
pak::pak("biomaRt")
# pak::pak("BSgenome")
pak::pak("org.Hs.eg.db")

# 单细胞分析   some dependencies ：data.table, patchwork, plotly, ggrepel ,reticulate,shiny
pak::pak("Seurat")
pak::pak("satijalab/seurat-wrappers")
remotes::install_github('satijalab/azimuth', ref = 'master')
pak::pak("SingleR")
pak::pak("celldex")

# 
# pak::pak("scater") # pheatmap,scuttle
# pak::pak("scran") # edgeR limma


# 免疫浸润分析
pak::pak("TCGAbiolinks")
pak::pak("GEOquery")
pak::pkg_install("IOBR/IOBR")
pak::pkg_install("omnideconv/immunedeconv")


# 下游分析包

pak::pak("BiocManager")
pak::pak("rtracklayer")
pak::pak("SingleCellExperiment") # deps ： SummarizedExperiment
# 差异表达分析
BiocManager::install("DESeq2")
# 富集分析
pak::pak("clusterProfiler") # deps: DOSE, enrichplot, fgsea,GO.db, qvalue


cat("---------------------创建环境快照---------------------\n")
renv::snapshot()

cat("\n=========== 环境设置完成, 关键包安装成功)===========\n")
