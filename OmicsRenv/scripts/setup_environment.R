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

# 会话信息报告
pak::pak("devtools")

# 生存分析    some dependencies :ggpubr，rstatix，ggsci，ggsignif，corrplot
pak::pak("survminer")


# 安装 Seurat 包    some dependencies ：data.table, patchwork, plotly, ggrepel ,reticulate,shiny    
pak::pak("Seurat")      

# 单细胞分析
pak::pak("scater")              # pheatmap,scuttle
pak::pak("scran")               # edgeR limma
# 下游分析包

pak::pak("BiocManager")
pak::pak("rtracklayer")
pak::pak("SingleCellExperiment")      # deps ： SummarizedExperiment
    # 注释数据库
pak::pak("AnnotationHub")
pak::pak("biomaRt")
pak::pak("BSgenome")
pak::pak("org.Hs.eg.db")

    # 差异表达分析           
BiocManager::install("DESeq2")     
    # 富集分析
pak::pak("clusterProfiler")           # deps: DOSE, enrichplot, fgsea,GO.db, qvalue

# HPO.db
# MPO.db
    

cat("---------------------验证安装---------------------\n")

cat("---------------------创建环境快照---------------------\n")
renv::snapshot()

cat("\n=========== 环境设置完成, 关键包安装成功)===========\n")