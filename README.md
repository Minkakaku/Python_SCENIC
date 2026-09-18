# R-Python-scenic-regulon-workflow-template

独立转录调控网络分析点。正式入口均位于仓库根目录：

- scenic_workflow.R：Seurat RDS → gene filtering → GENIE3 → RcisTarget → AUCell → RSS。
- pyscenic_workflow.py：AnnData h5ad → GRNBoost2 → cisTarget motif pruning → AUCell → RSS。

## R SCENIC（旧项目兼容入口）

    Rscript scenic_workflow.R --input annotated.rds --output results/scenic --organism hgnc --database-dir cisTarget --databases db1.feather,db2.feather

R SCENIC 官方仓库已标记 deprecated，新项目优先使用下方 pySCENIC；本入口用于复现既有 R 项目。默认使用 RNA assay 的 counts layer，也可通过 `--layer data` 显式选择非中心化 log-normalized 表达。organism、ranking 数据库和 motif annotations 必须属于同一物种及相容版本。

创建 Conda 环境后，按官方安装方式补装 R SCENIC：

    Rscript -e "remotes::install_github('aertslab/SCENIC')"

## pySCENIC

    python pyscenic_workflow.py --input annotated.h5ad --output results/pyscenic --tf-list allTFs.txt --databases db1.feather db2.feather --motif-annotations motifs.tbl

正式脚本按官方三步 CLI 运行，不使用 os.system，不包含固定服务器路径。输出包括 GRNBoost2 adjacency、motif-pruned regulon、AUCell、RSS、带 regulon 活性的 h5ad 和配对 UMAP 图。

## 解释边界

- 共表达边不是因果调控。
- motif 富集支持潜在 TF 结合，但不能替代 ChIP、扰动或其他直接实验验证。
- RSS 表示细胞类型特异性，不是 CON/EXP 生物学重复检验。
- 数据库版本、物种、基因 ID 和 TF 列表必须在结果说明中记录。

旧平台脚本和示例原样保存在 legacy/original/；source_catalog/ 仅用于来源追溯。

方法依据：SCENIC（PMID: 28991892）、pySCENIC（PMID: 32561888）及官方文档。