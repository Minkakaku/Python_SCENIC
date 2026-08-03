# R-Python-scenic-regulon-workflow-template

SCENIC、pySCENIC 和 regulon 活性分析代码素材模板。

## 主要代码

- `platform/scenic.py`：AnnData 预处理及 pySCENIC 输入准备。
- `platform/pyscenic.sh`、`platform/pyscenic_mouse.sh`：人和小鼠 pySCENIC 命令流程。
- `platform/scenic.r`：AUCell/SCENIC 结果读取与可视化。
- `platform/report_scenic.rmd`：结果报告模板。
- `examples/`：原项目运行示例。
- `platform_wrappers/`：原平台启动代码。

## 使用前注意

运行前需要准备匹配物种的 motif、ranking 数据库及 loom 文件，并重新配置数据库路径和计算资源。平台封装与项目脚本按原样保留，不作为新项目首选入口。

原始副本保存在 `../R-Python-bioinformatics-code-archive/legacy/original/07_scenic_regulon/`。

## 其他人员来源代码

`source_catalog/` 收录 62 份 SCENIC/pySCENIC/regulon 语义副本；`SOURCE_CATALOG.tsv` 可按方法、语言、owner 和场景筛选。来源脚本可能依赖不同版本的 motif/ranking 数据库，运行前必须核对物种和数据库版本。
