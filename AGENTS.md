# SCENIC 项目约定

- 本仓库是独立分析点，根目录 R/Python 脚本分别接收 RDS/H5AD。
- TF 列表、ranking 数据库、motif annotations、物种和路径通过命令行传入。
- 正式代码不得写固定服务器路径，不使用 os.system。
- 图保存 PNG/PDF，AUC、RSS、映射和中文说明放同一输出目录。
- regulon 是计算推断，不作为直接 TF 结合证据。
- legacy/original/ 和 source_catalog/ 保持原字节。