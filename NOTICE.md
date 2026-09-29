# 许可证与第三方来源

本项目原创软件代码及原创说明文档采用根目录 [MIT 许可证](LICENSE)，版权署名为 `songsiyi2006-chem`。

MIT 许可证不替代第三方数据、依赖和引用材料各自的许可条件。数据集原始快照、由其整理的数据库及训练模型使用的数据来源，详见 [data/SOURCES.md](data/SOURCES.md)。该文件保留了已核实的许可信息、引用、来源链接和尚未独立核实的范围；不能将整个数据目录宣称为 MIT 数据。

- AqSolDB 作者仓库的许可证原文保存在 [data/raw/AqSolDB-LICENSE.txt](data/raw/AqSolDB-LICENSE.txt)，其上游数据的权利和引用要求继续适用。
- ESOL、Lipophilicity、FreeSolv 的来源和相应条款见上述数据来源说明。
- QM9 原始数据登记的 CC0 信息及元数据保存在 `data/raw/qm9-figshare-metadata.json`，原论文和其他引用材料保留各自许可。
- FastAPI、Uvicorn、RDKit、NumPy、SciPy、scikit-learn 等依赖通过 `requirements.txt` 安装，未将本机虚拟环境打包进仓库；这些依赖保留其各自许可证。

模型文件和结果可用于复现本项目，但软件许可证、内部测试成绩或成功运行均不构成实验验证、准确性保证或任意用途的适用性保证。
