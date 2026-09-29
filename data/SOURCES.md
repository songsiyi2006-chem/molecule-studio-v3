# 本地理化性质数据来源与收录说明

核验日期：2026-09-28。本目录包含四个指定公开基准的数据快照，覆盖水溶解度、pH 7.4 的辛醇/水分配系数和水合自由能。它不是 PubChem 或 ChEMBL 的完整镜像。`raw/` 保留取得的原始 CSV；`sources.json` 记录来源、标签列、单位、字节数、SHA-256、原始取得时间及本次复用/校验信息。

## 实际收录数量

| 数据集 | 原始行数 | 本项目收录观测 | 排除 | 目标与单位 |
|---|---:|---:|---:|---|
| ESOL | 1128 | 1128 | 0 | `logS`，log10(mol/L) |
| AqSolDB | 9982 | 8718 | 1264 | `logS`，log10(mol/L) |
| Lipophilicity | 4200 | 4197 | 3 | `logD74`，log10(D)，pH 7.4 |
| FreeSolv / SAMPL | 642 | 642 | 0 | `hydration_free_energy`，kcal/mol |
| 合计 | 15952 | 14685 | 1267 | 三种性质分别保存 |

SQLite 中包含 **12843 个规范化结构、14685 条来源观测**。去除立体信息用于分组的 `identity_smiles` 共 12410 个；此标识不表示已经统一了互变异构体、质子化状态或全部实验条件。目标观测数分别为 logS 9846、logD74 4197、水合自由能 642。

这些“观测”是源表中的标签记录，不意味着 14685 次相互独立、由本项目执行的实验。AqSolDB 已整合了多个文献来源，其中包含 Delaney 数据，与独立提供的 ESOL 表存在重叠。数据库保留分开的来源记录；训练前仍须按结构身份去重和隔离，不能把不同来源名当成外部验证。

## 1. AqSolDB：作者仓库完整原始表

- [作者仓库](https://github.com/mcsorkun/AqSolDB)
- [本次使用的原始 CSV](https://raw.githubusercontent.com/mcsorkun/AqSolDB/master/results/data_curated.csv)
- [原始论文：AqSolDB, a curated reference set of aqueous solubility and 2D descriptors for a diverse set of compounds](https://doi.org/10.1038/s41597-019-0151-1)
- [作者登记的数据存储记录](https://doi.org/10.7910/DVN/OVHAW8)

真实文件名是 `results/data_curated.csv`，本地保存为 `raw/aqsoldb.csv`，原始行数 9982。使用 `Solubility` 列作为 logS；源表内已有的计算描述符不作为新的实验标签。每行保留原 `ID`，并在 `quality` 中记录 `Group`、`SD` 和 `Occurrences`。

作者的分组是来源重复与标签选择过程的质量标记，不是经过统一测量误差校准的置信度。G2/G3 的候选实验值选择曾参考 ALOGPS 预测值；本项目保留的仍是作者选择的实验来源 `Solubility` 标签，但不将其描述为完全不受模型辅助整理影响的数据。不能把 G1 简单解释为最高可信度。

作者 GitHub 仓库附有 MIT 许可证，已保存为 `raw/AqSolDB-LICENSE.txt`。原始子数据集的权利与引用要求仍需保留；本次未成功独立取得 Dataverse 的单独数据许可字段，因此不将仓库代码许可证泛化为对所有上游数据权利的确认。

## 2. ESOL：MoleculeNet 分发快照

- [DeepChem/MoleculeNet 官方说明](https://deepchem.readthedocs.io/en/latest/api_reference/moleculenet.html)
- [分发 CSV](https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/delaney-processed.csv)
- [Delaney 原始论文](https://doi.org/10.1021/ci034243x)

使用本机此前取得且经 SHA-256 核验一致的 1128 行快照。目标只读取 `measured log solubility in mols per litre`；`ESOL predicted log solubility in mols per litre` 是计算预测列，不作为训练标签。CSV 没有稳定的唯一行 ID，因此本项目以包含表头计数的 `csv-row-N` 定位记录，`Compound ID` 保留为来源名称。

此分发表没有附带经本次独立核实的单独数据许可证。保留 MoleculeNet 与 Delaney 引用及原始数据条款，不宣称所有源数据均采用本项目软件许可证。

## 3. Lipophilicity：pH 7.4 的 logD

- [官方任务说明](https://deepchem.readthedocs.io/en/latest/api_reference/moleculenet.html#lipo-datasets)
- [分发 CSV](https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/Lipophilicity.csv)
- [ChEMBL AZ 数据登记](https://doi.org/10.6019/chembl3301361)

复用本机经哈希核验的 4200 行快照。使用 `exp`，结构标识为 `CMPD_CHEMBLID`。该任务是实验辛醇/水分配系数 **logD at pH 7.4**，不等同于 RDKit 计算的 Crippen LogP。模型不能由这一单一 pH 任务直接推断任意 pH 下的 logD。

记录 ChEMBL 与 AZ 来源并保留上游条款。本项目不会给源表重新指定许可证，也未逐条重新核验其测定方案和实验条件。

## 4. FreeSolv / SAMPL：实验水合自由能

- [FreeSolv 作者仓库与许可说明](https://github.com/MobleyLab/FreeSolv)
- [本次使用的 MoleculeNet 分发 CSV](https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/SAMPL.csv)
- [作者原始数据库字段说明](https://raw.githubusercontent.com/MobleyLab/FreeSolv/master/database.txt)
- [FreeSolv 论文](https://doi.org/10.1007/s10822-014-9747-x)

复用本机经哈希核验的 642 行 MoleculeNet 快照，只读取 `expt`，单位 kcal/mol。`calc` 是计算水合自由能，未作为实验标签。此 CSV 没有原始数据库中所有逐条不确定度及参考文字，本项目不会补造这些信息；它也不声称是上游最新版本的完整镜像。

作者仓库对可由作者授权的数据指定 CC-BY-4.0，同时明确保留原始来源可能另有条款的限制；代码许可证另行规定。本项目保留该界限与作者引用。

## 导入规则与排除

导入只接受共享 RDKit 解析器支持的单片段闭壳层结构，元素范围为 H、B、C、N、O、F、P、S、Cl、Br、I，解析后最多 200 个原子。不会从盐中擅自剥离片段后继承原标签，也不会为超范围元素创建替代结构。

AqSolDB 排除：1097 条多片段、158 条不支持元素、6 条自由基、1 条超过原子数上限、2 条当前 RDKit 无法解析。Lipophilicity 排除：2 条不支持元素、1 条多片段。全部记录保存在 `../reports/data_audit.json`，含源 ID、CSV 行号、原始 SMILES 和原因。

输入在数据库阶段不做平均，也不删除跨来源重叠记录。`canonical_smiles` 保留可表示的立体信息，`identity_smiles` 用非异构规范 SMILES 进行保守分组；两者均通过共享解析流程去除原子映射号。骨架使用 Bemis–Murcko，空骨架字符串表示无环分子。结构标准化不等于实验条件标准化。

## 重建与核验

在项目环境执行：

```powershell
python scripts/prepare_data.py --offline
python -m unittest discover -s tests -p test_database.py -v
```

已有 `raw/` 时可完全离线重建。缺少原始文件时，不带 `--offline` 的准备脚本先检查本机已验证的化学数据缓存，再从记录的公开地址取得；只有与固定 SHA-256 一致的文件才会进入整理。脚本生成 `chemistry.sqlite`、`curated_observations.csv`、`sources.json` 与逐条数据审计。

本次 13 项数据库测试实际通过，检查全部收录标签与源 CSV 的指定实验列逐值一致、四份原始文件的哈希与行数、读写边界、分页和性质筛选、来源查询、已知中文别名精确匹配、结果排序、SQL 注入与通配符处理、表完整性及排除计数。它们验证数据处理与查询的一致性，不替代原始实验质量审查。


## V3 QM9 quantum reference catalog

- Dataset: Ramakrishnan R., Dral P. O., Rupp M., von Lilienfeld O. A. (2014), *Quantum chemistry structures and properties of 134 kilo molecules*. DOI: https://doi.org/10.1038/sdata.2014.22
- Original data: https://doi.org/10.6084/m9.figshare.c.978904 ; original data item https://springernature.figshare.com/articles/dataset/Data_for_133885_GDB-9_molecules/1057646 . Figshare API metadata saved under raw/qm9-figshare-metadata.json reports **CC0** (data license; article license is separate).
- Reused DeepChem CSV snapshot: https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/qm9.csv ; SHA-256 3e668f8c34e4bc392a90d417a50a5eed3b64b842a817a633024bdc054c68ccb4. Stored losslessly as raw/qm9.csv.gz.
- Official geometry-consistency exclusion list: https://ndownloader.figshare.com/files/3195404 . All 3054 listed source IDs are excluded. 130831 accepted reference rows represent 130744 distinct canonical structures; duplicate structure IDs remain separately traceable.
- The original README is saved as raw/qm9-readme.txt. Its eleven special convergence cautions are flagged on matching returned records (unless already excluded by the geometry list). Passing parsing/consistency checks is not an independent DFT validation.
- Returned properties: dipole moment in Debye; isotropic polarizability in Bohr³; HOMO, LUMO and gap converted from Hartree to eV using 27.211386245988. Method: B3LYP/6-31G(2df,p). Raw precision is retained; gap can differ slightly from subtracting rounded orbital energies.
- This catalog is a lookup of prior calculations. It supplies neither experimental measurements nor newly computed DFT. It is **not used to fit any of the three physicochemical prediction models**. It contains small H/C/N/O/F molecules and cannot answer arbitrary structures.
- Database construction: `python scripts/build_catalog.py`; row counts, exclusions and hashes: reports/catalog_audit.json.

## V3 runtime three-dimensional calculations

ETKDGv3 embeds conformers; MMFF94 supplies local molecular-mechanics optimization. These force-field energies are not DFT energies, aqueous free energies, or cross-molecule ranking scores. Only relative conformer energies within the same input molecule are shown.
RDKit documentation: https://www.rdkit.org/docs/RDKit_Book.html and https://www.rdkit.org/docs/source/rdkit.Chem.rdForceFieldHelpers.html .

实验参考库相似检索使用不含手性标记的 Morgan radius 2 / 2048 位 Tanimoto；因此相似度 1 不证明是同一立体异构体。精确结构检索使用规范异构 SMILES。模型适用域相似度使用独立训练集合与含手性指纹，不从整个参考数据库推定训练成员。
