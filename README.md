# 线束 3D 拓扑解析工具（本地运行说明）

从 CATIA 导出的 STEP(AP242) 中提取线束拓扑：分支中心线 → 拓扑节点 → 关联连接器。
只读 STEP 里的 3D 几何（装配树 + 实体 + 曲线），**不需要 CATIA、不需要线束模块授权**。

## 核心能力

- **分支判定不依赖命名**：纯结构化判定——实体有成组等面积平面端盖、长径比≥4、侧壁面很少
  → 判为管状分支；连接器（粗短）、卡扣（碎面）被自动排除。
- **中心线自动获取**：分支有 Flexible Curve 线框就用线框；**没有线框时自动从管状实体反推**
  （侧壁面逐站位轮廓圆拟合求轴线 → 半壳去重 → 拼接 → 端点对齐端盖中心）。
  反推精度已用 M1E-DRD 验证：与真实中心线对比，平均偏差 0.10mm、最大 0.34mm。

## 一、环境准备（只需做一次）

1. 安装 Anaconda（已有可跳过）。
2. 新建独立环境（不要装进旧的 Python 3.8 里）：
```bash
conda create -n harness python=3.12 -y
conda activate harness
pip install "cadquery-ocp==8.0.1.0.0" numpy matplotlib openpyxl
```
`cadquery-ocp` 是 OpenCASCADE 几何内核的 Python 绑定（约两三百 MB），STEP 的装配树读取和几何计算全靠它。

## 二、使用方法（图形界面，推荐）

1. 把本工具包解压到任意文件夹（比如 `D:\harness_tool\`）。
2. 双击 `启动线束拓扑解析.bat`（它会自动找到 harness 环境并打开界面）。
3. 点"选择 STEP 文件并分析"，选中你的 `.stp` / `.step` 文件。
4. 等待分析完成，输出文件自动生成在 **STEP 文件所在文件夹**：
   - `<模型名>_线束3D拓扑分析.xlsx` —— 分支清单 / 拓扑节点 / 连续走线 / 连接器卡扣清单 / 说明
   - `<模型名>_拓扑3D图.png` —— 3D 拓扑图（分支=彩色线，分支点=红菱形，终端=黑点，灰盒=连接器/卡扣）
   - `<模型名>_topology.json` —— 中间数据（可删）

界面上的"端点聚类容差"默认 3.0mm：分支端点在此距离内会被聚成同一个拓扑节点。

## 三、命令行方式（可选）

```bash
conda activate harness
cd D:\harness_tool
python harness_topology.py 你的模型.stp 3.0
python viz_topology.py
python report_xlsx.py
```

## 四、适用条件与局限

- STEP 需为 AP242（AP214 也可），保留装配结构；
- 分支产品最好带有 Flexible Curve 中心线（如 CATIA Electrical Harness 导出的 Multi-branchable）；
  若只有扫掠实体、没有中心线，本工具会提示无法识别分支，需要另做中心线反推；
- 连接器关联基于空间邻近启发式，结果需对照 3D 图复核；
- STEP 为哑几何，不含导线代号/回路等电气属性；电气拓扑需用 CATIA EHA 导出的 from-to/XML 另行对齐。

## 五、文件清单

| 文件 | 说明 |
|---|---|
| `启动线束拓扑解析.bat` | 双击启动图形界面 |
| `harness_gui.pyw` | 图形界面主程序 |
| `harness_topology.py` | 拓扑分析核心（可独立命令行运行） |
| `viz_topology.py` | 3D 图生成 |
| `report_xlsx.py` | Excel 报告生成 |
| `README.md` | 本说明 |
