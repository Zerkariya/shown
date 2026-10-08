[English](README.md) | **简体中文**

# shown

把**列数不固定**的患者 / 队列表格（CSV、TSV、Excel、VCF）交给 AI，生成一个**本地双击即可打开、完全离线**的 HTML 报告。报告的图按 [Delly](https://github.com/dellytools/delly)（结构变异、读深拷贝数）和 [ASCAT](https://github.com/VanLoo-lab/ascat)（等位基因特异拷贝数、纯度、倍性）的输出样式绘制。

这是一个 Claude Code / Codex **插件**，里面只有一个 skill：`csv-report`。分工如下：

- **AI（Claude 或 Codex）**：读表头和数据，判断每一列的含义（患者编号、染色体、位置、nMajor/nMinor、SV 类型、纯度、临床字段……），写出 `mapping.json`。
- **固定脚本**：按映射读数据、计算 ASCAT 指标、画图、输出单个 HTML。只依赖 `python3` 标准库，不用装任何包。

报告默认是英文界面，右上角的 **中文** 按钮可以一键把整页切换成中文（AI 写的标题、字段名和解读也会一起切换），另有浅色 / 深色主题切换。

示例报告见 [`examples/report.html`](examples/report.html)，下载后双击打开。数据是随机生成的，不是真实患者。

## 安装

同事都需要能访问本仓库；仓库是私有的话，需要先配好 git 凭据。

### Claude Code

```text
/plugin marketplace add zerkariya/shown
/plugin install shown@shown
```

也可以在终端里执行：

```bash
claude plugin marketplace add zerkariya/shown
claude plugin install shown@shown
```

装好后新开一个会话（或执行 `/reload-plugins`）就能用。更新插件：`claude plugin marketplace update shown`。

### Codex

```bash
codex plugin marketplace add zerkariya/shown
codex plugin add shown@shown
```

Codex 直接读取本仓库的 `.claude-plugin/marketplace.json` 和 `plugins/shown/.claude-plugin/plugin.json`，所以两边共用同一套文件。

### 不想用插件系统时（两边通用）

直接把 skill 文件夹复制到个人 skills 目录：

```bash
git clone https://github.com/zerkariya/shown.git
# Claude Code
cp -r shown/plugins/shown/skills/csv-report ~/.claude/skills/
# Codex
mkdir -p ~/.agents/skills && cp -r shown/plugins/shown/skills/csv-report ~/.agents/skills/
```

## 使用

在 Claude Code 或 Codex 里直接说需求即可，例如：

> 把 `data/患者信息.xlsx` 和 `data/sv.csv` 做成可视化报告

Claude Code 里也可以显式调用：`/shown:csv-report data/患者信息.xlsx`。Codex 里用 `$csv-report` 提及这个 skill。

AI 会按下面的顺序做：

1. 运行 `profile_table.py` 看每一列
2. 写出 `report/mapping.json`
3. 运行 `build_report.py --check`，把报错改到通过
4. **生成报告前先问你确认**：发一条简短的消息，说明它对每个文件的理解、要隐藏哪些列，再列出几个带建议答案的问题。你回复修改意见，或者回"OK"全部按默认；如果一开始就说"不用问"，会跳过这一步
5. 生成 `report/report.html`
6. 告诉你每张表是怎么映射的、哪些是你确认过的、哪些用了默认值

`mapping.json` 会保留下来，数据更新后可以直接复用。

### 提示词怎么写

一句话就能用，拿不准的地方 AI 会回头问。下面每多交代一项都能少一次来回，越靠前越有用：

```text
用 csv-report skill 把下面的数据做成可视化报告：
- 文件：data/clinical.xlsx、data/ascat/*.segments.txt、data/delly_sv.csv
- 这些是什么：临床信息表；ASCAT 拷贝数片段（每个样本一个文件）；Delly 的 SV（bcftools query 导出）
- 参考基因组：hg38
- 患者编号：临床表里是 "Sample"，ASCAT 文件里是 sample 列，Delly 表里是 "SAMPLE"
- 按 "Histology" 分组上色
- 不要放进报告：姓名、MRN、手机号
- 标题："XX 队列拷贝数报告"，做好中文版
- 输出到 report/，做完告诉我每个文件怎么映射的、做了哪些假设
```

### 什么样的输入可以用

列名、列的顺序、列的数量都随意。真正决定能画什么图的，是"一行代表什么"：

| 想看的图 | 一行是 | 至少要有 | 有更好 |
|---|---|---|---|
| 队列概览 | 一个患者 | 患者 ID | 任意临床或分析指标（每列一张图） |
| 纯度 vs 倍性 | 一个患者 | 患者 ID、纯度、倍性 | 拟合优度 |
| ASCAT 图谱、拷贝数热图 | 一个拷贝数片段 | 患者 ID、染色体、起点、终点、拷贝数 | nMajor/nMinor、未取整值、logR |
| 结构变异图 | 一个 SV | 染色体、位置（多患者时加患者 ID） | 终点或第二断点、类型、FILTER、支持 reads 数 |
| LogR / BAF / 读深散点 | 一个 bin 或 SNP | 染色体、位置 | logR、BAF、拷贝数估计 |

会影响结果的几条规则：
- **表头：** 只能有一行。上面有标题行没关系，合并单元格和多级表头不行。
- **患者 ID：** 各个文件里要一致，或者能从文件名里取出来。
- **坐标：** 所有坐标用同一个参考基因组。
- **文件格式：** CSV、TSV、TXT、`.xlsx`、文本 VCF，可以是 `.gz` 压缩的。
  - `.xls` 要先另存为 `.xlsx`。
  - BCF 要先用 `bcftools view` 转成 VCF。
- **其他形状**（一列一个样本、一列一个染色体臂、基因层面的表）：AI 会先写一个小转换脚本。
- **目前还没有专门的图：** 基因层面的结果（表达、MAF 突变、富集分析）和生存曲线。

### 不用 AI，手动跑

```bash
SKILL=plugins/shown/skills/csv-report
python3 $SKILL/scripts/profile_table.py data/*.csv data/*.xlsx      # 看每列的类型和提示
python3 $SKILL/scripts/build_report.py my_mapping.json --check      # 校验映射
python3 $SKILL/scripts/build_report.py my_mapping.json -o report.html
```

映射格式见 [`references/mapping.md`](plugins/shown/skills/csv-report/references/mapping.md)，Delly / ASCAT 原生输出的对应写法见 [`references/formats.md`](plugins/shown/skills/csv-report/references/formats.md)。

## 支持的数据

| 类型 (`kind`) | 每行是 | 例子 |
|---|---|---|
| `patients` | 一位患者 | 临床信息表。列数随意，没映射的列会自动做成图表；也可以是 ASCAT 的 metrics / summary 表 |
| `segments` | 一个拷贝数片段 | ASCAT `*.segments.txt` / `segments_raw.txt`、Delly `seg.bed`，或任何带 染色体 / 起止 / 拷贝数 的表 |
| `variants` | 一个 SV / CNV | Delly VCF（文本）、`bcftools query` 导出表，或任何 SV 列表 |
| `bins` | 一个 bin / SNP | Delly `cov.gz`、ASCAT `Tumor_LogR.txt` / `Tumor_BAF.txt` |

- **自动处理**：GBK / UTF-8 编码、Excel 表头上方的标题行和表格下方的备注行、Excel 日期、逗号小数（`0,55`）、`chr1` / `1` / `23` 这类染色体写法、`deletion` / `TRA` / `<DEL>` / `缺失` 这类 SV 类型写法、百分数形式的纯度、参考基因组版本推断（hg19 / hg38 / CHM13）
- **靠 AI 写的 `transforms` 处理**：以 Mb 为单位的位置、`"2+1"` 这种合在一格里的拷贝数等不规范写法

## 报告内容

- **概览**
  - 关键指标卡片
  - AI 解读（明确标注由 AI 生成）
  - ASCAT 纯度 vs 倍性散点图
  - 每个患者字段一张分布图（字段多少就有多少张）
  - 可排序、可搜索的患者表
- **拷贝数**
  - 队列扩增 / 缺失频率图
  - 患者 × 基因组热图
- **结构变异**
  - 每位患者的 SV 数量，按 Delly SVTYPE 堆叠
  - 大小分布、染色体分布
  - SV 列表，可筛选 PASS 和类型
- **患者详情**
  - ASCAT 拷贝数图谱：取整版红 = nMajor、蓝 = nMinor；未取整版紫 = 总拷贝数、绿 = nMinor；横轴为真实染色体坐标
  - Delly 读深拷贝数图：黑点 + 绿色分段线
  - LogR / BAF：红点 + 蓝色分段 / 拟合值
  - SV 弧线图
  - 拖动放大，双击恢复
- **数据与映射**
  - 每个字段来自哪个文件的哪一列
  - 未使用的列、提示信息
  - 指标定义（WGD / GI / LOH 按 `ascat.metrics`）
  - `mapping.json` 原文，方便核对 AI 的判断

每个页面都可以在 English / 中文、浅色 / 深色之间切换。

## 隐私提醒

HTML 里嵌入了用于画图的数据，分享报告等于分享这些数据。AI 默认不会把姓名、身份证号、电话、住址、病历号等直接标识符放进报告。必要时可以在 `mapping.json` 里用 `exclude` / `hide_attributes` 再确认一遍。

## 仓库结构

```
.claude-plugin/marketplace.json        # 插件市场清单（Claude Code 和 Codex 都读）
plugins/shown/.claude-plugin/plugin.json
plugins/shown/skills/csv-report/
  SKILL.md                             # 给 AI 的工作流说明
  scripts/tableio.py                   # 读 CSV/TSV/XLSX/VCF（纯标准库）
  scripts/profile_table.py             # 列画像 + 角色提示
  scripts/build_report.py              # 映射 → 数据 → HTML
  assets/report_template.html          # 报告模板（原生 JS，无外部依赖）
  references/mapping.md                # mapping.json 说明
  references/formats.md                # Delly / ASCAT 输出格式
examples/                              # 假数据生成器、示例映射、示例报告
tests/                                 # python3 -m unittest discover -s tests
```

重新生成示例：

```bash
python3 examples/make_fake_data.py
python3 plugins/shown/skills/csv-report/scripts/build_report.py examples/mapping.json -o examples/report.html
```
