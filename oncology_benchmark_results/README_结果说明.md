# oncology_test_benchmark 盲测结果

按 csv-report 技能的流程跑：先用 profile_table.py 看列，再写 mapping.json，`--check` 通过后生成 HTML。全程没有打开 ANSWER_KEY.md 和 make_test_data.py。5 份报告都生成成功，浏览器里打开没有 JS 报错。每个数据集的 mapping 放在 `mappings/`，每个页签的截图放在 `screenshots/`。

## 01 临床极简（case_notes.csv，8 人 × 2 列）
- 映射：`病例号` 映射为患者 ID，`生存(月)` 作为数值属性，单位为月。
- 页面：只有“概览 / 患者 / 数据与映射”三个页签。拷贝数和 SV 页签自动隐藏，患者页显示“无基因组数据”，没有凭空生成纯度、CNV 或 SV。
- 问题：页眉写着 “Genome hg38”，可这个数据集根本没有坐标，容易误导。

## 02 ASCAT 实验室流程（7 人，3 个文件）
- 映射：
  - `puirty`（拼写错误）映射为 purity，`psi_fit`（全部为 2）映射为 ploidy。
  - 分段文件里的 `chromsome / start.pos / end.pos / nMajor / nMinor / logR_mean / nSNP` 都映射上了。
  - 信号文件里 `PID` 为患者，`position_hg19` 为坐标，`L2R` 为 logR，`B_allele_frequency` 为 BAF，基因组按 hg19。
- 页面：ASCAT 红蓝 profile、LogR、BAF（带拟合线）三条轨道都正常显示。`cn.em` 和 nMajor+nMinor 逐行一致，没有矛盾。
- 问题：
  - 分析器把 `psi_fit` 提示成 goodness_of_fit，实际上它是倍性。如果照着提示映射就会出错。
  - 空列 `Unnamed: 9` 被标成“疑似直接身份标识”，属于误报。
  - `PID` 和 `position_hg19` 没有得到 patient / pos 提示。
  - 分组图例里 BG3 被归到 “Other”，虽然一共只有 4 组。
  - 没有性别列，所以 chrY 71%、chrX 57% 的“缺失”多半只是男女混合造成的，而 builder 对此没有任何警告。

## 03 Delly 断点（3 人）
- 映射：
  - 多样本 VCF 按 P001/P002/P003 各写一个条目，按 GT 过滤掉 0/0，`DV` 作为 PE 支持数。
  - `coverage.cov.gz` 中的 `depth_ratio` 经过 log2 转成 logR，起点 `from0` 是 0-based。
  - `seg.bed` 没有样本列，但它的边界和 P001 的覆盖度文件逐行吻合，所以按 P001 处理（这是推断）。
- 页面：每人 SV 堆叠图、大小分布、染色体分布、调用表都正常显示。TRA/BND 显示为 BND 弧线，没有当成 DEL。
- 问题：
  - 只有 logR 分段、没有 CN 时，患者页仍然画出一条空的 “ASCAT copy-number profile” 轨道。
  - 分析器把 VCF 的 `ID` 列（SV0001…）提示成 patient。

## 04 混乱多表（5 人）
- 映射：
  - `audit_export_38cols.tsv` 作为主临床表：`purity_pct` 自动除以 100，`psi` 映射为 ploidy，`xingbie` 映射为 sex，`zl_cd` 在这里是肿瘤类型并用来分组。`fake_name / tel_fake / hospital_fake / 空列 / col_7` 已隐藏。
  - Excel 的 Clinical notes 表：标题行被自动跳过；`P-` 前缀已去掉；尾注行 “End of sheet…” 会被当成患者，我用 row_filter 把它去掉了。
  - `片段_本地流程.csv` 是 GB18030 编码、分号分隔：`zl_cd` 在这里是患者 ID，`1:1-3187722` 用正则拆开，`2+2` 拆成 major/minor，逗号小数先替换成点再解析。
  - Arm-CN wide 表：写了一个小脚本，汇总成每人获得 / 缺失的臂数。
- 页面：概览、拷贝数热图、患者 profile 都能画出来。
- 发现：
  - **分析器把逗号小数列 `0,0093` 当成整数**（显示 min -23855、max 9000）。如果不手动加 replace，logR 会差 1 万倍。
  - Arm 表里每条染色体的 p 臂和 q 臂数值完全相同，就等于用分段算出的整条染色体长度加权均值，所以它不是独立的臂水平数据。
  - 审核表的 `bh` 列，P001 为 5.1，而 `psi` 为 4，分段推出的倍性约 4.0，两者对不上。我只把它作为属性展示，没有用来设倍性。
  - HRD（32–42 共 5 个值）被自动归成分类型，画成了条形图，应该是数值型。38 列全部展示，概览页很长。
  - `chr1` 最大坐标 248,956,422 与 hg38 吻合。
- unsafe_coordinates.tsv：没有 CN 值，不能当分段画图。我单独做了 check 探测：
  - chr1:249,250,616 超出 hg38 长度，有警告 ✓
  - chrM 和 chrUn 被排除，有警告 ✓
  - **P002 的 start=90000 > end=45000 没有任何警告，被静默接受了。**

## 05 CNVkit 规模（200 人）
- 映射：
  - `cohort.call.cns` 作为分段：log2、cn、cn1/cn2、probes。
  - `cohort.snp.tsv.gz`（30 万行）作为 bins：`value` 为 logR，`xb` 为 BAF。
  - `P001.cnr` 的样本名从文件名中提取。
- 页面：200 人的热图和频率图都能渲染，整个 check 约 6.5 秒。HTML 有 9.5 MB，仍能离线打开。所有样本都覆盖 1–22、X、Y。
- 问题：
  - 没有性别信息时，X/Y 一大片显示成“缺失”（X 62%，Y 49%），没有警告。
  - `genome: auto` 推断出 hg19，但没有像 SKILL.md 说的那样给出提示警告。
  - P001 同时有 SNP 和 .cnr 两个 bins 来源，会被混在一起画，也没有警告。
  - `P001.cnr` 只有 chr1–6 的 40 行。

## 值得考虑修的地方（需要你决定）
1. 逗号小数被解析成整数（最严重，会静默算错）。
2. 分段 start > end 不报警告。
3. 未映射 sex 时，X/Y 大面积“缺失”不报警告；genome 自动推断不报警告。
4. 分析器提示的问题：psi 被提示成 GoF、空列被误报为身份标识、VCF ID 被提示成 patient、HRD 被当成分类型。
5. 展示上的小问题：没有坐标时也显示 genome；没有 CN 时画出空轨道；4 个分组里有一组变成 “Other”。
