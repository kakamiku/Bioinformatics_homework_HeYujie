# EMP-Web 课程项目报告

会话 ID: K2ScCtzRaI78WfeTQn54VhoW
生成时间: 2026-09-22 15:54:29

## 视频测验进度
已通过 5 个步骤测验。

## 科学解读与假设

## 任务反思

- **course_transcriptomics / s1_import**
  先检查 count 矩阵行名：常见为 Ensembl、Symbol 或 Entrez。

若下游富集/可视化工具要求 Symbol，则需转换；注意数据库版本、一对多映射和基因丢失。

colData 的样本 ID 必须与 count 矩阵列名一一对应。可用 match() 检查，确保无重复、无缺失、无多余；顺序可不同，但必须按 ID 对齐。

- **course_transcriptomics / s2_prepare**
  DESeq2/edgeR 基于负二项分布，需要整数 count 和文库大小信息；TPM/FPKM 已标准化，会破坏模型假设。

低表达基因：用 filterByExpr（edgeR）或要求至少在一组中 count ≥ 10 等策略过滤，减少噪声并提高检验功效。

批次效应：设计阶段尽量平衡；分析时在模型中加批次协变量，如 ~ batch + condition。必要时可用 ComBat-seq、RUV 等，但优先模型校正，避免直接修改原始 count。

- **course_transcriptomics / s3_analysis**
  对比：如 condition: treated vs control。

协变量：如批次、性别、年龄、测序批次等，设计公式如 ~ batch + sex + condition。

显著性：常用 BH 校正 FDR < 0.05。

log2FC 阈值：常用 |log2FC| > 1（2 倍变化），但应依据生物学效应大小、数据变异和统计功效调整，并在报告中说明依据。

- **course_transcriptomics / s4_visualization**
  选择依据：padj 显著、|log2FC| 较大、表达丰度足够、与表型/通路相关、有文献支持。

示例格式：

基因 A：padj = …，log2FC = …，已知参与免疫应答；

基因 B：padj = …，log2FC = …，与细胞周期相关；

基因 C：padj = …，log2FC = …，可能是标志物。

警惕：统计显著不等于因果；功能注释可能来自其他细胞/条件；需 qPCR、蛋白或功能实验验证，避免仅凭单一基因过度解读。

- **course_transcriptomics / s5_interpretation**
  科学解读：结合 FDR、富集分数、基因覆盖度、方向一致性和通路冗余，联系实验背景，不把富集直接当作因果机制。

AI 解读人工修正：核对 GO/KEGG/Reactome 原始数据库和文献，纠正 AI 可能编造的文献、方向错误、过度概括或忽略基因方向性。

声明：AI 仅辅助整理和提示，最终生物学结论由研究者人工确认并负

## Learning Trace 摘要
共 47 条事件。

