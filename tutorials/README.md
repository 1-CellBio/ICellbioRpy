# ICellbioRpy 教程集合

欢迎使用ICellbioRpy包！这里提供了完整的教程集合，帮助您掌握单细胞数据格式转换的各种方法。

## 📚 教程概览

### 🔧 环境安装教程（必须先完成）

#### [教程0a: Windows环境完整安装指南](tutorial_0_environment_setup_windows.html)
**适用用户**: Windows系统的完全初学者

**主要内容**:
- 从零开始安装R和RStudio
- 下载和安装Anaconda
- 创建conda虚拟环境
- 安装ICellbioRpy包及依赖
- 配置Python-R交互环境
- 超详细的故障排除指南

**预计时间**: 1-2小时

---

#### [教程0b: macOS环境完整安装指南](tutorial_0_environment_setup_macos.html)
**适用用户**: macOS系统的完全初学者

**主要内容**:
- 安装Xcode Command Line Tools
- 安装R和RStudio（包含Apple Silicon优化）
- 下载和安装Anaconda
- 创建conda虚拟环境
- 处理macOS特有的编译问题
- Apple Silicon (M1/M2/M3) 特殊优化

**预计时间**: 1-2小时

---

### 📊 数据转换教程

### [教程1: 将1CellBio ZIP文件转换为H5AD格式](tutorial_1_1cellbio_to_h5ad.html)
**适用场景**: 您有1CellBio分析结果，需要在Python/scanpy中进行下游分析

**主要内容**:
- Python环境配置详解
- 1CellBio数据结构解析
- H5AD格式转换步骤
- 数据验证和质量检查
- 常见问题排查

**预计时间**: 30-45分钟

---

### [教程2: 将1CellBio ZIP文件转换为Seurat对象](tutorial_2_1cellbio_to_seurat.html)  
**适用场景**: 您有1CellBio分析结果，需要在R/Seurat中进行分析

**主要内容**:
- Seurat环境准备
- 1CellBioData对象解析
- Seurat对象创建和验证
- 数据可视化示例
- 进阶分析建议

**预计时间**: 45-60分钟

---

### [教程3: 将Seurat对象转换为H5AD格式](tutorial_3_seurat_to_h5ad.html)
**适用场景**: 您在R中分析了数据，现在需要转到Python环境

**主要内容**:
- 跨语言数据转换原理
- Seurat对象结构检查
- 多种转换选项详解
- 数据一致性验证
- Python中的后续分析

**预计时间**: 40-50分钟

---

### [教程4: 读取10X MTX格式数据并转换为H5AD](tutorial_4_10x_mtx_to_h5ad.html)
**适用场景**: 您有多个10X Cell Ranger输出，需要整合并转换格式

**主要内容**:
- 10X数据格式详解
- CSV样本信息配置
- 高速数据读取方法
- 质量控制和过滤
- 多样本整合策略

**预计时间**: 50-70分钟

## 🚀 快速开始

### 准备工作

1. **安装ICellbioRpy包**
```r
# 如果还没有安装
# install.packages("devtools")
# devtools::install_github("your_username/ICellbioRpy")
library(ICellbioRpy)
```

2. **配置Python环境**（所有教程都需要）
```r
# 基本配置
configure_python_env(verbose = TRUE)
check_anndata_available()
```

3. **选择合适的教程**
   - 根据您的数据源选择对应教程
   - 建议按顺序学习，每个教程都有独特的知识点

### 学习建议

- **完全初学者**: 
  1. **必须先完成教程0a（Windows）或教程0b（macOS）**
  2. 然后从教程1开始学习数据转换
  3. 遇到问题时查看对应教程的排查部分
- **有基础用户**: 
  - 如果已有R、Python环境，可以直接从教程1开始
  - 但建议先运行环境测试脚本确认配置正确
- **疑难排查**: 每个教程都有详细的问题排查部分

## 🔧 环境要求

### R环境
- R ≥ 4.0.0
- 必需包: ICellbioRpy, Matrix, reticulate
- 推荐包: Seurat, data.table

### Python环境
- Python ≥ 3.7
- 必需包: anndata
- 推荐包: scanpy, pandas, numpy

### 硬件建议
- **内存**: 至少8GB RAM（大数据集需要更多）
- **存储**: 确保有足够的磁盘空间
- **CPU**: 多核CPU可提高处理速度

## 📊 教程对比

| 教程 | 输入格式 | 输出格式 | 难度 | 用时 | 主要用途 |
|------|---------|---------|------|------|---------|
| 教程0a | - | 环境配置 | ⭐ | 1-2h | Windows环境安装 |
| 教程0b | - | 环境配置 | ⭐ | 1-2h | macOS环境安装 |
| 教程1 | 1CellBio ZIP | H5AD | ⭐⭐ | 30-45min | Python分析 |
| 教程2 | 1CellBio ZIP | Seurat | ⭐⭐ | 45-60min | R分析 |
| 教程3 | Seurat RDS | H5AD | ⭐⭐⭐ | 40-50min | 跨语言转换 |
| 教程4 | 10X MTX | H5AD | ⭐⭐⭐⭐ | 50-70min | 多样本整合 |

## 🆘 获取帮助

### 遇到问题时

1. **查看对应教程的排查部分**
2. **检查函数文档**: `?函数名`
3. **查看包的整体文档**: `vignette("ICellbioRpy")`
4. **提交GitHub Issue**: [项目地址]

### 常见问题快速索引

- **Python环境配置失败** → 所有教程都有详细说明
- **内存不足** → 教程3和4有优化建议
- **文件路径问题** → 教程4有详细的路径检查方法
- **数据格式不兼容** → 每个教程都有数据验证步骤

## 🔄 教程更新

这些教程会随着包的更新而持续改进。建议：

- 定期查看教程更新
- 使用最新版本的ICellbioRpy包
- 遇到新问题时及时反馈

## 📖 相关资源

### 官方文档
- [Seurat官方教程](https://satijalab.org/seurat/)
- [scanpy官方文档](https://scanpy.readthedocs.io/)
- [AnnData格式说明](https://anndata.readthedocs.io/)

---

**开始您的单细胞数据转换之旅吧！** 🧬✨

如有疑问，欢迎联系维护者或在GitHub上提交Issue。