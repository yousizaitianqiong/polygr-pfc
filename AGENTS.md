# Repository Guidelines

## 项目结构与模块组织

源码位于 `src/`。CPU 主求解器是 `polygr.c`，MPI 版本在
`polygr_mpi.c`，可选 CUDA 后端在 `polygr_cuda.cu`。密度场与 PNG
输出工具在 `field_image.*`，坐标提取逻辑在 `xyz.*`。`figure.c` 和
`figure_samples.c` 用于生成复合图和示例数据。

示例输入、生成的 `.dat` 快照、PNG 和 XYZ 输出位于 `example/`。
旧版 Java 参考工具位于 `java/`。集群运行脚本在 `slurm/`，设计和优化说明在
`docs/`。

## 构建、测试与开发命令

Windows + MSYS2 UCRT64 环境下使用：

```powershell
.\build-windows.ps1
```

该脚本会构建 `polygr.exe`、`polygr_figure.exe` 和
`polygr_figure_samples.exe`，并复制所需运行时 DLL。

Linux 环境下使用：

```bash
make          # 构建 CPU 求解器和绘图工具
make mpi      # 构建 MPI 求解器
make cuda     # 构建可选 CUDA 后端
make clean
```

快速本地运行示例：

```powershell
cd example
..\polygr.exe quick1 quick2 --xyz quick.xyz --png quick.png --threads 8
```

## 代码风格与命名约定

C 源码使用 C11，并保持 `-Wall -Wextra` 下无警告。遵循现有风格：四空格缩进、
`snake_case` 函数和变量名、文件内 `static` 小工具函数，以及简短明确的
`die()` 风格错误处理。源码和文档默认使用 ASCII，除非已有文件明确需要中文或
其他 Unicode 内容。

## 测试指南

当前没有正式单元测试框架。最低验证标准是构建成功并运行有针对性的示例。修改求
解器后，至少运行 quick 示例并确认原子数、XYZ 和 PNG 输出正常。修改绘图逻辑后，
在 `example/` 中重新生成代表性 PNG，并进行人工视觉检查。

## 提交与 Pull Request 指南

提交历史采用简洁的祈使式摘要，例如 `Refine CVD graphene figure generation`
和 `Enforce warning-free C builds`。每个提交应聚焦一个行为或子系统。

PR 应包含简短说明、已运行的命令、受影响平台，以及绘图或可视化变更的前后对比图。
如有关联 issue，应在 PR 中链接；如有意加入生成产物，也需明确说明。
