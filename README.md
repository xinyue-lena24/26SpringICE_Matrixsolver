# 小型 C 语言矩阵计算库基础版

本基础版本按照 `src/core`、`src/util` 与 `src/algorithm` 的工程结构组织代码。

## 工程结构

```text
project/
├── Makefile
├── README.md
├── main_basic.c
├── main_timing_ops.c
├── main_timing_mul.c
├── src/
│   ├── core/
│   │   ├── matrix_core.h
│   │   ├── matrix_core.c
│   │   ├── matrix_ops.h
│   │   └── matrix_ops.c
│   ├── util/
│   |   ├── matrix_rand.h
│   |   ├── matrix_rand.c
│   |   ├── timer.h
│   |   └── timer.c
|   └── algorithm/
|       ├── matrix_lu.c
|       ├── matrix_lu.h
|       ├── matrix_solve.c
|       └── matrix_solve.h
└── results/
    └── .gitkeep
```

## 模块说明

- `src/core/matrix_core.h` / `src/core/matrix_core.c`：矩阵结构体、错误码、创建、释放、元素访问、填充、拷贝和打印。
- `src/core/matrix_ops.h` / `src/core/matrix_ops.c`：矩阵加法、减法、标量乘法、转置、乘法、范数和误差度量。
- `src/util/matrix_rand.h` / `src/util/matrix_rand.c`：随机矩阵填充与随机矩阵生成。
- `src/util/timer.h` / `src/util/timer.c`：跨平台计时工具。
- `src/algorithm/matrix_lu.h` / `src/algorithm/matrix_lu.c`：LU分解、求行列式。
- `src/algorithm/matrix_solve.h` / `src/algorithm/matrix_solve.c`：求解线性方程组。
- `main_basic.c`：小规模矩阵正确性测试。
- `main_timing_ops.c`：基础矩阵运算运行时间测试。
- `main_timing_mul.c`：不同矩阵乘法循环顺序的效率比较。
- `results/`：保存运行时间测试生成的 CSV 文件，例如 `ops_results.csv` 和 `mul_results.csv`。

## 当前已包含的基础功能

**核心模块 (`matrix_core`)**

- `MatrixInit` — 初始化矩阵结构体
- `MatrixCreate` — 创建矩阵并分配内存
- `MatrixFree` — 释放矩阵内存
- `MatrixIsValid` — 检查矩阵是否有效（非空指针、data 非空）
- `MatrixHasShape` — 检查矩阵是否具有指定行列数
- `MatrixIndex` — 计算一维存储偏移量
- `MatrixSet` — 设置矩阵元素值
- `MatrixGet` — 获取矩阵元素值
- `MatrixFillZero` — 填零
- `MatrixFillConstant` — 填常数
- `MatrixFillIdentity` — 填单位矩阵
- `MatrixFillSequence` — 填等差序列
- `MatrixCopy` — 矩阵拷贝
- `MatrixPrint` — 打印矩阵
- `MatrixPrintIndexMap` — 打印矩阵索引映射
- `MatrixErrorMessage` — 返回错误码对应的提示字符串

**运算模块 (`matrix_ops`)**

- `MatrixAdd` — 矩阵加法
- `MatrixSub` — 矩阵减法
- `MatrixScale` — 标量乘法
- `MatrixTranspose` — 矩阵转置
- `MatrixMultiply` — 矩阵乘法（ijk 顺序）
- `MatrixMultiplyIKJ` — 矩阵乘法（ikj 顺序）
- `MatrixNormFrobenius` — Frobenius 范数
- `MatrixNormOne` — 1-范数
- `MatrixNormInfinity` — ∞-范数
- `MatrixMaxAbsDiff` — 最大绝对差
- `MatrixRelativeErrorFrobenius` — Frobenius 相对误差

**随机工具 (`matrix_rand`)**

- `MatrixFillRandom` — 随机填充已有矩阵
- `MatrixGenerateRandom` — 创建并填充随机矩阵

**计时工具 (`timer`)**

- `timer_start` — 开始计时
- `timer_elapsed_ms` — 获取经过毫秒数
- `timer_elapsed_sec` — 获取经过秒数

**LU 分解 (`matrix_lu`)**

- `LUDecomposeNoPivot` — 无选主元 LU 分解
- `ForwardSubstitution` — 前向代入
- `BackSubstitution` — 回代
- `LUDeterminant` — 由 U 计算行列式

**线性方程组求解 (`matrix_solve`)**

- `MatrixSwapRows` — 交换矩阵两行
- `GaussianSolvePartialPivot` — 部分选主元 Gauss 消元求解
- `GaussianSolveMultiple` — Gauss 消元求解多右端项
- `LUSolve` — 由已有 LU 分解求解单方程
- `LUDecomposeSolveMultiple` — LU 分解后求解多右端项

## 编译与运行

在工程根目录运行：

```bash
make
```

生成三个可执行文件：

```text
main_basic
main_timing_ops
main_timing_mul
```

分别运行：

```bash
make run-basic
make run-ops
make run-mul
```

一次性运行全部测试：

```bash
make run-all
```

运行 `make run-ops` 会生成 `results/ops_results.csv`，运行 `make run-mul` 会生成 `results/mul_results.csv`。

清理编译产物和 `results/*.csv`：

```bash
make clean
```

## 后续扩展建议