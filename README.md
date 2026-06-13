# 小型 C 语言矩阵计算库扩展版

本项目按照 `src/core`、`src/util` 与 `src/algorithm` 的工程结构组织代码。在基础矩阵运算的基础上，进一步加入 LU 分解、Gauss 消元法、带主元交换的 LU 分解、多右端项线性系统求解等数值线性代数功能。

项目默认使用纯 C 语言实现矩阵运算与线性方程组求解，不依赖外部 BLAS/LAPACK 库。同时，项目提供可选的 OpenBLAS/LAPACKE 对比测试，用于比较自实现求解器与 LAPACKE 求解器的运行时间。

---

## 工程结构

```text
project/
├── Makefile
├── README.md
├── environment.yml
├── main_basic.c
├── main_timing_ops.c
├── main_timing_mul.c
├── main_timing_solve.c
├── src/
│   ├── core/
│   │   ├── matrix_core.h
│   │   ├── matrix_core.c
│   │   ├── matrix_ops.h
│   │   └── matrix_ops.c
│   ├── util/
│   │   ├── matrix_rand.h
│   │   ├── matrix_rand.c
│   │   ├── timer.h
│   │   └── timer.c
│   └── algorithm/
│       ├── matrix_lu.h
│       ├── matrix_lu.c
│       ├── matrix_solve.h
│       └── matrix_solve.c
└── results/
    └── .gitkeep
```

---

## 模块说明

* `src/core/matrix_core.h` / `src/core/matrix_core.c`：矩阵结构体、错误码、创建、释放、元素访问、填充、拷贝和打印。
* `src/core/matrix_ops.h` / `src/core/matrix_ops.c`：矩阵加法、减法、标量乘法、转置、矩阵乘法、范数和误差度量。
* `src/util/matrix_rand.h` / `src/util/matrix_rand.c`：随机矩阵填充与随机矩阵生成。
* `src/util/timer.h` / `src/util/timer.c`：跨平台计时工具。
* `src/algorithm/matrix_lu.h` / `src/algorithm/matrix_lu.c`：LU 分解、带主元交换的 LU 分解、三角方程组求解和行列式计算。
* `src/algorithm/matrix_solve.h` / `src/algorithm/matrix_solve.c`：Gauss 消元法、批量 Gauss 消元法、LU 求解线性方程组。
* `main_basic.c`：小规模矩阵正确性测试。
* `main_timing_ops.c`：基础矩阵运算运行时间测试。
* `main_timing_mul.c`：不同矩阵乘法循环顺序的效率比较。
* `main_timing_solve.c`：线性方程组求解方法的运行时间测试，并可选与 OpenBLAS/LAPACKE 进行对比。
* `environment.yml`：可选 OpenBLAS/LAPACKE 对比环境的 Conda 配置文件。
* `results/`：保存运行时间测试生成的 CSV 文件。

---

## 当前已包含的功能

### 核心模块 `matrix_core`

* `MatrixInit` — 初始化矩阵结构体
* `MatrixCreate` — 创建矩阵并分配内存
* `MatrixFree` — 释放矩阵内存
* `MatrixIsValid` — 检查矩阵是否有效
* `MatrixHasShape` — 检查矩阵是否具有指定行列数
* `MatrixIndex` — 计算 row-major 存储下的一维数组偏移量
* `MatrixSet` — 设置矩阵元素值
* `MatrixGet` — 获取矩阵元素值
* `MatrixFillZero` — 填零
* `MatrixFillConstant` — 填常数
* `MatrixFillIdentity` — 填单位矩阵
* `MatrixFillSequence` — 填等差序列
* `MatrixCopy` — 矩阵拷贝
* `MatrixPrint` — 打印矩阵
* `MatrixPrintIndexMap` — 打印矩阵索引映射
* `MatrixErrorMessage` — 返回错误码对应的提示字符串

### 运算模块 `matrix_ops`

* `MatrixAdd` — 矩阵加法
* `MatrixSub` — 矩阵减法
* `MatrixScale` — 标量乘法
* `MatrixTranspose` — 矩阵转置
* `MatrixMultiply` — 矩阵乘法，采用 `ijk` 循环顺序
* `MatrixMultiplyIKJ` — 矩阵乘法，采用 `ikj` 循环顺序
* `MatrixNormFrobenius` — Frobenius 范数
* `MatrixNormOne` — 1-范数
* `MatrixNormInfinity` — 无穷范数
* `MatrixMaxAbsDiff` — 最大绝对差
* `MatrixRelativeErrorFrobenius` — Frobenius 相对误差

### 随机工具 `matrix_rand`

* `MatrixFillRandom` — 随机填充已有矩阵
* `MatrixGenerateRandom` — 创建并填充随机矩阵

### 计时工具 `timer`

* `timer_start` — 开始计时
* `timer_elapsed_ms` — 获取经过毫秒数
* `timer_elapsed_sec` — 获取经过秒数

### LU 分解模块 `matrix_lu`

* `LUDecomposeNoPivot` — 无主元 Doolittle LU 分解，计算 `A = LU`
* `LUDecomposeByEliminationNoPivot` — 通过记录 Gauss 消元乘子得到无主元 LU 分解
* `LUDecomposePartialPivot` — 带部分选主元的 LU 分解，计算 `PA = LU`
* `LUDecomposeByEliminationPartialPivot` — 通过记录 Gauss 消元乘子得到带主元交换的 LU 分解
* `ForwardSubstitution` — 前向代入，求解 `Ly = b`
* `ForwardSubstitutionMultiple` — 多右端项前向代入，求解 `LY = B`
* `BackSubstitution` — 回代，求解 `Ux = y`
* `BackSubstitutionMultiple` — 多右端项回代，求解 `UX = Y`
* `LUDeterminant` — 由无主元 LU 分解中的 `U` 计算行列式

### 线性方程组求解模块 `matrix_solve`

* `MatrixSwapRows` — 交换矩阵两行
* `GaussianSolvePartialPivot` — 部分选主元 Gauss 消元法求解 `Ax = b`
* `GaussianSolveMultiple` — 逐列调用 Gauss 消元法求解多右端项系统 `AX = B`
* `GaussianSolveMultipleBatch` — 对增广矩阵 `[A|B]` 进行一次 Gauss 消元，批量求解 `AX = B`
* `LUSolve` — 由已有无主元 LU 分解求解 `LUx = b`
* `LUSolveMultiple` — 由已有无主元 LU 分解求解 `LUX = B`
* `LUSolveWithPivot` — 由带主元 LU 分解求解 `LUx = Pb`
* `LUSolveMultipleWithPivot` — 由带主元 LU 分解求解 `LUX = PB`
* `LUDecomposeSolveMultiple` — 先计算无主元 LU 分解，再求解多右端项系统 `AX = B`

---

## 环境配置

本项目默认编译不依赖 OpenBLAS/LAPACK。基础矩阵运算、LU 分解、Gauss 消元法和线性方程组求解函数均由 C 语言显式循环实现。

默认情况下，只需要系统中已经安装：

```bash
gcc
make
```

默认编译只链接 C 标准数学库：

```text
-lm
```

### 默认编译

在工程根目录运行：

```bash
make clean
make all
```

该模式不会启用 OpenBLAS/LAPACKE 对比代码，因此不需要安装 OpenBLAS 或 LAPACKE。

### 可选：配置 OpenBLAS/LAPACKE 对比环境

本项目提供可选的 OpenBLAS/LAPACKE 对比测试，用于在 `main_timing_solve.c` 中比较自实现线性方程组求解函数与 LAPACKE 求解器的运行时间。该部分由宏 `COMPARE_TO_BLAS` 控制，默认不会启用。

项目提供 `environment.yml` 用于创建 OpenBLAS/LAPACKE 环境。

如果本机有权限创建普通 Conda 环境，可以运行：

```bash
conda env create -f environment.yml
conda activate c_matrix
```

如果在无管理员权限的服务器上，推荐把环境创建到用户目录：

```bash
mkdir -p ~/conda_envs
conda env create -p ~/conda_envs/c_matrix -f environment.yml
conda activate ~/conda_envs/c_matrix
```

环境创建并激活后，使用：

```bash
make clean
make all COMPARE_TO_BLAS=1
```

此时 Makefile 会定义：

```text
-DCOMPARE_TO_BLAS
```

并链接：

```text
-llapacke -llapack -lopenblas
```

因此，`main_timing_solve.c` 中由

```c
#ifdef COMPARE_TO_BLAS
```

保护的 LAPACKE 对比代码才会参与编译。

---

## `environment.yml` 内容

项目中的 `environment.yml` 内容如下：

```yaml
name: c_matrix

channels:
  - conda-forge

dependencies:
  - openblas
  - liblapack
  - liblapacke
```

该文件只用于配置 OpenBLAS/LAPACKE 对比环境。默认纯 C 编译不需要该环境。

---

## 编译与运行

### 编译全部程序

在工程根目录运行：

```bash
make
```

或：

```bash
make all
```

生成四个可执行文件：

```text
main_basic
main_timing_ops
main_timing_mul
main_timing_solve
```

### 分别运行测试

```bash
make run-basic
make run-ops
make run-mul
make run-solve
```

其中：

```bash
make run-basic
```

运行小规模正确性测试。

```bash
make run-ops
```

运行基础矩阵运算计时测试，并生成：

```text
results/ops_results.csv
```

```bash
make run-mul
```

运行不同矩阵乘法循环顺序的计时测试，并生成：

```text
results/mul_results.csv
```

```bash
make run-solve
```

运行线性方程组求解计时测试，并生成：

```text
results/lu_solve_matrix.csv
results/gaussian_solve_matrix.csv
results/lu_pivot_solve_matrix.csv
```

若编译时使用了 `COMPARE_TO_BLAS=1`，则 `make run-solve` 还会包含 OpenBLAS/LAPACKE 对比，并额外生成：

```text
results/blas_solve.csv
```

### 一次性运行全部测试

```bash
make run-all
```

### 清理编译产物

```bash
make clean
```

该命令会删除可执行文件以及 `results/*.csv`。

---

## OpenBLAS/LAPACKE 对比说明

`main_timing_solve.c` 中的 LAPACKE 对比代码默认不会参与编译。只有在编译时使用：

```bash
make all COMPARE_TO_BLAS=1
```

才会启用如下代码块：

```c
#ifdef COMPARE_TO_BLAS
...
#endif
```

其中调用的 LAPACKE 接口为：

```c
LAPACKE_dgesv
```

该函数用于求解一般稠密线性方程组。项目中的 `BlasSolve` 函数仅作为性能对比基准，不影响默认纯 C 版本的编译和运行。

---

## 结果文件

运行不同测试后，`results/` 文件夹中会生成 CSV 文件。常见输出包括：

```text
results/ops_results.csv
results/mul_results.csv
results/lu_solve_matrix.csv
results/gaussian_solve_matrix.csv
results/lu_pivot_solve_matrix.csv
results/blas_solve.csv
```

其中 `results/blas_solve.csv` 只会在启用 `COMPARE_TO_BLAS=1` 编译并运行 `make run-solve` 后生成。

---

## 清理与重新构建

如果修改了代码或 Makefile，建议重新构建：

```bash
make clean
make all
```

如果需要重新构建并启用 OpenBLAS/LAPACKE 对比：

```bash
make clean
make all COMPARE_TO_BLAS=1
```

---

## 备注

1. 默认版本不依赖 OpenBLAS/LAPACK，便于在普通 C 环境中直接编译。
2. OpenBLAS/LAPACKE 只用于性能对比，不是项目核心矩阵库的强制依赖。
3. 如果在服务器上没有写入系统 Conda 环境的权限，可以使用 `conda env create -p ~/conda_envs/c_matrix -f environment.yml` 在用户目录下创建环境。
4. 如果需要运行 OpenBLAS/LAPACKE 对比，请先激活对应 Conda 环境，再使用 `COMPARE_TO_BLAS=1` 编译。
