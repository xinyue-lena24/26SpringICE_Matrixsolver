# 小型 C 语言矩阵计算库基础版

本基础版本按照 `src/core` 与 `src/util` 的工程结构组织代码，当前只合并已经完成的基础矩阵功能、辅助工具和测试程序。LU 分解、线性方程组求解等后续扩展暂时不放入本版本，之后可以继续按同一结构添加。

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
│   └── util/
│       ├── matrix_rand.h
│       ├── matrix_rand.c
│       ├── timer.h
│       └── timer.c
└── results/
    └── .gitkeep
```

## 模块说明

- `src/core/matrix_core.h` / `src/core/matrix_core.c`：矩阵结构体、错误码、创建、释放、元素访问、填充、拷贝和打印。
- `src/core/matrix_ops.h` / `src/core/matrix_ops.c`：矩阵加法、减法、标量乘法、转置、乘法、范数和误差度量。
- `src/util/matrix_rand.h` / `src/util/matrix_rand.c`：随机矩阵填充与随机矩阵生成。
- `src/util/timer.h` / `src/util/timer.c`：跨平台计时工具。
- `main_basic.c`：小规模矩阵正确性测试。
- `main_timing_ops.c`：基础矩阵运算运行时间测试。
- `main_timing_mul.c`：不同矩阵乘法循环顺序的效率比较。
- `results/`：保存运行时间测试生成的 CSV 文件，例如 `ops_results.csv` 和 `mul_results.csv`。

## 当前已包含的基础功能

- `MatrixInit`
- `MatrixCreate`
- `MatrixFree`
- `MatrixSet`
- `MatrixGet`
- `MatrixFillZero`
- `MatrixFillConstant`
- `MatrixFillIdentity`
- `MatrixFillSequence`
- `MatrixCopy`
- `MatrixPrint`
- `MatrixAdd`
- `MatrixSub`
- `MatrixScale`
- `MatrixTranspose`
- `MatrixMultiply`
- `MatrixMultiplyIKJ`
- `MatrixNormFrobenius`
- `MatrixNormOne`
- `MatrixNormInfinity`
- `MatrixMaxAbsDiff`
- `MatrixRelativeErrorFrobenius`
- `MatrixFillRandom`
- `MatrixGenerateRandom`

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

之后可以继续添加：

```text
src/core/matrix_lu.h
src/core/matrix_lu.c
src/core/matrix_solve.h
src/core/matrix_solve.c
main_lu.c
```

并在 `Makefile` 中把新增源文件加入 `CORE_SRCS`，再添加对应的测试目标。
