# Matrix Operations Performance Benchmark

基于 C 语言实现的矩阵运算库，包含加法、减法、标量乘法、转置和矩阵乘法，并附带性能基准测试功能。

## 工程结构

```
code/
├── main.c                    # 主程序入口（正确性验证 + 性能基准测试）
├── timer.h / timer.c         # 高精度计时器（Windows: QueryPerformanceCounter / Linux: clock_gettime）
├── Makefile                  # 构建配置
├── src/
│   ├── core/
│   │   ├── matrix_core.h     # 矩阵核心定义（Matrix 结构体、MatrixError 枚举、创建/释放/索引等接口）
│   │   ├── matrix_core.c     # 矩阵核心实现（初始化、创建、释放、索引映射、填充、拷贝、打印）
│   │   ├── matrix_ops.h      # 矩阵运算接口（加法、减法、标量乘法、转置、乘法、Frobenius范数）
│   │   └── matrix_ops.c      # 矩阵运算实现
│   └── util/
│       ├── matrix_rand.h     # 随机矩阵生成接口
│       └── matrix_rand.c     # 随机矩阵生成实现
└── *_results.csv             # 性能测试结果（运行后生成）
```

### 核心模块说明

- **matrix_core**: 定义 `Matrix` 结构体（行优先一维数组存储）和 `MatrixError` 错误码体系，提供矩阵的创建、释放、元素访问、填充和打印等基础操作。
- **matrix_ops**: 实现五种矩阵运算——`MatrixAdd`、`MatrixSub`、`MatrixScale`、`MatrixTranspose`、`MatrixMultiply`，以及 Frobenius 范数计算 `MatrixNormFrobenius`。乘法采用转置优化策略（先将 B 转置为 BT，再按行计算点积）。
- **matrix_rand**: 提供随机矩阵生成功能，用于性能测试数据准备。
- **timer**: 跨平台高精度计时器，Windows 下使用 `QueryPerformanceCounter`，Linux 下使用 `clock_gettime(CLOCK_MONOTONIC)`。

## 编译方法

依赖：GCC 编译器

```bash
make          # 编译生成可执行文件 demo
make clean    # 清理编译产物（.o 文件和 demo 可执行文件）
```

Makefile 编译选项：`-std=c99 -Wall -Wextra -pedantic -O0 -g`，头文件搜索路径包含 `./src/core`、`./src/util` 和当前目录。

需注意，如需在linux环境下使用 `make clean` 命令，需对Makefile文件做一定修改，将 `powershell -Command "Remove-Item -Force -ErrorAction SilentlyContinue demo, *.o, src\core\*.o, src\util\*.o"` 命令替换为 `rm -f $(TARGET) $(OBJ) main.o` 。

## 运行方法

```bash
make run      # 编译并运行
./demo        # 直接运行已编译的可执行文件（Windows 下为 demo.exe）
```

程序运行流程：
1. **正确性验证** — 使用 3x3 矩阵验证四种运算（加法、标量乘法、转置、乘法）的结果正确性并打印。
2. **性能基准测试** — 对不同规模的矩阵（100x100、1000x1000、2000x2000）分别测试加法、标量乘法、转置和乘法的平均耗时，结果输出到 CSV 文件：
   - `addition_results.csv`
   - `scaling_results.csv`
   - `transpose_results.csv`
   - `multiply_results.csv`

CSV 文件格式：各列分别为矩阵维度参数、平均耗时（ms）、测试重复次数。