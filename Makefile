# 编译器
CC = gcc

# 编译选项
CFLAGS = -std=c99 -Wall -Wextra -pedantic -O0 -g -lm \
         -I./src/core \
		 -I./src/util \
		 -I./

# 目标
TARGET = demo

# 源文件
SRC = src/core/matrix_core.c \
	  src/core/matrix_ops.c \
	  src/util/matrix_rand.c \
	  timer.c

# 目标文件
OBJ = $(SRC:.c=.o)

# 默认目标
all: $(TARGET)

# 链接
$(TARGET): main.o $(OBJ)
	$(CC) $(CFLAGS) -o $(TARGET) main.o $(OBJ)

# 编译规则（支持子目录）
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# 运行
run: $(TARGET)
	./$(TARGET)

# 清理
clean:
	rm -f $(TARGET) $(OBJ) main.o