#ifndef TIMER_H
#define TIMER_H

typedef struct Timer {
#ifdef _WIN32
    long long start;
#else
    long long start_sec;
    long long start_nsec;
#endif
} Timer;

/**
 * 开始计时
 */
void timer_start(Timer* t);

/**
 * 获取经过的毫秒数
 */
double timer_elapsed_ms(const Timer* t);

/**
 * 获取经过的秒数
 */
double timer_elapsed_sec(const Timer* t);

#endif