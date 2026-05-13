#include "timer.h"

#ifdef _WIN32

#include <windows.h>

static LARGE_INTEGER g_freq;
static int g_initialized = 0;

static void timer_init() {
    if (!g_initialized) {
        QueryPerformanceFrequency(&g_freq);
        g_initialized = 1;
    }
}

void timer_start(Timer* t) {
    timer_init();

    LARGE_INTEGER now;
    QueryPerformanceCounter(&now);

    t->start = now.QuadPart;
}

double timer_elapsed_ms(const Timer* t) {
    LARGE_INTEGER now;
    QueryPerformanceCounter(&now);

    return
        (double)(now.QuadPart - t->start)
        * 1000.0
        / (double)g_freq.QuadPart;
}

double timer_elapsed_sec(const Timer* t) {
    return timer_elapsed_ms(t) / 1000.0;
}

#else

#include <time.h>

void timer_start(Timer* t) {
    struct timespec ts;

    clock_gettime(CLOCK_MONOTONIC, &ts);

    t->start_sec = ts.tv_sec;
    t->start_nsec = ts.tv_nsec;
}

double timer_elapsed_ms(const Timer* t) {
    struct timespec ts;

    clock_gettime(CLOCK_MONOTONIC, &ts);

    long long sec =
        ts.tv_sec - t->start_sec;

    long long nsec =
        ts.tv_nsec - t->start_nsec;

    return
        (double)sec * 1000.0
        + (double)nsec / 1000000.0;
}

double timer_elapsed_sec(const Timer* t) {
    return timer_elapsed_ms(t) / 1000.0;
}

#endif