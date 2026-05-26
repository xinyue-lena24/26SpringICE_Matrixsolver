#include "timer.h"

#ifdef _WIN32

#include <windows.h>

// Windows-specific implementation using QueryPerformanceCounter for high-resolution timing
static LARGE_INTEGER g_freq;
static int g_initialized = 0;

/**
 * @brief Initialize the timer by retrieving the frequency of the high-resolution performance counter. This function is called internally before starting the timer to ensure that the frequency is available for calculating elapsed time.
 *        The function checks if the timer has already been initialized to avoid redundant calls to QueryPerformanceFrequency, which can be costly. If the timer is not initialized, it retrieves the frequency and sets the initialized flag to true.
 */
static void timer_init() {
    if (!g_initialized) {
        QueryPerformanceFrequency(&g_freq);
        g_initialized = 1;
    }
}

/**
 * @brief Start the timer by recording the current value of the high-resolution performance counter. This function initializes the timer if it has not been initialized yet and then retrieves the current counter value, storing it in the Timer structure for later use in calculating elapsed time.
 * @param t Pointer to a Timer structure where the start time will be recorded.
 */
void timer_start(Timer* t) {
    timer_init();

    LARGE_INTEGER now;
    QueryPerformanceCounter(&now);

    t->start = now.QuadPart;
}

/**
 * @brief Calculate the elapsed time in milliseconds since the timer was started. This function retrieves the current value of the high-resolution performance counter, calculates the difference from the start time recorded in the Timer structure, and converts it to milliseconds using the frequency of the counter.
 * @param t Pointer to a Timer structure that contains the start time.
 * @return The elapsed time in milliseconds as a double.
 */
double timer_elapsed_ms(const Timer* t) {
    LARGE_INTEGER now;
    QueryPerformanceCounter(&now);

    return
        (double)(now.QuadPart - t->start)
        * 1000.0
        / (double)g_freq.QuadPart;
}

/**
 * @brief Calculate the elapsed time in seconds since the timer was started. This function retrieves the current value of the high-resolution performance counter, calculates the difference from the start time recorded in the Timer structure, and converts it to seconds using the frequency of the counter.
 * @param t Pointer to a Timer structure that contains the start time.
 * @return The elapsed time in seconds as a double.
 */
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