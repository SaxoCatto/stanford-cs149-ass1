#include <stdio.h>
#include <thread>
#include <iostream>
#include <vector>

using std::vector;
using std::cout;
using std::endl;

#include "CycleTimer.h"

typedef struct {
    float x0, x1;
    float y0, y1;
    unsigned int width;
    unsigned int height;
    int maxIterations;
    int* output;
    int threadId;
    int numThreads;
    int startRow;
    int numRows;
    double time;
} WorkerArgs;


extern void mandelbrotSerial(
    float x0, float y0, float x1, float y1,
    int width, int height,
    int startRow, int numRows,
    int maxIterations,
    int output[]);

// Determine whether a point belongs to the mandel set, and return the number of iterations. Re = real, im = imaginary
static inline int mandel(float c_re, float c_im, int count)
{
    float z_re = c_re, z_im = c_im;
    int i;
    for (i = 0; i < count; ++i) {
        // Determine whether the current complex number has tendency to infinity during the iteration of the Mandelbrot set.
        // Points in the Mandelbrot set either approach to infinity or remain bounded during the iteration 
        // (on the complex plane, points in the Mandelbrot set always remain within a circle with a radius of 2 during the iteration)
        if (z_re * z_re + z_im * z_im > 4.f)
            break;
        // iteration
        float new_re = z_re*z_re - z_im*z_im;
        float new_im = 2.f * z_re * z_im;
        z_re = c_re + new_re;
        z_im = c_im + new_im;
    }

    return i;
}

//
// workerThreadStart --
//
// Thread entrypoint.
void workerThreadStart(WorkerArgs * const args) {

    // TODO FOR CS149 STUDENTS: Implement the body of the worker
    // thread here. Each thread should make a call to mandelbrotSerial()
    // to compute a part of the output image.  For example, in a
    // program that uses two threads, thread 0 could compute the top
    // half of the image and thread 1 could compute the bottom half.

    printf("Hello world from thread %d\n", args->threadId);
    double startTime = CycleTimer::currentSeconds();

#ifdef PLAN_A
    float dx = (args->x1 - args->x0) / args->width;                 
    float dy = (args->y1 - args->y0) / args->height;                
    int endRow = args->startRow + args->numRows;
    for (int j = args->startRow; j < endRow; ++j) {
        for (unsigned int i = 0; i < args->width; ++i) {
            float x = args->x0 + i * dx;
            float y = args->y0 + j * dy;
            int index = j * args->width + i ;
            args->output[index] = mandel(x, y, args->maxIterations);
        }
    }
#else
    // Assign tasks according to the i-th grid and thread x, instead of according to the image area
    int step = args->numThreads;
    int id = args->threadId;
    float dx = (args->x1 - args->x0) / args->width;                 
    float dy = (args->y1 - args->y0) / args->height;                
    int endRow = args->startRow + args->numRows;
    // Thread x only processes pixels in the specific i-th grid of each row. i % numThreads = x.
    for (int j = args->startRow; j < endRow; ++j) {
        for (unsigned int i = id; i < args->width; i += step) {
            float x = args->x0 + i * dx;
            float y = args->y0 + j * dy;
            int index = j * args->width + i;
            if ((int) i % step != id) {
                cout<<"unexpected index and step!"<<index<<" "<<step<<" "<<" "<<id<<endl;
                exit(-1);
            }
            args->output[index] = mandel(x, y, args->maxIterations);
        }
    }
    // for (int j = id; j < endRow; j += step) {
    //     for (unsigned int i = 0; i < args->width; ++i) {
    //         float x = args->x0 + i * dx;
    //         float y = args->y0 + j * dy;
    //         int index = j * args->width + i ;
    //         args->output[index] = mandel(x, y, args->maxIterations);
    //     }
    // }

#endif

    double endTime = CycleTimer::currentSeconds();
    args->time = (endTime - startTime) * 1000;

}


//
// MandelbrotThread --
//
// Multi-threaded implementation of mandelbrot set image generation.
// Threads of execution are created by spawning std::threads.
vector<double> mandelbrotThread(
    int numThreads,
    float x0, float y0, float x1, float y1,
    int width, int height,
    int maxIterations, int output[])
{
    static constexpr int MAX_THREADS = 32;

    if (numThreads > MAX_THREADS)
    {
        fprintf(stderr, "Error: Max allowed threads is %d\n", MAX_THREADS);
        exit(1);
    }

    vector<double> vec(numThreads);

    // Creates thread objects that do not yet represent a thread.
    std::thread workers[MAX_THREADS];
    WorkerArgs args[MAX_THREADS];
    
    int y_interval = height / numThreads;    
    for (int i=0; i<numThreads; i++) {
        // TODO FOR CS149 STUDENTS: You may or may not wish to modify
        // the per-thread arguments here.  The code below copies the
        // same arguments for each thread
    #ifdef PLAN_A
        args[i].x0 = x0;
        args[i].y0 = y0;
        args[i].x1 = x1;
        args[i].y1 = y1;
        args[i].width = width;
        args[i].height = height;
        args[i].maxIterations = maxIterations;
        args[i].numThreads = numThreads;
        args[i].output = output;
        args[i].threadId = i;
        if (i == numThreads - 1) {                  
            args[i].numRows = height - y_interval * i;
        } else {
            args[i].numRows = y_interval;
        }
        args[i].startRow = i * y_interval;          
    #else
        args[i].x0 = x0;
        args[i].y0 = y0;
        args[i].x1 = x1;
        args[i].y1 = y1;
        args[i].width = width;
        args[i].height = height;
        args[i].maxIterations = maxIterations;
        args[i].numThreads = numThreads;
        args[i].output = output;
        args[i].threadId = i;
        args[i].numRows = height;
        args[i].startRow = 0;
    #endif
    }

    // Spawn the worker threads.  Note that only numThreads-1 std::threads
    // are created and the main application thread is used as a worker
    // as well.
    for (int i=1; i<numThreads; i++) {
        workers[i] = std::thread(workerThreadStart, &args[i]);
    }
    
    workerThreadStart(&args[0]);

    // join worker threads
    for (int i=1; i<numThreads; i++) {
        workers[i].join();
    }

    for (int i = 0; i < numThreads; ++i) {
        vec[i] = args[i].time;
    }

    return vec;
}
