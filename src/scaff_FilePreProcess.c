#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <x86intrin.h>
#include <stdatomic.h>
#include <zlib.h>
#include <pthread.h>
#include <unistd.h>

typedef int8_t s08;
typedef int16_t s16;
typedef int32_t s32;
typedef int64_t s64;

typedef uint8_t u08;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;

typedef double f64;

#define global_function static
#define global_variable static

#define Min(x, y) (x < y ? x : y)
#define Pow2(N) (1 << N)

#define ArrayCount(array) (sizeof(array) / sizeof(array[0]))
#define ForLoop(n) for (u32 index = 0; index < (n); ++index)
#define ForLoop2(n) for (u32 index2 = 0; index2 < (n); ++index2)
#define TraverseLinkedList(startNode, type) for (type *(node) = (startNode); node; node = node->next)

#define ArgCount argc
#define ArgBuffer argv
#define Main s32 main()
#define MainArgs s32 main(s32 ArgCount, const char *ArgBuffer[])
#define EndMain return(0)

#define ThreadFence __asm__ volatile("" ::: "memory")
#define FenceIn(x) ThreadFence; \
	x; \
	ThreadFence

typedef pthread_t thread;
typedef pthread_mutex_t mutex;
typedef pthread_cond_t cond;

typedef volatile u32 threadSig;

#define CreateThread(x) thread *x
#define CreateMutex(x) mutex *x
#define CreateCond(x) cond *x

#define InitialiseMutex(x) x = (mutex)PTHREAD_MUTEX_INITIALIZER
#define InitialiseCond(x) x = (cond)PTHREAD_COND_INITIALIZER

#define LaunchThread(thread, func, dataIn) pthread_create(&thread, NULL, func, dataIn);
#define WaitForThread(x) pthread_join(*x, NULL)
#define DetachThread(thread) pthread_detach(thread)
#define PauseThread(thread) pthread_kill(thread, SIGUSR1)

#define LockMutex(x) pthread_mutex_lock(&x)
#define TryToLockMutex(x) pthread_mutex_trylock(&x)
#define UnlockMutex(x) pthread_mutex_unlock(&x)
#define WaitOnCond(cond, mutex) pthread_cond_wait(&cond, &mutex) 
#define SignalCondition(x) pthread_cond_signal(&x)
#define BroadcastCondition(x) pthread_cond_broadcast(&x)

struct
binary_semaphore
{
    mutex mut;
    cond con;
    u64 v;
};

struct
thread_job
{
    struct thread_job *prev;
    void (*function)(void *arg);
    void *arg;
};

struct
job_queue
{
    mutex rwMutex;
    struct thread_job *front;
    struct thread_job *rear;
    struct binary_semaphore *hasJobs;
    u32 len;
    u32 nFree;
    struct binary_semaphore *hasFree;
    struct thread_job *freeFront;
    struct thread_job *freeRear;
};

struct thread_context;

struct
thread_pool
{
    struct thread_context **threads;
    threadSig numThreadsAlive;
    threadSig numThreadsWorking;
    mutex threadCountLock;
    cond threadsAllIdle;
    struct job_queue jobQueue;
};

struct
thread_context
{
    u64	id;
    thread th;
    struct thread_pool *pool;
};

#define KiloByte(x) 1024*x
#define MegaByte(x) 1024*KiloByte(x)
#define GigaByte(x) 1024*MegaByte(x)

#define Default_Memory_Alignment_Pow2 4

struct
memory_arena
{
	u08 *base;
	u64 currentSize;
	u64 maxSize;
};

struct
memory_arena_snapshot
{
    u64 size;
};

global_function
void
TakeMemoryArenaSnapshot(struct memory_arena *arena, struct memory_arena_snapshot *snapshot)
{
    snapshot->size = arena->currentSize;
}

global_function
void
RestoreMemoryArenaFromSnapshot(struct memory_arena *arena, struct memory_arena_snapshot *snapshot)
{
    arena->currentSize = snapshot->size;
}

global_function
void
CreateMemoryArena_(struct memory_arena *arena, u64 size)
{
	posix_memalign((void **)&arena->base, Pow2(Default_Memory_Alignment_Pow2), size);
	arena->currentSize = 0;
	arena->maxSize = size;
}

#define CreateMemoryArena(arena, size, ...) CreateMemoryArena_(&arena, size, ##__VA_ARGS__)
#define CreateMemoryArenaP(arena, size, ...) CreateMemoryArena_(arena, size, ##__VA_ARGS__)

global_function
void
ResetMemoryArena_(struct memory_arena *arena)
{
	arena->currentSize = 0;
}

#define ResetMemoryArena(arena) ResetMemoryArena_(&arena)
#define ResetMemoryArenaP(arena) ResetMemoryArena_(arena)

global_function
void
FreeMemoryArena_(struct memory_arena *arena)
{
	free(arena->base);
}

#define FreeMemoryArena(arena) FreeMemoryArena_(&arena)
#define FreeMemoryArenaP(arena) FreeMemoryArena_(arena)

global_function
u64
GetAlignmentPadding(u64 base, u32 alignment_pow2)
{
	u64 alignment = (u64)Pow2(alignment_pow2);
	u64 result = ((base + alignment - 1) & ~(alignment - 1)) - base;

	return(result);
}

global_function
u32
AlignUp(u32 x, u32 alignment_pow2)
{
	u32 alignment_m1 = Pow2(alignment_pow2) - 1;
	u32 result = (x + alignment_m1) & ~alignment_m1;

	return(result);
}

global_function
void *
PushSize_(struct memory_arena *arena, u64 size)
{
	u64 padding = GetAlignmentPadding((u64)(arena->base + arena->currentSize), Default_Memory_Alignment_Pow2);
	
	void *result;
	if((size + arena->currentSize + padding + sizeof(u64)) > arena->maxSize)
	{
		result = 0;
		fprintf(stderr, "Push of %lu bytes failed, out of memory.\n", size);
		exit(1);
	}
	else
	{
		result = arena->base + arena->currentSize + padding;
		arena->currentSize += (size + padding + sizeof(u64));
		*((u64 *)(arena->base + arena->currentSize - sizeof(u64))) = (size + padding);
	}
	
	return(result);
}

global_function
void
FreeLastPush_(struct memory_arena *arena)
{
	if (arena->currentSize)
	{
		u64 sizeToRemove = *((u64 *)(arena->base + arena->currentSize - sizeof(u64)));
		arena->currentSize -= (sizeToRemove + sizeof(u64));
	}
}

#define PushStruct(arena, type, ...) (type *)PushSize_(&arena, sizeof(type), ##__VA_ARGS__)
#define PushArray(arena, type, n, ...) (type *)PushSize_(&arena, sizeof(type) * n, ##__VA_ARGS__)
#define PushStructP(arena, type, ...) (type *)PushSize_(arena, sizeof(type), ##__VA_ARGS__)
#define PushArrayP(arena, type, n, ...) (type *)PushSize_(arena, sizeof(type) * n, ##__VA_ARGS__)

#define FreeLastPush(arena) FreeLastPush_(&arena)
#define FreeLastPushP(arena) FreeLastPush_(arena)

global_function
struct memory_arena *
PushSubArena_(struct memory_arena *mainArena, u64 size)
{
	struct memory_arena *subArena = PushStructP(mainArena, struct memory_arena);
	subArena->base = PushArrayP(mainArena, u08, size);
	subArena->currentSize = 0;
	subArena->maxSize = size;

	return(subArena);
}

#define PushSubArena(arena, size, ...) PushSubArena_(&arena, size, ##__VA_ARGS__)
#define PushSubArenaP(arena, size, ...) PushSubArena_(arena, size, ##__VA_ARGS__)

global_variable
threadSig
Threads_KeepAlive;

global_function
void
BinarySemaphoreInit(struct binary_semaphore *bsem, u32 value)
{
    InitialiseMutex(bsem->mut);
    InitialiseCond(bsem->con);
    bsem->v = value;
}

global_function
void
BinarySemaphoreWait(struct binary_semaphore *bsem)
{
    LockMutex(bsem->mut);
    
    while (bsem->v != 1)
    {
	WaitOnCond(bsem->con, bsem->mut);
    }
    
    bsem->v = 0;
    UnlockMutex(bsem->mut);
}

global_function
void
BinarySemaphorePost(struct binary_semaphore *bsem)
{
    LockMutex(bsem->mut);
    bsem->v = 1;
    SignalCondition(bsem->con);
    UnlockMutex(bsem->mut);
}

global_function
void
BinarySemaphorePostAll(struct binary_semaphore *bsem)
{
    LockMutex(bsem->mut);
    bsem->v = 1;
    BroadcastCondition(bsem->con);
    UnlockMutex(bsem->mut);
}

global_function
struct thread_job *
JobQueuePull(struct job_queue *jobQueue)
{
    LockMutex(jobQueue->rwMutex);
    struct thread_job *job = jobQueue->front;
    
    switch (jobQueue->len)
    {
	case 0:
	    break;

	case 1:
	    jobQueue->front = 0;
	    jobQueue->rear  = 0;
	    jobQueue->len = 0;
	    break;

	default:
	    jobQueue->front = job->prev;
	    --jobQueue->len;
	    BinarySemaphorePost(jobQueue->hasJobs);
    }

    UnlockMutex(jobQueue->rwMutex);
    
    return(job);
}

global_function
struct thread_job *
GetFreeThreadJob(struct job_queue *jobQueue)
{
    LockMutex(jobQueue->rwMutex);
    struct thread_job *job = jobQueue->freeFront;

    switch (jobQueue->nFree)
    {
	case 0:
	    break;
	
	case 1:
	    jobQueue->freeFront = 0;
	    jobQueue->freeRear = 0;
	    jobQueue->nFree = 0;
	    break;
	
	default:
	    jobQueue->freeFront = job->prev;
	    --jobQueue->nFree;
	    BinarySemaphorePost(jobQueue->hasFree);
    }

    UnlockMutex(jobQueue->rwMutex);

    return(job);
}

global_function
void
FreeThreadJob(struct job_queue *jobQueue, struct thread_job *job)
{
    LockMutex(jobQueue->rwMutex);
    job->prev = 0;

    switch (jobQueue->nFree)
    {
	case 0:
	    jobQueue->freeFront = job;
	    jobQueue->freeRear  = job;
	    break;

	default:
	    jobQueue->freeRear->prev = job;
	    jobQueue->freeRear = job;
    }
    ++jobQueue->nFree;

    BinarySemaphorePost(jobQueue->hasFree);	
    UnlockMutex(jobQueue->rwMutex);   
}

global_function
void *
ThreadFunc(void *in)
{
    struct thread_context *context = (struct thread_context *)in;
    
    struct thread_pool *pool = context->pool;

    LockMutex(pool->threadCountLock);
    pool->numThreadsAlive += 1;
    UnlockMutex(pool->threadCountLock);

    while (Threads_KeepAlive)
    {
	    BinarySemaphoreWait(pool->jobQueue.hasJobs);

        if (Threads_KeepAlive)
        {
            LockMutex(pool->threadCountLock);
            ++pool->numThreadsWorking;
            UnlockMutex(pool->threadCountLock);
                
            void (*funcBuff)(void*);
            void *argBuff;
            struct thread_job *job = JobQueuePull(&pool->jobQueue);

            if (job)
            {
            funcBuff = job->function;
            argBuff  = job->arg;
            funcBuff(argBuff);
            FreeThreadJob(&pool->jobQueue, job);
            }
                
            LockMutex(pool->threadCountLock);
            --pool->numThreadsWorking;

            if (!pool->numThreadsWorking)
            {
            SignalCondition(pool->threadsAllIdle);
            }

            UnlockMutex(pool->threadCountLock);
        }
    }
    
    LockMutex(pool->threadCountLock);
    --pool->numThreadsAlive;
    UnlockMutex(pool->threadCountLock);

    return(NULL);
}

global_function
void 
ThreadInit(struct memory_arena *arena, struct thread_pool *pool, struct thread_context **context, u32 id)
{
    *context = PushStructP(arena, struct thread_context);

    (*context)->pool = pool;
    (*context)->id = id;

    LaunchThread((*context)->th, ThreadFunc, *context);
    DetachThread((*context)->th);
}

#define Number_Thread_Jobs 1024
global_function
void
JobQueueInit(struct memory_arena *arena, struct job_queue *jobQueue)
{
    jobQueue->hasJobs = PushStructP(arena, struct binary_semaphore);
    jobQueue->hasFree = PushStructP(arena, struct binary_semaphore);

    InitialiseMutex(jobQueue->rwMutex);
    BinarySemaphoreInit(jobQueue->hasJobs, 0);
    BinarySemaphoreInit(jobQueue->hasFree, 0);

    jobQueue->len = 0;
    jobQueue->front = 0;
    jobQueue->rear = 0;

    jobQueue->nFree = 0;
    for (   u32 index = 0;
	    index < Number_Thread_Jobs;
	    ++index )
    {
	    struct thread_job *job = PushStructP(arena, struct thread_job);
	    FreeThreadJob(jobQueue, job);
    }
}

global_function
void
JobQueueClear(struct job_queue *jobQueue)
{
    while (jobQueue->len)
    {
	    struct thread_job *job = JobQueuePull(jobQueue);
	    FreeThreadJob(jobQueue, job);
    }

    jobQueue->front = 0;
    jobQueue->rear  = 0;
    BinarySemaphoreInit(jobQueue->hasJobs, 0);
    jobQueue->len = 0;
}

global_function
void
JobQueuePush(struct job_queue *jobQueue, struct thread_job *job)
{
    LockMutex(jobQueue->rwMutex);
    job->prev = 0;

    switch (jobQueue->len)
    {
        case 0:
            jobQueue->front = job;
            jobQueue->rear  = job;
            break;

        default:
            jobQueue->rear->prev = job;
            jobQueue->rear = job;
    }
    ++jobQueue->len;

    BinarySemaphorePost(jobQueue->hasJobs);	
    UnlockMutex(jobQueue->rwMutex);
}

global_function
struct thread_pool *
ThreadPoolInit(struct memory_arena *arena, u32 nThreads)
{
    Threads_KeepAlive = 1;

    struct thread_pool *threadPool = PushStructP(arena, struct thread_pool);
    threadPool->numThreadsAlive = 0;
    threadPool->numThreadsWorking = 0;

    JobQueueInit(arena, &threadPool->jobQueue);

    threadPool->threads = PushArrayP(arena, struct thread_context*, nThreads);
	
    InitialiseMutex(threadPool->threadCountLock);
    InitialiseCond(threadPool->threadsAllIdle);
	
    for (   u32 index = 0;
	    index < nThreads;
	    ++index )
    {
	    ThreadInit(arena, threadPool, threadPool->threads + index, index);
    }

    while (threadPool->numThreadsAlive != nThreads) {}

    return(threadPool);
}

#define ThreadPoolAddTask(pool, func, args) ThreadPoolAddWork(pool, (void (*)(void *))func, (void *)args)

global_function
void
ThreadPoolAddWork(struct thread_pool *threadPool, void (*function)(void*), void *arg)
{
    struct thread_job *job;
    
    BinarySemaphoreWait(threadPool->jobQueue.hasFree);
    while (!(job = GetFreeThreadJob(&threadPool->jobQueue)))
    {
	    sleep(1);
	    BinarySemaphoreWait(threadPool->jobQueue.hasFree);
    }
    
    job->function = function;
    job->arg = arg;

    JobQueuePush(&threadPool->jobQueue, job);
}

global_function
void
ThreadPoolWait(struct thread_pool *threadPool)
{
    LockMutex(threadPool->threadCountLock);
    
    while (threadPool->jobQueue.len || threadPool->numThreadsWorking)
    {
	    WaitOnCond(threadPool->threadsAllIdle, threadPool->threadCountLock);
    }
    
    UnlockMutex(threadPool->threadCountLock);
}

global_function
void
ThreadPoolDestroy(struct thread_pool *threadPool)
{
    if (threadPool)
    {
        Threads_KeepAlive = 0;

        f64 timeout = 1.0;
        time_t start, end;
        f64 tPassed = 0.0;
        time (&start);
        while (tPassed < timeout && threadPool->numThreadsAlive)
        {
            BinarySemaphorePostAll(threadPool->jobQueue.hasJobs);
            time (&end);
            tPassed = difftime(end, start);
        }

        while (threadPool->numThreadsAlive)
        {
            BinarySemaphorePostAll(threadPool->jobQueue.hasJobs);
            sleep(1);
        }

        JobQueueClear(&threadPool->jobQueue);
    }
}

global_function
u32
AreNullTerminatedStringsEqual(u08 *string1, u08 *string2)
{
	u32 result;
	do
	{
		result = (*string1 == *(string2++));
	} while(result && (*(string1++) != '\0'));

	return(result);
}

global_function
u32
CopyNullTerminatedString(u08 *source, u08 *dest)
{
	u32 stringLength = 0;

	while(*source != '\0')
	{
		*(dest++) = *(source++);
		++stringLength;
	}
	*dest = '\0';

	return(stringLength);
}

global_function
u32
StringLength(u08 *string)
{
    	u32 length = 0;
    
    	while(*string++ != '\0') ++length;
    	
    	return(length);
}

global_function
u32
StringToInt(u08 *stringEnd, u32 length)
{
    u32 result = 0;
    u32 pow = 1;

    while (--length > 0)
    {
	result += (u32)(*--stringEnd - '0') * pow;
	pow *= 10;
    }
    result += (u32)(*--stringEnd - '0') * pow;

    return(result);
}

global_function
u32
NullTerminatedStringToInt(u08 *string)
{
    u32 result = 0;
    
    u32 length = 0;
    while(*string++ != '\0') ++length;
    --string;

    u32 pow = 1;

    while (--length > 0)
    {
	result += (u32)(*--string - '0') * pow;
	pow *= 10;
    }
    result += (u32)(*--string - '0') * pow;

    return(result);
}

global_variable
struct memory_arena
Working_Set;

global_variable
struct thread_pool *
Thread_Pool;

#define Default_Number_of_Threads 8

#define Line_Buffer_Size 508
#define Number_Of_Line_Buffers_Per_Queue 32
#define Number_Of_Line_Buffer_Queues 32
#define Default_Number_Reads_Per_File 2

global_variable
u64
Total_Memory;

global_variable
u32
Number_of_Threads = Default_Number_of_Threads;

global_variable
u32
Number_Reads_Per_File = Default_Number_Reads_Per_File;

enum
output_mode
{
    file,
    std
};

global_variable
enum output_mode
Output_Mode = file;

global_variable
mutex
Std_Out_Lock;

struct
line_buffer
{
    u08 line[Line_Buffer_Size];
    u32 homeIndex;
    struct line_buffer *prev;
};

struct
single_line_buffer_queue
{
    u32 queueLength;
    u32 pad;
    mutex rwMutex;
    struct line_buffer *front;
    struct line_buffer *rear;
};

struct
line_buffer_queue
{
    struct single_line_buffer_queue **queues;
    threadSig index;
    u32 pad;
};

global_variable
struct line_buffer_queue *
Line_Buffer_Queue;

global_function
void
InitialiseSingleLineBufferQueue(struct single_line_buffer_queue *queue)
{
    InitialiseMutex(queue->rwMutex);
    queue->queueLength = 0;
}

global_function
void
AddSingleLineBufferToQueue(struct single_line_buffer_queue *queue, struct line_buffer *buffer);

global_function
void
InitialiseLineBufferQueue(struct memory_arena *arena, struct line_buffer_queue *queue)
{
    queue->queues = PushArrayP(arena, struct single_line_buffer_queue *, Number_Of_Line_Buffer_Queues);
    queue->index = 0;

    ForLoop(Number_Of_Line_Buffer_Queues)
    {
        queue->queues[index] = PushStructP(arena, struct single_line_buffer_queue);
        InitialiseSingleLineBufferQueue(queue->queues[index]);

        ForLoop2(Number_Of_Line_Buffers_Per_Queue)
        {
            struct line_buffer *buffer = PushStructP(arena, struct line_buffer);
            buffer->homeIndex = index;
            AddSingleLineBufferToQueue(queue->queues[index], buffer);
        }
    }
}

global_function
void
AddSingleLineBufferToQueue(struct single_line_buffer_queue *queue, struct line_buffer *buffer)
{
    LockMutex(queue->rwMutex);
    buffer->prev = 0;

    switch (queue->queueLength)
    {
        case 0:
            queue->front = buffer;
            queue->rear = buffer;
            break;

        default:
            queue->rear->prev = buffer;
            queue->rear = buffer;
    }

    ++queue->queueLength;
    UnlockMutex(queue->rwMutex);
}

global_function
void
AddLineBufferToQueue(struct line_buffer_queue *queue, struct line_buffer *buffer)
{
    struct single_line_buffer_queue *singleQueue = queue->queues[buffer->homeIndex];
    AddSingleLineBufferToQueue(singleQueue, buffer);
}

global_function
struct line_buffer *
TakeSingleLineBufferFromQueue(struct single_line_buffer_queue *queue)
{
    LockMutex(queue->rwMutex);
    struct line_buffer *buffer = queue->front;

    switch (queue->queueLength)
    {
        case 0:
            break;

        case 1:
            queue->front = 0;
            queue->rear = 0;
            queue->queueLength = 0;
            break;

        default:
            queue->front = buffer->prev;
            --queue->queueLength;
    }

    UnlockMutex(queue->rwMutex);

    return(buffer);
}

global_function
struct single_line_buffer_queue *
GetSingleLineBufferQueue(struct line_buffer_queue *queue)
{
    u32 index = __atomic_fetch_add(&queue->index, 1, 0) % Number_Of_Line_Buffer_Queues;
    return(queue->queues[index]);
}

global_function
struct line_buffer *
TakeLineBufferFromQueue(struct line_buffer_queue *queue)
{
    return(TakeSingleLineBufferFromQueue(GetSingleLineBufferQueue(queue)));
}

global_function
struct line_buffer *
TakeLineBufferFromQueue_Wait(struct line_buffer_queue *queue)
{
    struct line_buffer *buffer = 0;
    while (!buffer)
    {
        buffer = TakeLineBufferFromQueue(queue);
    }
    return(buffer);
}

global_variable
mutex
Working_Set_Lock;

#define Default_ProcessingPairMemorySize MegaByte(8)

global_variable
u32
ProcessingPairMemorySize = Default_ProcessingPairMemorySize;

enum
file_type
{
    text,
    gzip
};

#define Gzip_Buffer_Size KiloByte(256)
#define Gzip_Decompress_Memory_Size KiloByte(64)
#define Gzip_Compress_Memory_Size KiloByte(512)

struct
gzip_data
{
    u08 *gzipInBuffer;
    u08 *gzipOutBuffer;
    u32 gzipBufferPtr;
    u32 gzipBufferPtrEnd;
    z_stream *zstream;
    struct memory_arena *zarena;
};

struct
file_stream
{
    FILE *file;
    enum file_type type;
    u08 name[20];
    struct line_buffer *buffer;
    struct gzip_data *gzip;
};

struct
output_file_set
{
    struct file_stream *out;
    mutex lock;
};

struct
file_io_set
{
    struct file_stream *in1;
    struct file_stream *in2;
    threadSig inDoneCount;
    u32 groupsToSkip;
    struct output_file_set *out1;
    struct output_file_set *out2;
};

struct
buffer_save
{
    struct line_buffer *in;
    struct line_buffer *out;
    struct output_file_set *outSet;
};

global_function
void *
ZAllocCallBack(void *arena, u32 items, u32 size)
{
    return (void *)(PushArrayP((struct memory_arena *)arena, u08, items * size));
}

global_function
void
ZFreeCallBack(void *arena, void *ptr)
{
    (void)arena;
    (void)ptr;
}

struct
gzip_file_node
{
    struct file_stream *file;
    struct gzip_file_node *next;
};

global_variable
struct gzip_file_node *
First_GZip_Output_File = 0;

global_variable
struct gzip_file_node *
Current_GZip_Output_File = 0;

global_function
void
InitialiseGzipRead(struct gzip_data *gzip, struct memory_arena *arena)
{
    gzip->gzipOutBuffer = PushArrayP(arena, u08, Gzip_Buffer_Size);
    gzip->gzipInBuffer = PushArrayP(arena, u08, Gzip_Buffer_Size);
    gzip->gzipBufferPtr = Gzip_Buffer_Size;
    gzip->gzipBufferPtrEnd = Gzip_Buffer_Size + 1;
    gzip->zstream = PushStructP(arena, z_stream);
    gzip->zarena = PushSubArenaP(arena, Gzip_Decompress_Memory_Size);
    
    gzip->zstream->zalloc = ZAllocCallBack;
    gzip->zstream->zfree = ZFreeCallBack;
    gzip->zstream->opaque = (void *)gzip->zarena;

    gzip->zstream->avail_in = 0;
    gzip->zstream->next_in = Z_NULL;
    inflateInit2(gzip->zstream, 16+MAX_WBITS);
}

global_function
void
InitialiseGzipWrite(struct gzip_data *gzip, struct memory_arena *arena)
{
    gzip->gzipOutBuffer = PushArrayP(arena, u08, Gzip_Buffer_Size);
    gzip->gzipInBuffer = PushArrayP(arena, u08, Gzip_Buffer_Size);
    gzip->gzipBufferPtr = 0;
    gzip->zstream = PushStructP(arena, z_stream);
    gzip->zarena = PushSubArenaP(arena, Gzip_Compress_Memory_Size);

    gzip->zstream->zalloc = ZAllocCallBack;
    gzip->zstream->zfree = ZFreeCallBack;
    gzip->zstream->opaque = (void *)gzip->zarena;

    deflateInit2 (  gzip->zstream, Z_BEST_COMPRESSION, Z_DEFLATED,
                    31,
                    8,
                    Z_DEFAULT_STRATEGY);
}

global_function
u32
OpenFileForReading(u08 *fileName, struct file_stream *file, struct memory_arena *arena)
{
    u32 result = 0;
    FILE *f = fopen((const char *)fileName, "rb");
    if (f)
    {
        u08 buff[2];
        fread(buff, 2, 1, f);
        fclose(f);
        enum file_type type;
        
        if (buff[0] == 0x1f && buff[1] == 0x8b)
        {
            type = gzip;
            f = fopen((const char *)fileName, "rb");
        }
        else
        {
            type = text;
            f = fopen((const char *)fileName, "r");
        }
        
        result = f ? 1 : 0;
        if (result)
        {
            file->file = f;
            file->type = type;

            if (type == gzip)
            {
                file->gzip = PushStructP(arena, struct gzip_data);
                InitialiseGzipRead(file->gzip, arena);
            }
        }
    }
    
    return(result);
}

global_function
u32
OpenFileForWriting(u08 *fileName, struct file_stream *file, struct memory_arena *arena)
{
    u32 result = 0;
    CopyNullTerminatedString(fileName, file->name); 
    FILE *f = fopen((const char *)fileName, "wb");

    if (f)
    {
        file->file = f;
        file->type = gzip;
        file->gzip = PushStructP(arena, struct gzip_data);
        InitialiseGzipWrite(file->gzip, arena);
        result = 1;
    }

    return(result);
}

global_variable
threadSig
Output_Number = 0;

global_function
void
WriteToGzipFile(struct file_stream *file, u08 *in, u32 nIn);

global_function
void
WriteFile(void *in)
{
    struct line_buffer *bufferIn = (struct line_buffer *)in;
    struct buffer_save *mem = (struct buffer_save *)bufferIn->line;
    struct output_file_set *file = mem->outSet;
    
    if (TryToLockMutex(file->lock) == 0)
    {
        
        WriteToGzipFile(file->out, mem->out->line + 4, *((u32 *)mem->out->line));

        UnlockMutex(file->lock);
        AddLineBufferToQueue(Line_Buffer_Queue, mem->out);
        AddLineBufferToQueue(Line_Buffer_Queue, bufferIn);
    }
    else
    {
        ThreadPoolAddTask(Thread_Pool, WriteFile, in);
    }
}

global_variable
u32
StdOut_Queue_Length = 4;

global_variable
threadSig
Number_Waiting_For_StdOut = 0;

global_function
void
WriteToStdOut(void *in)
{
    u32 waitInQueue = 0;
    u32 nWaiting;
    __atomic_load(&Number_Waiting_For_StdOut, &nWaiting, 0);
    if (nWaiting < StdOut_Queue_Length)
    {
	__atomic_add_fetch(&Number_Waiting_For_StdOut, 1, 0);
	waitInQueue = 1;
	LockMutex(Std_Out_Lock);
    }
    
    if (waitInQueue || TryToLockMutex(Std_Out_Lock) == 0)
    {
        struct line_buffer *bufferIn = (struct line_buffer *)in;
        struct buffer_save *mem = (struct buffer_save *)bufferIn->line;
        struct line_buffer *out1 = mem[0].out;
        struct line_buffer *out2 = mem[1].out;

        u08 buff1[512];
        u08 buff2[512];

        CopyNullTerminatedString(out1->line, buff1);
	CopyNullTerminatedString(out2->line, buff2);

        AddLineBufferToQueue(Line_Buffer_Queue, out1);
        AddLineBufferToQueue(Line_Buffer_Queue, out2);
        
        printf("%s", buff1);
        printf("%s", buff2);
        
	if (waitInQueue)
	{
	    __atomic_sub_fetch(&Number_Waiting_For_StdOut, 1, 0);
	}
	
	UnlockMutex(Std_Out_Lock);
        AddLineBufferToQueue(Line_Buffer_Queue, bufferIn);
    }
    else
    {
        ThreadPoolAddTask(Thread_Pool, WriteToStdOut, in);
    }
}

global_function
void
FinishRead_StdOut(void *in)
{
    struct line_buffer *bufferIn = (struct line_buffer *)in;
    struct buffer_save *mem = (struct buffer_save *)bufferIn->line;
    struct line_buffer *in1 = mem[0].in;
    struct line_buffer *out1 = mem[0].out;
    struct line_buffer *in2 = mem[1].in;
    struct line_buffer *out2 = mem[1].out;

    u08 *readOut = out1->line;
    u08 *readIn = in1->line;

    while (*readOut++ != '\n') {}
    while (*readIn++ != '\n') {}
    readIn += 23;

    u32 nNewLines = 2;
    while (nNewLines)
    {
        if (*readIn == '\n')
        {
            --nNewLines;
        }
        *readOut++ = *readIn++;
    }
    readIn += 23;
    
    nNewLines = 1;
    while (nNewLines)
    {
        if (*readIn == '\n')
        {
            --nNewLines;
        }
        *readOut++ = *readIn++;
    }
    *readOut = '\0';

    AddLineBufferToQueue(Line_Buffer_Queue, in1);
    
    readOut = out2->line;
    readIn = in2->line;

    while (*readOut++ != '\n') {}
    while (*readIn++ != '\n') {}

    nNewLines = 2;
    while (nNewLines)
    {
        if (*readIn == '\n')
        {
            --nNewLines;
        }
        *readOut++ = *readIn++;
    }
    
    nNewLines = 1;
    while (nNewLines)
    {
        if (*readIn == '\n')
        {
            --nNewLines;
        }
        *readOut++ = *readIn++;
    }
    *readOut = '\0';

    AddLineBufferToQueue(Line_Buffer_Queue, in2);
   
    ThreadPoolAddTask(Thread_Pool, WriteToStdOut, in);
}

global_function
void
FinishRead1(void *in)
{
    struct line_buffer *bufferIn = (struct line_buffer *)in;
    struct buffer_save *mem = (struct buffer_save *)bufferIn->line;

    u08 *readOut = mem->out->line + 4;
    u08 *readIn = mem->in->line;
    u32 nBytesOut = 1;

    while (*readOut++ != '\n')
    {
        ++nBytesOut;
    }
    while (*readIn++ != '\n') {}
    ForLoop(23)
    {
        ++readIn;
    }

    u32 nNewLines = 2;
    while (nNewLines)
    {
        if (*readIn == '\n')
        {
            --nNewLines;
        }
        *readOut++ = *readIn++;
        ++nBytesOut;
    }
    ForLoop(23)
    {
        ++readIn;
    }
    nNewLines = 1;
    while (nNewLines)
    {
        if (*readIn == '\n')
        {
            --nNewLines;
        }
        *readOut++ = *readIn++;
        ++nBytesOut;
    }

    AddLineBufferToQueue(Line_Buffer_Queue, mem->in);

    *((u32 *)mem->out->line) = nBytesOut;

    ThreadPoolAddTask(Thread_Pool, WriteFile, in);
}

global_function
void
FinishRead2(void *in)
{
    struct line_buffer *bufferIn = (struct line_buffer *)in;
    struct buffer_save *mem = (struct buffer_save *)bufferIn->line;

    u08 *readOut = mem->out->line + 4;
    u08 *readIn = mem->in->line;
    u32 nBytesOut = 1;

    while (*readOut++ != '\n')
    {
        ++nBytesOut;
    }
    while (*readIn++ != '\n') {}

    u32 nNewLines = 2;
    while (nNewLines)
    {
        if (*readIn == '\n')
        {
            --nNewLines;
        }
        *readOut++ = *readIn++;
        ++nBytesOut;
    }
    
    nNewLines = 1;
    while (nNewLines)
    {
        if (*readIn == '\n')
        {
            --nNewLines;
        }
        *readOut++ = *readIn++;
        ++nBytesOut;
    }

    AddLineBufferToQueue(Line_Buffer_Queue, mem->in);

    *((u32 *)mem->out->line) = nBytesOut;

    ThreadPoolAddTask(Thread_Pool, WriteFile, in);
}

global_function
u32
FileNotEmpty(struct file_stream *file)
{
    u32 result = 0;

    switch (file->type)
    {
        case text:
            {
                result = !feof(file->file);
            } break;
        case gzip:
            {
                result = file->gzip->gzipBufferPtrEnd != file->gzip->gzipBufferPtr;
            } break;
    }

    return(result);
}

global_function
void
ProcessReadPair_Step1(void *in);

global_function
void
ProcessReadPair_Step2(void *in)
{
    struct file_io_set *fileSet = (struct file_io_set *)in;
    fileSet->out1->out->buffer = TakeLineBufferFromQueue_Wait(Line_Buffer_Queue);
    fileSet->out2->out->buffer = TakeLineBufferFromQueue_Wait(Line_Buffer_Queue);

    u08 *out1;
    u08 *out2;
    
    if (Output_Mode == file)
    {
        out1 = fileSet->out1->out->buffer->line + 4;
        out2 = fileSet->out2->out->buffer->line + 4;
    }
    else
    {
        out1 = fileSet->out1->out->buffer->line;
        out2 = fileSet->out2->out->buffer->line;
    }
    
    u08 *read1 = fileSet->in1->buffer->line;

    while (*read1 != 0x20)
    {
        *out1++ = *read1;
        *out2++ = *read1++;
    }
    *out1++ = '_';
    *out2++ = '_';
    while (*read1++ != '\n') {}
    ForLoop(16)
    {
        *out1++ = *read1;
        *out2++ = *read1++;
    }
    *out1 = '\n';
    *out2 = '\n';

    if (Output_Mode == file)
    {
        struct line_buffer *mem1 = TakeLineBufferFromQueue_Wait(Line_Buffer_Queue); 
        struct buffer_save *save = (struct buffer_save *)mem1->line;
        save->outSet = fileSet->out1;
        save->in = fileSet->in1->buffer;
        save->out = fileSet->out1->out->buffer;

        struct line_buffer *mem2 = TakeLineBufferFromQueue_Wait(Line_Buffer_Queue); 
        save = (struct buffer_save *)mem2->line;
        save->outSet = fileSet->out2;
        save->in = fileSet->in2->buffer;
        save->out = fileSet->out2->out->buffer;

        ThreadPoolAddTask(Thread_Pool, FinishRead1, (void *)mem1);
        ThreadPoolAddTask(Thread_Pool, FinishRead2, (void *)mem2);
    }
    else
    {
        struct line_buffer *mem = TakeLineBufferFromQueue_Wait(Line_Buffer_Queue);
        struct buffer_save *save = (struct buffer_save *)mem->line;
        save[0].in = fileSet->in1->buffer;
        save[0].out = fileSet->out1->out->buffer;
        save[1].in = fileSet->in2->buffer;
        save[1].out = fileSet->out2->out->buffer;

        ThreadPoolAddTask(Thread_Pool, FinishRead_StdOut, (void *)mem);
    }

    if (FileNotEmpty(fileSet->in1) && FileNotEmpty(fileSet->in2))
    {
        ThreadPoolAddTask(Thread_Pool, ProcessReadPair_Step1, in);
    }
    else
    {
        fclose(fileSet->in1->file);
        fclose(fileSet->in2->file);
    }
}

global_function
void
WriteToGzipFile(struct file_stream *file, u08 *in, u32 nIn)
{
    while (nIn)
    {
        u32 nSpaceLeft = Gzip_Buffer_Size - file->gzip->gzipBufferPtr;
        u32 nBytesToCopy = Min(nIn, nSpaceLeft);
        ForLoop(nBytesToCopy)
        {
            file->gzip->gzipInBuffer[file->gzip->gzipBufferPtr++] = in[index];
        }
        nIn -= nBytesToCopy;

        if (file->gzip->gzipBufferPtr == Gzip_Buffer_Size)
        {
            file->gzip->gzipBufferPtr = 0;
            
            z_stream *zstream = file->gzip->zstream;
            zstream->next_in = file->gzip->gzipInBuffer;
            zstream->avail_in = Gzip_Buffer_Size;

            FILE *outFile = file->file;

            do
            {
                zstream->avail_out = Gzip_Buffer_Size;
                zstream->next_out = file->gzip->gzipOutBuffer;

                deflate(zstream, Z_NO_FLUSH);
                u32 nToWrite = Gzip_Buffer_Size - zstream->avail_out;

                fwrite(file->gzip->gzipOutBuffer, 1, nToWrite, outFile);

            } while (zstream->avail_out == 0);
        }
    }
}

global_function
void
CloseGzipFile(struct file_stream *file)
{
    z_stream *zstream = file->gzip->zstream;
    zstream->next_in = file->gzip->gzipInBuffer;
    zstream->avail_in = file->gzip->gzipBufferPtr;
    
    do
    {
        zstream->next_out = file->gzip->gzipOutBuffer;
        zstream->avail_out = Gzip_Buffer_Size;
        
        deflate(zstream, Z_FINISH);
        fwrite(file->gzip->gzipOutBuffer, 1, Gzip_Buffer_Size - zstream->avail_out,
                file->file);
    
    } while (zstream->avail_out == 0);
    
    (void)deflateEnd(zstream);
    fclose(file->file);
}

global_function
void
CloseGzipFileJob(void *in)
{
    CloseGzipFile((struct file_stream *)in);
}

global_function
void
DecompressIntoGzipBuffer(struct file_stream *file)
{
    z_stream *zstream = file->gzip->zstream;
    u32 bufferFilled = 0;
    s32 inflateResult = 0;
    u32 noFill = 0;
    u32 eof = 0;

    while (bufferFilled != Gzip_Buffer_Size)
    {
        u32 oldIn = zstream->avail_in;
	u32 oldOut = zstream->avail_out;
	
	if (zstream->avail_in == 0)
        {
            zstream->avail_in = (u32)fread(file->gzip->gzipInBuffer, 1, Gzip_Buffer_Size, file->file);
            zstream->next_in = file->gzip->gzipInBuffer;
        }
        if (zstream->avail_in == 0)
        {
            eof = 1;
	    break;
        }

        zstream->avail_out = Gzip_Buffer_Size - bufferFilled;
        zstream->next_out = file->gzip->gzipOutBuffer + bufferFilled;
        inflateResult = inflate(zstream, Z_NO_FLUSH);
        
	bufferFilled = Gzip_Buffer_Size - zstream->avail_out;

	if (inflateResult == Z_STREAM_ERROR)
	{
	    fprintf(stderr, "Z_STREAM_ERROR with file %s\n", file->name);
	    exit(1);
	}
	else if (inflateResult == Z_STREAM_END)
	{
	    (void)inflateEnd(zstream);
	    ResetMemoryArenaP(file->gzip->zarena);
	    inflateInit2(zstream, 16+MAX_WBITS);
	}

	noFill = oldIn == zstream->avail_in && oldOut == zstream->avail_out;

	if (noFill)
	{
	    break;
	}
    }

    if (noFill || eof)
    {
	file->gzip->gzipBufferPtrEnd = bufferFilled;
    }
}

global_function
s32
GzipGetChar(struct file_stream *file)
{
    u08 *gzipBuffer = file->gzip->gzipOutBuffer;
    u32 gzipBufferPtr = file->gzip->gzipBufferPtr;

    if (gzipBufferPtr == Gzip_Buffer_Size)
    {
        DecompressIntoGzipBuffer(file);
        gzipBufferPtr = 0;
    }

    u32 gzipBufferPtrEnd = file->gzip->gzipBufferPtrEnd;

    s32 character = gzipBufferPtr == gzipBufferPtrEnd ? EOF : (s32)gzipBuffer[gzipBufferPtr++];
    file->gzip->gzipBufferPtr = gzipBufferPtr;

    return(character);
}

global_function
void
SkipReads_Text(struct file_stream *file, u32 groupsToSkip)
{
    u32 nLines = 0;
    u32 nLinesToSkip = 4 * groupsToSkip;

    s32 character;
    while (nLines < nLinesToSkip && (character = fgetc(file->file)) != EOF)
    {
        if (character == 10)
        {
            ++nLines;
        }
    }
}

global_function
void
SkipReads_Gzip(struct file_stream *file, u32 groupsToSkip)
{
    u32 nLines = 0;
    u32 nLinesToSkip = 4 * groupsToSkip;

    s32 character;
    while (nLines < nLinesToSkip && (character = GzipGetChar(file)) != EOF)
    {
        if (character == 10)
        {
            ++nLines;
        }
    }
}


global_function
u32
FillInputBuffer_Text(struct file_stream *file)
{
    u32 bufferPtr = 0;
    u32 nLines = 0;
    s32 character;
    while (nLines < 4 && (character = fgetc(file->file)) != EOF)
    {
        file->buffer->line[bufferPtr++] = (u08)character;

        if (character == 10)
        {
            ++nLines;
        }
    }

    return(bufferPtr);
}

global_function
u32
FillInputBuffer_Gzip(struct file_stream *file)
{
    u32 bufferPtr = 0;
    u32 nLines = 0;
    s32 character;
    while (nLines < 4 && (character = GzipGetChar(file)) != EOF)
    {
        file->buffer->line[bufferPtr++] = (u08)character;

        if (character == 10)
        {
            ++nLines;
        }
    }

    return(bufferPtr);
}

global_function
void
SkipReads(struct file_stream *file, u32 groupsToSkip)
{
    switch (file->type)
    {
        case text:
            {
                SkipReads_Text(file, groupsToSkip);
            } break;

        case gzip:
            {
                SkipReads_Gzip(file, groupsToSkip);
            } break;
    }
}

global_function
u32
FillInputBuffer(struct file_stream *file)
{
    u32 result = 0;
    
    switch (file->type)
    {
        case text:
            {
                result = FillInputBuffer_Text(file);
            } break;

        case gzip:
            {
                result = FillInputBuffer_Gzip(file);
            } break;
    }

    return(result);
}

global_function
void
ReadInputLines1(void *in)
{
    struct file_io_set *fileSet = (struct file_io_set *)in;
    struct file_stream *file = fileSet->in1;
    
    file->buffer = TakeLineBufferFromQueue_Wait(Line_Buffer_Queue);
   
    if (FillInputBuffer(file))
    {
        SkipReads(file, fileSet->groupsToSkip);
        
        if (__atomic_add_fetch(&fileSet->inDoneCount, 1, 0) == 2)
        {
            ThreadPoolAddTask(Thread_Pool, ProcessReadPair_Step2, in);
        }
    }
}

global_function
void
ReadInputLines2(void *in)
{
    struct file_io_set *fileSet = (struct file_io_set *)in;
    struct file_stream *file = fileSet->in2;

    file->buffer = TakeLineBufferFromQueue_Wait(Line_Buffer_Queue);
  
    if (FillInputBuffer(file))
    {
        SkipReads(file, fileSet->groupsToSkip);
        
        if (__atomic_add_fetch(&fileSet->inDoneCount, 1, 0) == 2)
        {
            ThreadPoolAddTask(Thread_Pool, ProcessReadPair_Step2, in);
        }
    }
}

global_function
void
ProcessReadPair_Step1(void *in)
{
    struct file_io_set *fileSet = (struct file_io_set *)in;
    fileSet->inDoneCount = 0;

    ThreadPoolAddTask(Thread_Pool, ReadInputLines1, in);
    ThreadPoolAddTask(Thread_Pool, ReadInputLines2, in);
}

global_function
void
StartProcessingReadPair(void *in)
{
    struct line_buffer *buffer = (struct line_buffer *)in;
    u08 *line = buffer->line;

    u08 *read1 = line;
    u08 *read2 = line;
    while (*read2++ != '\0') {}

    if (read1[0] == 'q' && read1[1] == '1' && read1[2] == '=')
    {
	read1 += 3;
    }
    
    if (read2[0] == 'q' && read2[1] == '2' && read2[2] == '=')
    {
	read2 += 3;
    }
    
    LockMutex(Working_Set_Lock);
    struct memory_arena *arena = PushSubArena(Working_Set, ProcessingPairMemorySize);
    UnlockMutex(Working_Set_Lock);
    
    ForLoop(Number_Reads_Per_File)
    {
        struct file_io_set *fileSet = PushStructP(arena, struct file_io_set); 
        fileSet->in1 = PushStructP(arena, struct file_stream);
        fileSet->in2 = PushStructP(arena, struct file_stream);
        fileSet->out1 = PushStructP(arena, struct output_file_set);
        fileSet->out1->out = PushStructP(arena, struct file_stream);
        InitialiseMutex(fileSet->out1->lock);
        fileSet->out2 = PushStructP(arena, struct output_file_set);
        fileSet->out2->out = PushStructP(arena, struct file_stream);
        InitialiseMutex(fileSet->out2->lock);

        fileSet->groupsToSkip = Number_Reads_Per_File - 1;

        if (Output_Mode == file)
        {
            LockMutex(Working_Set_Lock);
            
            struct gzip_file_node *node1 = PushStruct(Working_Set, struct gzip_file_node);
            node1->file = fileSet->out1->out;
            struct gzip_file_node *node2 = PushStruct(Working_Set, struct gzip_file_node);
            node2->file = fileSet->out2->out;

            node1->next = node2;
            node2->next = 0;

            if (First_GZip_Output_File == 0)
            {
                First_GZip_Output_File = node1;
            }
            else
            {
                Current_GZip_Output_File->next = node1;
            }
            Current_GZip_Output_File = node2;
            
            UnlockMutex(Working_Set_Lock);
        }

        u32 success = 0;
        if (OpenFileForReading(read1, fileSet->in1, arena))
        {
            if (OpenFileForReading(read2, fileSet->in2, arena))
            {
                if (Output_Mode == file)
                {
                    u32 n = __atomic_fetch_add(&Output_Number, 1, 0);
                    u08 buff[32];
                    sprintf((char *)buff, "R_%d_0.fqgz", n);
                    
                    if (OpenFileForWriting(buff, fileSet->out1->out, arena))
                    {
                        sprintf((char *)buff, "R_%d_1.fqgz", n);
                        
                        if (OpenFileForWriting(buff, fileSet->out2->out, arena))
                        {
                            success = 1;
                        }
                    }
                }
                else
                {
                    success = 1;
                }
            }
        }

        if (!success)
        {
            fprintf(stderr, "Could not open read pair: %s :: %s\n", read1, read2);
        }
        else
        {
            SkipReads(fileSet->in1, index);
            SkipReads(fileSet->in2, index);

            ThreadPoolAddTask(Thread_Pool, ProcessReadPair_Step1, (void *)fileSet);
        }
    }

    AddLineBufferToQueue(Line_Buffer_Queue, buffer);
}

global_function
void
ProcessReadFiles(void *in)
{
    u08 *readFileListName = (u08 *)in;
    FILE *readFileList = fopen((const char *)readFileListName, "r");
    
    struct line_buffer *buffer = 0;
    u32 bufferPtr = 0;
    u32 pairFlag = 0;
    s32 character;
    while ((character = fgetc(readFileList)) != EOF)
    {
        while (!buffer)
        {
            buffer = TakeLineBufferFromQueue(Line_Buffer_Queue);
        }
        buffer->line[bufferPtr++] = (u08)character;

        if (character == 10 && bufferPtr > 1)
        {
            buffer->line[bufferPtr-1] = '\0';
            if (pairFlag)
            {
                ThreadPoolAddTask(Thread_Pool, StartProcessingReadPair, (void *)buffer);
                buffer = 0;
                bufferPtr = 0;
                pairFlag = 0;
            }
            else
            {
                pairFlag = 1;
            }
        }
        else if (character == 10 && bufferPtr == 1)
	{
	    bufferPtr = 0;
	}
	else if (bufferPtr == Line_Buffer_Size)
        {
            fprintf(stderr, "Error: line found larger than %d characters.", Line_Buffer_Size);
            break;
        }
    }

    if (bufferPtr)
    {
        if (pairFlag)
        {
            ThreadPoolAddTask(Thread_Pool, StartProcessingReadPair, (void *)buffer);
            pairFlag = 0;
        }
        else
        {
            pairFlag = 1;
        } 
    }

    if (pairFlag)
    {
        fprintf(stderr, "Error: input file: %s does not have an even number of file names\n", readFileListName);
    }

    fclose(readFileList);
}

global_function
void
ConcatenateReads(u08 *outname, u32 skip)
{
    FILE *outFile = fopen((const char *)outname, "wb");
    if (outFile)
    {
        struct gzip_file_node *firstNode = skip ? First_GZip_Output_File->next : First_GZip_Output_File;

        TraverseLinkedList(firstNode, struct gzip_file_node)
        {
            struct file_stream *inFileStream = node->file;
            FILE *inFile = fopen((const char *)inFileStream->name, "rb");

            if (inFile)
            {
                u32 nRead;
                do
                {
                    nRead = (u32)fread(inFileStream->gzip->gzipOutBuffer, 1, Gzip_Buffer_Size, inFile);
                    fwrite(inFileStream->gzip->gzipOutBuffer, 1, nRead, outFile);
                } while (nRead == Gzip_Buffer_Size);

                fclose(inFile);
                remove((const char *)inFileStream->name);
            }
            else
            {
                fprintf(stderr, "Could not open %s for reading!\n", inFileStream->name);
            }

            node = node->next;
            if (!node)
            {
                break;
            }
        }

        fclose(outFile);
    }
    else
    {
        fprintf(stderr, "Could not open %s for writing!\n", outname);
    }
}

global_function
void
ConcatenateRead1(void *in)
{
    ConcatenateReads((u08 *)in, 0);
}

global_function
void
ConcatenateRead2(void *in)
{
    ConcatenateReads((u08 *)in, 1);
}

MainArgs
{
    u32 memSet = 0;
    u32 printHelp = ArgCount == 1;
    u32 index = 1;

    for (   ;
            index < ((u32)ArgCount - 1);
            ++index )
    {
        u08 *arg = (u08 *)ArgBuffer[index];
        
        if (AreNullTerminatedStringsEqual(arg, (u08 *)"-t"))
        {
            Number_of_Threads = NullTerminatedStringToInt((u08 *)ArgBuffer[index + 1]);
            ++index;
        }
        else if (AreNullTerminatedStringsEqual(arg, (u08 *)"-n"))
        {
            Number_Reads_Per_File = NullTerminatedStringToInt((u08 *)ArgBuffer[index + 1]);
            ++index;
        }
        else if (AreNullTerminatedStringsEqual(arg, (u08 *)"-m"))
        {
            memSet = 1;
            
            u08 *value = (u08 *)ArgBuffer[index + 1];
            while (*value != '\0')
            {
                ++value;
            }
            --value;
            if (*value == 'k' || *value == 'K')
            {
                *value = '\0';
                Total_Memory = KiloByte(NullTerminatedStringToInt((u08 *)ArgBuffer[index + 1]));
            }
            else if (*value == 'm' || *value == 'M')
            {
                *value = '\0';
                Total_Memory = MegaByte(NullTerminatedStringToInt((u08 *)ArgBuffer[index + 1]));
            }
            else if (*value == 'g' || *value == 'G')
            {
                *value = '\0';
                Total_Memory = GigaByte(NullTerminatedStringToInt((u08 *)ArgBuffer[index + 1]));
            }
            else
            {
                Total_Memory = NullTerminatedStringToInt((u08 *)ArgBuffer[index + 1]);
            }
            ++index;
        }
        else if (AreNullTerminatedStringsEqual(arg, (u08 *)"-h"))
        {
            printHelp = 1;
            break;
        }
        else
        {
            break;
        }
    }

    if (printHelp)
    {
        printf("\nScaff10x_FilePreProcess {options} input.fofn output1 {output2}\n\ninput.fofn:\ta list of full paths to 10x fastq files (text or gzipped) ordered by read pairs\n\t\t(i.e. read1,read2,read1,read2)\n\noutput1:\tthe name of the concatenated gzipped first read pair fastq reads,\n\t\tif '-' is passed then both read pairs are printed in pair order to stdout\n\t\t(and output2 is not required)\n\noutput2:\tthe name of the concatenated gzipped second read pair fastq reads\n\nOptions:\t-t number of threads {%d}\n\t\t-n number of simultaneous reads per file {%d}\n\t\t-m total heap memory to use (use k/m/g suffixes)\n\t\t   {if not set, memory usage will be determined by number of threads,\n\t\t   reads per file and output type}\n\t\t-h print this help", Default_Number_of_Threads, Default_Number_Reads_Per_File);
        return(0);
    }

    if (((u32)ArgCount - index) < 2)
    {
        fprintf(stderr, "At least two arguments required\n");
        return(1);
    }

    u08 *inputFile = (u08 *)ArgBuffer[index];
    u08 *out1 = (u08 *)ArgBuffer[index + 1];
    u08 *out2 = 0;

    if (*out1 != '-' && ((u32)ArgCount - index) < 3)
    {
        fprintf(stderr, "Two output arguments required if writing to files\n");
        return(1);
    }
    else if (*out1 == '-')
    {
        Output_Mode = std;
        InitialiseMutex(Std_Out_Lock);
    }
    else
    {
        out2 = (u08 *)ArgBuffer[index + 2];
        Output_Mode = file;
    }
    
    if (!memSet)
    {
        Total_Memory = Default_ProcessingPairMemorySize * Number_Reads_Per_File * Number_of_Threads * 16;
    }
    
    if (Output_Mode == std)
    {
	StdOut_Queue_Length = Number_of_Threads >> 1;
    }

    CreateMemoryArena(Working_Set, Total_Memory);
    Thread_Pool = ThreadPoolInit(&Working_Set, Number_of_Threads);
    InitialiseMutex(Working_Set_Lock);

    Line_Buffer_Queue = PushStruct(Working_Set, struct line_buffer_queue);
    InitialiseLineBufferQueue(&Working_Set, Line_Buffer_Queue);

    ThreadPoolAddTask(Thread_Pool, ProcessReadFiles, (void *)inputFile);
    
    if (Output_Mode == file)
    {
        ThreadPoolWait(Thread_Pool);

        TraverseLinkedList(First_GZip_Output_File, struct gzip_file_node)
        {
            ThreadPoolAddTask(Thread_Pool, CloseGzipFileJob, (void *)node->file);
        }
        ThreadPoolWait(Thread_Pool);
        
        ThreadPoolAddTask(Thread_Pool, ConcatenateRead1, (void *)out1);
        ThreadPoolAddTask(Thread_Pool, ConcatenateRead2, (void *)out2);
    }
    
    ThreadPoolWait(Thread_Pool);
    ThreadPoolDestroy(Thread_Pool);
    EndMain;
}
