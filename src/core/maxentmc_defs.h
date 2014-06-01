#ifndef MAXENTMC_DEFS_H_INCLUDED
#define MAXENTMC_DEFS_H_INCLUDED

typedef unsigned char maxentmc_index_t;

#ifdef MAXENTMC_SINGLE_PRECISION
typedef float maxentmc_float_t;
#else
typedef double maxentmc_float_t;
#endif

#include <stdlib.h>
#include <errno.h>

#ifndef MAXENTMC_CACHE_LINE_SIZE
#warning MAXENTMC_CACHE_LINE_SIZE not defined, assuming 64 bytes
#define MAXENTMC_CACHE_LINE_SIZE 64
#endif

#define MAXENTMC_FREAD(_pt_,_size_,_stream_)                \
    if(fread((_pt_),1,(_size_),(_stream_)) != (_size_)){    \
        MAXENTMC_MESSAGE(stderr,"stream read error");       \
        return -1;                                          \
    }

#define MAXENTMC_FREAD_PT(_pt_,_size_,_stream_)             \
    if(fread((_pt_),1,(_size_),(_stream_)) != (_size_)){    \
        MAXENTMC_MESSAGE(stderr,"stream read error");       \
        return NULL;                                        \
    }

#define MAXENTMC_FWRITE(_pt_,_size_,_stream_)               \
    if(fwrite((_pt_),1,(_size_),(_stream_)) != (_size_)){   \
        MAXENTMC_MESSAGE(stderr,"stream write error");      \
        return -1;                                          \
    }

#define MAXENTMC_FWRITE_PT(_pt_,_size_,_stream_)            \
    if(fwrite((_pt_),1,(_size_),(_stream_)) != (_size_)){   \
        MAXENTMC_MESSAGE(stderr,"stream write error");      \
        return NULL;                                        \
    }

#define MAXENTMC_MESSAGE(stream,message)                \
{                                                       \
    fprintf((stream)," %s: " message "\n",__func__);    \
    fflush((stream));                                   \
}

#define MAXENTMC_MESSAGE_VARARG(stream,message,...)                 \
{                                                                   \
    fprintf((stream)," %s: " message "\n",__func__,__VA_ARGS__);    \
    fflush((stream));                                               \
}

#define MAXENTMC_CHECK_NULL(pointer)                        \
  if((pointer)==NULL){                                      \
    MAXENTMC_MESSAGE(stderr,"error: data pointer is NULL"); \
    return -1;                                              \
  }

#define MAXENTMC_CHECK_NULL_PT(pointer)                     \
  if((pointer)==NULL){                                      \
    MAXENTMC_MESSAGE(stderr,"error: data pointer is NULL"); \
    return NULL;                                            \
  }

#define MAXENTMC_ALLOC(pt,size,status)                                              \
(status) = posix_memalign((void **)(&(pt)),MAXENTMC_CACHE_LINE_SIZE,(size));        \
switch((status)){                                                                   \
    case 0:                                                                         \
        break;                                                                      \
    case ENOMEM:                                                                    \
        MAXENTMC_MESSAGE(stderr,"error: insufficient memory");                      \
        break;                                                                      \
    case EINVAL:                                                                    \
        MAXENTMC_MESSAGE(stderr,"error: incorrect boundary in posix_memalign");     \
        break;                                                                      \
    default:                                                                        \
        MAXENTMC_MESSAGE(stderr,"error: unknown error in posix_memalign (?!?!?!)"); \
}

#define MAXENTMC_INCREMENT_POINTER(pt,offset) ((void *)(((char *)(pt))+(offset)))

#ifdef __AVX__
#define MAXENTMC_FLOAT_ALIGNMENT 32
#else
#define MAXENTMC_FLOAT_ALIGNMENT 16
#endif

#define MAXENTMC_ALIGNED_SIZE(boundary,size)  ((size)+((boundary)-((size)%(boundary)))%(boundary))

#endif // MAXENTMC_DEFS_H_INCLUDED
