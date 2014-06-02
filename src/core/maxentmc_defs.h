/** This file is part of MaxEntMC, a maximum entropy algorithm with moment constraints. **/
/** Copyright (C) 2014 Rafail V. Abramov.                                               **/
/**                                                                                     **/
/** This program is free software: you can redistribute it and/or modify it under the   **/
/** terms of the GNU General Public License as published by the Free Software           **/
/** Foundation, either version 3 of the License, or (at your option) any later version. **/
/**                                                                                     **/
/** This program is distributed in the hope that it will be useful, but WITHOUT ANY     **/
/** WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A     **/
/** PARTICULAR PURPOSE.  See the GNU General Public License for more details.           **/
/**                                                                                     **/
/** You should have received a copy of the GNU General Public License along with this   **/
/** program.  If not, see <http://www.gnu.org/licenses/>.                               **/

#ifndef MAXENTMC_DEFS_H_INCLUDED
#define MAXENTMC_DEFS_H_INCLUDED

#include <stdlib.h>
#include <errno.h>

#ifndef MAXENTMC_CACHE_LINE_SIZE
#warning MAXENTMC_CACHE_LINE_SIZE not defined, assuming 64 bytes
#define MAXENTMC_CACHE_LINE_SIZE 64
#endif

#ifdef __AVX__
#define MAXENTMC_FLOAT_ALIGNMENT 32
#else
#define MAXENTMC_FLOAT_ALIGNMENT 16
#endif

#define MAXENTMC_ALIGNED_SIZE(boundary,size)  ((size)+((boundary)-((size)%(boundary)))%(boundary))

#define MAXENTMC_INCREMENT_POINTER(pt,offset) ((void *)(((char *)(pt))+(offset)))

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

#endif // MAXENTMC_DEFS_H_INCLUDED
