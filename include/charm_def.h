//
// Created by zhrv on 05.01.2020.
//

#ifndef CHARM_3D_CHARM_DEF_H
#define CHARM_3D_CHARM_DEF_H

#include "charm_version.h"


#define CHARM_DIM           P4EST_DIM
#define CHARM_ALLOC         P4EST_ALLOC
#define CHARM_FREE          P4EST_FREE
#define CHARM_HALF          P4EST_HALF
#define CHARM_CHILDREN      P4EST_CHILDREN
#define CHARM_FACES         P4EST_FACES
#define CHARM_CONNECT_FULL  P4EST_CONNECT_FULL
#define CHARM_CONNECT_FACE  P4EST_CONNECT_FACE
#define CHARM_REALLOC       P4EST_REALLOC
#define CHARM_QUADRANT_LEN  P4EST_QUADRANT_LEN

#define CHARM_CONFIG_YAML

#ifdef POGGI
#warning "POGGI!!!!!!!!!!!"
#endif


#ifdef CHARM_DEBUG

#define CHARM_LOG_LEVEL SC_LP_ESSENTIAL
#define DBG_CH(R) {printf("Rank: %d. File: %s. Line: %d\n", (R), __FILE__, __LINE__);fflush(stdout);}
#define CHARM_ASSERT P4EST_ASSERT

#else

#define CHARM_LOG_LEVEL SC_LP_ESSENTIAL
#define DBG_CH(R) ((void)0)
#define CHARM_ASSERT(R) ((void)0)

#endif



/* log helper macros */
#define CHARM_GLOBAL_LOG(p,s)                           \
  SC_GEN_LOG (charm_package_id, SC_LC_GLOBAL, (p), (s))
#define CHARM_LOG(p,s)                                  \
  SC_GEN_LOG (charm_package_id, SC_LC_NORMAL, (p), (s))
void                CHARM_GLOBAL_LOGF (int priority, const char *fmt, ...)
__attribute__ ((format (printf, 2, 3)));
void                CHARM_LOGF (int priority, const char *fmt, ...)
__attribute__ ((format (printf, 2, 3)));
//#ifndef __cplusplus
#define CHARM_GLOBAL_LOGF(p,f,...)                                      \
  SC_GEN_LOGF (charm_package_id, SC_LC_GLOBAL, (p), (f), __VA_ARGS__)
#define CHARM_LOGF(p,f,...)                                             \
  SC_GEN_LOGF (charm_package_id, SC_LC_NORMAL, (p), (f), __VA_ARGS__)
//#endif

/* convenience global log macros will only print if identifier <= 0 */
#define CHARM_GLOBAL_TRACE(s) CHARM_GLOBAL_LOG (SC_LP_TRACE, (s))
#define CHARM_GLOBAL_LDEBUG(s) CHARM_GLOBAL_LOG (SC_LP_DEBUG, (s))
#define CHARM_GLOBAL_VERBOSE(s) CHARM_GLOBAL_LOG (SC_LP_VERBOSE, (s))
#define CHARM_GLOBAL_INFO(s) CHARM_GLOBAL_LOG (SC_LP_INFO, (s))
#define CHARM_GLOBAL_STATISTICS(s) CHARM_GLOBAL_LOG (SC_LP_STATISTICS, (s))
#define CHARM_GLOBAL_PRODUCTION(s) CHARM_GLOBAL_LOG (SC_LP_PRODUCTION, (s))
#define CHARM_GLOBAL_ESSENTIAL(s) CHARM_GLOBAL_LOG (SC_LP_ESSENTIAL, (s))
#define CHARM_GLOBAL_LERROR(s) CHARM_GLOBAL_LOG (SC_LP_ERROR, (s))
void                CHARM_GLOBAL_TRACEF (const char *fmt, ...)
__attribute__ ((format (printf, 1, 2)));
void                CHARM_GLOBAL_LDEBUGF (const char *fmt, ...)
__attribute__ ((format (printf, 1, 2)));
void                CHARM_GLOBAL_VERBOSEF (const char *fmt, ...)
__attribute__ ((format (printf, 1, 2)));
void                CHARM_GLOBAL_INFOF (const char *fmt, ...)
__attribute__ ((format (printf, 1, 2)));
void                CHARM_GLOBAL_STATISTICSF (const char *fmt, ...)
__attribute__ ((format (printf, 1, 2)));
void                CHARM_GLOBAL_PRODUCTIONF (const char *fmt, ...)
__attribute__ ((format (printf, 1, 2)));
void                CHARM_GLOBAL_ESSENTIALF (const char *fmt, ...)
__attribute__ ((format (printf, 1, 2)));
void                CHARM_GLOBAL_LERRORF (const char *fmt, ...)
__attribute__ ((format (printf, 1, 2)));
//#ifndef __cplusplus
#define CHARM_GLOBAL_TRACEF(f,...)                      \
  CHARM_GLOBAL_LOGF (SC_LP_TRACE, (f), __VA_ARGS__)
#define CHARM_GLOBAL_LDEBUGF(f,...)                     \
  CHARM_GLOBAL_LOGF (SC_LP_DEBUG, (f), __VA_ARGS__)
#define CHARM_GLOBAL_VERBOSEF(f,...)                    \
  CHARM_GLOBAL_LOGF (SC_LP_VERBOSE, (f), __VA_ARGS__)
#define CHARM_GLOBAL_INFOF(f,...)                       \
  CHARM_GLOBAL_LOGF (SC_LP_INFO, (f), __VA_ARGS__)
#define CHARM_GLOBAL_STATISTICSF(f,...)                         \
  CHARM_GLOBAL_LOGF (SC_LP_STATISTICS, (f), __VA_ARGS__)
#define CHARM_GLOBAL_PRODUCTIONF(f,...)                         \
  CHARM_GLOBAL_LOGF (SC_LP_PRODUCTION, (f), __VA_ARGS__)
#define CHARM_GLOBAL_ESSENTIALF(f,...)                          \
  CHARM_GLOBAL_LOGF (SC_LP_ESSENTIAL, (f), __VA_ARGS__)
#define CHARM_GLOBAL_LERRORF(f,...)                     \
  CHARM_GLOBAL_LOGF (SC_LP_ERROR, (f), __VA_ARGS__)
//#endif
#define CHARM_GLOBAL_NOTICE     CHARM_GLOBAL_STATISTICS
#define CHARM_GLOBAL_NOTICEF    CHARM_GLOBAL_STATISTICSF

/* convenience log macros that are active on every processor */
#define CHARM_TRACE(s) CHARM_LOG (SC_LP_TRACE, (s))
#define CHARM_LDEBUG(s) CHARM_LOG (SC_LP_DEBUG, (s))
#define CHARM_VERBOSE(s) CHARM_LOG (SC_LP_VERBOSE, (s))
#define CHARM_INFO(s) CHARM_LOG (SC_LP_INFO, (s))
#define CHARM_STATISTICS(s) CHARM_LOG (SC_LP_STATISTICS, (s))
#define CHARM_PRODUCTION(s) CHARM_LOG (SC_LP_PRODUCTION, (s))
#define CHARM_ESSENTIAL(s) CHARM_LOG (SC_LP_ESSENTIAL, (s))
#define CHARM_LERROR(s) CHARM_LOG (SC_LP_ERROR, (s))
void                CHARM_TRACEF (const char *fmt, ...)
__attribute__ ((format (printf, 1, 2)));
void                CHARM_LDEBUGF (const char *fmt, ...)
__attribute__ ((format (printf, 1, 2)));
void                CHARM_VERBOSEF (const char *fmt, ...)
__attribute__ ((format (printf, 1, 2)));
void                CHARM_INFOF (const char *fmt, ...)
__attribute__ ((format (printf, 1, 2)));
void                CHARM_STATISTICSF (const char *fmt, ...)
__attribute__ ((format (printf, 1, 2)));
void                CHARM_PRODUCTIONF (const char *fmt, ...)
__attribute__ ((format (printf, 1, 2)));
void                CHARM_ESSENTIALF (const char *fmt, ...)
__attribute__ ((format (printf, 1, 2)));
void                CHARM_LERRORF (const char *fmt, ...)
__attribute__ ((format (printf, 1, 2)));
//#ifndef __cplusplus
#define CHARM_TRACEF(f,...)                     \
  CHARM_LOGF (SC_LP_TRACE, (f), __VA_ARGS__)
#define CHARM_LDEBUGF(f,...)                    \
  CHARM_LOGF (SC_LP_DEBUG, (f), __VA_ARGS__)
#define CHARM_VERBOSEF(f,...)                   \
  CHARM_LOGF (SC_LP_VERBOSE, (f), __VA_ARGS__)
#define CHARM_INFOF(f,...)                      \
  CHARM_LOGF (SC_LP_INFO, (f), __VA_ARGS__)
#define CHARM_STATISTICSF(f,...)                        \
  CHARM_LOGF (SC_LP_STATISTICS, (f), __VA_ARGS__)
#define CHARM_PRODUCTIONF(f,...)                        \
  CHARM_LOGF (SC_LP_PRODUCTION, (f), __VA_ARGS__)
#define CHARM_ESSENTIALF(f,...)                         \
  CHARM_LOGF (SC_LP_ESSENTIAL, (f), __VA_ARGS__)
#define CHARM_LERRORF(f,...)                    \
  CHARM_LOGF (SC_LP_ERROR, (f), __VA_ARGS__)
//#endif
#define CHARM_NOTICE            CHARM_STATISTICS
#define CHARM_NOTICEF           CHARM_STATISTICSF

#define CHARM_STRING "charm_3d"

#define CHARM_RIM_NEWTON_STEPS 5000
#define CHARM_RIM_EPS 1.e-5

#define CHARM_EPS 1.e-12

#define CHARM_NO_NORM_ZEROS

#ifndef CHARM_NO_NORM_ZEROS
#define _NORM_(X) ( (fabs(X) <= CHARM_EPS) ? 0. : (X) )
#else
#define _NORM_(X) ( (X) )
#endif

#define _MAX_(X,Y) ((X)>(Y) ? (X) : (Y))
#define _MIN_(X,Y) ((X)<(Y) ? (X) : (Y))
#define _SQR_(X) ((X)*(X))
#define _MAG_(X,Y,Z) (_SQR_(X)+_SQR_(Y)+_SQR_(Z))


#define CHARM_FACE_TYPE_INNER 0
#define CHARM_BND_MAX 128

#define CHARM_BASE_FN_COUNT 4
#define CHARM_FACE_GP_COUNT 6
#define CHARM_QUAD_GP_COUNT 8
#define CHARM_UNIQUE_COMPONENT_COUNT 7 ////i_max

#define CHARM_MAX_COMPONETS_COUNT 128

#define CHARM_ARR_SET_ZERO(A) {int i; for (i = 0; i < CHARM_BASE_FN_COUNT; i++) A[i] = 0.; }


#endif //CHARM_3D_CHARM_DEF_H
