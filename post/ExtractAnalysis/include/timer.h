/**
 * \file    timer.h
 * \ingroup dg4est_group
 * \author  akirby
 * Created on August 16, 2019, 10:11 PM
 */

#ifndef TIMER_H
#define TIMER_H

/* header files */
#ifndef DOXYGEN_IGNORE
#  include <mpi.h>
#endif

#define TIMER(x,time) \
    Real t1_loc = MPI_Wtime();  \
        (x);                    \
    Real t2_loc = MPI_Wtime();  \
    (time) += (t2_loc - t1_loc);

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __cplusplus
}
#endif
#endif /* TIMER_H */
