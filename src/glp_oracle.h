/** glp_oracle.h  -- facet separation oracle **/

/***********************************************************************
 * This code is part of MAXE, a helper program for maximum entropy method.
 *
 * Copyright (C) 2025 Laszlo Csirmaz, https://github.com/lcsirmaz/MAXE
 *
 * This program is free, open-source software. You may redistribute it
 * and/or modify under the terms of the GNU General Public License (GPL).
 *
 * There is ABSOLUTELY NO WARRANTY, use at your own risk.
 ***********************************************************************/
         
/**********************************************************************
* The *facet separation oracle* has a (hidden) polytope. A question is
*  a point 'q' (presumably outside of the polytope, possibly ideal). The
*  answer is a supporting hyperplane 'f' (presumably a facet) of the
*  polytope which separates 'q' and the polytope. This module realizes
*  this oracle using an internal point which is connected to 'q'. Where
*  the segment intersects the boundary, the dual solution provides the
*  coordinates of the supporting hyperplane.
*
* OracleData: struct containing the question to be asked from the
*    oracle, and the answer it returns. The structure is allocated
*    by load_vlp(), and never released.
*
* int load_vlp(void)
*  Read the polytope description from the vlp file; store problem 
*    sizes in PARAMS(); create the glpk LP instance.
*  Return value:
*    0: file read, memory allocated, glpk LP object initialized
*    1: some error (syntax error, out of bound values, no memory);
*      errors are reported as R_fatal. The vlp file might be left
*      open. The program should abort, no way to recover.
*/


typedef struct {
    double *overtex;		/* overtex[0..objs] the request */
    double *ofacet;		/* ofacet[0..objs] the response */
} Oracle_t;

extern Oracle_t OracleData;
int load_vlp(void);

/**********************************************************************
* Oracle manipulation
*
* int initialize_oracle()
*  Check the consistency of the loaded vlp problem; add the ideal points 
*  if they were missing. Return values:
*    ORACLE_OK     success
*    ORACLE_EMPTY  there is no feasible solution
*    ORACLE_FAIL   other error condition (LP solver failed for some reason)
*
* int get initial_vertex()
*  Create an outside point to OracleData.overtex for the first approximation
*    ORACLE_OK     success, point is in OracleData.overtex
*    ORACLE_UNBND  the problem is unbounded from below
*    ORACLE_FAIL   other error condition (LP solver failed for some reason)
*
* int ask_oracle()
*  The question and the answer are provided in OracleData.
*    ORACLE_OK    separating facet is stored in OracleData.ofacet
*                 coordinates are rounded to the nearest rational value
*                 when "RoundFacets" is set.
*    ORACLE_UNBND the vertex is inside or at the boundary
*    ORACLE_FAIL  the LP solver failed to solve the problem
*
* Errors are reported in R_fatal
*/
#define ORACLE_OK	0
#define ORACLE_UNBND	1	/* the projection is unbounded */
#define ORACLE_EMPTY	2	/* the polytope is empty or x is not inside */
#define ORACLE_LIMIT	3	/* iteration or time limit reached */
#define ORACLE_FAIL	4	/* the oracle failed */

int initialize_oracle(void);
int ask_oracle(void);

/**********************************************************************
* Get oracle statistics
*
* void get_oracle_stat(int *no, int *it, unsigned long *t)
*    the number of LP instances, total number of iterations (rounds),
*    total time spent by glpk in 0.01 seconds, and version string
*/
void get_oracle_stat(int *no, int *it, unsigned long *t, const char **ver);

/* EOF */

