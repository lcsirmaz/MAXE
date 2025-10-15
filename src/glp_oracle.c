/** glp_oracle.c  --  facet separation oracle **/

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
         
/*
* This part of OUTER realizes the facet separation oracle using glpk
* (Gnu Linear Program Kit).
* OracleData structure is used to communicate with the calling routine.
*/

#include "glp_oracle.h"
#include "report.h"
#include "params.h"

Oracle_t OracleData;

/**********************************************************************
* The problem dimensions, question and answer space
*   int vcols   number of columns in the constraint matrix
*   int vrows   number of rows in the constraint matrix
*   int vobjs   dimension of the objective space: number of objectives
*   double vvertex[vobjs+1]
*               the question direction
*   double vfacet[vobjs+1]
*               the answer vertex which minimizes the facet direction
*/

#define vcols	PARAMS(ProblemColumns)
#define vrows	PARAMS(ProblemRows)
#define vobjs	PARAMS(ProblemObjects)
#define vvertex OracleData.overtex
#define vfacet	OracleData.ofacet

/*==================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h> 	// sqrt
#include "glpk.h"

/*********************************************************************
* glpk routines used
*    P->valid:  the base is valid
*  P=glp_create_prob()    - create empty glpk problem P
*  glp_add_cols(P,ncols)  - add additional "ncols"
*  glp_add_rows(P,nrows)  - add additional "nrows"
*  glp_set_col_bnds(P,idx,type,lwb,upb) 
*  glp_set_row_bnds(P,idx,type,lwb,upb)
*  glp_set_mat_col(P,idx,n,idxarr[1..n],valarr[1..n])
*                          - set or replace the matrix column "idx"
*  glp_set_obj_coef(P,idx,v) 
*  glp_set_obj_dir(P,GLP_MIN/GLP_MAX)
*  glp_sort_matrix(P)      - rebuild row and column linked lists (ascending order)
*  glp_term_out(GLP_OFF/GLP_ON)   - no messages
*  glp_scale_prob(P,GLP_SF_AUTO)  - scale rows and columd
*  glp_adv_basis(P,0)
*  glp_init_smcp(&parm)    - get default values of parameters
*  glp_simplex(P,&parm)    - solve the LP problem
*  glp_get_status(P)       - get the status 
*  glp_get_obj_val(P)      - retrieve  object value
*  glp_get_row_dual(P,idx) - retrieve dual solution
*  glp_get_it_cnt(P)       - retrieve iteration count
**********************************************************************
* The glpk LP object and parameter structure
*   glp_prob *P    the glpk problem object
*   glp_smcp parm  the glpk parameter structure
*/

static glp_prob *P=NULL;
static glp_smcp parm; /* glp parameter structure */

/**********************************************************************
* Read the constraint matrix and objectives from a vlp file
*
* int load_vlp()
*   open and read the problem from the vlp file. Return value:
*     0: file read, memory allocated, dimensions stored,
*        glpk LP object P is initialized and ready to use.
*     1: some error; errors are reported as R_fatal. The vlp file
*        is not necessarily closed; the program should abort.
*
*/

/*---------------------------------------------------------------------
* void perm_array(int len, int array[1..len])
*   make a random permutation of the array
*/
static inline int mrandom(int v)
{ return v<=1?0 : (random()&0x3fffffff)%v; }
static inline void perm_array(int len, int arr[/* 1:len */])
{int i,j,t;
    for(i=1;i<len;i++){
        j=i+mrandom(len+1-i);
        t=arr[i];arr[i]=arr[j];arr[j]=t;
    }
}
/*---------------------------------------------------------------------
* Storage for the constraint matrix, the objectives and the shuffle
*   arrays. The matrix and the shuffle arrays will be freed after
*   loading the problem into 'P'.
*
* double M(row,col)             temporary storage for the constraint
*                               matrix; indices go from 1
* int vlp_rowidx[1:rows+objs]   shuffling rows, temporary
* int vlp_colidx[1:cols+1]      shuffling columns, temporary
* int vlp_objidx[1:objs]        objective indices in rows
* int lambdaidx                 lambda column index
* double vlp_lambda[1:objs]     obj coeffs in the lambda column
* double vlp_init[1:objs]       initial internal point
*
* int allocate_vlp(row,cols,objs)
*   allocate memory for vfacet, vvertex, and temporal storage.
*   Return value: -1: out of memory (should not happen)
*
*    x (cols)        lambda        RHS
*
*   AAAAAAAAAAA         0         ? <c>
*   AAAAAAAAAAA         0         ? <c>
*
*   PPPPPPPPPPP        -w         = <e>
*   PPPPPPPPPPP        -w         = <e>
*  -------------------------------------
*   00000000000         1         maximize
*/

static double *vlp_M;		/* temporary storage for M */
static int *vlp_rowidx;		/* temporary row permutation */
static int *vlp_colidx;		/* temporary column permutation */
static int *vlp_objidx;		/* object indices */
static double *vlp_lambda;	/* lambda column */
static double *vlp_init;	/* internal point from x lines */
static int lambda_idx;		/* lambda column index */

/* M has size row+objs times cols+1, 1<=r<=row+objs, 1<=c<=cols+1 */
#define M(r,c)		vlp_M[(r)+((c)-1)*(vrows+vobjs)-1]

#define xalloc(ptr,type,size)	\
    (ptr=(type*)calloc(size,sizeof(type)))==NULL

static int allocate_vlp(int rows, int cols, int objs)
{ int i;
    vrows=rows; vcols=cols; vobjs=objs; // store these values in PARAMS()
    if(xalloc(vfacet,double,objs+2) ||
       xalloc(vvertex,double,objs+2) ||
       xalloc(vlp_objidx,int,objs+1) ||
       xalloc(vlp_lambda,double,objs+1) ||
       xalloc(vlp_init,double,objs+1) ||
       xalloc(vlp_M,double,(rows+objs)*(cols+1)) ||
       xalloc(vlp_rowidx,int,rows+objs+1) ||
       xalloc(vlp_colidx,int,cols+2))
         return -1; // out of memory
    // indirect indices to rows and columns
    for(i=0;i<=rows+objs;i++) vlp_rowidx[i]=i;
    for(i=0;i<=cols+1;i++) vlp_colidx[i]=i;
    if(PARAMS(ShuffleMatrix)){ // permute randomly
        perm_array(rows+objs,vlp_rowidx);
        perm_array(cols+1,vlp_colidx);
    }
    // objective rows and the lambda column
    for(i=1;i<=objs;i++) vlp_objidx[i]=vlp_rowidx[i+rows];
    lambda_idx=vlp_colidx[cols+1];
    return 0;
}

/*---------------------------------------------------------------------
* Character input from the vlp file
* int MAX_LINELEN = 80
*   maximum line length expected in a vlp file
* char inpline[MAX_LINELEN]
*   the next line read form the vlp file
*
* int nextline(FILE *f)
*   read the next line to inpline[] from the open stream 'f'. Ignore
*       leading spaces and merge spaces. Convert A-Z to a-z.
*   Return value:
*        1: next line is in inpline[]
*        0: EOF
*
* int vlp_type_ok(char ctrl,int parno)
*   checks the ctrl char against the supplied number of args:
*    'f' (free) parno==0
*    'u','l','s' (upper,lower, fixed) parno==1
*    'd' (both) parno==2
*
* int glp_type(char ctrl)
*   returns the glpk version of the ctrl char, or -1 if error.
*
*/

/* read a single line from a VLP file */
#define MAX_LINELEN	80 /* maximal length of a file */

static char inpline[MAX_LINELEN+1]; /* contains the next vlp line */

/* read the next vlp line to inpline */
static int nextline(FILE *f)
{int i,sp,ch;
    i=0;sp=0; memset(inpline,0,MAX_LINELEN+1);
    while((ch=getc(f))>=0){
        if(ch=='\n'){ if(i==0){sp=0; continue;} return 1; }
        if(ch==' '||ch=='\t'){sp=1; continue;}
        if(ch<=0x20 || ch>126) continue; /* ignore these characters */
        if(sp && i>0){if (i<MAX_LINELEN){inpline[i]=' '; i++;} }
        if('A'<=ch && ch<='Z') ch +='a'-'A'; /* upper case => lower case */
        sp=0; if(i<MAX_LINELEN){inpline[i]=ch; i++; }
    }
    /** EOF **/
    return i>0?1:0;
}

/* direction character and the number of following real numbers */
static int vlp_type_ok(char ctrl, int parno)
{   switch(ctrl){
  case 'f': return parno==2;
  case 'u': case 'l': case 's': return parno==3;
  case 'd': return parno==4;
    }
    return 0;
}

/* glpk version of directions */
static int glp_type(char ctrl)
{   switch(ctrl){
  case 'f': return GLP_FR;
  case 'u': return GLP_UP;
  case 'l': return GLP_LO;
  case 's': return GLP_FX;
  case 'd': return GLP_DB;
    }
    return -1;
}

/* read a vlp problem from a file as an LP instance
 * set lambda >=0  */
int load_vlp(void)
{FILE *f; int rows,cols,objs; int i,j,cnt; double p,b1,b2; char ctrl;
 double dir=1.0;
    f=fopen(PARAMS(VlpFile),"r");
    if(!f){
        report(R_fatal,"Cannot open vlp file %s for reading\n",PARAMS(VlpFile));
        return 1;
    }
    rows=cols=objs=0;
    P = glp_create_prob();
    while(nextline(f)) switch(inpline[0]){ /* read next input line */
       case 'c': // comment line, print out nonempty comment lines before the first p line
                 if(rows==0 && inpline[1]){
                     report(R_warn,"C%s\n",inpline+1);
                 }
                 continue;
       case 'e': break; // end
       case 'p': if(rows>0){
                    report(R_fatal,"read_vlp: second p line in %s:\n   %s\n",
                              PARAMS(VlpFile),inpline); return 1; 
                 }
                 dir=+1.0; PARAMS(Direction)=0;
                                    /**   <rows> <cols> <> <objs> <> **/
                 cnt=sscanf(inpline,"p vlp min %d %d %*d %d %*d",&rows,&cols,&objs);
                 if(cnt==0){
                    dir=-1.0; PARAMS(Direction)=1;
                    cnt=sscanf(inpline,"p vlp max %d %d %*d %d %*d",&rows,&cols,&objs);
                 }
                 if(cnt!=3 || rows<=1 || cols<=1 || objs<1 ){
                    report(R_fatal,"read_vlp: wrong p line in %s\n   %s\n",
                               PARAMS(VlpFile),inpline); return 1;
                 }
                 if(allocate_vlp(rows,cols,objs)){
                     report(R_fatal,"read_vlp: out of memory for %s:\n   %s\n",
                               PARAMS(VlpFile),inpline); return 1;
                 }
		 glp_add_cols(P,cols+1); glp_add_rows(P,rows+objs);
                 glp_set_col_bnds(P,lambda_idx,GLP_LO,0.0,0.0); // lambda is >=0
                 continue;
       case 'j': // j <col> [ f || l <val> | u <val> | d <val1> <val2> | s <val> ]
                 if(rows==0){
                    report(R_fatal,"read_vlp: j line before p in %s\n  %s\n",
                               PARAMS(VlpFile),inpline); return  1;
                 }
                 b1=b2=0.0;
                 cnt=sscanf(inpline,"j %d %c %lg %lg",&j,&ctrl,&b1,&b2);
                 if(cnt<2 || cols<j || j<1 || !vlp_type_ok(ctrl,cnt)){
                    report(R_fatal,"read_vlp: wrong j line in %s\n   %s\n",
                                PARAMS(VlpFile),inpline); return 1;
                 }
                 if(cnt<4) b2=b1; // GLP_UP uses b2 as the value
                 glp_set_col_bnds(P,vlp_colidx[j],glp_type(ctrl),b1,b2);
                 continue;
       case 'i': if(rows==0){
                    report(R_fatal,"read_vlp: i line before p in %s\n   %s\n",
                                PARAMS(VlpFile),inpline); return 1;
                 }
                 // i <row> [ f | l <val> | u <val> | d <val1> <val2> | s <val> ]
                 b1=b2=0.0;
                 cnt=sscanf(inpline,"i %d %c %lg %lg",&i,&ctrl,&b1,&b2);
                 if(cnt<2 || rows<i || i<1 || !vlp_type_ok(ctrl,cnt)){
                    report(R_fatal,"read_vlp: wrong i line in %s\n   %s\n",
                                PARAMS(VlpFile),inpline); return 1;
                 }
                 if(cnt<4) b2=b1; // GLP_UP uses b2 as the value
                 glp_set_row_bnds(P,vlp_rowidx[i],glp_type(ctrl),b1,b2);
                 continue;
       case 'a': if(rows==0){
                    report(R_fatal,"read_vlp: a line before p in %s\n   %s\n",
                                PARAMS(VlpFile),inpline); return 1;
                 }
                 cnt=sscanf(inpline,"a %d %d %lg",&i,&j,&p);
                 if(cnt!=3 || rows<i || i<1 || cols<j || j<1){
                    report(R_fatal,"read_vlp: wrong a line in %s\n   %s\n",
                                PARAMS(VlpFile),inpline); return 1;
                 }
                 M(vlp_rowidx[i],vlp_colidx[j])=p; // store it
                 continue;
       case 'o': if(rows==0){
                    report(R_fatal,"read_vlp: o line before p in %s\n   %s\n",
                                PARAMS(VlpFile),inpline); return 1;
                 }
                 cnt=sscanf(inpline,"o %d %d %lg",&i,&j,&p);
                 if(cnt!=3 || objs<i || i<1 || cols<j || j<1){
                    report(R_fatal,"read_vlp: wrong o line in %s\n   %s\n",
                                PARAMS(VlpFile),inpline); return 1;
                 }
                 M(vlp_objidx[i],vlp_colidx[j])=dir*p; // store it
                 continue;
       case 'x': if(rows==0){
                     report(R_fatal,"read_vlp: x line before p in %s\n   %s\n",
                                 PARAMS(VlpFile),inpline); return 1;
                 }
                 cnt=sscanf(inpline,"x %d %lg",&i,&p);
                 if(cnt!=2 || objs<i||i<1){
                    report(R_fatal,"read_vlp: wrong x line in %s\n   %s\n",
                                 PARAMS(VlpFile),inpline); return 1;
                 }
                 vlp_init[i]=p; // store it
                 continue;
       default: report(R_fatal,"read_vlp: unknown line in %s\n  %s\n",
                           PARAMS(VlpFile),inpline); return 1;
    }
    fclose(f);
    if(rows==0){
       report(R_fatal,"read_vlp: no 'p' line in %s\n",PARAMS(VlpFile)); return 1; 
    }
    /* the vlp file has been read; set the glpk LP instance */
    // check if vlp_init[] is all positive
    for(i=1;i<=objs;i++){
        if(vlp_init[i]<PARAMS(PolytopeEps)){
           report(R_fatal,"read_vlp: initial value[%d]=%lg not positive\n",
                           i,vlp_init[i]); return 1;
        }
    }
    // upload constraints into P
    for(i=0;i<=rows+objs;i++) vlp_rowidx[i]=i; // index file
    for(j=1;j<=cols+1;j++){
        if(j!=lambda_idx) glp_set_mat_col(P,j,rows+objs,vlp_rowidx,&M(1,j)-1);
    }
    free(vlp_colidx); free(vlp_rowidx); free(vlp_M);
    // set objective lines to =vlp_init[]
    for(i=1;i<=objs;i++) glp_set_row_bnds(P,vlp_objidx[i],GLP_FX,vlp_init[i],vlp_init[i]);
    // LP objective: maximize lambda
    glp_set_obj_coef(P,lambda_idx,1.0);
    glp_set_obj_dir(P,GLP_MAX);
    return 0;
}

/**********************************************************************
* void set_oracle_parameters(void)
*   Set the LP solver parameters from the configuration:
*    verbosity: 0: no, 1: error; 2: normal, 3: verbose
*    output frequency: indicate that the LP solver is working
*    method: primal, dual
*    pricing: standard or steepest edge
*    ratio test: standard of Harris
*    iteration limit
*    time limit (in seconds)
*/

static void set_oracle_parameters(void)
{   glp_init_smcp(&parm);
    // verbosity
    switch(PARAMS(OracleMessage)){
      case  0: parm.msg_lev=GLP_MSG_OFF; break;
      case  1: parm.msg_lev=GLP_MSG_ERR; break;
      case  2: parm.msg_lev=GLP_MSG_ON; break;
      default: parm.msg_lev=GLP_MSG_ALL; break;
    }
    // method, pricing, ratio test
    parm.meth = GLP_PRIMAL;	// PRIMAL,DUAL,DUALP
    if(PARAMS(OracleMethod)) parm.meth=GLP_DUAL;
    parm.pricing = GLP_PT_STD;
    if(PARAMS(OraclePricing)) parm.pricing=GLP_PT_PSE;
    parm.r_test = GLP_RT_STD;	// HAR
    if(PARAMS(OracleRatioTest)) parm.r_test=GLP_RT_HAR;
    // iteration and time limit
    parm.it_lim = 100000;	// iteration limit
    if(PARAMS(OracleItLimit)>=1000) parm.it_lim=PARAMS(OracleItLimit);
    if(PARAMS(OracleItLimit)==0) parm.it_lim=0; // no limit
    parm.tm_lim = 10000;	// time limit in milliseconds
    if(PARAMS(OracleTimeLimit)>=5) parm.tm_lim=1000*PARAMS(OracleTimeLimit);
    if(PARAMS(OracleTimeLimit)==0) parm.tm_lim=0; // no limit
}

/*********************************************************************&
* char *glp_status_msg(int code)
* char *glp_return_msg(int code)
*  return the verbatim error message corresponding to the glpk
*  code (to be printed out).
*/

static char *glp_status_msg(int stat)
{static char *statmsg[] = {
"the problem is undefined",		// GLP_UNDEF
"solution is feasible",			// GLP_FEAS
"solution is infeasible",		// GLP_INFEAS
"the problem has no feasible solution",	// GLP_NOFEAS
"solution is optimal",			// GLP_OPT
"the problem is unbounded",		// GLP_UNBND
};
    if(1<=stat && stat<=6) return statmsg[stat-1];
    return "unknown solution status";
}

static char *glp_return_msg(int retval)
{static char *retmsg[] = {
"invalid basis",			// GLP_EBADB	 *
"singular matrix",			// GLP_ESING	 *
"ill-conditioned matrix",		// GLP_ECOND	 *
"invalid bounds",			// GLP_EBOUND
"solver failed",			// GLP_EFAIL	 *
"objective lower limit reached",	// GLP_EOBJLL
"objective upper limit reached",	// GLP_EOBJUL
"iteration limit exceeded",		// GLP_EITLIM	 *
"time limit exceeded",			// GLP_ETMLIM	 *
"no primal feasible solution",		// GLP_ENOPFS
"no dual feasible solution",		// GLP_ENODFS
"root LP optimum not provided",		// GLP_EROOT
"search terminated by application",	// GLP_ESTOP
"relative mip gap tolerance reached",	// GLP_EMIPGAP
"no primal/dual feasible solution",	// GLP_ENOFEAS
"no convergence",			// GLP_ENOCVG
"numerical instability",		// GLP_EINSTAB
"invalid data",				// GLP_EDATA
"result out of range",			// GLP_ERANGE
};
    if(1<=retval && retval<=0x13) return retmsg[retval-1];
    return "unknown error";
}

/**********************************************************************
* int initialize_oracle(void)
*  Check the consistency of the loaded vlp problem
*
* int ask_oracle(void)
*  The question and the answer is provided in OracleData. Return values:
*/

#include "round.h" /* round_to */

/* count the number and measure time of LP calls */
static int oracle_calls=0;
static unsigned long oracle_time=0ul;

/* call the glpk simplex solver twice if necessary */
#include <sys/time.h> 
static int call_glp(void)
{int ret; struct timeval tv; unsigned long starttime;
    oracle_calls++;
    if(gettimeofday(&tv,NULL)==0){
        starttime=tv.tv_sec*1000 + (tv.tv_usec+500u)/1000u;
    } else {starttime=0ul;}
    if(PARAMS(OracleMessage)<2) glp_term_out(GLP_OFF);
    glp_sort_matrix(P);
    if(PARAMS(OracleScale)) glp_scale_prob(P,GLP_SF_AUTO);
    glp_adv_basis(P,0);  // make this optimization
    glp_term_out(GLP_ON); // enable messages form glpk
    ret=glp_simplex(P,&parm);
    if(ret==GLP_EBADB || ret==GLP_ESING){ // invalid base
       if(PARAMS(OracleMessage)<2) glp_term_out(GLP_OFF);
       if(PARAMS(OracleScale)) glp_scale_prob(P,GLP_SF_AUTO);
       glp_adv_basis(P,0);
       glp_term_out(GLP_ON);
       oracle_calls++;
       ret=glp_simplex(P,&parm);
    }
    if(ret==GLP_EFAIL){ // give it a second chance
        if(PARAMS(OracleMessage)<2) glp_term_out(GLP_OFF);
        glp_adv_basis(P,0);
        glp_term_out(GLP_ON);
        oracle_calls++;
        ret=glp_simplex(P,&parm);
    }
    if(gettimeofday(&tv,NULL)==0)
        oracle_time += (tv.tv_sec*1000 + (tv.tv_usec+500u)/1000u)-starttime;
    return ret;
}

int initialize_oracle(void)
{int ret;
    set_oracle_parameters();
    // check if E is an internal point
    // the lambda column is all zero
    glp_set_obj_dir(P,GLP_MIN);
    ret=call_glp();
    if(ret){
       report(R_fatal,"Internal point: the oracle says: %s\n",glp_return_msg(ret));
       return ORACLE_FAIL;
    }
    ret=glp_get_status(P);
    if(ret!=GLP_OPT){
       report(R_fatal,"Internal point, the oracle says: %s\n",glp_status_msg(ret));
       return ret==GLP_NOFEAS ? ORACLE_EMPTY : ORACLE_FAIL;
    } // otherwise OK
    glp_set_obj_dir(P,GLP_MAX);
    return ORACLE_OK;
}

/* int ask_oracle() 
*   ask oracle about vvertex[0:vobjs], return vfacet[0:vobjs] as the
*   separating supporting hyperplane; 
*   vvertex*vfacet<=0; vlp_init*vfacet>0
*/
int ask_oracle(void)
{int i,ret; double lambda,d;
    if(vvertex[vobjs]==0.0){ // ideal point
       for(i=1;i<=vobjs;i++){vlp_lambda[i]=0.0-vvertex[i-1]; }
    } else { // not an ideal point
       for(i=1;i<=vobjs;i++){vlp_lambda[i]=vlp_init[i]-vvertex[i-1]; }
    }
    glp_set_mat_col(P,lambda_idx,vobjs,vlp_objidx,vlp_lambda);
    ret=call_glp();
    if(ret){
        report(R_fatal,"The oracle says: %s (%d)\n",glp_return_msg(ret),ret);
        // one can continue if  ret==GLP_EITLIM || ret==GLP_ETMLIM
        return ORACLE_FAIL;
    }
    ret=glp_get_status(P);
    if(ret == GLP_UNBND){
        if(vvertex[vobjs]==0.0){ return ORACLE_UNBND;} // inside
        report(R_fatal,"The oracle says: problem unbounded\n");
        return ORACLE_FAIL;
    }
    if(ret != GLP_OPT){
        report(R_fatal,"The oracle says: %s (%d)\n",glp_status_msg(ret),ret);
        return ORACLE_FAIL; 
    }
    lambda=glp_get_obj_val(P); // optimal value of lambda
    /* if lambda ==1 the vertex is inside the polytope
     * the boundary point: vlp_init[i]-lambda*vlp_lambda[i]
     * the facet equation is the dual solution */
    if(lambda<10.0*PARAMS(PolytopeEps)){
      report(R_fatal,"Initial point is on the boundary\n");
      return ORACLE_FAIL;
    }
    if(vvertex[vobjs]!=0.0 && lambda > 1.0-PARAMS(PolytopeEps)){
        if(lambda>1.0+PARAMS(PolytopeEps)){
           report(R_fatal,"Numerical problem, lambda=%lg > 1.0\n",lambda);
           return ORACLE_FAIL;
        }
        return ORACLE_UNBND;
    }
    for(i=1;i<=vobjs;i++) vfacet[i-1]=glp_get_row_dual(P,vlp_objidx[i]);
    // normalize the equation to sum up to 1.0
    d=0.0;
    for(i=0;i<vobjs;i++) d += vfacet[i]<0 ? -vfacet[i]:vfacet[i];
    if(d<PARAMS(PolytopeEps)){
       report(R_fatal,"Numerical problem, facet all zero\n");
       return ORACLE_FAIL;
    }
    for(i=0;i<vobjs;i++) vfacet[i] /= d;
    if(PARAMS(RoundFacets)){
       for(i=0;i<vobjs;i++){
           d=vfacet[i]; round_to(&d);
           vfacet[i]=d; }
    }
    // the optimal solution is on the supporting hyperplane
    d=0.0;
    for(i=1;i<=vobjs;i++){
        d-=vfacet[i-1]*(vlp_init[i]-lambda*vlp_lambda[i]);
    }
    if(PARAMS(RoundFacets)) round_to(&d); 
    vfacet[vobjs]=d;
    // check that vvertex is on the negative side, vlp_init is on the positive side
    d=0.0; for(i=0;i<=vobjs;i++) d+=vvertex[i]*vfacet[i];
    if(d>0.0){
        report(R_fatal,"Numerical error: vertex is on the negative side (%lg)\n",d);
        return ORACLE_FAIL;
    }
    d=vfacet[vobjs]; for(i=1;i<=vobjs;i++) d+=vlp_init[i]*vfacet[i-1];
    if(d<PARAMS(PolytopeEps)){
        report(R_fatal,"Initial point is on the negative side (%lg) of the next facet\n",d);
        return ORACLE_FAIL;
    }
    return ORACLE_OK;
}

/**********************************************************************
* Get oracle statistics
*
* void get_oracle_stat(int *no, init *it, unsigned long *time)
*   return the number of LP calls and time. Use the 
*   undocumented glpk call to get the number of iterations.
*/

void get_oracle_stat(int *no, int *it, unsigned long *time)
{   *no=oracle_calls; *it=glp_get_it_cnt(P); 
    *time=(oracle_time+5ul)/10ul; // in 0.01 seconds
}

/* EOF */

