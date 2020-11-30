// This code is
// Copyright (C) 1998-2020 T. W. Yee, University of Auckland. All rights reserved.


#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<R.h>
#include<Rmath.h>


void vdecccc(int *row_index, int *col_index, int *dimm);
void m2accc(double *m, double *a, int *dimm, int *row_index,
       int *col_index, int *n, int *M, int *upper);
void a2mccc(double *a, double *m, int *dimm, int *row_index,
       int *col_index, int *n, int *M);
void mux2ccc(double *cc, double *ymat,
          double *ans, int *p, int *n, int *M);
void mux22ccc(double *cc, double *ymat, double *ans, int *dimm,
       int *row_index, int *col_index, 
       int *n, int *M, double *wk, int *upper);
void mux5ccc(double *cc, double *x,
          double *ans, int *M, int *n, int *r,
          int *dimm,
          int *dimr,
          int *matrix,
          double *wk, double *wk2,
          int *row_index_M, int *col_index_M,
          int *row_index_r, int *col_index_r);
void mux55ccc(double *evects, double *evals, double *ans, double *wk, 
           double *wk2, int *row_index, int *col_index,
           int *M, int *n);
void mux7ccc(double *cc, double *x,
          double *ans, int *M, int *q, int *n, int *r);
void mux111ccc(double *cc, double *txmat, int *M, int *R, int *n,
        double *wkcc, double *wk2, int *row_index, int *col_index,
	    int *dimm, int *upper);
void mux15ccc(double *cc, double *x,
           double *ans, int *M, int *n);
void vcholccc(double *cc, int *M, int *n, int *ok, double *wk,
           int *row_index, int *col_index, int *dimm);
void vforsubccc(double *cc, double *b, int *M, int *n,
             double *wk, int *row_index,
             int *col_index, int *dimm);
void vbacksubccc(double *cc, double *b, int *M, int *n,
              double *wk, int *row_index,
              int *col_index, int *dimm);
void tapply_mat1(double *mat, int *nr, int *nc, int *type);
















void vdecccc(int *row_index, int *col_index, int *dimm) {
  int ilocal;

  for(ilocal = 0; ilocal < *dimm; ilocal++) {
    row_index[ilocal] -= 1;
    col_index[ilocal] -= 1;
  }
}


void m2accc(double *m, double *a, int *dimm, int *row_index,
         int *col_index, int *n, int *M, int *upper) {
  int ilocal, klocal, MM = *M * *M, MMn = *M * *M * *n;

  if(*upper == 1 || *dimm != *M * (*M + 1) / 2)
    for(klocal = 0; klocal < MMn; klocal++)
        a[klocal] = 0.0;

  for(klocal = 0; klocal < *n; klocal++) {
    for(ilocal = 0; ilocal < *dimm; ilocal++) {
      a[row_index[ilocal] + col_index[ilocal] * *M] = m[ilocal];
      if(*upper == 0)
        a[col_index[ilocal] + row_index[ilocal] * *M] = m[ilocal];
    }
    a += MM;
    m += *dimm;
  }
}


void a2mccc(double *a, double *m, int *dimm, int *row_index,
         int *col_index, int *n, int *M) {
  int ilocal, klocal, MM= *M * *M;

  for(klocal = 0; klocal < *n; klocal++) {
    for(ilocal = 0; ilocal < *dimm; ilocal++)
      m[ilocal] = a[row_index[ilocal] + col_index[ilocal] * *M];
    a += MM;
    m += *dimm;
  }
}




void mux2ccc(double *cc, double *ymat,
          double *ans, int *p, int *n, int *M) {
  double slocal;
  int ilocal, jlocal, tlocal, Mp = *M * *p;

  for(ilocal = 0; ilocal < *n; ilocal++) {
    for(jlocal = 0; jlocal < *M; jlocal++) {
      slocal = 0.0;
      for(tlocal = 0; tlocal < *p; tlocal++)
        slocal += cc[jlocal + tlocal * *M] * ymat[tlocal];
      *ans++ = slocal;
    }
    ymat += *p;
    cc += Mp;
  }
}



void mux22ccc(double *cc, double *ymat, double *ans, int *dimm,
       int *row_index, int *col_index, 
       int *n, int *M, double *wk, int *upper) {
  double slocal;
  int jlocal, tlocal, klocal, one = 1, lower;

  vdecccc(row_index, col_index, dimm);
  for(klocal = 0; klocal < *n; klocal++) {
    m2accc(cc, wk, dimm, row_index, col_index, &one, M, upper);

    for(jlocal = 0; jlocal < *M; jlocal++) {
      slocal = 0.0;
      lower = *upper == 0 ? 0 : jlocal; 
      for(tlocal = lower; tlocal < *M; tlocal++)
          slocal += wk[jlocal + tlocal * *M] * ymat[tlocal];
      *ans++ = slocal;
    }
    ymat += *M;
    cc += *dimm;
  }
}


void mux5ccc(double *cc, double *x,
          double *ans, int *M, int *n, int *r,
          int *dimm,
          int *dimr,
          int *matrix,
          double *wk, double *wk2,
          int *row_index_M, int *col_index_M,
          int *row_index_r, int *col_index_r) {
  double slocal, *pd, *pd2;
  int ilocal, jlocal, klocal, tlocal, Mr = *M * *r,
      rr = *r * *r, MM = *M * *M, ulocal,
      jM, jr, kM, kr, one=1, upper=0;

  if(*matrix == 1) {
    vdecccc(row_index_M, col_index_M, dimm);
    vdecccc(row_index_r, col_index_r, dimr);
    pd = wk;
    pd2 = wk2;
  } else {
// Commented out on 2/5/06. Need to fix this up more cleanly.
//      Rprintf("Error: can only handle matrix.arg == 1\n");
//      exit(-1); 
//

// 20070926:
// The following line was added only to avoid a warning message from the compiler
//
    pd = pd2 = wk;

  }

  for(ilocal = 0; ilocal < *n; ilocal++) {
    if(*matrix == 1)
      m2accc(cc, pd, dimm, row_index_M, col_index_M, &one, M, &upper);
    else {
      pd = cc;
      pd2 = ans;
    }

    for(jlocal = 0; jlocal < *r; jlocal++) {
      jM = jlocal * *M;
      jr = jlocal * *r;
      for(klocal = jlocal; klocal < *r; klocal++) {
        kM = klocal * *M;
        kr = klocal * *r;
        slocal = 0.0;
        for(tlocal = 0; tlocal < *M; tlocal++)
          for(ulocal = 0; ulocal < *M; ulocal++)
            slocal +=  x[tlocal + jM] * pd[tlocal + ulocal * *M] *
                       x[ulocal + kM];
        pd2[jlocal + kr] =
        pd2[klocal + jr] = slocal;
      }
    }

    if(*matrix == 1)
      a2mccc(pd2, ans, dimr, row_index_r, col_index_r, &one, r);

    cc += (*matrix == 1 ? *dimm : MM);
    x += Mr;
    ans += (*matrix == 1 ? *dimr : rr);
  }
}



void mux55ccc(double *evects, double *evals, double *ans, double *wk, 
           double *wk2, int *row_index, int *col_index,
           int *M, int *n) {
  double *pd, *pd2, tlocal;
  int ilocal, jlocal, klocal, slocal, MM = *M * *M, one = 1,
      MM12 = *M * (*M + 1)/2;

  vdecccc(row_index, col_index, &MM12);

  for(ilocal = 0; ilocal < *n; ilocal++) {
    pd = evects;
    pd2 = wk2;
    for(jlocal = 0; jlocal < *M; jlocal++)
      for(klocal = 0; klocal < *M; klocal++)
        *pd2++ = *pd++ * evals[jlocal];

    for(jlocal = 0; jlocal < *M; jlocal++)
      for(klocal = jlocal; klocal < *M; klocal++) {
        tlocal = 0.0; 
        for(slocal = 0; slocal < *M; slocal++)
            tlocal +=    wk2[jlocal + slocal * *M] *
                      evects[klocal + slocal * *M];
        wk[jlocal + klocal * *M] =
        wk[klocal + jlocal * *M] = tlocal;
      }

    a2mccc(wk, ans, &MM12, row_index, col_index, &one, M);

    ans += MM12;
    evals += *M;
    evects += MM;
  }
}





void mux7ccc(double *cc, double *x,
          double *ans, int *M, int *q, int *n, int *r) {
  double slocal;
  int ilocal, jlocal, klocal, tlocal,
      Mq = *M * *q, qr = *q * *r, Mr = *M * *r,
      kq, kM;

  for(ilocal = 0; ilocal < *n; ilocal++) {
    for(jlocal = 0; jlocal < *M; jlocal++) {
      for(klocal = 0; klocal < *r; klocal++) {
        kq = klocal * *q;
        kM = klocal * *M;
        slocal = 0.0;
        for(tlocal = 0; tlocal < *q; tlocal++)
          slocal += cc[jlocal + tlocal * *M] * x[tlocal + kq];
        ans[jlocal + kM] = slocal;
      }
    }
    cc += Mq;
    ans += Mr;
    x += qr;
  }
}







/* 20170605; seems to crash for hdeff(), so looking at this again. */
/* 20170605; checked okay. */
void mux111ccc(double *cc, double *txmat,
            int *M, int *R, int *n,
            double *wkcc, double *wk2,
            int *row_index, int *col_index,
            int *dimm, int *upper) {
  double slocal, *pd2, tempdouble;
  int ilocal, jlocal, klocal, tlocal,
      MM = *M * *M, MR = *M * *R,
      lowlim;

/* dimm is all that is needed, although could use M*(M+1)/2 really */
  vdecccc(row_index, col_index, dimm);

  for(ilocal = 0; ilocal < MM; ilocal++)
    wkcc[ilocal] = 0.0;

  for(tlocal = 0; tlocal < *n; tlocal++) {
/* Copy a matrix from cc into wkcc. */
/* wkcc <- cc */
    for(ilocal = 0; ilocal < *dimm; ilocal++) {
      if(*upper == 0) {
/* The matrix is NOT upper triangular; it is symmetric. */
        tempdouble = *cc++;
        wkcc[row_index[ilocal] + col_index[ilocal] * *M] =
        wkcc[col_index[ilocal] + row_index[ilocal] * *M] = tempdouble;
      } else {
/* The matrix is upper triangular. */
        wkcc[row_index[ilocal] + col_index[ilocal] * *M] = *cc++;
      }
    }  /* ilocal */

/* Copy part of xmat (via txmat) into wk2. */
/* wk2 <- xmat  # M x R */
    pd2 = txmat;
    for(ilocal = 0; ilocal < *M; ilocal++)
      for(jlocal = 0; jlocal < *R; jlocal++)
        wk2[ilocal + jlocal * *M] = *pd2++;

/* Assign the product into the correct place in txmat. */
/* That is, a block of txmat is overwritten by the answer. */
    for(ilocal = 0; ilocal < *M; ilocal++) {
      lowlim = *upper == 0 ? 0 : ilocal;
      for(jlocal = 0; jlocal < *R; jlocal++) {
        slocal = 0.0;
        for(klocal = lowlim; klocal < *M; klocal++)
          slocal +=  wk2[klocal + jlocal * *M] *
                    wkcc[ilocal + klocal * *M];
        txmat[jlocal + ilocal * *R] = slocal;
      } /* jlocal */
    } /* ilocal */
    txmat += MR;
  }  /* tlocal */
}

void mux15ccc(double *cc, double *x,
           double *ans, int *M, int *n) {
  double *pd, *pd2;
  int ilocal, jlocal, klocal, MM = *M * *M;

  for(ilocal = 0; ilocal < *n; ilocal++) {
    pd = cc;
    pd2 = ans;
    for(jlocal = 0; jlocal < *M; jlocal++)
      for(klocal = 0; klocal < *M; klocal++)
        *pd2++ = *pd++ * x[jlocal];

    pd2 = ans;
    for(jlocal = 0; jlocal < *M; jlocal++)
      for(klocal = 0; klocal < *M; klocal++) {
        *pd2 *= x[klocal];
        pd2++;
      }

    ans += MM;
    x += *M;
  }
}




void vcholccc(double *cc, int *M, int *n, int *ok, double *wk,
           int *row_index, int *col_index, int *dimm) {
  double slocal, *pd;
  int tlocal, ilocal, jlocal, klocal, iM, iiM, upper = 0, one = 1;

  vdecccc(row_index, col_index, dimm);
  pd = wk;

  for(tlocal = 0; tlocal < *n; tlocal++) {
    *ok = 1; 

    m2accc(cc, wk, dimm, row_index, col_index, &one, M, &upper);

    for(ilocal = 0; ilocal < *M; ilocal++) {
      slocal = 0.0;
      iM = ilocal * *M;
      iiM = ilocal + iM;
      for(klocal = 0; klocal < ilocal; klocal++)
        slocal += pd[klocal + iM] * pd[klocal + iM];

      pd[iiM] -= slocal;
      if(pd[iiM] < 0.0) {
        *ok = 0;
        break;
      }
      pd[iiM] = sqrt(pd[iiM]);

      for(jlocal = ilocal+1; jlocal < *M; jlocal++) {
        slocal = 0.0;
        for(klocal = 0; klocal < ilocal; klocal++)
          slocal += pd[klocal + iM] * pd[klocal + jlocal * *M];
        pd[ilocal + jlocal * *M] = (pd[ilocal + jlocal * *M] -
                                    slocal) / pd[iiM];
      }
    }

    a2mccc(wk, cc, dimm, row_index, col_index, &one, M);

    cc += *dimm;
    ok++;
  }
}



void vforsubccc(double *cc, double *b, int *M, int *n,
             double *wk, int *row_index,
             int *col_index, int *dimm) {
  double slocal, *pd;
  int jlocal, klocal, tlocal, upper = 1, one = 1;

  pd = wk;
  vdecccc(row_index, col_index, dimm);

  for(tlocal = 0; tlocal < *n; tlocal++) {
    m2accc(cc, wk, dimm, row_index, col_index, &one, M, &upper);

    for(jlocal = 0; jlocal < *M; jlocal++) {
      slocal = b[jlocal];
      for(klocal = 0; klocal < jlocal; klocal++)
        slocal -= pd[klocal + jlocal * *M] * b[klocal];
      b[jlocal] = slocal / pd[jlocal + jlocal * *M];
    }
    cc += *dimm;
    b += *M;
  }
}




void vbacksubccc(double *cc, double *b, int *M, int *n,
              double *wk, int *row_index,
              int *col_index, int *dimm) {
  double slocal, *pd;
  int jlocal, klocal, tlocal, upper = 1, one = 1;

  pd = wk;
  vdecccc(row_index, col_index, dimm);

  for(tlocal = 0; tlocal < *n; tlocal++) {
    m2accc(cc, wk, dimm, row_index, col_index, &one, M, &upper);

    for(jlocal = *M - 1; jlocal >= 0; jlocal--) {
      slocal = b[jlocal];
      for(klocal = jlocal + 1; klocal < *M; klocal++)
        slocal -= pd[jlocal + klocal * *M] * b[klocal];
      b[jlocal] = slocal / pd[jlocal + jlocal * *M];
    }
    cc += *dimm;
    b += *M;
  }
}



void tapply_mat1(double *mat, int *nr, int *nc, int *type) {
  double *pd = mat, *pd2 = mat + *nr;
  int ilocal, jlocal;

  if(*type == 1)
    for(jlocal = 2; jlocal <= *nc; jlocal++)
      for(ilocal = 0; ilocal < *nr; ilocal++, pd2++)
        *pd2 += *pd++;

  if(*type == 2) {
    pd2 = mat + *nr * *nc - 1;
    pd = pd2 - *nr;
    for(jlocal = *nc; jlocal >= 2; jlocal--)
      for(ilocal = 0; ilocal < *nr; ilocal++, pd2--)
        *pd2 -= *pd--;
  }

  if(*type == 3)
    for(jlocal = 2; jlocal <= *nc; jlocal++)
      for(ilocal = 0; ilocal < *nr; ilocal++, pd2++)
        *pd2 *= *pd++;

  if(*type < 1 || *type > 3)
    Rprintf("Error: *type not matched\n");
}




