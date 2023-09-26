#ifndef _MOLEC_H
#define _MOLEC_H

#include <stdio.h>
#include "atypes.h"
#include "pdbconst.h"
#include "lapack.h"

#define ATOMNAMESIZ (PDB_ATOM_name_LEN+1)
#define RESNAMESIZ (PDB_ATOM_resName_LEN+1)

typedef struct {
  int id;
  Real *p;
  int *occ;
//  int disordered;
  char name[ATOMNAMESIZ];
  char resname[RESNAMESIZ];
  int resid;
//  int conformid;
  Real biso;
} Atom;

typedef struct {
  /*
  int id;
  */
  int index;
  // pdbChain *pdbchp; // int arche
  Atom *head;
  int natoms;
  Real *x, *y;                  // vec(X), vec(Y)
  Real *r, *t;			// elements of R and t. NEVER SET in MSFIT
  int *o;			// binary occupancy factor
  Real *g;
  int oty;                      // occupancy pattern. type of Sigma
  // Real *s;                      // vec(Sig)
  // Real *diags;                  // diag(S)
  struct covar *covp;
  Real rmsd;                    // weighted rmsd
  Real rmsd_uw;                 // unweighted rmsd
} Conformer;


typedef struct archetype {              // chain / template / model
  int id;
  Real *x, *y, *g;
  int *o;
  int nconfms;
  Conformer *confmbp;
  Real *r;              // elements of R
  Real *t;
} Archetype;

#define JV_prepUC	1
#define JV_prepInv	2

typedef struct jvar {
  Real *u;	// Ui = S (S+miV)^{-1} V. symmetric since Ui^{-1} symmetric
  Real *sv;	// (S+miV)^{1/2} (Cholesky-factorized S+miV)
  Real *svi;	// (S+miV)^{-1}. symmetric
  // Real *sviv;	// [(S+miV)^{-1} V]^T.
  Real *b;	// [(S+miV)^{-1} V]^T.
  Real *c;	// S^{-1} [(S+miV)^{-1} V]^T. symmetric (== S^{-1} Ui^T S^{-1})
  Real *sc;	// S^{-1} - Ci
  Real logdet;
  int mi;	// == nconfms[i]
  Real *h;
  int nens;	// num of ensembles with this mi
  struct jvar *next;
} jVar;

typedef struct ensemble {		// chain / template / model / archetype
  int id;
  Real *x, *y, *z; //, *r, *t;
  int *o;
  int nconfms;
  Conformer *confmbp;
  struct jvar *jvp;
} Ensemble;

#define NUMJVMATRICES	7		// u,sv,svi,b,c,sc,h

#endif
