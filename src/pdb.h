#ifndef _PDB_H
#define _PDB_H

#include "atypes.h"
#include "pdbconst.h"

/* #define PDB_NCHAINS_MAX 128	// max # of chains per entry		*/
#define PDB_NALTLOCS_MAX 8	// max # of locs per atom
#define PDB_NPARTNS_MAX 32	// max # of partitionings of atoms
#define PDB_OCCFRAC_BUFSIZ (PDB_NPARTNS_MAX*4)

/* markpdb */

typedef struct intrange {
  int lower, upper;
} intrange;

/* markpdb_utils.c でも ATOMNAMELEN を PDB_ATOM_name_LEN で置き換えること */
typedef struct landmark {
  // name[] must be large enough to hold `PDB_ATOM_name_LEN` chars
  char name[PDB_ATOM_name_LEN+1];		// atom name (CA,C,N,O,...)
  int nranges;
  intrange *range;
} landmark;

// atom type -> hasharray index (>=0)
struct atomtype_hashent {
  unsigned int atomtype;
  int index;
  struct atomtype_hashent *next;
};


/* [M=markpdb,C=mmfitc]---------- */

struct pdbSite {
  char altloc;					// [C]
  XYZ p;					// [C]
  Real occ;					// [C]
  Real Biso;					// [C]
  struct pdbSite *next;				// [C]
};
typedef struct pdbSite pdbSite;

typedef struct pdbAtom {
  int id;					// [C]
  int asnum;			/* Atom serial number, only in genconfms */
  char name[PDB_ATOM_name_LEN+1];		// [C]
  char resname[PDB_ATOM_resName_LEN+1];		// [C]
  int resseq;					// [C]
  //  //  char chainid;
  struct pdbSite *sp;				// [C]
  int naltlocs;					// [C]
} pdbAtom;

typedef struct {
  // int id;		/* [M] global chain id */
  int model;		/* [CM] MODEL id */
  int entry;		/* [CM] serial num for a PDB entry */

  int chainid;		/* [CM] a char in the chain ID field of a PDB ATOM record */
  int index;		/* [C] 上記 id で置き換え可能か？ index of pdbchainpool[] */
  char altlocs[PDB_NALTLOCS_MAX+1];	/* [M] */
  int resseq_max;			/* [M] */
  int resseq_min;			/* [M] */
  int occ_nparts;			/* [C] # of ways of partitioning occs */
  Real occ_frac[PDB_OCCFRAC_BUFSIZ];	/* [C] sequences of occupancies */
  int  occ_nfracs[PDB_NPARTNS_MAX];	/* [C] # of occs in each partitioning */
  int  occ_firstid[PDB_NPARTNS_MAX];
  int natoms;				/* [C] */
} pdbChain;

typedef struct {
  int id;		// redundant? [C] seq num of the entry (begins with 0)
  char *file;			// [CM] file path name
  int nranges;			// [M]
  intrange *range;		// [M] structures (chain/model ids)
  int nchains;			// [CM] num of structures read in
  //  int resseq_max;
  pdbChain *chainptr;		// [CM] ptr to the first chain of this pdbEntry
  //  /* pdbChain *chain[PDB_NCHAINS_MAX]; */
  int nchains_block;		// [C] num of chains in a MODEL block
} pdbEntry;


extern int
read_pdbm(FILE *fp, int id2ind[], int id_lim,
	  pdbAtom atom[], pdbSite site[],
	  int *id_minp, int *id_maxp, int *natomsp, int *nsitesp, int *modelp);

extern int
read_trans_write_pdbm(FILE *fp, pdbEntry *pe, Real *Rbuf, Real *tbuf,
		      int coordsonly);

#endif
