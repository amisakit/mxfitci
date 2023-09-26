#include <limits.h>
#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
       #include <sys/types.h>
       #include <sys/stat.h>
       #include <fcntl.h>
       #include <unistd.h>
#include "atypes.h"
#include "pdb.h"
#include "matrix.h"
#include "memory.h"
#include "msfit.h"
//#include "msfit.h"
#include "rand.h"


/*
extern int
centerM(Real *mp, struct covar *covp, int natoms, int dodiag, Real *t,
        Real *work);
extern int
alignM(Real *mp, int natoms, Real *vout);
*/

/*
 * TBD:
 */

#ifdef MKLROOT
int UseMKL = 1;
#else
int UseMKL = 0;
#endif


#define WORKAREASIZ 5000000

char *progname = 0;

int *LUTbuf = 0;

void abnormalexit(char *message, int code)
{
  if (message && code)
    fprintf(stderr, "(%d) %s\n", code, message);
  else if (message)
    fprintf(stderr, "%s\n", message);
  if (progname) free(progname);
  exit(code);
}

int quick_count_atoms(FILE *fp);
int read_structure_data(FILE *fp, Real *cor, pdbAtom *atom, pdbSite *site);
void print_conformer(pdbAtom *templ, Real *xp, Real *sd, Real *op, int natoms, int iconfm, int biso, int extraxyz, int nomit, landmark *alist, Real prob, void *work);

extern int lists_of_intranges(char *str, landmark *L, intrange *R);
#define LM_POOLSIZ 32
#define IR_POOLSIZ 1024
static intrange intrangepool[IR_POOLSIZ];
static landmark atomlist[LM_POOLSIZ];


char usagemessage[] =
  "  -s sig: scale factor of displacement [default: 1.0]\n"
  "  -n nconfms: num of conformers to be generated\n"
  "  -r: random rotate & shift [default: stay in the orig orientation/position]\n"
  "  -B: Biso in the tempFacotor field [default: SD of displacement in tempFac]\n"	// changed from -b [2021-05-13 Thu 09:06]
  "  -a: place the mean structure along its PAs before generating confms\n"
  "  -m: place the records in MODEL...ENDMODEL\n"
  "  -g seed: a seed for the random number generator\n"
  "  -x: write high precision xyz in the extra field\n"
  "  -o \"list\" p : omit atoms given in the list with probability p\n"
  "      list is a semicolon seperated groups of atod id. ex. \"1-3,5;10-20\"\n"
  "      group is a list of comma seperated numbers or ranges\n"
  "      of atom ids.  Membs of a group are omitted together with prob p\n"
  "      atom ids are renumbered (from zero) at input and are not those\n"
  "      given in the 'serial' field\n"
  "  -b v|a|m|x fname: binary I/O. see following comments.\n"
  "";

/*
 * binary I/O
 * -b m file: replace the contents of M[] with the coordinates in file,
 *            and **OVERWRITE file** with the centerM'ed coordinates.
 * -b v file: read nxn cov from file. sig^2 is applied.
 * -b a file: read 3nx3n cov from file. sig^2 is applied.
 * -b x file: write xyz to file in addition to stdout.
 *   When a character % is contained in 'file', it is regarded as a c-format
 *   string e.g., "crd%05d.bin", and the filenames are selected
 *   by `sprintf(filename, x, iconfm)`.
 */
char *binaryIO_v = 0;
char *binaryIO_m = 0;
char *binaryIO_x = 0;
int anisotropic = 0;

void usage()
{
  printf("usage: %s [OPTIONS] -s sig -n nconfms infile\n", progname);
  printf("%s", usagemessage);
  exit(-1);
}

int parse_args(int argc, char *argv[], char aglist[], Real *pomit,
	       char *infile, int *nconfms, Real *sig, int *extraxyz,
	       int *rotate, int *prerotate, int *biso, int *mdlout)
{
  int i, j, c, seed;

  aglist[0] = 0;			/* atom-group list */
  *sig = 1.0;
  *nconfms = 0;
  *rotate = *biso = *mdlout = *extraxyz = 0;
  *prerotate = 0;
  for (i = 1; i < argc; i++) {
    if (argv[i][0] != '-') {
      if (i == argc-1) {
	strcpy(infile, argv[argc-1]);
	return 0;
      } else
	usage();
    }
    for (j = 1; argv[i][j]; j++) {
      c = argv[i][j];
      switch (c) {
      case 'r': 
	*rotate = 1;  break;
      case 'B':
	*biso = 1;  break;      
      case 'a':
	*prerotate = 1;  break;      
      case 'm':
	*mdlout = 1;  break;      
      case 'x':
	*extraxyz = 1;  break;      
      case 'g':
	seed = atoi(argv[++i]);
	rand_seed(seed);
	srand48((long int)seed); break;
      case 's':
	*sig = atof(argv[++i]);  break;      
      case 'n':
	*nconfms = atoi(argv[++i]);  break;      
      case 'o':
	strcpy(aglist, argv[++i]);
	*pomit = atof(argv[++i]); break;      
      case 'b':
	switch (*argv[++i]) {
	case 'v':
	  binaryIO_v = argv[++i]; break;
	case 'a':
	  binaryIO_v = argv[++i]; anisotropic = 1; break;
	case 'm':
	  binaryIO_m = argv[++i]; break;
	case 'x':
	  binaryIO_x = argv[++i]; break;
	default:
	  usage();
	}
	break;
      default:
	usage();
      }
      if (c == 'g' || c == 's' || c == 'n' || c == 'o' || c == 'b') break;
    }
  }
  usage();
  return 0;
}

int
main(int argc, char *argv[])
{
  int nconfms, rotate, prerotate, biso, mdlout, natoms, iconfm, i, j, k, n;
  int extxyz, nomit;
  int fd, infixseq;
  Real sig, *cor, *cov, *sd, *mp, *xp, *occ, rp[9], tp[3], u[3], rmax, s;
  Real pomit;
  char infile[BUFSIZ], aglist[BUFSIZ], *work, outfile[BUFSIZ];
  FILE *fp;
  pdbAtom *atom;
  pdbSite *site;
  Matrix3 M, X, R;

  if (init_workarea(WORKAREASIZ) != 0)
    abnormalexit("cannot alloc workarea", 0);

  progname = strdup(argv[0]);
  parse_args(argc, argv, aglist, &pomit, infile, &nconfms, &sig,
	     &extxyz, &rotate, &prerotate, &biso, &mdlout);
  nomit = 0;
  if (aglist[0])
    nomit = lists_of_intranges(aglist, atomlist, intrangepool);
    
  if ((fp = fopen(infile, "r")) == 0)
    abnormalexit("cannot open infile", 0);
  natoms = quick_count_atoms(fp);
  rewind(fp);
  atom = (pdbAtom *) malloc(sizeof(pdbAtom) * natoms);
  site = (pdbSite *) malloc(sizeof(pdbSite) * natoms);
  cor = (Real *) malloc(sizeof(Real) * natoms * natoms);
  if (anisotropic)
    cov = (Real *) malloc(sizeof(Real) * natoms * natoms * 9);
  else
    cov = (Real *) malloc(sizeof(Real) * natoms * natoms);
  sd = (Real *) malloc(sizeof(Real) * natoms);
  occ = (Real *) malloc(sizeof(Real) * natoms);
  mp = (Real *) malloc(sizeof(Real) * natoms * 3);
  xp = (Real *) malloc(sizeof(Real) * natoms * 3);
  if (!atom || !site || !cor || !cov || !sd || !occ || !mp || !xp)
    abnormalexit("cannot alloc workarea", 0);

  identityvectmatrix(cor, natoms);
  n = read_structure_data(fp, cor, atom, site);
  fclose(fp);
  natoms = n;
  M = (Matrix3) mp;
  for (i = 0; i < natoms; i++) {
    M[i][0] = atom[i].sp->p.x;
    M[i][1] = atom[i].sp->p.y;
    M[i][2] = atom[i].sp->p.z;
    occ[i] = atom[i].sp->occ;
    if (biso)
      sd[i] = sqrt(atom[i].sp->Biso / 2.0) / (2.0 * M_PI);
    else
      sd[i] = atom[i].sp->Biso;
  }

  if (binaryIO_m && *binaryIO_m) {
    /* overwrite M with the contents of the file m */
    fd = open(binaryIO_m, O_RDONLY);
    if (fd < 0) {
      fprintf(stderr, "cannot open M file %s\n", binaryIO_m);
      exit(-1);
    }
    i = sizeof(Real) * natoms * 3;
    if (read(fd, mp, i) < i) {
      fprintf(stderr, "cannot read 3*natoms values\n");
      exit(-1);
    }
    close(fd);
  }


  centerM(mp, 0, natoms, 0);
  if (prerotate)
    alignM(mp, natoms, 0);
  if (binaryIO_m && *binaryIO_m) {
    fd = open(binaryIO_m, O_WRONLY);
    i = sizeof(Real) * natoms * 3;
    if (write(fd, mp, i) < i) {
      fprintf(stderr, "cannot write 3*natoms values in M[].\n");
      exit(-1);
    }
    close(fd);
  }

  rmax = 0.0;
  for (i = 0; i < natoms; i++) {
    if (fabs(M[i][0]) > rmax) rmax = fabs(M[i][0]);
    if (fabs(M[i][1]) > rmax) rmax = fabs(M[i][1]);
    if (fabs(M[i][2]) > rmax) rmax = fabs(M[i][2]);
  }

  if (binaryIO_v && *binaryIO_v) {
    /* overwrite the cov[] with the contents of the file v */
    fd = open(binaryIO_v, O_RDONLY);
    if (fd < 0) {
      fprintf(stderr, "cannot open covfile %s\n", binaryIO_v);
      exit(-1);
    }
    if (anisotropic)
      i = sizeof(Real) * natoms*natoms*9;
    else
      i = sizeof(Real) * natoms*natoms;
    if (read(fd, cov, i) < i) {
      fprintf(stderr, "cannot read natoms^2 values\n");
      exit(-1);
    }
    close(fd);
  } else {		/* use corrs and sds given in the pdb file */
    for (i = 0; i < natoms; i++)
      for (j = 0; j < natoms; j++)
	cov[i*natoms + j] = cor[i*natoms + j] * sd[i] * sd[j];
  }

  n = anisotropic ? natoms * 3 : natoms;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      cov[i*n + j] *= sig * sig;

  /* factorize at advance */
  i = cholesky_decomp(cov, n);
  if (i > 0) {
    fprintf(stderr, "Error: cholesky_decomp returns %d\n.  Aborting...", i);
    exit(-1);
  }

  if (binaryIO_x && *binaryIO_x) {
    /* binary xyz to x file */
    infixseq = 1;
    if (strchr(binaryIO_x, '%') == NULL) infixseq = 0;

    if (!infixseq) {
      fd = open(binaryIO_x, O_WRONLY | O_CREAT, 0644);
      if (fd < 0) {
	fprintf(stderr, "cannot open crdfile %s for writing\n", binaryIO_x);
	exit(-1);
      }
    }
  }

  R = (Matrix3) rp;
  X = (Matrix3) xp;
  work = (char *) cor;
  for (iconfm = 0; iconfm < nconfms; iconfm++) {
    if (anisotropic)
      i = rand_mvnorm(mp, cov, natoms*3, xp, 1);
    else
      i = rand_m3norm(mp, cov, natoms, xp, 1);
    if (i > 0) {
      fprintf(stderr, "Error: rand_mv/m3norm returns %d\n.  Aborting...", i);
      exit(-1);
    }
    if (rotate) {
      rand_rotmat(rp);
      tp[0] = (rand_uniformo() - 0.5) * 2.0 * rmax;
      tp[1] = (rand_uniformo() - 0.5) * 2.0 * rmax;
      tp[2] = (rand_uniformo() - 0.5) * 2.0 * rmax;
      for (i = 0; i < natoms; i++) {
	for (j = 0; j < 3; j++) {
	  s = 0.0;
	  for (k = 0; k < 3; k++)
	    s += R[j][k] * X[i][k];
	  u[j] = s + tp[j];
	  u[j] = s;
	}
	X[i][0] = u[0];  X[i][1] = u[1];  X[i][2] = u[2];
      }
    }
    if (iconfm == 0)
      print_conformer(atom, xp, sd, occ, natoms, (iconfm+1)*mdlout, biso,
		      extxyz,     0,        0,   0.0, work);
    else
      print_conformer(atom, xp, sd, occ, natoms, (iconfm+1)*mdlout, biso,
		      extxyz, nomit, atomlist, pomit, work);

    if (binaryIO_x && *binaryIO_x) {
      if (infixseq) {
	sprintf(outfile, binaryIO_x, iconfm);
	fd = open(outfile, O_WRONLY | O_CREAT, 0644);
	if (fd < 0) {
	  fprintf(stderr, "cannot open crdfile %s for writing\n", outfile);
	  exit(-1);
	}
      }
      i = sizeof(Real) * natoms * 3;
      if (write(fd, xp, i) < i) {
	fprintf(stderr, "cannot write natoms*3 values\n");
	exit(-1);
      }
      if (infixseq)
	close(fd);
    }
  }

  if (binaryIO_x && *binaryIO_x && !infixseq)
    close(fd);

  free(atom);
  free(site);
  free(cor);
  free(cov);
  free(sd);
  free(occ);
  free(mp);
  free(xp);
  finalize_workarea();
}


static double read_double(char rec[], int len)
{
  double d;
  char buf[BUFSIZ];
  strncpy(buf, rec, len); buf[len] = 0;
  sscanf(buf, "%lf", &d);
  return d;
}

static int read_int(char rec[], int len)
{
  int i;
  char buf[BUFSIZ];
  strncpy(buf, rec, len); buf[len] = 0;
  sscanf(buf, "%d", &i);
  return i;
}

static void read_string(char *rec, int len, char s[])
{
  for (; len > 0; rec++,len--)
    if (!isspace(*rec)) break;
  strncpy(s, rec, len);
  for (; len > 0; len--)
    if (!isspace(s[len-1])) break;
  s[len] = 0;
}

static int pdbrec(char *keyword, char *rec)
/* returns 1 if `rec` is the type of `keyword`, otherwise 0. */
{
  int i;
  char kbuf[7], rbuf[7];

  for (i = 0; i < 6; i++) {
    if (keyword[i] == 0 || isspace(keyword[i]))
      break;
    kbuf[i] = keyword[i];
  }
  for ( ; i < 6; i++)
    kbuf[i] = ' ';

  for (i = 0; i < 6; i++) {
    if (rec[i] == 0 || isspace(rec[i]))
      break;
    rbuf[i] = rec[i];
  }
  for ( ; i < 6; i++)
    rbuf[i] = ' ';

  for (i = 0; i < 6; i++)
    if (kbuf[i] != rbuf[i])
      return 0;
  return 1;
}

int
quick_count_atoms(FILE *fp)
{
  int natoms;
  char rec[BUFSIZ];

  natoms = 0;
  while (fgets(rec, BUFSIZ, fp))
    if (pdbrec("ATOM", rec)) 
      natoms++;
  return natoms;
}

/*
 * RESTRICTIONS:
 * altlocs (frac occs) are not supported
 */
int
read_structure_data(FILE *fp, Real *cor, pdbAtom *atom, pdbSite *site)
{
  int in_section_corr;
  int natoms, i, j, k, aseqi, aseqj;
  Real val;
  pdbSite *sp;
  char rec[BUFSIZ];

  in_section_corr = 0;
  natoms = 0;
  while (fgets(rec, BUFSIZ, fp)) {
    if (pdbrec("CORR", rec)) {
      sscanf(rec+6, "%d%d%lf", &aseqi, &aseqj, &val);
      i = j = -1;
      for (k = 0; k < natoms; k++) {
        if (atom[k].asnum == aseqi) i = k;
        if (atom[k].asnum == aseqj) j = k;
      }
      if (i < 0 || j < 0 || i == j)
        abnormalexit("CORR record contains an illegal pair.\n", 0);
      cor[i*natoms+j] = cor[i+j*natoms] = val;
      in_section_corr = 1;
    } else if (pdbrec("ATOM", rec)) {
      if (in_section_corr) {
	fclose(fp);
	abnormalexit("CORR record in ATOM section", 0);
      }
      k = natoms++;
      atom[k].id = k;
      atom[k].asnum = read_int(rec+PDB_ATOM_serial_LOC, PDB_ATOM_serial_LEN);
      atom[k].sp = site + k;
      atom[k].sp->next = 0;
      sp = atom[k].sp;
      atom[k].naltlocs = 0;
      if (strlen(rec) > PDBM_EXTRA_z_LOC) {
	sp->p.x = read_double(rec+PDBM_EXTRA_x_LOC, PDBM_EXTRA_x_LEN);
	sp->p.y = read_double(rec+PDBM_EXTRA_y_LOC, PDBM_EXTRA_y_LEN);
	sp->p.z = read_double(rec+PDBM_EXTRA_z_LOC, PDBM_EXTRA_z_LEN);
      } else {
	sp->p.x = read_double(rec+PDB_ATOM_x_LOC, PDB_ATOM_x_LEN);
	sp->p.y = read_double(rec+PDB_ATOM_y_LOC, PDB_ATOM_y_LEN);
	sp->p.z = read_double(rec+PDB_ATOM_z_LOC, PDB_ATOM_z_LEN);
      }
      sp->occ = read_double(rec+PDB_ATOM_occ_LOC, PDB_ATOM_occ_LEN);
      sp->Biso = read_double(rec+PDB_ATOM_tfac_LOC, PDB_ATOM_tfac_LEN);
      read_string(rec+PDB_ATOM_name_LOC, PDB_ATOM_name_LEN, atom[k].name);
      read_string(rec+PDB_ATOM_resName_LOC, PDB_ATOM_resName_LEN,
		  atom[k].resname);
      atom[k].resseq = read_int(rec+PDB_ATOM_resSeq_LOC, PDB_ATOM_resSeq_LEN);
    }
  }    
  return natoms;
}

void
print_conformer(pdbAtom *templ, Real *xp, Real *sd, Real *op, int natoms,
		int iconfm, int biso, int extraxyz,
		int nomit, landmark *alist, Real prob, void *work)
{
  int i, j, k, iatom, fwidth, fprec;
  Matrix3 X;
  pdbAtom *ap;
  Real s;
  char *omit, buf[BUFSIZ], buf2[BUFSIZ];

  if (iconfm > 0) printf("MODEL " "    " "%4d\n", iconfm);
  if (biso) s = 8.0*M_PI*M_PI;
  
  if (extraxyz) {
    fwidth = PDBM_EXTRA_x_LEN - 1;
    fprec = fwidth - 7;
  }

  omit = work;
  memset(omit, 0, natoms);
  if (nomit > 0) {
    for (k = 0; k < nomit; k++)
      if (drand48() < prob)
	for (i = 0; i < alist[k].nranges; i++)
	  for (j = alist[k].range[i].lower-1; j < alist[k].range[i].upper; j++)
	    omit[j] = 1;
  }

  ap = templ;
  X = (Matrix3) xp;

  if (!extraxyz) {
    for (iatom = 0; iatom < natoms; iatom++, ap++)
      if (!omit[ap->id])
	printf("ATOM  %5d %-4s %-3s A%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
	       ap->asnum, ap->name, ap->resname, ap->resseq,
	       X[iatom][0], X[iatom][1], X[iatom][2],
	       op[iatom], biso ? sd[iatom]*sd[iatom]*s : sd[iatom]);
  } else {
    for (iatom = 0; iatom < natoms; iatom++, ap++)
      if (!omit[ap->id]) {
	sprintf(buf, "ATOM  %5d %-4s %-3s A%4d    %8.3f%8.3f%8.3f%6.2f%6.2f",
		ap->asnum, ap->name, ap->resname, ap->resseq,
		X[iatom][0], X[iatom][1], X[iatom][2],
		op[iatom], biso ? sd[iatom]*sd[iatom]*s : sd[iatom]);
	sprintf(buf2, "@%*s %*.*e %*.*e %*.*e",
		PDBM_EXTRA_ID_LEN, "",
		fwidth, fprec, X[iatom][0],
		fwidth, fprec, X[iatom][1],
		fwidth, fprec, X[iatom][2]);
	printf("%-80s%s\n", buf, buf2);
      }
  }
  printf("TER\n");
  if (iconfm > 0) printf("ENDMDL\n");
  return;
}
