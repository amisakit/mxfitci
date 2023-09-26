#ifndef _PDBCONST_H
#define _PDBCONST_H

#define PDB_HEADER_IDcode_LOC 62
#define PDB_HEADER_IDcode_LEN 4
#define PDB_MODEL_serial_LOC 10
/* #define PDB_MODEL_serial_LEN 4 */
#define PDB_MODEL_serial_LEN 10
#define PDB_ATOM_serial_LOC 6
#define PDB_ATOM_serial_LEN 5
#define PDB_ATOM_name_LOC 12
#define PDB_ATOM_name_LEN 4
#define PDB_ATOM_altLoc_LOC 16
#define PDB_ATOM_resName_LOC 17
#define PDB_ATOM_resName_LEN 3
#define PDB_ATOM_chainID_LOC 21
#define PDB_ATOM_resSeq_LOC 22
#define PDB_ATOM_resSeq_LEN 4
#define PDB_ATOM_x_LOC 30
#define PDB_ATOM_x_LEN 8
#define PDB_ATOM_y_LOC 38
#define PDB_ATOM_y_LEN 8
#define PDB_ATOM_z_LOC 46
#define PDB_ATOM_z_LEN 8
#define PDB_ATOM_occ_LOC 54
#define PDB_ATOM_occ_LEN 6
#define PDB_ATOM_tfac_LOC 60
#define PDB_ATOM_tfac_LEN 6

#define PDBM_EXTRA_FLAG_LOC 80
#define PDBM_EXTRA_ID_LOC 81
#define PDBM_EXTRA_ID_LEN 5
#define PDBM_EXTRA_x_LOC 86
#define PDBM_EXTRA_x_LEN 25
#define PDBM_EXTRA_y_LOC 111
#define PDBM_EXTRA_y_LEN 25
#define PDBM_EXTRA_z_LOC 136
#define PDBM_EXTRA_z_LEN 25


/*
#define PDB_REC_ATOM	0x41544f4dU
#define PDB_REC_HETATM	0x48455441544dU
#define PDB_REC_HEADER	0x484541444552U
#define PDB_REC_END	0x454e44U
#define PDB_REC_ENDMDL	0x454e444d444cU
#define PDB_REC_MODEL	0x4d4f44454cU
#define PDB_REC_TER	0x544552U
*/

#endif
