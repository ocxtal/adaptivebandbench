/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#ifndef _PARASAIL_MATRIX_LOOKUP_H_
#define _PARASAIL_MATRIX_LOOKUP_H_

#include "matrices/blosum100.h"
#include "matrices/blosum30.h"
#include "matrices/blosum35.h"
#include "matrices/blosum40.h"
#include "matrices/blosum45.h"
#include "matrices/blosum50.h"
#include "matrices/blosum55.h"
#include "matrices/blosum60.h"
#include "matrices/blosum62.h"
#include "matrices/blosum65.h"
#include "matrices/blosum70.h"
#include "matrices/blosum75.h"
#include "matrices/blosum80.h"
#include "matrices/blosum85.h"
#include "matrices/blosum90.h"
#include "matrices/blosumn.h"
#include "matrices/pam10.h"
#include "matrices/pam100.h"
#include "matrices/pam110.h"
#include "matrices/pam120.h"
#include "matrices/pam130.h"
#include "matrices/pam140.h"
#include "matrices/pam150.h"
#include "matrices/pam160.h"
#include "matrices/pam170.h"
#include "matrices/pam180.h"
#include "matrices/pam190.h"
#include "matrices/pam20.h"
#include "matrices/pam200.h"
#include "matrices/pam210.h"
#include "matrices/pam220.h"
#include "matrices/pam230.h"
#include "matrices/pam240.h"
#include "matrices/pam250.h"
#include "matrices/pam260.h"
#include "matrices/pam270.h"
#include "matrices/pam280.h"
#include "matrices/pam290.h"
#include "matrices/pam30.h"
#include "matrices/pam300.h"
#include "matrices/pam310.h"
#include "matrices/pam320.h"
#include "matrices/pam330.h"
#include "matrices/pam340.h"
#include "matrices/pam350.h"
#include "matrices/pam360.h"
#include "matrices/pam370.h"
#include "matrices/pam380.h"
#include "matrices/pam390.h"
#include "matrices/pam40.h"
#include "matrices/pam400.h"
#include "matrices/pam410.h"
#include "matrices/pam420.h"
#include "matrices/pam430.h"
#include "matrices/pam440.h"
#include "matrices/pam450.h"
#include "matrices/pam460.h"
#include "matrices/pam470.h"
#include "matrices/pam480.h"
#include "matrices/pam490.h"
#include "matrices/pam50.h"
#include "matrices/pam500.h"
#include "matrices/pam60.h"
#include "matrices/pam70.h"
#include "matrices/pam80.h"
#include "matrices/pam90.h"
#include "matrices/nuc44.h"
#include "matrices/dnafull.h"
#include "matrices/blosum_map.h"
#include "matrices/pam_map.h"

#ifdef __cplusplus
extern "C" {
#endif

static const parasail_matrix_t * parasail_matrices[] = {
    &parasail_blosum100,
    &parasail_blosum30,
    &parasail_blosum35,
    &parasail_blosum40,
    &parasail_blosum45,
    &parasail_blosum50,
    &parasail_blosum55,
    &parasail_blosum60,
    &parasail_blosum62,
    &parasail_blosum65,
    &parasail_blosum70,
    &parasail_blosum75,
    &parasail_blosum80,
    &parasail_blosum85,
    &parasail_blosum90,
    &parasail_blosumn,
    &parasail_pam10,
    &parasail_pam100,
    &parasail_pam110,
    &parasail_pam120,
    &parasail_pam130,
    &parasail_pam140,
    &parasail_pam150,
    &parasail_pam160,
    &parasail_pam170,
    &parasail_pam180,
    &parasail_pam190,
    &parasail_pam20,
    &parasail_pam200,
    &parasail_pam210,
    &parasail_pam220,
    &parasail_pam230,
    &parasail_pam240,
    &parasail_pam250,
    &parasail_pam260,
    &parasail_pam270,
    &parasail_pam280,
    &parasail_pam290,
    &parasail_pam30,
    &parasail_pam300,
    &parasail_pam310,
    &parasail_pam320,
    &parasail_pam330,
    &parasail_pam340,
    &parasail_pam350,
    &parasail_pam360,
    &parasail_pam370,
    &parasail_pam380,
    &parasail_pam390,
    &parasail_pam40,
    &parasail_pam400,
    &parasail_pam410,
    &parasail_pam420,
    &parasail_pam430,
    &parasail_pam440,
    &parasail_pam450,
    &parasail_pam460,
    &parasail_pam470,
    &parasail_pam480,
    &parasail_pam490,
    &parasail_pam50,
    &parasail_pam500,
    &parasail_pam60,
    &parasail_pam70,
    &parasail_pam80,
    &parasail_pam90,
    &parasail_nuc44,
    &parasail_dnafull,
    NULL
};

#ifdef __cplusplus
}
#endif

#endif /* _PARASAIL_MATRIX_LOOKUP_H_ */

