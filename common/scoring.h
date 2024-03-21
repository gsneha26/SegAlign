#ifndef scoring_H                             // (prevent multiple inclusion)
#define scoring_H

#ifdef __cplusplus
extern "C" {
#endif

#include "dna_utilities.h"
#include "parameters.h"

void load_scoring_matrix(int scoring_matrix[][L_NT], char* scoreFilename);


#ifdef __cplusplus
}
#endif

#endif // scoring_H
