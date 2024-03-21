#include "scoring.h"

void load_scoring_matrix(int scoring_matrix[][L_NT], char* scoreFilename) {
	exscoreset*     xss = NULL;
	u8*             rowChars = (u8*)"ACGT";
	u8*             colChars = (u8*)"ACGT";
	u8*             r, *c;
	int             x, y;

	xss = read_score_set_by_name (scoreFilename);

	for (r = rowChars, x = 0; *r != 0; r++, x++) {
		for (c = colChars, y = 0; *c != 0; c++, y++) {
			scoring_matrix[x][y]  = (int) ((scoreset*) xss)->sub[*r][*c];
		}
	}

	if (xss != NULL) {
		free_if_valid("score set", xss);
		xss = NULL;
	}
}
