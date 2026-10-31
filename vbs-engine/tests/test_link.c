/* Toolchain smoke test: prove we link MRTKLIB and can call into it. */
#include <stdio.h>
#include "rtklib.h"

int main(void) {
	gtime_t t = epoch2time((double[]){2026, 7, 29, 0, 0, 0});
	char s[64];
	time2str(t, s, 0);
	printf("epoch2time/time2str: %s\n", s);
	printf("satno(SYS_GAL,12) = %d\n", satno(SYS_GAL, 12));
	printf("obs2code(\"7Q\") = %d, code2obs(%d) = %s\n",
	    obs2code("7Q"), obs2code("7Q"), code2obs(obs2code("7Q")));
	return 0;
}
