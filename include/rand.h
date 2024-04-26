#ifndef _RAND_H_
#define _RAND_H_

#include <stdint.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
Rand Public API

	rand_bytes

*/

int rand_bytes(uint8_t *buf, size_t buflen);


#ifdef __cplusplus
}
#endif
#endif

