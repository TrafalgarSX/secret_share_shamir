
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <stdint.h>
#include <windows.h>
#include <wincrypt.h>


int rand_bytes(uint8_t *buf, size_t len)
{
	HCRYPTPROV hCryptProv;
	int ret = -1;

	if (!buf) {
		return -1;
	}
	if (len > INT_MAX) {
		return -1;
	}
	if (CryptAcquireContextA(&hCryptProv, NULL, NULL, PROV_RSA_FULL,
		CRYPT_VERIFYCONTEXT|CRYPT_SILENT) != TRUE) {
		return -1;
	}
	if (CryptGenRandom(hCryptProv, (DWORD)len, buf) != TRUE) {
		goto end;
	}
	ret = 1;
end:
	if (CryptReleaseContext(hCryptProv, 0) != TRUE) {
		ret = -1;
	}
	return ret;
}
