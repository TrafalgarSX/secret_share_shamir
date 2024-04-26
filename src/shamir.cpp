#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <iostream>
#include <vector>
#include <string>
#include <rand.h>
#include <unordered_set>

#include <gmp.h>

#define MAXDEGREE 1024
#define MAXLINELEN (1 + 10 + 1 + MAXDEGREE / 4 + 10)
#define MAXSHARELEN (MAXDEGREE / 4 + 1 + MAXDEGREE / 4 + 1)

/* coefficients of some irreducible polynomials over GF(2) */
static const unsigned char irred_coeff[] = {
    4,  3,  1,  5,  3,  1,  4,  3,  1,  7,  3,  2,  5,  4,  3,  5,  3,  2,  7,
    4,  2,  4,  3,  1,  10, 9,  3,  9,  4,  2,  7,  6,  2,  10, 9,  6,  4,  3,
    1,  5,  4,  3,  4,  3,  1,  7,  2,  1,  5,  3,  2,  7,  4,  2,  6,  3,  2,
    5,  3,  2,  15, 3,  2,  11, 3,  2,  9,  8,  7,  7,  2,  1,  5,  3,  2,  9,
    3,  1,  7,  3,  1,  9,  8,  3,  9,  4,  2,  8,  5,  3,  15, 14, 10, 10, 5,
    2,  9,  6,  2,  9,  3,  2,  9,  5,  2,  11, 10, 1,  7,  3,  2,  11, 2,  1,
    9,  7,  4,  4,  3,  1,  8,  3,  1,  7,  4,  1,  7,  2,  1,  13, 11, 6,  5,
    3,  2,  7,  3,  2,  8,  7,  5,  12, 3,  2,  13, 10, 6,  5,  3,  2,  5,  3,
    2,  9,  5,  2,  9,  7,  2,  13, 4,  3,  4,  3,  1,  11, 6,  4,  18, 9,  6,
    19, 18, 13, 11, 3,  2,  15, 9,  6,  4,  3,  1,  16, 5,  2,  15, 14, 6,  8,
    5,  2,  15, 11, 2,  11, 6,  2,  7,  5,  3,  8,  3,  1,  19, 16, 9,  11, 9,
    6,  15, 7,  6,  13, 4,  3,  14, 13, 3,  13, 6,  3,  9,  5,  2,  19, 13, 6,
    19, 10, 3,  11, 6,  5,  9,  2,  1,  14, 3,  2,  13, 3,  1,  7,  5,  4,  11,
    9,  8,  11, 6,  5,  23, 16, 9,  19, 14, 6,  23, 10, 2,  8,  3,  2,  5,  4,
    3,  9,  6,  4,  4,  3,  2,  13, 8,  6,  13, 11, 1,  13, 10, 3,  11, 6,  5,
    19, 17, 4,  15, 14, 7,  13, 9,  6,  9,  7,  3,  9,  7,  1,  14, 3,  2,  11,
    8,  2,  11, 6,  4,  13, 5,  2,  11, 5,  1,  11, 4,  1,  19, 10, 3,  21, 10,
    6,  13, 3,  1,  15, 7,  5,  19, 18, 10, 7,  5,  3,  12, 7,  2,  7,  5,  1,
    14, 9,  6,  10, 3,  2,  15, 13, 12, 12, 11, 9,  16, 9,  7,  12, 9,  3,  9,
    5,  2,  17, 10, 6,  24, 9,  3,  17, 15, 13, 5,  4,  3,  19, 17, 8,  15, 6,
    3,  19, 6,  1};

#define mpz_lshift(A, B, l) mpz_mul_2exp(A, B, l)
#define mpz_sizeinbits(A) (mpz_cmp_ui(A, 0) ? mpz_sizeinbase(A, 2) : 0)

int field_size_valid(int deg) {
  return (deg >= 8) && (deg <= MAXDEGREE) && (deg % 8 == 0);
}

/* initialize 'poly' to a bitfield representing the coefficients of an
   irreducible polynomial of degree 'deg' */
void field_init(mpz_t *poly, int deg) {
  assert(field_size_valid(deg));
  mpz_init_set_ui(*poly, 0);
  mpz_setbit(*poly, deg);
  mpz_setbit(*poly, irred_coeff[3 * (deg / 8 - 1) + 0]);
  mpz_setbit(*poly, irred_coeff[3 * (deg / 8 - 1) + 1]);
  mpz_setbit(*poly, irred_coeff[3 * (deg / 8 - 1) + 2]);
  mpz_setbit(*poly, 0);
}

void field_deinit(mpz_t* poly) { mpz_clear(*poly); }

/* I/O routines for GF(2^deg) field elements, hexmode only*/
int field_import(mpz_t x, const char *s, int degree) {

    if (strlen(s) > degree / 4){
    //   printf("%ld, degree %d\n", strlen(s), degree / 4);
    //   fatal("input string too long");
      return -1;
    }
    if (strlen(s) < degree / 4)
      return -1;
    //   warning("input string too short, adding null padding on the left");

    if (mpz_set_str(x, s, 16) || (mpz_cmp_ui(x, 0) < 0))
      return -1;
    //   fatal("invalid syntax");
    return 0;
}

std::string field_out_hex_str(const mpz_t num, int degree) {
  std::string out_str;

  for (int i = degree / 4 - mpz_sizeinbase(num, 16); i; i--){
    out_str += "0";
  }
  char *temp = mpz_get_str(NULL, 16, num);
  out_str += temp;
  free(temp);
  return out_str;
}

/* basic field arithmetic in GF(2^deg) */
void field_add(mpz_t z, const mpz_t x, const mpz_t y) { mpz_xor(z, x, y); }

void field_mult(mpz_t *poly, mpz_t z, const mpz_t x, const mpz_t y, int degree) {
  mpz_t b;
  unsigned int i;
  assert(z != y);
  mpz_init_set(b, x);
  if (mpz_tstbit(y, 0))
    mpz_set(z, b);
  else
    mpz_set_ui(z, 0);
  for (i = 1; i < degree; i++) {
    mpz_lshift(b, b, 1);
    if (mpz_tstbit(b, degree))
      mpz_xor(b, b, *poly);
    if (mpz_tstbit(y, i))
      mpz_xor(z, z, b);
  }
  mpz_clear(b);
}

void field_invert(mpz_t *poly, mpz_t z, const mpz_t x) {
  mpz_t u, v, g, h;
  int i;
  assert(mpz_cmp_ui(x, 0));
  mpz_init_set(u, x);
  mpz_init_set(v, *poly);
  mpz_init_set_ui(g, 0);
  mpz_set_ui(z, 1);
  mpz_init(h);
  while (mpz_cmp_ui(u, 1)) {
    i = mpz_sizeinbits(u) - mpz_sizeinbits(v);
    if (i < 0) {
      mpz_swap(u, v);
      mpz_swap(z, g);
      i = -i;
    }
    mpz_lshift(h, v, i);
    mpz_xor(u, u, h);
    mpz_lshift(h, g, i);
    mpz_xor(z, z, h);
  }
  mpz_clear(u);
  mpz_clear(v);
  mpz_clear(g);
  mpz_clear(h);
}

/* routines for the random number generator */
int cprng_read(mpz_t x, int degree) {
  uint8_t buf[MAXDEGREE / 8];
  unsigned int count;
  int i;
  rand_bytes(buf, degree / 8);
  mpz_import(x, degree / 8, 1, 1, 0, 0, buf);
  return 0;
}

int xcprng_read(int degree, mpz_t *x, int opt_number) {
  std::unordered_set<std::string> unique_numbers;
  uint8_t buf[MAXDEGREE / 8];

  for (int j = 0; j < opt_number; ++j) {
    rand_bytes(buf, degree / 8);
    mpz_init(x[j]);
    mpz_import(x[j], degree / 8, 1, 1, 0, 0, buf);

    // Convert the number to a string
    char* number_str = mpz_get_str(NULL, 10, x[j]);
    std::string number_string(number_str);
    free(number_str);

    // Check if the number is unique
    if (unique_numbers.find(number_string) != unique_numbers.end()) {
        // The number is not unique, decrement j to try again
        --j;
    } else {
        // The number is unique, add it to the set
        unique_numbers.insert(number_string);
    }
  }
  return 0;
}

/* a 64 bit pseudo random permutation (based on the XTEA cipher) */
void encipher_block(uint32_t *v) {
  uint32_t sum = 0, delta = 0x9E3779B9;
  int i;
  for (i = 0; i < 32; i++) {
    v[0] += (((v[1] << 4) ^ (v[1] >> 5)) + v[1]) ^ sum;
    sum += delta;
    v[1] += (((v[0] << 4) ^ (v[0] >> 5)) + v[0]) ^ sum;
  }
}

void decipher_block(uint32_t *v) {
  uint32_t sum = 0xC6EF3720, delta = 0x9E3779B9;
  int i;
  for (i = 0; i < 32; i++) {
    v[1] -= ((v[0] << 4 ^ v[0] >> 5) + v[0]) ^ sum;
    sum -= delta;
    v[0] -= ((v[1] << 4 ^ v[1] >> 5) + v[1]) ^ sum;
  }
}

void encode_slice(uint8_t *data, int idx, int len,
                  void (*process_block)(uint32_t *)) {
  uint32_t v[2];
  int i;
  for (i = 0; i < 2; i++)
    v[i] = data[(idx + 4 * i) % len] << 24 |
           data[(idx + 4 * i + 1) % len] << 16 |
           data[(idx + 4 * i + 2) % len] << 8 | data[(idx + 4 * i + 3) % len];
  process_block(v);
  for (i = 0; i < 2; i++) {
    data[(idx + 4 * i + 0) % len] = v[i] >> 24;
    data[(idx + 4 * i + 1) % len] = (v[i] >> 16) & 0xff;
    data[(idx + 4 * i + 2) % len] = (v[i] >> 8) & 0xff;
    data[(idx + 4 * i + 3) % len] = v[i] & 0xff;
  }
}

enum encdec { ENCODE, DECODE };

void encode_mpz(mpz_t x, enum encdec encdecmode, int degree) {
  uint8_t v[(MAXDEGREE + 8) / 16 * 2];
  size_t t;
  int i;
  memset(v, 0, (degree + 8) / 16 * 2);
  mpz_export(v, &t, -1, 2, 1, 0, x);
  if (degree % 16 == 8)
    v[degree / 8 - 1] = v[degree / 8];
  if (encdecmode == ENCODE) /* 40 rounds are more than enough!*/
    for (i = 0; i < 40 * ((int)degree / 8); i += 2)
      encode_slice(v, i, degree / 8, encipher_block);
  else
    for (i = 40 * (degree / 8) - 2; i >= 0; i -= 2)
      encode_slice(v, i, degree / 8, decipher_block);
  if (degree % 16 == 8) {
    v[degree / 8] = v[degree / 8 - 1];
    v[degree / 8 - 1] = 0;
  }
  mpz_import(x, (degree + 8) / 16, -1, 2, 1, 0, v);
  assert(mpz_sizeinbits(x) <= degree);
}

/* evaluate polynomials efficiently */
void horner(mpz_t *poly, int n, mpz_t y, const mpz_t x, const mpz_t coeff[], int degree) {
  int i;
  mpz_set(y, x);
  for (i = n - 1; i; i--) {
    field_add(y, y, coeff[i]);
    field_mult(poly, y, y, x, degree);
  }
  field_add(y, y, coeff[0]);
}

/* calculate the secret from a set of shares solving a linear equation system */

#define MPZ_SWAP(A, B)                                                         \
  do {                                                                         \
    mpz_set(h, A);                                                             \
    mpz_set(A, B);                                                             \
    mpz_set(B, h);                                                             \
  } while (0)

// int restore_secret(mpz_t *poly, int n, mpz_t (*A)[n], mpz_t b[]) {
int restore_secret(mpz_t *poly, int n, int degree, mpz_t** AA, mpz_t b[]) {
  int i, j, k, found;
  mpz_t h;
  mpz_init(h);
  for (i = 0; i < n; i++) {
    if (!mpz_cmp_ui(AA[i][i], 0)) {
      for (found = 0, j = i + 1; j < n; j++)
        if (mpz_cmp_ui(AA[i][j], 0)) {
          found = 1;
          break;
        }
      if (!found)
        return -1;
      for (k = i; k < n; k++)
        MPZ_SWAP(AA[k][i], AA[k][j]);
      MPZ_SWAP(b[i], b[j]);
    }
    for (j = i + 1; j < n; j++) {
      if (mpz_cmp_ui(AA[i][j], 0)) {
        for (k = i + 1; k < n; k++) {
          field_mult(poly, h, AA[k][i], AA[i][j], degree);
          field_mult(poly, AA[k][j], AA[k][j], AA[i][i], degree);
          field_add(AA[k][j], AA[k][j], h);
        }
        field_mult(poly, h, b[i], AA[i][j], degree);
        field_mult(poly, b[j], b[j], AA[i][i], degree);
        field_add(b[j], b[j], h);
      }
    }
  }
  field_invert(poly, h, AA[n - 1][n - 1]);
  field_mult(poly, b[n - 1], b[n - 1], h, degree);
  mpz_clear(h);
  return 0;
}

int split(int opt_threshold, int opt_number, const std::string &secret, std::vector<std::string> &shares) {
  unsigned int fmt_len;
  mpz_t y;
  mpz_t *x = NULL;
  int degree, i;
  int opt_security = 0;
  mpz_t poly;
  mpz_t *coeff = NULL;

  coeff = (mpz_t*)malloc(opt_threshold * sizeof(mpz_t));
  x = (mpz_t*)malloc(opt_number * sizeof(mpz_t));
  opt_security = 4 * ((secret.size() + 1) & ~1) ;
  if (!field_size_valid(opt_security))
    return -1;
  //   fatal("security level invalid (secret too long?)");

  degree = opt_security;
  field_init(&poly, opt_security);

  mpz_init(coeff[0]);
  if(field_import(coeff[0], secret.c_str(), degree) != 0){
    return -1;
  }

  if (degree >= 64)
    encode_mpz(coeff[0], ENCODE, degree);

  for (i = 1; i < opt_threshold; i++) {
    mpz_init(coeff[i]);
    cprng_read(coeff[i], degree);
  }

  mpz_init(y);
  xcprng_read(degree, x, opt_number);
  for (i = 0; i < opt_number; i++) {
    horner(&poly, opt_threshold, y, x[i], (const mpz_t *)coeff, degree);
    shares.emplace_back(field_out_hex_str(x[i], degree) + "-" + field_out_hex_str(y, degree));
    mpz_clear(x[i]);
  }
  mpz_clear(y);

  for (i = 0; i < opt_threshold; i++)
    mpz_clear(coeff[i]);
  field_deinit(&poly);
  free(coeff);
  free(x);
  return 0;
}

int combine(const std::vector<std::string> secret_shares, std::string &secret) {
  int opt_threshold = secret_shares.size();
  mpz_t x;
  mpz_t *y = NULL;
  mpz_t **A = NULL;
  std::string a, b;
  int i, j;
  int degree;
  unsigned int s = 0;
  mpz_t poly;

  y = (mpz_t*)malloc(opt_threshold * sizeof(mpz_t));
  A = (mpz_t**)malloc(opt_threshold * sizeof(mpz_t*));
  for (int i = 0; i < opt_threshold; i++) {
    A[i] = (mpz_t*)malloc(opt_threshold * sizeof(mpz_t));
  }
  mpz_init(x);
  for (i = 0; i < opt_threshold; i++) {
    if(secret_shares[i].length() > MAXSHARELEN)
      return -1;

    if(secret_shares[i].find("-") == std::string::npos)
      return -1;

    a = secret_shares[i].substr(0, secret_shares[i].find('-'));
    b = secret_shares[i].substr(secret_shares[i].find('-') + 1);
    
    if (!s) {
      s = 4 * b.size();
      if (!field_size_valid(s))
        return -1;
        // fatal("share has illegal length");
      field_init(&poly, s);
      degree = s;
    } else if (s != 4 * b.size())
      return -1;
    //   fatal("shares have different security levels");

    mpz_init(x);
    field_import(x, a.c_str(), degree);

    mpz_init_set_ui(A[opt_threshold - 1][i], 1);
    for (j = opt_threshold - 2; j >= 0; j--) {
      mpz_init(A[j][i]);
      field_mult(&poly, A[j][i], A[j + 1][i], x, degree);
    }
    mpz_init(y[i]);
    field_import(y[i], b.c_str(), degree);
    field_mult(&poly, x, x, A[0][i], degree);
    field_add(y[i], y[i], x);
  }
  mpz_clear(x);
  if (restore_secret(&poly, opt_threshold, degree, A, y))
    return -1;
    // fatal("shares inconsistent. Perhaps a single share was used twice");

  if (degree >= 64)
      encode_mpz(y[opt_threshold - 1], DECODE, degree);
  else
      return -1;
    //   warning("security level too small for the diffusion layer");

  secret = field_out_hex_str(y[opt_threshold - 1], degree);

  for (i = 0; i < opt_threshold; i++) {
    for (j = 0; j < opt_threshold; j++)
      mpz_clear(A[i][j]);
    mpz_clear(y[i]);
  }

  for (int i = 0; i < opt_threshold; i++) {
    free(A[i]);
  }
  free(A);
  free(y);
  field_deinit(&poly);
  return 0;
}