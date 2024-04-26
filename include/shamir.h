#ifndef _SHAMIR_H_
#define _SHAMIR_H_

#include <vector>
#include <string>

int split(int opt_threshold, int opt_number, const std::string &secret, std::vector<std::string> &shares);

int combine(const std::vector<std::string> secret_shares, std::string &secret);

/**
u32 ETCALL EtCreateShares(IN u32 u32Threshold, IN u32 u32ShareNum, IN const u8* pu8Secret, u32 u32SecretLen, OUT char* pszShares, OUT u32* pu32SharesLen);
u32 ETCALL EtGetSecret(IN char* pszShares, OUT u8* pu8Secret, OUT u32* pu32SecretLen);
*/


#endif