#ifndef PROVER_H
#define PROVER_H

#include "util.h"

using namespace std;
using namespace NTL;

class Prover {
    private:

    public:
    void Query(query_result& QR, public_key& pk, string& pattern, authentication_information& auth_info);

    void reset_isMatch(authentication_information& auth_info);

    void polynomialization(vector<ZZ_pX>& G_poly, map<int, string>& rR, public_key& pk);

    void make_prefixes(vector<string>& prefixes, string& s);

    void ProveEBT(proof& proof, public_key& pk, ZZ_pX& poly, map<int, string>& rR, authentication_information& auth_info, ZZ_p& hash_pattern);

    void ProveEBT_for_mismatches(proof& proof, public_key& pk, ZZ_pX& poly, map<int, string>& rR, authentication_information& auth_info, ZZ_p& hash_pattern);

    void evaluate_polynomial_g1(element_t& result, ZZ_pX& poly, public_key& pk);

    void evaluate_polynomial_g2(element_t& result, ZZ_pX& poly, public_key& pk);

    void element_from_ZZ_p(element_t& e, ZZ_p& zp);

    void find_covers(vector<int>& covers, authentication_information& auth_info);
};
#endif