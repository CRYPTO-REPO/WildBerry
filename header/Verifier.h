#ifndef VERIFIER_H
#define VERIFIER_H

#include "util.h"

using namespace std;
using namespace NTL;

class Verifier {
    private:
    string pattern;

    public:
    void SetPattern(string& pattern);

    bool Verify(public_key& pk, query_result& QR, element_t& root);

    bool VerifyEBT(public_key& pk, query_value& qv, query_proof& qp, element_t& root);

    void polynomialization(vector<ZZ_pX>& G_poly, vector<string>& tuple, public_key& pk);

    void make_prefixes(vector<string>& prefixes, string& s);

    void evaluate_polynomial_g1(element_t& result, ZZ_pX& poly, public_key& pk);

    void evaluate_polynomial_g2(element_t& result, ZZ_pX& poly, public_key& pk);

    void element_from_ZZ_p(element_t& e, ZZ_p& zp);
};
#endif