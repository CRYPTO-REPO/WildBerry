#include "Verifier.h"

void Verifier::SetPattern(string& pattern){
    try {
        if (this->pattern.length() > length_of_bytes) throw out_of_range("pattern length exceeded 32 bytes");
    } catch(const exception& error) {
        cout << error.what() << endl;
        abort();
    }
    this->pattern = pattern;
}

bool Verifier::Verify(public_key& pk, query_result& QR, element_t& root){
        element_t tmp_p1[2], tmp_p2[2], tmp_gt, prod;
        element_t hash_pattern_pbc, pattern_g1;
        ZZ_p hash_pattern;
        bool b = false;

        element_init_Zr(hash_pattern_pbc, pk.bps.e);
        element_init_G1(pattern_g1, pk.bps.e);

        element_init_G1(tmp_p1[0], pk.bps.e);
        element_init_G1(tmp_p1[1], pk.bps.e);
        element_init_G2(tmp_p2[0], pk.bps.e);
        element_init_G2(tmp_p2[1], pk.bps.e);
        element_init_GT(tmp_gt, pk.bps.e);
        element_init_GT(prod, pk.bps.e);

        pk.hash(hash_pattern, this->pattern);
        hash_pattern = -(hash_pattern);
        element_from_ZZ_p(hash_pattern_pbc, hash_pattern);
        element_pow_zn(pattern_g1, pk.bps.g1, hash_pattern_pbc);
        element_mul(pattern_g1, pk.crs_g1_pbc[1], pattern_g1);      // g2^(alpha-p)

        // match verification
        b = VerifyEBT(pk, QR.qv, QR.qp, root);
        if (b != true){
            cout << "VerifyEBT failed" << endl;
            return false;
        }

        // mismatch verification
        element_pairing(tmp_gt, pk.bps.g1, pk.bps.g2);
        for(int itr = 0; itr < QR.qp.proof_mismatches.size(); itr++){
            element_set(tmp_p1[0], QR.qp.proof_mismatches[itr].s);
            element_set(tmp_p1[1], pattern_g1);
            element_set(tmp_p2[0], QR.qp.proof_mismatches[itr].nu);
            element_set(tmp_p2[1], QR.qp.proof_mismatches[itr].t);

            element_prod_pairing(prod, tmp_p1, tmp_p2, 2);

            if(element_cmp(prod, tmp_gt)){
                cout << "Verify mismatch failed" << endl;
                return false;
            }
        }
        return true;
    }

bool Verifier::VerifyEBT(public_key& pk, query_value& qv, query_proof& qp, element_t& root){
    vector<ZZ_pX> G_poly;
    element_t G_eval_pbc, eval_root;
    element_init_G1(G_eval_pbc, pk.bps.e);
    element_init_G1(eval_root, pk.bps.e);

    // Test if proof for matches 'nu' is equal to evaluation of g1^G(alpha)
    polynomialization(G_poly, qv.tuple, pk);

    for(int itr = 0; itr < G_poly.size(); itr++){
        element_set0(G_eval_pbc);
        evaluate_polynomial_g1(G_eval_pbc, G_poly[itr], pk);
        
        if(element_cmp(qp.proof_matches[itr].nu, G_eval_pbc)){
            cout << "proof_matches is not equal to poly evaluation" << endl;
            return false;
        }
    }

    // Build Subtree T & Compute and Compare root
    for(int itr = 0; itr < qp.proof_matches.size(); itr++){
        element_mul(eval_root, eval_root, pk.crs_g1_pbc[1]);
        element_mul(eval_root, eval_root, qp.proof_matches[itr].nu);
    }
    for(int itr = 0; itr < qp.w.size(); itr++){
        element_mul(eval_root, eval_root, qp.w[itr]);
    }

    if(element_cmp(eval_root, root) != 0){      // not same
        cout << "verifyEBT is not equal to root" << endl;
        return false;
    }
    return true;
}

void Verifier::polynomialization(vector<ZZ_pX>& G_poly, vector<string>& tuple, public_key& pk){
    int N = tuple.size();
    string s = "";
    vector<string> prefixes;

    for(vector<string>::iterator itr = tuple.begin(); itr != tuple.end(); ++itr){
        s = *itr;

        make_prefixes(prefixes, s);         // Slice a tuple to all suffixes && Slice suffixes to all prefixes

        vec_ZZ_p hash_value;                // Hash prefixes using openssl sha3-256 and convert to ZZ_p value
        hash_value.SetLength(prefixes.size());
        pk.hash(hash_value, prefixes);
        
        ZZ_pX F = ZZ_pX(INIT_MONO, 0, 1);   // Caculate polynomial G(X)
        for(auto value : hash_value){
            ZZ_pX tmp = ZZ_pX(INIT_MONO, 1, 1);
            SetCoeff(tmp, 0, -value);
            mul(F, F, tmp);
        }

        G_poly.push_back(F);

        prefixes.clear();
    }
}

void Verifier::make_prefixes(vector<string>& prefixes, string& s){
    int l = s.size();
    for(int itr_suf = l-1; itr_suf >= 0; --itr_suf){                // suffixes of a tuple
        for(int itr_pre = 1; itr_pre <= l - itr_suf; ++itr_pre){    // prefixes of a suffix
            prefixes.push_back(s.substr(itr_suf, itr_pre));
        }
    }

    sort(prefixes.begin(), prefixes.end());                         // Remove duplicate
    prefixes.erase(unique(prefixes.begin(), prefixes.end()), prefixes.end());
}

void Verifier::evaluate_polynomial_g1(element_t& result, ZZ_pX& poly, public_key& pk){
    element_t coeff, eval;
    element_init_Zr(coeff, pk.bps.e);
    element_init_G1(eval, pk.bps.e);
    
    for(int itr = 0; itr < poly.rep.length(); itr++){
        element_from_ZZ_p(coeff, poly[itr]);
        element_pow_zn(eval, pk.crs_g1_pbc[itr], coeff);

        element_mul(result, result, eval);
    }
}

void Verifier::evaluate_polynomial_g2(element_t& result, ZZ_pX& poly, public_key& pk){
    element_t coeff, eval;
    element_init_Zr(coeff, pk.bps.e);
    element_init_G2(eval, pk.bps.e);
    
    for(int itr = 0; itr < poly.rep.length(); itr++){
        element_from_ZZ_p(coeff, poly[itr]);
        element_pow_zn(eval, pk.crs_g2_pbc[itr], coeff);

        element_mul(result, result, eval);
    }
}

void Verifier::element_from_ZZ_p(element_t& e, ZZ_p& zp){ // convert from ZZ to element
    unsigned char bytes[length_of_bytes] = {0,};
    ZZ z = conv<ZZ>(zp);
    BytesFromZZ(bytes, z, length_of_bytes);
    reverse(bytes, bytes + length_of_bytes);
    element_from_bytes(e, bytes);
}
