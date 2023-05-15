#include "Prover.h"

void Prover::Query(query_result& QR, public_key& pk, string& pattern, authentication_information& auth_info){
    vector<ZZ_pX> G_poly;
    ZZ_p hash_pattern;
    ZZ_p G_eval;
    
    // Find all matches
    polynomialization(G_poly, auth_info.rR, pk);

    pk.hash(hash_pattern, pattern);

    for(int itr = 0; itr < auth_info.rR.size(); itr++){
        eval(G_eval, G_poly[itr], hash_pattern);
        if(G_eval == 0){
            QR.qv.matches.push_back(itr);
            QR.qv.tuple.push_back(auth_info.rR[itr]);
            auth_info.Tree[itr].isMatch = true;
        }
    }

    // Find mismatch nodes
    for(int itr = auth_info.Tree.size()/2 + 1; itr < auth_info.Tree.size(); itr++){
        auth_info.Tree[itr].isMatch = auth_info.Tree[itr].left->isMatch || auth_info.Tree[itr].right->isMatch;
    }
    
    // Find cover nodes
    find_covers(QR.qv.covers, auth_info);
    
    // Proofs for matches
    QR.qp.proof_matches.resize(QR.qv.matches.size());
    for(int itr = 0; itr < QR.qv.matches.size(); itr++){
        QR.qp.proof_matches[itr].index = QR.qv.matches[itr];
        ProveEBT(QR.qp.proof_matches[itr], pk, G_poly[QR.qv.matches[itr]], auth_info.rR, auth_info, hash_pattern);
    }
    
    vector<element_t> tmp(QR.qv.covers.size());
    for(int itr = 0; itr < tmp.size(); itr++){
        element_init_G1(tmp[itr], pk.bps.e);
        element_set(tmp[itr], auth_info.Tree[QR.qv.covers[itr]].digest);
    }
    QR.qp.w.swap(tmp);

    // Proofs for mismatches
    deque<Node*> dq;
    ZZ_pX G_poly_covers = ZZ_pX(INIT_MONO, 0, 0);
    QR.qp.proof_mismatches.resize(QR.qv.covers.size());
    for(int itr = 0; itr < QR.qv.covers.size(); itr++){
        QR.qp.proof_mismatches[itr].index = QR.qv.covers[itr];
        
        if(auth_info.Tree[QR.qv.covers[itr]].height > 0){
            Node* cur = &auth_info.Tree[QR.qv.covers[itr]];
            dq.push_back(cur);
            while(cur->left != NULL && cur->right != NULL){
                dq.push_back(cur->left);
                dq.push_back(cur->right);
                dq.pop_front();
                cur = dq.front();
            }
            for(int itr = 0; itr < dq.size(); itr++){
                G_poly_covers += G_poly[dq[itr]->index];
            }
        }
        else{
            G_poly_covers = G_poly[QR.qv.covers[itr]];
        }
        ProveEBT_for_mismatches(QR.qp.proof_mismatches[itr], pk, G_poly_covers, auth_info.rR, auth_info, hash_pattern);
        dq.clear();
    }

    // Reset isMatch
    reset_isMatch(auth_info);
}

void Prover::reset_isMatch(authentication_information& auth_info){
    for(int itr = 0; itr < auth_info.Tree.size(); itr++){
        auth_info.Tree[itr].isMatch = false;
    }
}

void Prover::polynomialization(vector<ZZ_pX>& G_poly, map<int, string>& rR, public_key& pk){
    int N = rR.size();
    string s = "";
    vector<string> prefixes;

    for(map<int, string>::iterator itr = rR.begin(); itr != rR.end(); ++itr){
        s = itr->second;

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

void Prover::make_prefixes(vector<string>& prefixes, string& s){
    int l = s.size();
    for(int itr_suf = l-1; itr_suf >= 0; --itr_suf){                // suffixes of a tuple
        for(int itr_pre = 1; itr_pre <= l - itr_suf; ++itr_pre){    // prefixes of a suffix
            prefixes.push_back(s.substr(itr_suf, itr_pre));
        }
    }

    sort(prefixes.begin(), prefixes.end());                         // Remove duplicate
    prefixes.erase(unique(prefixes.begin(), prefixes.end()), prefixes.end());
}

void Prover::ProveEBT(proof& proof, public_key& pk, ZZ_pX& poly, map<int, string>& rR, authentication_information& auth_info, ZZ_p& hash_pattern){

    element_init_G1(proof.nu, pk.bps.e);
    evaluate_polynomial_g1(proof.nu, poly, pk);

    ZZ_pX S, T;
    ZZ_pX poly_pattern = ZZ_pX(INIT_MONO, 1, 1);    // make X-p
    SetCoeff(poly_pattern, 0, -hash_pattern); 

    XGCD(poly_pattern, S, T, poly, poly_pattern);   // find S,T

    element_init_G1(proof.s, pk.bps.e);
    element_init_G2(proof.t, pk.bps.e);

    evaluate_polynomial_g1(proof.s, S, pk);
    evaluate_polynomial_g2(proof.t, T, pk);
}

void Prover::ProveEBT_for_mismatches(proof& proof, public_key& pk, ZZ_pX& poly, map<int, string>& rR, authentication_information& auth_info, ZZ_p& hash_pattern){

    element_init_G2(proof.nu, pk.bps.e);
    evaluate_polynomial_g2(proof.nu, poly, pk);

    ZZ_pX S, T;
    ZZ_pX poly_pattern = ZZ_pX(INIT_MONO, 1, 1);    // make X-p
    ZZ_pX poly_1 = ZZ_pX(INIT_MONO, 0, 1);
    SetCoeff(poly_pattern, 0, -(hash_pattern));

    XGCD(poly_1, S, T, poly, poly_pattern);   // find S,T

    element_init_G1(proof.s, pk.bps.e);
    element_init_G2(proof.t, pk.bps.e);

    evaluate_polynomial_g1(proof.s, S, pk);
    evaluate_polynomial_g2(proof.t, T, pk);
}

void Prover::evaluate_polynomial_g1(element_t& result, ZZ_pX& poly, public_key& pk){
    element_t coeff, eval;
    element_init_Zr(coeff, pk.bps.e);
    element_init_G1(eval, pk.bps.e);
    
    for(int itr = 0; itr < poly.rep.length(); itr++){
        element_from_ZZ_p(coeff, poly[itr]);
        element_pow_zn(eval, pk.crs_g1_pbc[itr], coeff);

        element_mul(result, result, eval);
    }
}

void Prover::evaluate_polynomial_g2(element_t& result, ZZ_pX& poly, public_key& pk){
    element_t coeff, eval;
    element_init_Zr(coeff, pk.bps.e);
    element_init_G2(eval, pk.bps.e);
    
    for(int itr = 0; itr < poly.rep.length(); itr++){
        element_from_ZZ_p(coeff, poly[itr]);
        element_pow_zn(eval, pk.crs_g2_pbc[itr], coeff);

        element_mul(result, result, eval);
    }
}

void Prover::element_from_ZZ_p(element_t& e, ZZ_p& zp){ // convert from ZZ to element
    unsigned char bytes[length_of_bytes] = {0,};
    ZZ z = conv<ZZ>(zp);
    BytesFromZZ(bytes, z, length_of_bytes);
    reverse(bytes, bytes + length_of_bytes);
    element_from_bytes(e, bytes);
}

void Prover::find_covers(vector<int>& covers, authentication_information& auth_info){
    // level-order traversal
    deque<Node*> dq;
    vector<int> tmp;

    Node* cur = &(auth_info.Tree[auth_info.Tree.size()-1]);
    dq.push_back(cur);

    if(cur->isMatch == false){  // no match => root's isMatch = false
        covers.push_back(auth_info.Tree.size()-1);
        return;
    }

    for(int itr =  auth_info.Tree.size()-1; !dq.empty(); itr--){
        if(cur->left == NULL && cur->right == NULL){    // leaf
            if(cur->isMatch == 0){      // cover node
                tmp.push_back(itr);
            }
        }
        else{                                           // internal
            if(cur->isMatch == 0){      // cover node
                tmp.push_back(itr);
            }
            dq.push_back(cur->right);
            dq.push_back(cur->left);
        }
        dq.pop_front();
        if(!dq.empty()){
            cur = dq.front();
        }
    }
    
    for(int itr = 0; itr < tmp.size(); itr++){
        if(auth_info.Tree[tmp[itr]].parent->isMatch != 0){
            covers.push_back(tmp[itr]);
        }
    }
    reverse(covers.begin(), covers.end());
}
