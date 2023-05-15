#ifndef DATAOWNER_H
#define DATAOWNER_H

#include "util.h"

class DataOwner{
    private:
    string filename;
    ZZ_p sk;

    int N;                                  // the number of tuples in rR
    int d;                                  // d = length_of_bytes(length_of_bytes+1)/2
    int attribute_index;                        // 0 ~ m-1 : arbitrary designation

    public:
    map<int, vector<string>> relation_schema;
    map<int, string> rR;                        // instansiation of relation schema
    deque<int> update_index_list;

    void load_data();

    void SetFileName(string& name);

    void SetAttribute(int& num);

    string readFileIntoString(const string& path);

    void instantiate_relation_schema();

    void Kgen(public_key& pk);

    void init_bps(bilinear_pairings_parameters& bps);

    void compute_crs(public_key& pk, const ZZ_p& alpha);

    void ZZ_from_element(ZZ& z, element_t& e);

    void element_from_ZZ(element_t& e, ZZ& z);

    void convert_to_alpha_pbc(public_key& pk, element_t& alpha_pbc);

    void Setup(public_key& pk, authentication_information& auth_info, element_t& root);

    void expand_N();

    void expand_rR();

    void polynomialization(vector<ZZ_pX>& G_poly, map<int, string>& rR, public_key& pk);

    void polynomialization(ZZ_pX& G_poly, string& tuple, public_key& pk);

    void make_prefixes(vector<string>& prefixes, string& s);

    void BuildEBT(authentication_information& auth_info, public_key& pk, vector<ZZ_pX>& G_poly);

    void compute_evaluation_binding_value(vector<element_t>& ebv, public_key& pk, vec_ZZ_p& G_eval);

    void make_Tree(vector<Node>& Tree, vector<element_t>& ebv, public_key& pk);

    bool check_Tree(vector<Node>& Tree, public_key& pk);

    void Update_Insertion(string& tuple, public_key& pk, authentication_information& auth_info, element_t& root);

    void Update_Deletion(string& tuple, public_key& pk, authentication_information& auth_info, element_t& root);

    void element_from_ZZ_p(element_t& e, ZZ_p& zp);
};

#endif