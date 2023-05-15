#ifndef UTIL_H
#define UTIL_H

#include <iostream>
#include <chrono>
#include <vector>
#include <string>
#include <algorithm>
#include <string.h>
#include <stdio.h>
#include <sstream>
#include <map>
#include <fstream>
#include <cmath>
#include <deque>

#include <pbc/pbc.h>

#include <iomanip>
#include <openssl/evp.h>
#include <openssl/sha.h>

#include<NTL/ZZ_p.h>
#include<NTL/ZZ_pX.h>
#include<NTL/vec_ZZ_p.h>

using namespace std;
using namespace NTL;

#define q_size 256          // order of q in NTL order of r in PBC. 256bit prime. 256bit hash
#define length_of_bytes 32  // byte length used to ZZ & element

struct proof{
    int index;
    element_t nu;
    element_t s;
    element_t t;
};

typedef struct bilinear_pairings_parameters{
    ZZ q;
    pairing_t e;
    element_t g1, g2;
} bp;

struct public_key{
    int d;
    bp bps;
    vector<element_t> crs_g1_pbc;
    vector<element_t> crs_g2_pbc;
    void hash(vec_ZZ_p& hash_value, vector<string>& prefixes);
    void hash(ZZ_p& hash_value, string& pattern);
    string bytes_to_hex_string(const vector<uint8_t>& bytes);
    void HexToDec(ZZ_p& ctxDec, string hex);
};

struct Node {
    element_t digest;
    bool isMatch;
    int height;
    int index;
    struct Node *parent;
    struct Node *left, *right;
};

typedef struct authentication_information{
    map<int, string> rR;
    vector<element_t> evaluation_binding_value;
    vector<Node> Tree;
} authentication_information;

struct query_value{
    vector<int> matches;
    vector<int> covers;
    vector<string> tuple;
};

struct query_proof{
    vector<proof> proof_matches;
    vector<proof> proof_mismatches;
    vector<element_t> w;
};

struct query_result{
    struct query_value qv;
    struct query_proof qp;
};

template <typename S, typename V>
void print_map(map<S, V>& m);

void print_vector_pbc(vector<element_t>& v, int size);

template <typename T>
void print_vector(vector<T>& v);

#endif