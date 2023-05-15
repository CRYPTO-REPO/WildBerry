#include "DataOwner.h"

void DataOwner::load_data(){ // Load data set to relation_schema
    string dataset_string;
    char delimiter = ',';
    dataset_string = readFileIntoString(this->filename); // Save contents of filename's data set to dataset_string as string

    istringstream sstream(dataset_string);
    vector<string> row;
    string record;

    int counter = 0;
    while (getline(sstream, record, '\n')) {        // Import a row in dataset_string
        istringstream line(record);

        while (getline(line, record, delimiter)) {  // Import a word from a row & Push it to vector
            if (record[record.size()-1] == '\r'){
                record = record.substr(0, record.size()-1);
            }

            try {
                if (record.size() > length_of_bytes) throw out_of_range("tuple length exceeded 32 bytes");
            } catch(const exception& error) {
                cout << error.what() << endl;
                abort();
            }

            row.push_back(record);
        }

        this->relation_schema[counter] = row;

        row.clear();
        line.clear();
        counter += 1;
    }
}

void DataOwner::SetFileName(string& name){
    this->filename = name;
}

void DataOwner::SetAttribute(int& num){
    this->attribute_index = num;
}

string DataOwner::readFileIntoString(const string& path) { // Import all file contents of csv as string
    auto ss = ostringstream{};
    ifstream input_file(path);
    if (!input_file.is_open()) {
        cerr << "Could not open the file - '" << path << "'" << endl;
        exit(EXIT_FAILURE);
    }
    ss << input_file.rdbuf();
    return ss.str();
}

void DataOwner::instantiate_relation_schema(){ // Construct rR
    int counter = 0;
    for (map<int, vector<string>>::iterator itr = this->relation_schema.begin(); itr != this->relation_schema.end(); ++itr){
        vector<string> row = itr->second;
        this->rR[counter] = row[this->attribute_index];

        counter += 1;
    }
    this->N = counter;
    this->d = length_of_bytes * (length_of_bytes + 1) / 2;
}

void DataOwner::Kgen(public_key& pk) { // Key generation
    pk.bps.q = conv<ZZ>("16283262548997601220198008118239886026907663399064043451383740756301306087801"); // f.param
    ZZ_p::init(pk.bps.q);
    random(this->sk);   // Sample alpha in prime q

    pk.d = this->d;     // Set d

    init_bps(pk.bps);

    compute_crs(pk, this->sk); //Calculate crs
}

void DataOwner::init_bps(bilinear_pairings_parameters& bps){ // Initialize g1, g2 & Sample random g1, g2 on Zr
    char param[3072];
    FILE* f = NULL;
    f = fopen("f.param", "rt");
    size_t count = fread(param, 1, 3072, f);
    if (!count) pbc_die("input error");
    pairing_init_set_buf(bps.e, param, count);

    element_init_G1(bps.g1, bps.e);
    element_init_G2(bps.g2, bps.e);
    element_random(bps.g1);
    element_random(bps.g2);

    fclose(f);
}

void DataOwner::compute_crs(public_key& pk, const ZZ_p& alpha){
    element_pp_t g1_pp, g2_pp;
    element_pp_init(g1_pp, pk.bps.g1);
    element_pp_init(g2_pp, pk.bps.g2);

    element_t alpha_pbc;
    element_init_Zr(alpha_pbc, pk.bps.e);
    convert_to_alpha_pbc(pk, alpha_pbc);

    vector<element_t> tmp1(pk.d+1), tmp2(pk.d+1);
    element_t alpha_mul_g1, alpha_mul_g2;
    element_init_Zr(alpha_mul_g1, pk.bps.e);
    element_init_Zr(alpha_mul_g2, pk.bps.e);
    element_set1(alpha_mul_g1);
    element_set1(alpha_mul_g2);
    
    // compute over g1
    // alpha_mul_g1 : alpha_pbc^0 ~ alpha_pbc^d
    element_init_G1(tmp1[0], pk.bps.e);
    element_pp_pow_zn(tmp1[0], alpha_mul_g1, g1_pp);    // alpha_mul_g1 = 1 = alpha_pbc^0
    for(int itr = 1; itr < pk.d+1; ++itr){              // alpha_mul_g1^1 ~ alpha_mul_g1^d
        element_init_G1(tmp1[itr], pk.bps.e);

        element_mul(alpha_mul_g1, alpha_mul_g1, alpha_pbc); // alpha_mul_g1 = alpha_pbc^1 ~ alpha_pbc^d
        element_pp_pow_zn(tmp1[itr], alpha_mul_g1, g1_pp);
    }

    // compute over g2
    // alpha_mul_g2 : alpha_pbc^0 ~ alpha_pbc^d
    element_init_G2(tmp2[0], pk.bps.e);
    element_pp_pow_zn(tmp2[0], alpha_mul_g2, g2_pp);    // alpha_mul_g2 = 1 = alpha_pbc^0
    for(int itr = 1; itr < pk.d+1; ++itr){              // alpha_mul_g2^1 ~ alpha_mul_g2^d
        element_init_G2(tmp2[itr], pk.bps.e);

        element_mul(alpha_mul_g2, alpha_mul_g2, alpha_pbc); // alpha_mul_g2 = alpha_pbc^1 ~ alpha_pbc^d
        element_pp_pow_zn(tmp2[itr], alpha_mul_g2, g2_pp);
    }

    pk.crs_g1_pbc.swap(tmp1);
    pk.crs_g2_pbc.swap(tmp2);
    
    element_pp_clear(g1_pp);
    element_pp_clear(g2_pp);
}

void DataOwner::ZZ_from_element(ZZ& z, element_t& e){ // convert from element to ZZ
    unsigned char bytes[length_of_bytes] = {0,};
    element_to_bytes(bytes, e);
    reverse(bytes, bytes + length_of_bytes);
    ZZFromBytes(z, bytes, length_of_bytes);
}

void DataOwner::element_from_ZZ(element_t& e, ZZ& z){ // convert from ZZ to element
    unsigned char bytes[length_of_bytes] = {0,};
    BytesFromZZ(bytes, z, length_of_bytes);
    reverse(bytes, bytes + length_of_bytes);
    element_from_bytes(e, bytes);
}

void DataOwner::convert_to_alpha_pbc(public_key& pk, element_t& alpha_pbc){ // Convert data type of sk from ZZ_p to element_t
    ZZ alpha_ZZ = conv<ZZ>(this->sk);
    element_from_ZZ(alpha_pbc, alpha_ZZ);   // ZZ -> element
}

void DataOwner::Setup(public_key& pk, authentication_information& auth_info, element_t& root){
    vector<ZZ_pX> G_poly;

    expand_N();
    for(int itr = (this->rR).size(); itr < this->N; itr++){
        update_index_list.push_back(itr);
    }
    expand_rR();

    polynomialization(G_poly, this->rR, pk);
    for(int itr = (this->N) - update_index_list.size(); itr < this->N; itr++){
        clear(G_poly[itr]);
    }
    BuildEBT(auth_info, pk, G_poly);

    auth_info.rR = this->rR;

    element_init_G1(root, pk.bps.e);
    element_set(root, auth_info.Tree[2*this->N-2].digest);
}

void DataOwner::expand_N(){
    if(log2(this->N) - (int)(log2(this->N)) == 0){  // N = 2^n -> expand 2 times
        this->N = 2*(this->N);
    }else{                              // 2^(n-1) < N < 2^n -> expand more 2 times
        this->N = pow(2, (int)(log2(2*(this->N)))+1);
    }
}

void DataOwner::expand_rR(){
    string s = "";
    for(int itr = (this->N)/2; itr < this->N; itr++){
        (this->rR).insert(pair<int, string>(itr, ""));
    }
}

void DataOwner::polynomialization(vector<ZZ_pX>& G_poly, map<int, string>& rR, public_key& pk){
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

void DataOwner::polynomialization(ZZ_pX& G_poly, string& tuple, public_key& pk){
    string s = "";
    vector<string> prefixes;

    s = tuple;
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

    G_poly = F;
    prefixes.clear();
}

void DataOwner::make_prefixes(vector<string>& prefixes, string& s){
    int l = s.size();                                               // length of tuple
    for(int itr_suf = l-1; itr_suf >= 0; --itr_suf){                // suffixes of a tuple
        for(int itr_pre = 1; itr_pre <= l - itr_suf; ++itr_pre){    // prefixes of a suffix
            prefixes.push_back(s.substr(itr_suf, itr_pre));
        }
    }

    sort(prefixes.begin(), prefixes.end());                         // Remove duplicate
    prefixes.erase(unique(prefixes.begin(), prefixes.end()), prefixes.end());
}

void DataOwner::BuildEBT(authentication_information& auth_info, public_key& pk, vector<ZZ_pX>& G_poly){
    vec_ZZ_p G_eval;
    G_eval.SetLength(this->N);
    
    for(int itr = 0; itr < this->N; itr++){     // Evaluate G(X) at point alpha = G(sk)
        eval(G_eval[itr], G_poly[itr], this->sk);
    }

    element_t alpha_pbc;
    element_init_Zr(alpha_pbc, pk.bps.e);
    convert_to_alpha_pbc(pk, alpha_pbc);        // Convert data type of sk from ZZ_p to element_t

    compute_evaluation_binding_value(auth_info.evaluation_binding_value, pk, G_eval);    // Compute ebv

    make_Tree(auth_info.Tree, auth_info.evaluation_binding_value, pk);                   // Make Tree
    
    if(check_Tree(auth_info.Tree, pk) == false){
        abort();
    }
}

void DataOwner::compute_evaluation_binding_value(vector<element_t>& ebv, public_key& pk, vec_ZZ_p& G_eval){ // Compute g^G(sk)
    vector<element_t> tmp_ebv(this->N);
    element_t tmp_G_eval_pbc, alpha_pbc;
    element_pp_t g1_pp;

    element_init_Zr(tmp_G_eval_pbc, pk.bps.e);
    element_init_Zr(alpha_pbc, pk.bps.e);
    element_pp_init(g1_pp, pk.bps.g1);
    convert_to_alpha_pbc(pk, alpha_pbc);            // Convert data type of sk from ZZ_p to element_t

    for(int itr = 0; itr < this->N; ++itr){         // Compute evaluation_binding_value = g1^G(sk)
        element_init_G1(tmp_ebv[itr], pk.bps.e);
        ZZ tmp_ZZ = conv<ZZ>(G_eval[itr]);
        element_from_ZZ(tmp_G_eval_pbc, tmp_ZZ);
        element_pp_pow_zn(tmp_ebv[itr], tmp_G_eval_pbc, g1_pp); // tmp_ebv[itr] = g1^tmp_G_eval_pbc = g1^G(sk)
    }
    
    ebv.swap(tmp_ebv);   // save to evaluation_binding_value
    element_pp_clear(g1_pp);
}

void DataOwner::make_Tree(vector<Node>& Tree, vector<element_t>& ebv, public_key& pk){
    Tree.resize(2*(this->N)-1);

    element_t g1_alpha, alpha_pbc;                      // init params
    element_t identity;
    element_pp_t g1_pp;

    element_init_G1(g1_alpha, pk.bps.e);
    element_init_G1(identity, pk.bps.e);
    element_init_Zr(alpha_pbc, pk.bps.e);
    element_pp_init(g1_pp, pk.bps.g1);
    convert_to_alpha_pbc(pk, alpha_pbc);                // Convert data type of sk from ZZ_p to element_t
    element_pp_pow_zn(g1_alpha, alpha_pbc, g1_pp);      // compute g1^sk

    for(int h = 0, counter = 0; h <= log2(this->N); h++){
        if(h == 0){                                                                         // leaf node
            for(int itr = counter; itr < this->N; itr++){
                element_init_G1(Tree[itr].digest, pk.bps.e);
                if(element_cmp(identity, ebv[itr]) == 0){   // same
                    element_set(Tree[itr].digest, ebv[itr]);
                } else{                                     // not same
                    element_mul(Tree[itr].digest, ebv[itr], g1_alpha);
                }

                Tree[itr].height = 0;
                Tree[itr].index = itr;
                Tree[itr].isMatch = false;
            }
            counter = this->N;
        }
        else{                                                                               // internal node
            for(int itr = counter, j = 0; itr < counter + this->N / pow(2, h); itr++, j++){
                element_init_G1(Tree[itr].digest, pk.bps.e);

                Tree[itr].left = &Tree[itr - this->N/pow(2, h-1) + j];
                Tree[itr].right = &Tree[itr - this->N/pow(2, h-1) + j + 1];

                Tree[itr].left->parent = &Tree[itr];
                Tree[itr].right->parent = &Tree[itr];

                element_mul(Tree[itr].digest, Tree[itr].left->digest, Tree[itr].right->digest);

                Tree[itr].height = h;
                Tree[itr].index = itr;
                Tree[itr].isMatch = false;
            }
            counter += this->N / pow(2,h);
        }
    }
    element_pp_clear(g1_pp);
}

bool DataOwner::check_Tree(vector<Node>& Tree, public_key& pk){
    element_t check;
    element_init_G1(check, pk.bps.e);

    for(int i = 0; i < this->N; ++i){   // make Tree's root
        element_mul(check, check, Tree[i].digest);
    }
    
    if(!element_cmp(check, Tree[2*this->N-2].digest)){
        std::cout << "EBT is successfully made" << endl;
        return true;    // equal
    }
    else{
        std::cout << "EBT is not successfully made" << endl;
        return false;   // not equal
    }
}

void DataOwner::Update_Insertion(string& tuple, public_key& pk, authentication_information& auth_info, element_t& root){
    try {
        if (tuple.length() > length_of_bytes) throw out_of_range("insert tuple length exceeded 32 bytes");
    } catch(const exception& error) {
        cout << error.what() << endl;
        abort();
    }

    int update_index = (this->update_index_list).front();

    (this->rR)[update_index] = tuple;
    auth_info.rR[update_index] = tuple;

    deque<Node*> path;
    Node* cur = &auth_info.Tree[update_index];

    ZZ_pX G_poly;
    ZZ_p G_eval;

    polynomialization(G_poly, tuple, pk);
    G_eval = eval(G_poly, this->sk);

    element_t G_eval_pbc, g1_G_eval;
    element_t g1_alpha, alpha_pbc;
    element_init_Zr(G_eval_pbc, pk.bps.e);
    element_init_Zr(alpha_pbc, pk.bps.e);
    element_init_G1(g1_G_eval, pk.bps.e);
    element_init_G1(g1_alpha, pk.bps.e);
    element_from_ZZ_p(G_eval_pbc, G_eval);
    convert_to_alpha_pbc(pk, alpha_pbc);
    element_pow_zn(g1_alpha, pk.bps.g1, alpha_pbc);     // compute g1^sk
    element_pow_zn(g1_G_eval, pk.bps.g1, G_eval_pbc);     // compute g1^G(sk)
    
    element_set(auth_info.evaluation_binding_value[update_index], g1_G_eval);
    element_mul(auth_info.Tree[update_index].digest, g1_G_eval, g1_alpha);

    // find path to root
    while(cur->parent != NULL){
        cur = cur->parent;
        path.push_back(cur);
    }

    // compute path's element from children
    cur = NULL;
    while(!path.empty()){
        cur = path.front();
        element_mul(cur->digest, cur->left->digest, cur->right->digest);
        path.pop_front();
    }

    element_set(root, cur->digest); // set auth_info.root

    (this->update_index_list).pop_front();
}

void DataOwner::Update_Deletion(string& tuple, public_key& pk, authentication_information& auth_info, element_t& root){
    try {
        if (tuple.length() > length_of_bytes) throw out_of_range("delete tuple length exceeded 32 bytes");
    } catch(const exception& error) {
        cout << error.what() << endl;
        abort();
    }
    
    int update_index = -1;
    for(auto& itr : auth_info.rR){  // find first matching tuple's index
        if(itr.second == tuple){
            update_index = itr.first;
            break;
        }
    }
    if(update_index == -1){
        return;
    }

    (this->rR)[update_index] = "";
    auth_info.rR[update_index] = "";

    deque<Node*> path;
    Node* cur = &auth_info.Tree[update_index];

    element_t identity;
    element_init_G1(identity, pk.bps.e);
    element_set0(identity);
    
    element_set(auth_info.evaluation_binding_value[update_index], identity);
    element_set(auth_info.Tree[update_index].digest, identity);

    // find path to root
    while(cur->parent != NULL){
        cur = cur->parent;
        path.push_back(cur);
    }

    // compute path's element from children
    cur = NULL;
    while(!path.empty()){
        cur = path.front();
        element_mul(cur->digest, cur->left->digest, cur->right->digest);
        path.pop_front();
    }

    element_set(root, cur->digest); // set auth_info.root

    (this->update_index_list).push_front(update_index);
}

void DataOwner::element_from_ZZ_p(element_t& e, ZZ_p& zp){ // convert from ZZ_p to element
    unsigned char bytes[length_of_bytes] = {0,};
    ZZ z = conv<ZZ>(zp);
    BytesFromZZ(bytes, z, length_of_bytes);
    reverse(bytes, bytes + length_of_bytes);
    element_from_bytes(e, bytes);
}
