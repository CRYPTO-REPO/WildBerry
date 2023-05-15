#include "util.h"

string public_key::bytes_to_hex_string(const vector<uint8_t>& bytes)
{
    ostringstream stream;
    for (uint8_t b : bytes)
    {
        stream << setw(2) << setfill('0') << hex << static_cast<int>(b);
    }
    return stream.str();
};

void public_key::hash(vec_ZZ_p& hash_value, vector<string>& prefixes){
    ZZ_p ctxDec;
    for(int itr = 0; itr < prefixes.size(); itr++){ // sha256 to all prefixes elements
            
        uint32_t digest_length = SHA256_DIGEST_LENGTH;
        const EVP_MD* algorithm = EVP_sha3_256();
        uint8_t* digest = static_cast<uint8_t*>(OPENSSL_malloc(digest_length));
        EVP_MD_CTX* context = EVP_MD_CTX_new();
        EVP_DigestInit_ex(context, algorithm, nullptr);
        EVP_DigestUpdate(context, prefixes[itr].c_str(), prefixes[itr].size());
        EVP_DigestFinal_ex(context, digest, &digest_length);
        EVP_MD_CTX_destroy(context);
        
        string output = bytes_to_hex_string(vector<uint8_t>(digest, digest + digest_length));
        OPENSSL_free(digest);

        HexToDec(ctxDec, output);                 // convert sha256 hex value to ZZ_p value

        hash_value[itr] = ctxDec;
    }
};

void public_key::hash(ZZ_p& hash_value, string& pattern){ // hash 1 string
    ZZ_p ctxDec;
    uint32_t digest_length = SHA256_DIGEST_LENGTH;
    const EVP_MD* algorithm = EVP_sha3_256();
    uint8_t* digest = static_cast<uint8_t*>(OPENSSL_malloc(digest_length));
    EVP_MD_CTX* context = EVP_MD_CTX_new();
    EVP_DigestInit_ex(context, algorithm, nullptr);
    EVP_DigestUpdate(context, pattern.c_str(), pattern.size());
    EVP_DigestFinal_ex(context, digest, &digest_length);
    EVP_MD_CTX_destroy(context);
        
    string output = bytes_to_hex_string(vector<uint8_t>(digest, digest + digest_length));
    OPENSSL_free(digest);

    HexToDec(ctxDec, output);                 // convert sha256 hex value to ZZ_p value

    hash_value = ctxDec;
};

void public_key::HexToDec(ZZ_p& ctxDec, string hex){  // convert string hex value to ZZ_p value
	ZZ_p base;
    base = 1;
    ctxDec = 0;

	for (int itr = hex.size() - 1; itr >= 0; itr--) {
		if (hex[itr] >= '0' && hex[itr] <= '9') {
			ctxDec += (int(hex[itr]) - 48) * base;
			base = base * 16;
		}
		else if (hex[itr] >= 'A' && hex[itr] <= 'F') {
			ctxDec += (int(hex[itr]) - 55) * base;
			base = base * 16;
		}
	}
};

template <typename S, typename V>
void print_map(map<S, V>& m) {
    for (typename map<S, V>::iterator itr = m.begin(); itr != m.end(); ++itr)
    {
        cout << itr->first << " : ";
        for (const auto &tuple : itr->second) {
        cout << tuple << " ";
        }
        cout << endl;
    }
}

void print_vector_pbc(vector<element_t>& v, int size){
    for(int itr = 0; itr < size; ++itr){
        element_printf("%d : %B\n", itr, v[itr]);
    }
};

template <typename T>
void print_vector(vector<T>& v){
    for(auto& i : v){
        cout << i << ' ';
    }
    cout << endl;
};
