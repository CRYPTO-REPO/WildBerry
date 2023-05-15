#include "util.h"
#include "DataOwner.h"
#include "Prover.h"
#include "Verifier.h"

using namespace std;
using namespace NTL;

int main(){
    typedef chrono::high_resolution_clock elapse_time;
    auto v_duration_total = 0, p_duration_total = 0, do_duration_total = 0;

    // Initialize parties
    DataOwner DO;
    Prover P;
    Verifier V;

    // Initialize public parameters
    public_key pk;
    authentication_information auth_info;
    element_t root;

    // Initialize return value
    query_result QR;
    bool verification_result = false;

    // Load dataset & Instantiate relation schema
    string filename = "dataset_1024x6.csv";
    int attribute_index = 3;
    DO.SetFileName(filename);
    DO.SetAttribute(attribute_index);
    DO.load_data();
    DO.instantiate_relation_schema();

    // Key generation & Set up
    auto do_start = elapse_time::now();
    DO.Kgen(pk);
    DO.Setup(pk, auth_info, root);
    auto do_end = elapse_time::now();

    string patterns[10] = {"es", "e", "pot", "o", "take", "jingo", "pro", "ci", "te", "shak"};
    int iteration = 100;

    // 100 iteration
    for(auto& s : patterns){
         for(int itr = 0; itr < iteration; itr++){
            cout << itr << endl;
            // Set pattern for query
            string pattern = s;
            V.SetPattern(pattern);

            // Query
            auto p_start = elapse_time::now();
            P.Query(QR, pk, pattern, auth_info);
            auto p_end = elapse_time::now();

            // Check the number of matches & covers
            if(itr == 0){
                if (QR.qv.matches.empty()){
                    cout << "no matches" << endl;
                } else{
                    cout << "mathces : " << QR.qv.matches.size() << endl;
                }
                if (QR.qv.covers.empty()){
                    cout << "no covers" << endl;
                } else{
                    cout << "covers : " << QR.qv.covers.size() << endl;
                }
            }

            auto v_start = elapse_time::now();
            verification_result = V.Verify(pk, QR, root);
            auto v_end = elapse_time::now();

            if(verification_result){
                // cout << "Verification is successed" << endl;
            } else{
                cout << "Verification is failed" << endl;
                cout << "s : " << s << "\tindex : " << itr << endl;
                abort();
            }

            int p_elapsed_in_msec = chrono::duration_cast<std::chrono::milliseconds>(p_end - p_start).count();
            p_duration_total += p_elapsed_in_msec;
            int v_elapsed_in_msec = chrono::duration_cast<std::chrono::milliseconds>(v_end - v_start).count();
            v_duration_total += v_elapsed_in_msec;

            QR.qp.proof_matches.clear();
            QR.qp.proof_mismatches.clear();
            QR.qp.w.clear();
            QR.qv.covers.clear();
            QR.qv.matches.clear();
            QR.qv.tuple.clear();
        }
        cout << "Average 100 times of string s :" << s << endl;
        cout << "prover duration average time(ms) : " << p_duration_total/iteration << endl;
        cout << "verifier duration average time(ms) : " << v_duration_total/iteration << endl;
        cout << endl;
        p_duration_total = 0;
        v_duration_total = 0;
    }

    int do_elapsed_in_msec = chrono::duration_cast<std::chrono::milliseconds>(do_end - do_start).count();
    do_duration_total += do_elapsed_in_msec;
    cout << "data owner duration total time(ms) : " << do_duration_total << endl;

    element_clear(pk.bps.g1);
    element_clear(pk.bps.g2);
    pairing_clear(pk.bps.e);

    cout << "main finished" << endl;
}