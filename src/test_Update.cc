#include "util.h"
#include "DataOwner.h"
#include "Prover.h"
#include "Verifier.h"

using namespace std;
using namespace NTL;

int main(){
    typedef chrono::high_resolution_clock elapse_time;
    auto insertion_duration_total = 0, deletion_duration_total = 0, do_duration_total = 0;

    // Initialize parties
    DataOwner DO;

    // Initialize public parameters
    public_key pk;
    authentication_information auth_info;
    element_t root;

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

    string update_patterns[5] = {"tuple", "string", "update insertion", "update deletion", "test string for update"};
    int iteration = 100;

    // 100 iteration
    for(auto& s : update_patterns){
        // insertion
        for(int itr = 0; itr < iteration; itr++){
            cout << itr << endl;

            auto insertion_start = elapse_time::now();
            DO.Update_Insertion(s, pk, auth_info, root);
            auto insertion_end = elapse_time::now();

            int insert_elapsed_in_msec = chrono::duration_cast<std::chrono::microseconds>(insertion_end - insertion_start).count();
            insertion_duration_total += insert_elapsed_in_msec;
        }
        // deletion
        for(int itr = 0; itr < iteration; itr++){
            cout << itr << endl;

            auto deletion_start = elapse_time::now();
            DO.Update_Deletion(s, pk, auth_info, root);
            auto deletion_end = elapse_time::now();

            int delete_elapsed_in_msec = chrono::duration_cast<std::chrono::microseconds>(deletion_end - deletion_start).count();
            deletion_duration_total += delete_elapsed_in_msec;
        }
        cout << "Average 100 times of string s : " << s << endl;
        cout << "Insertion duration average time(us) : " << insertion_duration_total / iteration << endl;
        cout << "Deletion duration average time(us) : " << deletion_duration_total / iteration << endl;
        cout << endl;

        insertion_duration_total = 0;
        deletion_duration_total = 0;
    }

    int do_elapsed_in_msec = chrono::duration_cast<std::chrono::milliseconds>(do_end - do_start).count();
    do_duration_total += do_elapsed_in_msec;
    cout << "data owner duration total time(ms) : " << do_duration_total << endl;

    element_clear(pk.bps.g1);
    element_clear(pk.bps.g2);
    pairing_clear(pk.bps.e);

    cout << "main finished" << endl;
}