#include <atomic>
#include "ntcoding.h"
#include "graph.h"

std::mutex io_lock;
int print_header = 1;

size_t maf_printer_body::operator()(printer_input input){
    auto &payload = get<0>(input); 

    auto &reads = get<0>(payload);

    auto &data = get<1>(payload);

    size_t token = get<1>(input);

    int mat_offset[] = { 0, 1, 3, 6 };

    for (auto e: data) {
        int score = 0;
        int num_r = 0;
        int num_q = 0;

        int length = e.aligned_reference_str.length();
        int open = 0;
        int index = 0;
        for (int i = 0; i < length; i++) {
            char r = e.aligned_reference_str[i];
            char q = e.aligned_query_str[i];
            if (r == '-') {
                if (open == 1) {
                    score += cfg.gap_extend;
                }
                else {
                    score += cfg.gap_open;
                }
                num_q++;
                open = 1;
            }
            else if (q == '-') {
                if (open == 1) {
                    score += cfg.gap_extend;
                }
                else {
                    score += cfg.gap_open;
                }
                num_r++;
                open = 1;
            }
            else {
                int r_nt = NtChar2IntCaseInsensitive(r);
                int q_nt = NtChar2IntCaseInsensitive(q);
                assert (r_nt <= 4);
                assert (q_nt <= 4);
                if ((r_nt <= 3) && (q_nt <= 3)) {
                    if (r_nt > q_nt) {
                        index = q_nt * 4 + r_nt - mat_offset[q_nt];
                    }
                    else {
                        index = r_nt * 4 + q_nt - mat_offset[r_nt];
                    }
                    score += cfg.gact_sub_mat[index];
                }
                else {
                    score += cfg.gact_sub_mat[10];
                }
                num_r++;
                num_q++;
                open = 0;
            }
        }

        std::string r_chrom = r_chr_id[e.chr_id]; 
        std::string q_chrom = reads.description;
        uint32_t r_len = e.reference_length;
        uint32_t q_len = e.query_length;

        if (print_header == 1) {
            io_lock.lock();
            fprintf(mafFile, "##maf version=1\n");
            print_header = 0;
            io_lock.unlock();
        }

        assert(score >= cfg.first_tile_score_threshold);

        printf("%d\n", score);

        if (score >= cfg.extension_threshold) {
            io_lock.lock();
            //fprintf(mafFile, "a\tscore=%d\n", score);
            //fprintf(mafFile, "s\t%s\t%u\t%d\t+\t%u\t%s\n", r_chrom.c_str(), 1+e.reference_start_offset, num_r, r_len, e.aligned_reference_str.c_str());
            //fprintf(mafFile, "s\t%s\t%u\t%d\t%c\t%u\t%s\n", q_chrom.c_str(), 1+e.query_start_offset, num_q, e.strand, q_len, e.aligned_query_str.c_str());
            //fprintf(mafFile, "\n");
            io_lock.unlock();
        }
    }

    return token;
};

