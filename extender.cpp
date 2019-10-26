#include <atomic>
#include <mutex>
#include "graph.h"
#include <unordered_map>

#define TB_MASK (1<<2)-1

std::atomic<uint64_t> extender_body::num_extend_tiles(0);

void extender_body::operator()(extender_input input, extender_node::output_ports_type &op)
{
    auto &payload = get<0>(input);

    auto &read = get<0>(payload);
    auto &fwData = get<1>(payload);
    auto &rcData = get<2>(payload);

    size_t token = get<1>(input);

    extender_output output;

    std::unordered_map<uint64_t, bool> fwAnchors; 
    std::unordered_map<uint64_t, bool> rcAnchors; 

    fwAnchors.clear();
    rcAnchors.clear();
    

    for (auto anc: fwData) {
        uint64_t key = anc.reference_offset;
        key = (key << 32) + anc.query_offset;
        fwAnchors[key] = true;
    }

    char* query = (char*) read.seq.data();
    for (auto anc: fwData) {
        uint64_t key = anc.reference_offset;
        key = (key << 32) + anc.query_offset;
        int anc_score = anc.score;

        if (fwAnchors[key]) {
            Alignment e = makeAlignment(read, anc, '+');
            bool right_ext_done = false;
            bool left_ext_done = false;
            uint8_t align_fields = 0;
            while (!right_ext_done) {
                uint32_t r_start = e.curr_reference_offset;
                uint32_t r_end = std::min(e.curr_reference_offset + cfg.tile_size, e.reference_length);
                uint32_t q_start = e.curr_query_offset;
                uint32_t q_end = std::min(e.curr_query_offset + cfg.tile_size, e.query_length);
                extend_tile tile(e.reference_start_addr + r_start, q_start, r_end-r_start, q_end-q_start);

                extend_output op = g_GACTXRequest(tile, align_fields);
                num_extend_tiles++;
                e.num_right_tiles++;

                int tb_size = op.tb_pointers.size();
                char* ref_buf = (char*) malloc(16 * tb_size);
                char* query_buf = (char*) malloc(16 * tb_size);

                uint32_t ref_pos = e.curr_reference_offset + op.max_ref_offset;
                uint32_t query_pos = e.curr_query_offset + op.max_query_offset;

                uint32_t rp = e.reference_start_addr + ref_pos;
                uint32_t qp = query_pos;
                
                int num_r_bases = 0, num_q_bases = 0;
                uint32_t r_tile_begin = e.reference_start_addr + r_start + cfg.tile_size - cfg.tile_overlap;
                uint32_t q_tile_begin = q_start + cfg.tile_size - cfg.tile_overlap;

                int tb_pos = 0;

                bool begin_tb = false;
                //TODO: add condition

                for (int i = 0; i < tb_size; i++) {
                    uint32_t tb_ptr = op.tb_pointers[i];
                    for(int j = 0; j < 16; j++){
                        int dir = ((tb_ptr >> (2*j)) & TB_MASK);
                        switch(dir) {
                            case Z:
                                break;
                            case D: 
                                if (begin_tb) {
                                    ref_buf[tb_pos] = g_DRAM->buffer[rp];
                                    query_buf[tb_pos] = '-';
                                    tb_pos++;
                                    num_r_bases++;
                                }
                                rp--;
                                break;
                            case I:
                                if (begin_tb) {
                                    ref_buf[tb_pos] = '-';
                                    query_buf[tb_pos] = query[qp];
                                    tb_pos++;
                                    num_q_bases++;
                                }
                                qp--;
                                break;
                            case M:
                                if ((rp < r_tile_begin) && (qp < q_tile_begin)) {
                                    begin_tb = true;
                                }
                                if (begin_tb) {
                                    uint64_t k = rp;
                                    k = (k << 32) + qp;
                                    if (fwAnchors.find(k) != fwAnchors.end()) {
                                        fwAnchors[k] = false;
                                    }
                                    ref_buf[tb_pos] = g_DRAM->buffer[rp];
                                    query_buf[tb_pos] = query[qp];
                                    tb_pos++;
                                    num_r_bases++;
                                    num_q_bases++;
                                }
                                rp--;
                                qp--;
                                break;
                        }
                    }
                }

                e.curr_reference_offset += num_r_bases;
                e.curr_query_offset += num_q_bases;

                std::string aligned_reference_str (ref_buf, tb_pos);
                std::string aligned_query_str (query_buf, tb_pos);

                std::reverse(aligned_reference_str.begin(), aligned_reference_str.end());
                std::reverse(aligned_query_str.begin(), aligned_query_str.end());

                e.aligned_reference_str = e.aligned_reference_str + aligned_reference_str; 
                e.aligned_query_str = e.aligned_query_str + aligned_query_str; 

                std::string tile_ref(g_DRAM->buffer + e.reference_start_addr + r_start, r_end-r_start);
                std::string tile_query(g_DRAM->buffer + g_DRAM->referenceSize + q_start, q_end-q_start);

                free(ref_buf);
                free(query_buf);


                if ((tb_size == 0) || (r_tile_begin >= e.reference_start_addr + e.reference_length) || (q_tile_begin >= e.query_length)) {
                    right_ext_done = true;
                    e.reference_end_offset = e.curr_reference_offset;
                    e.query_end_offset = e.curr_query_offset;
                    e.curr_reference_offset = e.reference_start_offset;
                    e.curr_query_offset = e.query_start_offset;
                }
            }

            if ((e.curr_reference_offset == 0) || (e.curr_query_offset == 0)) {
                left_ext_done = true;
            }
            align_fields = reverse_ref + reverse_query;
            while (!left_ext_done) {
                uint32_t r_end = e.curr_reference_offset;
                uint32_t r_start = std::max(r_end, (uint32_t) cfg.tile_size) - cfg.tile_size;
                uint32_t q_end = e.curr_query_offset;
                uint32_t q_start = std::max(q_end, (uint32_t) cfg.tile_size) - cfg.tile_size;

                std::string t_r(g_DRAM->buffer + e.reference_start_addr + r_start, r_end - r_start);
                std::string t_q(g_DRAM->buffer + g_DRAM->referenceSize + q_start, q_end - q_start);
                std::reverse(t_r.begin(), t_r.end());
                std::reverse(t_q.begin(), t_q.end());
                extend_tile tile(e.reference_start_addr + r_start, q_start, r_end-r_start, q_end-q_start);

                extend_output op = g_GACTXRequest(tile, align_fields);
                num_extend_tiles++;
                e.num_left_tiles++;

                int tb_size = op.tb_pointers.size();
                char* ref_buf = (char*) malloc(16 * tb_size);
                char* query_buf = (char*) malloc(16 * tb_size);

                uint32_t ref_pos = e.curr_reference_offset - op.max_ref_offset - 1;
                uint32_t query_pos = e.curr_query_offset - op.max_query_offset - 1;

                uint32_t rp = e.reference_start_addr + ref_pos;
                uint32_t qp = query_pos;

                int num_r_bases = 0, num_q_bases = 0;
                uint32_t r_tile_begin = (e.reference_start_addr + std::max(r_end, (uint32_t) (cfg.tile_size-cfg.tile_overlap))) - (cfg.tile_size - cfg.tile_overlap);
                uint32_t q_tile_begin = (std::max(q_end, (uint32_t) (cfg.tile_size-cfg.tile_overlap))) - (cfg.tile_size - cfg.tile_overlap);

                int tb_pos = 0;

                bool begin_tb = false;

                for (int i = 0; i < tb_size; i++) {
                    uint32_t tb_ptr = op.tb_pointers[i];
                    for(int j = 0; j < 16; j++){
                        int dir = ((tb_ptr >> (2*j)) & TB_MASK);
                        switch(dir) {
                            case Z:
                                break;
                            case D: 
                                if (begin_tb) {
                                    ref_buf[tb_pos] = g_DRAM->buffer[rp];
                                    query_buf[tb_pos] = '-';
                                    tb_pos++;
                                    num_r_bases++;
                                }
                                rp++;
                                break;
                            case I:
                                if (begin_tb) {
                                    ref_buf[tb_pos] = '-';
                                    query_buf[tb_pos] = query[qp];
                                    tb_pos++;
                                    num_q_bases++;
                                }
                                qp++;
                                break;
                            case M:
                                if ((rp >= r_tile_begin) && (qp >= q_tile_begin)) {
                                    begin_tb = true;
                                }
                                if (begin_tb) {
                                    uint64_t k = rp;
                                    k = (k << 32) + qp;
                                    if (fwAnchors.find(k) != fwAnchors.end()) {
                                        fwAnchors[k] = false;
                                    }
                                    ref_buf[tb_pos] = g_DRAM->buffer[rp];
                                    query_buf[tb_pos] = query[qp];
                                    tb_pos++;
                                    num_r_bases++;
                                    num_q_bases++;
                                }
                                rp++;
                                qp++;
                                break;
                        }
                    }
                }
                e.curr_reference_offset -= num_r_bases;
                e.curr_query_offset -= num_q_bases;
                std::string aligned_reference_str (ref_buf, tb_pos);
                std::string aligned_query_str (query_buf, tb_pos);

                e.aligned_reference_str = aligned_reference_str + e.aligned_reference_str; 
                e.aligned_query_str = aligned_query_str + e.aligned_query_str; 

                free(ref_buf);
                free(query_buf);

                if ((tb_size == 0) || (r_tile_begin == e.reference_start_addr) || (q_tile_begin == 0)) {
                    left_ext_done = true;
                    e.reference_start_offset = e.curr_reference_offset;
                    e.query_start_offset = e.curr_query_offset;
                    if (tb_size == 0) {
                    }
                }
            }

            output.push_back(e);
        }
    }
    
    for (auto anc: rcData) {
        uint64_t key = anc.reference_offset;
        key = (key << 32) + anc.query_offset;
        rcAnchors[key] = true;
    }
    
    char* rc_query = (char*) read.rc_seq.data();
    for (auto anc: rcData) {
        uint64_t key = anc.reference_offset;
        key = (key << 32) + anc.query_offset;

        if (rcAnchors[key]) {
            Alignment e = makeAlignment(read, anc, '-');
            bool right_ext_done = false;
            bool left_ext_done = false;
            uint8_t align_fields = reverse_query + complement_query;
            
            while (!right_ext_done) {
                uint32_t r_start = e.curr_reference_offset;
                uint32_t r_end = std::min(e.curr_reference_offset + cfg.tile_size, e.reference_length);
                uint32_t q_start = e.curr_query_offset;
                uint32_t q_end = std::min(e.curr_query_offset + cfg.tile_size, e.query_length);

                uint32_t q_tile_start = e.query_length - q_end;
                extend_tile tile(e.reference_start_addr + r_start, q_tile_start, r_end-r_start, q_end-q_start);

                extend_output op = g_GACTXRequest(tile, align_fields);
                num_extend_tiles++;
                e.num_right_tiles++;

                int tb_size = op.tb_pointers.size();
                char* ref_buf = (char*) malloc(16 * tb_size);
                char* query_buf = (char*) malloc(16 * tb_size);

                uint32_t ref_pos = e.curr_reference_offset + op.max_ref_offset;
                uint32_t query_pos = e.curr_query_offset + op.max_query_offset;

                uint32_t rp = e.reference_start_addr + ref_pos;
                uint32_t qp = query_pos;

                int num_r_bases = 0, num_q_bases = 0;

                uint32_t r_tile_begin = e.reference_start_addr + r_start + cfg.tile_size - cfg.tile_overlap;
                uint32_t q_tile_begin = q_start + cfg.tile_size - cfg.tile_overlap;

                int tb_pos = 0;

                bool begin_tb = false;

                for (int i = 0; i < tb_size; i++) {
                    uint32_t tb_ptr = op.tb_pointers[i];
                    for(int j = 0; j < 16; j++){
                        int dir = ((tb_ptr >> (2*j)) & TB_MASK);
                        switch(dir) {
                            case Z:
                                break;
                            case D: 
                                if (begin_tb) {
                                    ref_buf[tb_pos] = g_DRAM->buffer[rp];
                                    query_buf[tb_pos] = '-';
                                    tb_pos++;
                                    num_r_bases++;
                                }
                                rp--;
                                break;
                            case I:
                                if (begin_tb) {
                                    ref_buf[tb_pos] = '-';
                                    query_buf[tb_pos] = rc_query[qp];
                                    tb_pos++;
                                    num_q_bases++;
                                }
                                qp--;
                                break;
                            case M:
                                if ((rp < r_tile_begin) && (qp < q_tile_begin)) {
                                    begin_tb = true;
                                }
                                if (begin_tb) {
                                    uint64_t k = rp;
                                    k = (k << 32) + qp;
                                    if (rcAnchors.find(k) != rcAnchors.end()) {
                                        rcAnchors[k] = false;
                                    }
                                    ref_buf[tb_pos] = g_DRAM->buffer[rp];
                                    query_buf[tb_pos] = rc_query[qp];
                                    tb_pos++;
                                    num_r_bases++;
                                    num_q_bases++;
                                }
                                rp--;
                                qp--;
                                break;
                        }
                    }
                }
                e.curr_reference_offset += num_r_bases;
                e.curr_query_offset += num_q_bases;

                std::string aligned_reference_str (ref_buf, tb_pos);
                std::string aligned_query_str (query_buf, tb_pos);

                std::reverse(aligned_reference_str.begin(), aligned_reference_str.end());
                std::reverse(aligned_query_str.begin(), aligned_query_str.end());


                e.aligned_reference_str = e.aligned_reference_str + aligned_reference_str; 
                e.aligned_query_str = e.aligned_query_str + aligned_query_str; 

                free(ref_buf);
                free(query_buf);

                if ((tb_size == 0) || (r_tile_begin >= e.reference_start_addr + e.reference_length) || (q_tile_begin >= e.query_length)) {
                    right_ext_done = true;
                    e.reference_end_offset = e.curr_reference_offset;
                    e.query_end_offset = e.curr_query_offset;
                    e.curr_reference_offset = e.reference_start_offset;
                    e.curr_query_offset = e.query_start_offset;
                }
            }
            
            if ((e.curr_reference_offset == 0) || (e.curr_query_offset == 0)) {
                left_ext_done = true;
            }

            align_fields = reverse_ref + complement_query;
            while (!left_ext_done) {
                uint32_t r_end = e.curr_reference_offset;
                uint32_t r_start = std::max(r_end, (uint32_t) cfg.tile_size) - cfg.tile_size;
                uint32_t q_end = e.curr_query_offset;
                uint32_t q_start = std::max(q_end, (uint32_t) cfg.tile_size) - cfg.tile_size;

                uint32_t q_tile_start = e.query_length - q_end;
                extend_tile tile(e.reference_start_addr + r_start, q_tile_start, r_end-r_start, q_end-q_start);

                extend_output op = g_GACTXRequest(tile, align_fields);
                num_extend_tiles++;
                e.num_left_tiles++;

                int tb_size = op.tb_pointers.size();
                char* ref_buf = (char*) malloc(16 * tb_size);
                char* query_buf = (char*) malloc(16 * tb_size);

                uint32_t ref_pos = e.curr_reference_offset - op.max_ref_offset - 1;
                uint32_t query_pos = e.curr_query_offset - op.max_query_offset - 1;

                uint32_t rp = e.reference_start_addr + ref_pos;
                uint32_t qp = query_pos;

                int num_r_bases = 0, num_q_bases = 0;

                uint32_t r_tile_begin = (e.reference_start_addr + std::max(r_end, (uint32_t) (cfg.tile_size-cfg.tile_overlap))) - (cfg.tile_size - cfg.tile_overlap);
                uint32_t q_tile_begin = std::max(q_end, (uint32_t) (cfg.tile_size-cfg.tile_overlap)) - (cfg.tile_size - cfg.tile_overlap);

                int tb_pos = 0;

                bool begin_tb = false;

                for (int i = 0; i < tb_size; i++) {
                    uint32_t tb_ptr = op.tb_pointers[i];
                    for(int j = 0; j < 16; j++){
                        int dir = ((tb_ptr >> (2*j)) & TB_MASK);
                        switch(dir) {
                            case Z:
                                break;
                            case D: 
                                if (begin_tb) {
                                    ref_buf[tb_pos] = g_DRAM->buffer[rp];
                                    query_buf[tb_pos] = '-';
                                    tb_pos++;
                                    num_r_bases++;
                                }
                                rp++;
                                break;
                            case I:
                                if (begin_tb) {
                                    ref_buf[tb_pos] = '-';
                                    query_buf[tb_pos] = rc_query[qp];
                                    tb_pos++;
                                    num_q_bases++;
                                }
                                qp++;
                                break;
                            case M:
                                if ((rp >= r_tile_begin) && (qp >= q_tile_begin)) {
                                    begin_tb = true;
                                }
                                if (begin_tb) {
                                    uint64_t k = rp;
                                    k = (k << 32) + qp;
                                    if (rcAnchors.find(k) != rcAnchors.end()) {
                                        rcAnchors[k] = false;
                                    }
                                    ref_buf[tb_pos] = g_DRAM->buffer[rp];
                                    query_buf[tb_pos] = rc_query[qp];
                                    tb_pos++;
                                    num_r_bases++;
                                    num_q_bases++;
                                }
                                rp++;
                                qp++;
                                break;
                        }
                    }
                }
                e.curr_reference_offset -= num_r_bases;
                e.curr_query_offset -= num_q_bases;
                
                std::string aligned_reference_str (ref_buf, tb_pos);
                std::string aligned_query_str (query_buf, tb_pos);

                e.aligned_reference_str = aligned_reference_str + e.aligned_reference_str; 
                e.aligned_query_str = aligned_query_str + e.aligned_query_str; 


                free(ref_buf);
                free(query_buf);

                if ((tb_size == 0) || (r_tile_begin == e.reference_start_addr) || (q_tile_begin == 0)) {
                    left_ext_done = true;
                    e.reference_start_offset = e.curr_reference_offset;
                    e.query_start_offset = e.curr_query_offset;
                }
            }
            
            output.push_back(e);
        }
    }
    
    get<1>(op).try_put(token);
    get<0>(op).try_put(printer_input(printer_payload(read, output), token));
}


Alignment extender_body::makeAlignment(reader_output read, anchor anc, char strand)
{
    Alignment extend_alignment;
    
    const size_t read_len = read.seq.size();
    
    int chr_id = std::upper_bound(r_chr_coord.cbegin(), r_chr_coord.cend(), anc.reference_offset) - r_chr_coord.cbegin() - 1;
    uint32_t chr_start = r_chr_coord[chr_id];

    extend_alignment.chr_id = chr_id;

    extend_alignment.curr_reference_offset = anc.reference_offset - chr_start;
    extend_alignment.curr_query_offset = anc.query_offset;

    extend_alignment.reference_start_offset = anc.reference_offset - chr_start;
    extend_alignment.query_start_offset = anc.query_offset;
    
    extend_alignment.reference_end_offset = anc.reference_offset - chr_start;
    extend_alignment.query_end_offset = anc.query_offset;

    extend_alignment.reference_start_addr = chr_start;
    extend_alignment.query_start_addr = g_DRAM->referenceSize;

    extend_alignment.reference_length = r_chr_len[chr_id];
    extend_alignment.query_length = read_len;

    extend_alignment.aligned_reference_str = "";
    extend_alignment.aligned_query_str = "";

    extend_alignment.score = 0;

    extend_alignment.strand = strand;


    extend_alignment.num_left_tiles = 0;
    extend_alignment.num_right_tiles = 0;

    return extend_alignment;
}

