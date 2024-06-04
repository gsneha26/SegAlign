#include <condition_variable>
#include <mutex>
#include <vector>

extern std::mutex mu;
extern std::condition_variable cv;
extern std::vector<int> available_gpus;

extern int NUM_DEVICES;

extern char** d_ref_seq;
extern uint32_t ref_len;
 
extern uint32_t** d_index_table;
extern uint32_t** d_pos_table;
