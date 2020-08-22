#include <condition_variable>
#include <vector>

extern std::mutex mu;
extern std::condition_variable cv;
extern std::vector<int> available_gpus;

extern int NUM_DEVICES;

extern uint32_t** d_index_table;
extern uint32_t** d_pos_table;
