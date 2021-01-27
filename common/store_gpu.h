#include <condition_variable>
#include <vector>
#include <claraparabricks/genomeworks/utils/device_buffer.hpp>
#include <claraparabricks/genomeworks/utils/cudautils.hpp>
#include <claraparabricks/genomeworks/cudaextender/utils.hpp>

using namespace claraparabricks::genomeworks;
using namespace cudautils;
using namespace cudaextender;

extern std::mutex mu;
extern std::condition_variable cv;
extern std::vector<int> available_gpus;

extern int NUM_DEVICES;

extern std::vector<device_buffer<int8_t>> d_ref_seq;
extern int32_t ref_len;
 
extern std::vector<device_buffer<uint32_t>> d_index_table;
extern std::vector<device_buffer<uint32_t>> d_pos_table;

extern std::vector<DefaultDeviceAllocator> allocator_;
extern std::vector<CudaStream> stream_;
