#include <stdint.h>

typedef int(*InitializeInterface_ptr)(int num_gpu);
typedef void(*SendRefWriteRequest_ptr)(char* seq, size_t addr, uint32_t len);
typedef void(*ClearRef_ptr)();
typedef void(*ShutdownProcessor_ptr)();

extern InitializeInterface_ptr g_InitializeInterface;
extern SendRefWriteRequest_ptr g_SendRefWriteRequest;
extern ClearRef_ptr g_ClearRef;
extern ShutdownProcessor_ptr g_ShutdownProcessor;
