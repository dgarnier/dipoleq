#if MULTITASK
#define MULTI doMultiTask(1)
#else
#define MULTI {}
#endif

#if MULTITASK
#ifdef __cplusplus
extern "C" {
#endif

   void doMultiTask(long);

#ifdef __cplusplus
}
#endif

#endif
