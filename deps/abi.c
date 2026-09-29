#include <stddef.h>
#include <dace/dacebase.h>
#include <dace/daceaux.h>

/* Introspection only: all mathematical calls bind directly to dacebase.h. */
DACE_API size_t daceDirectSizeofDA(void) { return sizeof(DACEDA); }
DACE_API size_t daceDirectOffsetMem(void) { return offsetof(DACEDA, mem); }
DACE_API int daceDirectMemoryModel(void) { return DACE_MEMORY_MODEL; }
DACE_API size_t daceDirectSizeofMonomial(void) { return sizeof(monomial); }
DACE_API size_t daceDirectStringLength(void) { return DACE_STRLEN; }
DACE_API size_t daceDirectAlignMonomial(void) {
    struct aligned_monomial { char prefix; monomial value; };
    return offsetof(struct aligned_monomial, value);
}
