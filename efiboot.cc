extern "C" {
#if defined(_USERLAND_)
#include <stdio.h>
#include <assert.h>
  static void cnputc(int c) { putchar(c); fflush(stdout); }
#else
#include <sys/types.h>
#include "efiboot.c.punched"
  void cnputc(int c) { efi_cons_putc(cn_tab->cn_dev, c); }
  int cnischar(void) { return (*cn_tab->cn_getc)(cn_tab->cn_dev|0x80); }
  static void addbootarg(int t, size_t l, void* p) { ; }
  int atexit(void (*function)(void)) { return 0; }
  static void puts(const char *fmt) {
    for(int i = 0; fmt[i]; i ++) cnputc(fmt[i]);
  }
#endif
  static void puthex(const int& x0) {
    const int x(x0 < 0 ? - x0 : x0);
    if(x0 < 0) cnputc('-');
    for(int i = sizeof(int) * 8 - 4; 0 <= i; i -= 4)
      cnputc("0123456789abcdef"[(i ? (x >> i) : x) & 0x0f]);
  }
#if !defined(_USERLAND_)
    __dead void panic(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    // vprintf(fmt, ap);
    puts(fmt);
    va_end(ap);
    _rtt();
    /*NOTREACHED*/
  }
  static void assert(const uint64_t f) {
    if(! f) {
      uint64_t* stack;
      puts("assertion failed\n");
      puts("f: "); puthex((size_t)(void*)assert); cnputc('\n');
      for(int i = 0; i < 16; i ++) {
        puthex(stack[i]); cnputc('\n');
      }
      for(; efi_cons_getc(0) != ' ';);
      _rtt();
    }
    return;
  }
  static struct consdev constab[] = {
    { efi_cons_probe, efi_cons_init, efi_cons_getc, efi_cons_putc },
    // { efi_com_probe, efi_com_init, efi_com_getc, efi_com_putc },
    { NULL }
  };
  struct consdev *cn_tab = constab;
#endif
  const bool prngdata[] = {
//#include "prng.h"
  };
  int pdata[] = {
#include "pdata.h"
  };
  static unsigned long long last;
  static unsigned long long lastptr;
  static unsigned long long sam_upper;
  static unsigned long long *v_alloc;
  static int *in_use;
  // N.B. thanks to musl-1.2.3/src/prng/rand.c
  static unsigned long long prng_seed(1234);
  static int prng_idx(0);
  static int myrandom() {
    if(prng_idx < sizeof(prngdata) / sizeof(bool))
      return prngdata[prng_idx ++] ? 0xffffff : 0;
    cnputc('!');
    prng_seed = 6364136223846793005ULL * prng_seed + 1;
    return int(prng_seed >> 33);
  }
  static unsigned long long pprng_seed(8765);
  static int myprandom() {
    pprng_seed = 6364136223846793005ULL * pprng_seed + 1;
    return int(pprng_seed >> 33);
  }
#define M_ALLOC (32 * 1024 * 1024)
#if !defined(_USERLAND_)
  EFI_STATUS calc() ;
  static unsigned long long simplealloc_init() {
    unsigned long long heap = 0;
    EFI_STATUS status;
    status = BS->AllocatePages(AllocateAnyPages, EfiLoaderData,
      1200 * 1024 * 1024 / (4 * 1024), &heap);
    if (status != EFI_SUCCESS)
      panic("BS->AllocatePages()");
    unsigned long long stack(heap);
    v_alloc = reinterpret_cast<unsigned long long *>(heap + 64 * 1024 * 1024);
    in_use  = reinterpret_cast<int*>(heap + M_ALLOC * sizeof(unsigned long long) + 64 * 1024 * 1024);
    last    = heap + M_ALLOC * sizeof(unsigned long long) * 2 + 64 * 1024 * 1024;
    sam_upper = heap + 1200 * 1024 * 1024;
    lastptr = 0;
    return stack;
  }
  
  EFI_STATUS efi_main(EFI_HANDLE image, EFI_SYSTEM_TABLE *systab) {
    EFI_STATUS status;
    ST = systab;
    BS = ST->BootServices;
    RS = ST->RuntimeServices;
    IH = image;
    BS->SetWatchdogTimer(0, 0, 0, NULL);
    efi_video_init();
    unsigned long long stack(simplealloc_init());
#if defined(IA32)
    asm("movl %0, %%esp; call calc;" : : "r"(stack + 64 * 1024 * 1024 - 64) );
#else
    asm("movq %0, %%rsp; call calc;" : : "r"(stack + 64 * 1024 * 1024 - 64) );
#endif
    return EFI_SUCCESS;
  }
#endif
}

#if defined(_WITH_CLANGPP_)
#if defined(IA32)
void* operator new(unsigned long size);
void* operator new[](unsigned long size);
void* operator new(unsigned long size, void* ptr) noexcept;
void* operator new[](unsigned long size, void* ptr) noexcept;
#else
void* operator new(size_t size);
void* operator new[](size_t size);
void* operator new(size_t size, void* ptr) noexcept;
void* operator new[](size_t size, void* ptr) noexcept;
#endif
#elif defined(_WITH_GPP_)
#if defined(IA32)
inline void* operator new(long unsigned int, void* __p) noexcept { return __p; }
inline void* operator new[](long unsigned int, void* __p) noexcept { return __p; }
#else
inline void* operator new(size_t, void* __p) noexcept { return __p; }
inline void* operator new[](size_t, void* __p) noexcept { return __p; }
#endif
#else
#include <new>
#endif

#define _BURN_ 11
#define _P_RAW_ 1
#define _FLOAT_BITS_ 32
#define _SIMPLEALLOC_ 64
namespace {
#include "cppimport.hh"
#include "lieonn.hh"
typedef myfloat num_t;

static void show3(const num_t& x) {
  cnputc(x == num_t(int(0)) ? '0' : (x < num_t(int(0)) ? '-' : '+'));
}

static void showOut(const num_t& sum, int& ctr, int& tctr) {
  if(sum != num_t(int(0))) tctr ++;
  if(num_t(int(0)) < sum) ctr ++;
  const num_t rr(num_t(2 * ctr) / num_t(max(int(tctr), int(1))) - num_t(int(1)) );
  puthex(tctr);
  cnputc(':');
  puthex(int(rr  * num_t(int(0x1000))) );
  cnputc('\n');
}

}

#if defined(_USERLAND_)
int main(int argc, char* argv[]) {
  unsigned long long heap = (unsigned long long)malloc(1200 * 1024 * 1024);
  unsigned long long stack(heap);
  v_alloc = reinterpret_cast<unsigned long long *>(heap + 32 * 1024 * 1024);
  in_use  = reinterpret_cast<int*>(heap + M_ALLOC * sizeof(unsigned long long) + 32 * 1024 * 1024);
  last    = heap + M_ALLOC * sizeof(unsigned long long) * 2 + 32 * 1024 * 1024;
  sam_upper = heap + 1200 * 1024 * 1024;
  lastptr = 0;
#else
EFI_STATUS calc() {
#endif
  puts("preparing...");
  const int cols(1);
  const int length(loop22<num_t>() + 2 + (1 << cols));
  const num_t bmsqpi(binMargin<num_t>(sqrt(num_t().pi()) / num_t(int(4))));
  puthex(bmsqpi.m); puthex(bmsqpi.e);
#if defined(_USERLAND_)
  const int m('r');
  for(int i = length; 0 < i; i --) pnextcacher<num_t>(i, 1);
#else
 mode:
  int m('r');
  puts("mode? (n for number | r for prng | k for keyboard [a-z])\n");
  while(true) switch(m = efi_cons_getc(0)) {
    case 'n': case 'r': case 'k': goto bbreak; }
 bbreak:
  for(int i = length; 0 < i; i --) pnextcacher<num_t>(i, 1);
#endif
  const int Mlen(max(length + 15, int(sizeof(pdata) / sizeof(int)) ));
  SimpleVector<num_t> M0(_BURN_ * cols);
  M0.O();
  idFeeder<SimpleVector<num_t> > D0(length);
  idFeeder<SimpleVector<num_t> > D1(0);
  D1.res.entity.reserve(Mlen);
#if !defined(_USERLAND_)
  pncr_cp->reserve(Mlen);
#endif
#if 0
    int ctr(0);
    int tctr(0);
    SimpleVector<num_t> M1(M0.size());
    SimpleVector<num_t> M2(M0.size());
    SimpleVector<num_t> M3(M0.size());
    M1.O(); M2.O(); M3.O();
    idFeeder<SimpleVector<num_t> > D2(0);
    idFeeder<SimpleVector<num_t> > D3(0);
    D2.res.entity.reserve(Mlen);
    D3.res.entity.reserve(Mlen);
#endif
    int ctr(0);
    int tctr(0);
  for(int lc = 0; 0 <= lc && lc < Mlen; lc ++) {
    SimpleVector<num_t> vbuf(M0.size() * 3);
    vbuf.O();
    switch(m) {
    case 'n': {
      char buf[0x200];
      int i;
      for(i = 0; i < sizeof(buf) - 1; i ++)
#if defined(_USERLAND_)
        if((buf[i] = getchar()) == '\n') break;
#else
        if((buf[i] = efi_cons_getc(0)) == '\n') break;
#endif
      buf[i] = '\0';
      num_t nx;
      vbuf[0] = nx;
      break;
    } case 'r': {
      for(int i = 0; i < cols; i ++)
        vbuf[i] = myrandom() & 1 ? - num_t(int(1)) : num_t(int(1));
      break;
    } case 'd': {
      if(sizeof(pdata) / sizeof(int) <= lc) goto bbbreak;
      for(int i = 0; i < cols; i ++) {
        vbuf[i] = (num_t(pdata[lc]) - num_t(int(999)) / num_t(int(2))) /
          num_t(int(1000));
        pdata[lc] <<= 1;
        pdata[lc]  |= pdata[lc] >> 10;
        pdata[lc]  &= (1LL << 10) - 1;
      }
      break;
#if !defined(_USERLAND_)
    } case 'k': {
      const int x(efi_cons_getc(0));
      if('a' <= x && x <= 'z') {
        vbuf[0] = sgn<num_t>(num_t(x - 'a') / num_t('z' - 'a') - num_t(int(1)) /
            num_t(int(2)));
      } else goto lnext;
      break;
#endif
    } }
    for(int i = 1; i < _BURN_; i ++) for(int j = 0; j < cols; j ++)
      vbuf[i * cols + j] = myprandom() & 1 ? - vbuf[j] : vbuf[j];
    vbuf.setVector(0, offsetHalf<num_t>(vbuf.subVector(0, M0.size() )) );
    vbuf.setVector(M0.size(), mother<num_t>(vbuf.subVector(0, M0.size() )));
    vbuf.setVector(M0.size() * 2, mother<num_t>(vbuf.subVector(M0.size(),
      M0.size() )));
    D0.next(vbuf);
    vbuf.setVector(0, unOffsetHalf<num_t>(vbuf.subVector(0, M0.size() )) );
    if(D0.full) {
      D1.next(vbuf.subVector(0, M0.size()).tdot(M0));
#if 0
      for(int i = 0; i < cols; i ++) {
        num_t vpsum(int(0));
        num_t rvpsum(int(0));
        for(int j = 0; j < _BURN_; j ++) {
          const int ii(j * cols + i);
          // vpsum  += M2[ii] + M3[ii];
          rvpsum += D1.res[D1.res.size() - 1][ii];
        }
        show3(rvpsum);
        showOut(rvpsum, ctr, tctr);
//        show3(vpsum *= rvpsum);
//        showOut(vpsum, ctr, tctr);
      }
      D2.next(D1.res[D1.res.size() - 1] - M1);
      D3.next(M1);
      M1 = p0maxNext<num_t>(D1.res);
      M2 = p0maxNext<num_t>(D2.res);
      M3 = p0maxNext<num_t>(D3.res);
#endif
      pair<SimpleVector<SimpleVector<num_t> >, num_t> d0(normalizeS<num_t>(
        D0.res) );
      M0 = unOffsetHalf<num_t>(predMasp<num_t, 20>(move(d0.first), 2,
        string("")).subVector(0, M0.size() ) ) * d0.second;
    }
#if !defined(_USERLAND_)
    pncr_cp->resize(length);
#endif
   lnext: ;
  }
#if 1
  {
    SimpleVector<num_t> M1(M0.size());
    SimpleVector<num_t> M2(M0.size());
    SimpleVector<num_t> M3(M0.size());
    M1.O(); M2.O(); M3.O();
    int ctr(0);
    int tctr(0);
    idFeeder<SimpleVector<num_t> > D2(0);
    idFeeder<SimpleVector<num_t> > D3(0);
    D2.res.entity.reserve(Mlen);
    D3.res.entity.reserve(Mlen);
    for(int i0 = 0; i0 < D1.res.size(); i0 ++) {
      for(int i = 0; i < cols; i ++) {
        num_t vpsum(int(0));
        num_t rvpsum(int(0));
        for(int j = 0; j < _BURN_; j ++) {
          const int ii(j * cols + i);
          vpsum  += M2[ii] + M3[ii];
          rvpsum += D1.res[i0][ii];
        }
        show3(vpsum *= rvpsum);
        showOut(vpsum, ctr, tctr);
      }
      D2.next(D1.res[i0] - M1);
      D3.next(M1);
      M1 = p0maxNext<num_t>(D1.res.subVector(0, i0 + 1));
      M2 = p0maxNext<num_t>(D2.res);
      M3 = p0maxNext<num_t>(D3.res);
    }
  }
#endif
 bbbreak:
#if defined(_USERLAND_)
  return 0;
#else
  goto mode;
  return EFI_SUCCESS;
#endif
}

// google search (2026/05/11) results:
/*
// cl.exe /c /Zl /W3 main.c # /GS-
// link.exe /nodefaultlib /subsystem:10 /entry:EfiMain /machine:IX86 /out:boootia32.efi main.obj
// EFI System Table
typedef struct _EFI_SIMPLE_TEXT_OUTPUT_INTERFACE {
    void* Reset;
    long (*OutputString)(struct _EFI_SIMPLE_TEXT_OUTPUT_INTERFACE*, unsigned short*);
} EFI_SIMPLE_TEXT_OUTPUT_INTERFACE;

typedef struct {
    char Header[44];
    EFI_SIMPLE_TEXT_OUTPUT_INTERFACE *ConOut;
} EFI_SYSTEM_TABLE;

long __stdcall EfiMain(void* ImageHandle, EFI_SYSTEM_TABLE *SystemTable) {
  SystemTable->ConOut->OutputString(SystemTable->ConOut, L"Hello VC6 EFI!\r\n");
  while(1);
  return 0;
}

// we abandoned this vc++ 6.0 compatible code outputs bootia32.efi because
// we don't have 2006-2007 imac/macbook nor 2006-2007 UMPC cf. sony which
// boots bootia32.efi in bare metal.
// we can test bootia32.efi with qemu with pflash=... options but this is on
// base system qemu things.
//
// so we test only very easy electro magnetic champer like alminium foil wrapped
// battery powered laptop pc on the calculation also gets same result, we stop
// to implement low level implementations other than this efiboot.cc.
//
// there's another selections for bare metals s.t. some consumer game console
// with un*x systems nor some chips we didn't test : ppc64, riscv64, loongson,
// esp, sh4, etc etc... but test with them is only to ease our mind because
// the terminal we implement and compile MUST not be infected if this is false
// positive also from the infections. in the case, we should calculate raw
// data stream with analogue dice and analogue abacuses with pen and papers,
// which can be exchanged by some of the another actors. 
// --- if we start to doubt anything, there's no end to it. ---
*/

