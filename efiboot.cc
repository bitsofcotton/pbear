/* p-patch parts */
extern "C" {
#include <sys/types.h>
void addbootarg(int, size_t, void *);
#include "efiboot.c.punched"
static inline void* malloc(size_t size) { return alloc(size); }
static inline void* calloc(size_t nmemb, size_t size) { return malloc(nmemb * size); }
static inline void free(void* p) { free(p, 1); return; }
int atexit(void (*function)(void)) { return 0; }
}

#define assert (void)

#define _P_BIT_  3
#define _P_MLEN_ 21
#define _P_PRNG_ 11
#define _SIMPLEALLOC_ 64
#define M_ALLOC 65536
#include "cppimport.hh"
#include "lieonn.hh"
typedef myfloat num_t;

num_t next(const num_t& x, idFeeder<SimpleVector<num_t> >& f) {
  // lieonn p2/p chain.
  SimpleVector<num_t> v(1);
  v[0] = x;
  f.next(v);
  if(! f.full) return num_t(int(0));
  // XXX: we need this:
  for(int i = 0; i < f.res.size(); i ++) printf("%d,", f.res[i].size());
  printf("\n");
  // XXX: 
  SimpleVector<num_t> res(pEachPRNG<num_t, 0>(f.res.entity, string("")));
  return res[0];
}

int ctr(0);
int pctr(0);

num_t nextprng(idFeeder<SimpleVector<num_t> >& f, const num_t& b) {
  num_t nb(num_t(random() & 0x1fff) / num_t(0x1fff) - num_t(int(1)) / num_t(int(2)) );
  num_t n(next(nb, f));
  if(num_t(int(0)) < b * nb) pctr ++;
  if(n != num_t(int(0))) ctr ++;
  return n;
}

num_t nextstr(const char* x, idFeeder<SimpleVector<num_t> >& f, const num_t& b) {
  num_t nx;
  num_t n(next(nx, f));
  if(num_t(int(0)) < b * nx) pctr ++;
  if(n != num_t(int(0))) ctr ++;
  return n;
}

num_t nextkey(const char* x, idFeeder<SimpleVector<num_t> >& f, const num_t& b) {
  num_t nb;
  while(*(x ++) != '\0') {
    if(! ('a' <= *x && *x <= 'z')) continue;
    num_t nx(num_t(*x - 'a') / num_t('z' - 'a') - num_t(int(1)) / num_t(int(2)));
    if(num_t(int(0)) < nb * nx) pctr ++;
    num_t n(next(nx, f));
    if(n != num_t(int(0))) ctr ++;
    nb = n;
  }
  return nb;
}

void simplealloc_init() {
  last    = reinterpret_cast<size_t>(efi_loadaddr);
  lastptr = 0;
  return;
}

extern "C" {
  EFI_STATUS calc() ;
  
  EFI_STATUS efi_main(EFI_HANDLE image, EFI_SYSTEM_TABLE *systab) {
    EFI_DEVICE_PATH  *dp0 = NULL;
    EFI_LOADED_IMAGE *imgp;
    EFI_STATUS status;
    ST = systab;
    BS = ST->BootServices;
    RS = ST->RuntimeServices;
    IH = image;
    BS->SetWatchdogTimer(0, 0, 0, NULL);
    efi_video_init();
    efi_heap_init();
    status = BS->HandleProtocol(image, &imgp_guid, (void **)&imgp);
    if (status == EFI_SUCCESS)
            status = BS->HandleProtocol(imgp->DeviceHandle, &devp_guid,
                (void **)&dp0);
    efi_memprobe();
    simplealloc_init();
    asm("movq %0, %%rsp; call calc;" :: "r"(efi_loadaddr + KERN_LOADSPACE_SIZE - 8) );
    return EFI_SUCCESS;
  }
  
  EFI_STATUS calc() {
    printf("mode? (n for number | r for prng | k for keyboard [a-z]\r\n\0");
    //int m(efi_cons_getc(0));
    int m('r');
    int lc;
    idFeeder<SimpleVector<num_t> > f(_P_MLEN_ / 2 + 1);
    num_t b;
    for(lc = 0; ; lc ++) {
      char buf[0x200];
      int i;
      switch(m) {
      case 'n':
/*
        for(i = 0; i < 0x1ff; i ++)
          if('\n' != (buf[i] = efi_cons_getc(dev) )) break;
*/
        buf[i] = '\0';
        b = nextstr(buf, f, b);
        break;
      case 'r':
        b = nextprng(f, b);
        break;
      case 'k':
/*
        for(i = 0; i < 0x1ff; i ++)
          if('\n' != (buf[i] = efi_cons_getc(dev) )) break;
*/
        b = nextkey(buf, f, b);
        break;
      }
      int per10000(num_t(pctr) / num_t(ctr) * num_t(int(10000)));
      printf("%d: %d%c%d, %d\r\n\0", lc, (per10000 / 100 ) % 100, '.', per10000 % 100, ctr);
    }
    return EFI_SUCCESS;
  }
}

