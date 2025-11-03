/* p-patch parts */
extern "C" {
#include <sys/types.h>
void addbootarg(int, size_t, void *);
#include "efiboot.c.punched"
int atexit(void (*function)(void)) { return 0; }
}

#define M_ALLOC (512 * 1024)
#define assert (void)

#define _P_BIT_  3
#define _P_MLEN_ 21
#define _P_PRNG_ 11
#define _SIMPLEALLOC_ 64
#include "cppimport.hh"
#include "lieonn.hh"
typedef myfloat num_t;

const int pdata[] = {
#include "pdata.h"
};

#define _VDIM_ 11

num_t next(const num_t& x, idFeeder<SimpleVector<num_t> >& f) {
  SimpleVector<num_t> v(1);
  v[0] = x;
  if(! f.full) {
    f.next(v);
    return num_t(int(0));
  }
  const num_t M(unOffsetHalf<num_t>(pEachPRNG<num_t, 0>(offsetHalf<num_t>(f.res.entity), string(""))[0]));
  f.next(v);
  return x * M;
}

extern "C" {
  EFI_STATUS calc() ;
  
  void simplealloc_init() {
    unsigned long long heap = 0;
    EFI_STATUS status;
    status = BS->AllocatePages(AllocateAnyPages, EfiLoaderData,
      768 * 1024 * 1024 / (4 * 1024), &heap);
    if (status != EFI_SUCCESS)
      panic("BS->AllocatePages()");
    v_alloc = reinterpret_cast<unsigned long long *>(heap);
    in_use  = reinterpret_cast<int*>(heap + M_ALLOC * sizeof(unsigned long long));
    last    = heap + M_ALLOC * sizeof(unsigned long long) * 2;
    sam_upper = heap + 768 * 1024 * 1024;
    lastptr = 0;
    return;
  }
  
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
}

EFI_STATUS calc() {
  printf("preparing...\n");
  char buf[0x200];
  int m('r');
  const num_t sq2(sqrt(num_t(int(2)) ));
  const num_t bmqpi(binMargin<num_t>(num_t().quatpi()));
  for(int i = 0; i < _P_MLEN_; i ++)
    printf("%d:", pnextcacher<num_t>(i + 1, 1).size());
  printf("mem usage temporal efficiency nil: %d, %d, %d, %d\n", sq2.m, sq2.e, bmqpi.m, bmqpi.e);
  printf("mode? (n for number | r for prng | k for keyboard [a-z] | d for pdata.h)\n");
  idFeeder<SimpleVector<num_t> > f((_P_MLEN_ + 1) / 2);
  int pctr(0);
  while(true)
    switch(m = efi_cons_getc(0)) { case 'n': case 'r': case'k': case'd': goto bbreak; }
 bbreak:
  for(int lc = 0; ; lc ++) {
    int i;
    switch(m) {
    case 'n': {
      for(i = 0; i < sizeof(buf) - 1; i ++)
        if((buf[i] = efi_cons_getc(0)) == '\n') break;
      buf[i] = '\0';
      num_t nx;
      if(num_t(int(0)) < next(nx, f)) pctr ++;
      break;
    } case 'r': {
      if(num_t(int(0)) <
        next(num_t(random() & 0x1fff) / num_t(0x1fff) - num_t(int(1)) /
          num_t(int(2)), f) ) pctr ++;
      break;
    } case 'd': {
      if(sizeof(pdata) / sizeof(int) <= lc) goto bbbreak;
      if(num_t(int(0)) <
        next(num_t(pdata[lc]) / num_t(999) - num_t(int(1)) /
          num_t(int(2)), f) ) pctr ++;
      break;
    } case 'k': {
      const int x(efi_cons_getc(0));
      if('a' <= x && x <= 'z') {
        if(num_t(int(0)) <
          next(num_t(x - 'a') / num_t('z' - 'a') - num_t(int(1)) /
            num_t(int(2)), f) ) pctr ++;
      } else lc --;
      break;
    } }
    int per10000(num_t(pctr) / num_t((_P_MLEN_ + 3) / 2 < lc ? lc - (_P_MLEN_ + 3) / 2 : 1) * num_t(int(10000)));
    printf("%c%d: %d%c%d\r\n\0", m, lc, (per10000 / 100 ) % 100, '.', per10000 % 100);
  }
 bbbreak:
  return EFI_SUCCESS;
}

