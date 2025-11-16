/* p-patch parts */
extern "C" {
#include <sys/types.h>
void addbootarg(int, size_t, void *);
#include "efiboot.c.punched"
int atexit(void (*function)(void)) { return 0; }
}

#define M_ALLOC (1024 * 1024)
#define assert (void)

#define _P_PRNG_ 11
#define _SIMPLEALLOC_ 64
#include "cppimport.hh"
#include "lieonn.hh"
typedef myfloat num_t;

const int pdata[] = {
#include "pdata.h"
};

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
    unsigned long long stack = 0;
    status = BS->AllocatePages(AllocateAnyPages, EfiLoaderData,
      8 * 1024 * 1024 / (4 * 1024), &stack);
    if (status != EFI_SUCCESS)
      panic("BS->AllocatePages()");
    asm("movq %0, %%rsp; call calc;" :: "r"(stack + 8 * 1024 * 1024 - 32) );
    return EFI_SUCCESS;
  }
}

EFI_STATUS calc() {
  printf("preparing...\n");
  char buf[0x200];
  int ctr(0);
  int tctr(0);
  int length(0);
  int m('r');
  idFeeder<SimpleVector<num_t> > b;
  const num_t sq2(sqrt(num_t(int(2)) ));
  const num_t bmqpi(binMargin<num_t>(num_t().quatpi()));
  printf("mem usage temporal efficiency nil: %d, %d, %d, %d\n", sq2.m, sq2.e, bmqpi.m, bmqpi.e);
  printf("mode? (n for number | r for prng | k for keyboard [a-z] | d for pdata.h)\n");
  while(true)
    switch(m = efi_cons_getc(0)) { case 'n': case 'r': case'k': case'd': goto bbreak; }
 bbreak:
  printf("mode: %c, length? \n", m);
  while(true) {
    int buf(efi_cons_getc(0) - '0');
    if(0 <= buf && buf < 10) length = length * 10 + buf;
    else break;
  }
  printf("length: %d\n", length);
  pncr_cp = SimpleAllocator<vector<vector<SimpleVector<myfloat> > > >().allocate(1);
  ::new ((void*)pncr_cp) vector<vector<SimpleVector<myfloat> > >();
  pncr_cp->resize(length + 1);
  for(int i = 1; i <= length; i ++) {
    pnextcacher<num_t>(i, 1);
    printf("%d:", i);
  }
  b = idFeeder<SimpleVector<num_t> >(length);
  for(int lc = 0; 0 <= lc; lc ++) {
    SimpleVector<num_t> vbuf(1);
    vbuf[0] = num_t(int(0));
    int i;
    switch(m) {
    case 'n': {
      for(i = 0; i < sizeof(buf) - 1; i ++)
        if((buf[i] = efi_cons_getc(0)) == '\n') break;
      buf[i] = '\0';
      num_t nx;
      vbuf[0] = nx;
      break;
    } case 'r': {
      vbuf[0] = num_t(random() & 0x1fff) / num_t(0x1fff) - num_t(int(1)) /
          num_t(int(2));
      break;
    } case 'd': {
      if(sizeof(pdata) / sizeof(int) <= lc) goto bbbreak;
      vbuf[0] = num_t(pdata[lc]) / num_t(999) - num_t(int(1)) /
          num_t(int(2));
      break;
    } case 'k': {
      const int x(efi_cons_getc(0));
      if('a' <= x && x <= 'z') {
        vbuf[0] = num_t(x - 'a') / num_t('z' - 'a') - num_t(int(1)) /
            num_t(int(2));
      } else goto lnext;
      break;
    } }
    b.next(vbuf);
   lnext:
    if(b.full) {
      pair<SimpleVector<SimpleVector<num_t> >, num_t> work(normalizeS(b.res));
      work.first = delta<SimpleVector<num_t> >(work.first);
      for(int i = 0; i < work.first.size(); i ++) work.first[i] /= num_t(int(2));
      // XXX:
      b.res = work.first;
      for(int i = 0; i < work.first.size(); i += 2)
        work.first[i] = - work.first[i];
      SimpleVector<SimpleVector<num_t> > p(pPRNG0<num_t, 0, false>(
        offsetHalf<num_t>(work.first), 10, string("") ));
      for(int i = 0; i < p.size(); i ++)
        p[i] = unOffsetHalf<num_t>(p[i] * work.second);
      SimpleVector<SimpleVector<num_t> > w(p.size() - 1);
      for(int i = 0; i < w.size(); i ++) {
        w[i].resize(p[i].size() * 2);
        w[i].setVector(0, b.res[i - w.size() + b.res.size()] - p[i]);
        w[i].setVector(p[i].size(), p[i]);
      }
      for(int i = 1; i < w.size(); i ++) w[i] += w[i - 1];
      SimpleVector<SimpleVector<num_t> > x(w.size());
      for(int i = 0; i < w.size(); i ++) {
        x[i].resize(w[i].size() / 2);
        for(int j = 0; j < x[i].size(); j ++)
          x[i][j] = (w[i][j] + w[i][w[i].size() / 2 + j]) * w[i][w[i].size() / 2 + j];
      }
      for(int i = x.size() - 1; 0 < i; i --)
        for(int j = 0; j < x[i].size(); j ++) x[i][j] *= x[i - 1][j];
      x[0].O();
      for(int i = 0; i < x.size(); i ++) {
        const num_t j(x[i][0]);
        if(j == num_t(int(0))) continue;
        else if(num_t(int(0)) < j) ctr ++;
        tctr ++;
        const int per10000(num_t(ctr) / num_t(max(int(tctr), int(1))) * num_t(int(10000)));
        printf("%c%d: %d%c%d\r\n\0", m, lc, per10000 / 100, '.', per10000 % 100);
      }
      b.t = 0;
      b.full = 0;
    }
  }
 bbbreak:
  return EFI_SUCCESS;
}

