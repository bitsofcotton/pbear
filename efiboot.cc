/* p-patch parts */
extern "C" {
#include <sys/types.h>
#include "efiboot.c.punched"
  int atexit(void (*function)(void)) { return 0; }
  __dead void panic(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    vprintf(fmt, ap);
    printf("\n");
    va_end(ap);
    _rtt();
    /*NOTREACHED*/
  }
  void assert(const uint64_t f) {
    if(! f) {
      uint64_t* stack;
      __asm__ __volatile__ ("movq %%rsp, %0\n" : "=r" (stack) :: "memory");
      printf("assertion failed\n");
      printf("f: %x\n", (size_t)(void*)assert);
      for(int i = 0; i < 16; i ++)
        printf("s: %x\n", stack[i]);
      for(; efi_cons_getc(0) != ' ';);
      _rtt();
    }
    return;
  }
  struct consdev constab[] = {
    { efi_cons_probe, efi_cons_init, efi_cons_getc, efi_cons_putc },
    { efi_com_probe, efi_com_init, efi_com_getc, efi_com_putc },
    { NULL }
  };
  struct consdev *cn_tab = constab;
}

const bool prngdata[] = {
#include "prng.h"
};

#define _FLOAT_BITS_ 64
#define M_ALLOC (32 * 1024 * 1024)
//#define _BURN_ 3
#define _SIMPLEALLOC_ 64
#include "cppimport.hh"
#include "lieonn.hh"
typedef myfloat num_t;

const int pdata[] = {
#include "pdata.h"
};
#include "cg.h"

extern "C" {
  EFI_STATUS calc() ;
  
  unsigned long long simplealloc_init() {
    unsigned long long heap = 0;
    EFI_STATUS status;
    status = BS->AllocatePages(AllocateAnyPages, EfiLoaderData,
      1700 * 1024 * 1024 / (4 * 1024), &heap);
    if (status != EFI_SUCCESS)
      panic("BS->AllocatePages()");
    unsigned long long stack(heap);
    v_alloc = reinterpret_cast<unsigned long long *>(heap + 64 * 1024 * 1024);
    in_use  = reinterpret_cast<int*>(heap + M_ALLOC * sizeof(unsigned long long) + 64 * 1024 * 1024);
    last    = heap + M_ALLOC * sizeof(unsigned long long) * 2 + 64 * 1024 * 1024;
    sam_upper = heap + 1700 * 1024 * 1024;
    lastptr = 0;
    return stack;
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
    unsigned long long stack(simplealloc_init());
    asm("movq %0, %%rsp; call calc;" :: "r"(stack + 64 * 1024 * 1024 - 64) );
    return EFI_SUCCESS;
  }
}

char n2byte(const num_t c) {
  return max(char(0), min(char(255), char(int(c * num_t(int(255)) )) ));
}

void gm() {
  int gopmode(0);
  UINTN sz, gopsiz, bestsiz = 0;
  EFI_GRAPHICS_OUTPUT_MODE_INFORMATION *gopi;
  for (int i = 0; i < gop->Mode->MaxMode; i++) {
    EFI_STATUS status = gop->QueryMode(gop, i, &sz, &gopi);
    if (EFI_ERROR(status))
      continue;
    gopsiz = gopi->HorizontalResolution *
      gopi->VerticalResolution;
    if (gopsiz > bestsiz) {
      gopmode = i;
      bestsiz = gopsiz;
    }
  }
  if (gopmode >= 0 && gopmode != gop->Mode->Mode) {
    int curmode = gop->Mode->Mode;
    if (efi_gop_setmode(gopmode) != EFI_SUCCESS)
      (void)efi_gop_setmode(curmode);
  }
  return;
}

void setPixel(const int& y, const int& x, const num_t& r, const num_t& g, const num_t& b) {
  const int& height(gop->Mode->Info->VerticalResolution);
  const int& width(gop->Mode->Info->HorizontalResolution);
  EFI_PHYSICAL_ADDRESS fb(gop->Mode->FrameBufferBase);
  reinterpret_cast<unsigned long*>(fb)[(y % height) * width + (x % width)] = n2byte(r) | (n2byte(g) << 8) | (n2byte(b) << 16);
}

EFI_STATUS calc() {
  printf("preparing...\n");
  int ctr(0);
  int tctr(0);
  int length(0);
  int m('r');
  idFeeder<SimpleVector<num_t> > b;
  const num_t bmsqpi(binMargin<num_t>(sqrt(num_t().pi()) / num_t(int(4))));
  printf("mem usage temporal efficiency nil: %d, %d\n", bmsqpi.m, bmsqpi.e);
  pncr_cp = SimpleAllocator<vector<vector<SimpleVector<myfloat> > > >().allocate(1);
  ::new ((void*)pncr_cp) vector<vector<SimpleVector<myfloat> > >();
 mode:
  printf("mode? (n for number | r for prng | k for keyboard [a-z] | d for pdata.h | g for cg.h)\n");
  while(true)
    switch(m = efi_cons_getc(0)) { case 'g': case 'n': case 'r': case'k': case'd': goto bbreak; }
 bbreak:
  if(m == 'g') {
    pncr_cp->resize(cg_n + 1);
    for(int i = 1; i <= cg_n; i ++) {
      (*pncr_cp)[i].resize(2);
      (*pncr_cp)[i][1].resize(i);
      const int& size(i);
      const int  step(1);
      const SimpleVector<num_t> w((dft<num_t>(- size) * (dft<num_t>(size * 2).subMatrix(0, 0, size, size * 2) * taylorc<num_t>(size * 2, num_t(step < 0 ? step * 2 : (size + step) * 2 - 1), num_t(step < 0 ? step * 2 + 2 : (size + step) * 2 - 3)) )).template real<num_t>());
      for(int j = 0; j < w.size(); j ++) (*pncr_cp)[i][1][j] = w[j];
      printf("%d:", (*pncr_cp)[i][1].size());
    }
    SimpleVector<SimpleVector<num_t> > d(cg_n);
    for(int i = 0; i < d.size(); i ++) {
      d[i].resize(cg_height * cg_width);
      for(int j = 0; j < d[i].size(); j ++) d[i][j] = num_t(cg[i][j]);
    }
    SimpleVector<num_t> p(normalize<num_t>(pPRNG<num_t, 0>(normalize<num_t>(d), 8, string("") )));
    printf("done\n");
    for(int i = 0; i < p.size(); i ++) printf("%d ", int(p[i] * num_t(int(255)) ));
    printf("\n");
    gm();
    for(int i = 0; i < p.size() / cg_width; i ++)
      for(int j = 0; j < cg_width; j ++)
        setPixel(i, j, p[i * cg_width + j], p[i * cg_width + j], p[i * cg_width + j]);
    while(efi_cons_getc(0) != ' ') ;
    goto mode;
  }
  printf("mode: %c, length? \n", m);
  while(true) {
    int buf(efi_cons_getc(0) - '0');
    if(0 <= buf && buf < 10) length = length * 10 + buf;
    else break;
  }
  printf("length: %d\n", length);
  pncr_cp->resize(length + 1);
  for(int i = 1; i <= length; i ++) {
    (*pncr_cp)[i].resize(2);
    (*pncr_cp)[i][1].resize(i);
    const int& size(i);
    const int  step(1);
    const SimpleVector<num_t> w((dft<num_t>(- size) * (dft<num_t>(size * 2).subMatrix(0, 0, size, size * 2) * taylorc<num_t>(size * 2, num_t(step < 0 ? step * 2 : (size + step) * 2 - 1), num_t(step < 0 ? step * 2 + 2 : (size + step) * 2 - 3)) )).template real<num_t>());
    for(int j = 0; j < w.size(); j ++) (*pncr_cp)[i][1][j] = w[j];
    printf("%d:", (*pncr_cp)[i][1].size());
  }
  b = idFeeder<SimpleVector<num_t> >(length);
  for(int lc = 0; 0 <= lc; lc ++) {
    SimpleVector<num_t> vbuf(8);
    vbuf.O();
    switch(m) {
    case 'n': {
      char buf[0x200];
      int i;
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
      vbuf[0] = (num_t(pdata[lc]) - num_t(int(999)) / num_t(int(2))) / num_t(int(1000));
      break;
    } case 'k': {
      const int x(efi_cons_getc(0));
      if('a' <= x && x <= 'z') {
        vbuf[0] = num_t(x - 'a') / num_t('z' - 'a') - num_t(int(1)) /
            num_t(int(2));
      } else goto lnext;
      break;
    } }
    for(int j = 1; j < vbuf.size(); j ++)
      vbuf[j] = random() & 1 ? - vbuf[0] : vbuf[0];
    b.next(vbuf);
   lnext:
    if(b.full) {
      SimpleVector<SimpleVector<num_t> > db(delta<SimpleVector<num_t> >(b.res));
      MFENCE();
      SimpleVector<SimpleVector<num_t> > p(
        pPRNGM<num_t, 0>(offsetHalf<num_t>(db), 8, string("") ));
      MFENCE();
      for(int i = 0; i < db.size() - (p.size() - 1); i ++) p[0] += db[i];
      MFENCE();
      for(int i = 1; i < p.size(); i ++) p[i] += p[i - 1];
      MFENCE();
      for(int i = 0; i < p.size() - 1; i ++) {
        num_t j(int(0));
        for(int k = 0; k < p[i].size(); k ++)
          j += p[i][k] * b.res[i - (p.size() - 1) + b.res.size()][k];
        tctr ++;
        if(j == num_t(int(0))) tctr --;
        if(num_t(int(0)) < j) ctr ++;
        const num_t rr(num_t(2 * ctr) / num_t(max(int(tctr), int(1))) - num_t(int(1)) );
        printf("%c%d: %c0.%04d\r\n\0", m, lc, rr < num_t(int(0)) ? '-' : ' ',
          int(abs(rr) * num_t(int(10000)) ) );
      }
      b.t = 0;
      b.full = 0;
    }
  }
 bbbreak:
  goto mode;
  return EFI_SUCCESS;
}

