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

idFeeder<SimpleVector<num_t> > *fb, *ff;
SimpleVector<num_t> *bb;
SimpleVector<num_t> *s0, *s1;
SimpleMatrix<num_t> *pm;
num_t next(const num_t& x) {
  SimpleVector<num_t> v(_P_PRNG_);
  for(int i = 0; i < v.size(); i ++) v[i] = random() & 1 ? - x : x;
  if(! fb) {
    fb = SimpleAllocator<idFeeder<SimpleVector<num_t> > >().allocate(1);
    ff = SimpleAllocator<idFeeder<SimpleVector<num_t> > >().allocate(1);
    bb = SimpleAllocator<SimpleVector<num_t> >().allocate(1);
    s0 = SimpleAllocator<SimpleVector<num_t> >().allocate(1);
    s1 = SimpleAllocator<SimpleVector<num_t> >().allocate(1);
    pm = SimpleAllocator<SimpleMatrix<num_t> >().allocate(1);
    ::new ((void*)fb) idFeeder<SimpleVector<num_t> >(_P_MLEN_);
    ::new ((void*)ff) idFeeder<SimpleVector<num_t> >(_P_MLEN_);
    ::new ((void*)bb) SimpleVector<num_t>(_P_PRNG_);
    ::new ((void*)s0) SimpleVector<num_t>(_P_PRNG_);
    ::new ((void*)s1) SimpleVector<num_t>(_P_PRNG_);
    ::new ((void*)pm) SimpleMatrix<num_t>(3, _P_PRNG_);
    bb->O();
    s0->O();
    s1->O();
    pm->O();
    fb->next(*bb);
    ff->next(*bb);
  }
  (*bb) = v * num_t(int(2)) - (*bb);
  if(! fb->full) {
    fb->next(v);
    ff->next(v);
    fb->next(   *bb );
    ff->next(- (*bb));
    return num_t(int(0));
  }
  num_t M(int(0));
  {
    SimpleMatrix<num_t> off(pm->rows(), pm->cols());
    off.O();
    off.row(off.rows() - 2) = ff->res[ff->res.size() - 2];
    SimpleVector<num_t> buf0(pm->cols());
    SimpleVector<num_t> buf1(pm->cols());
    buf0.O();
    buf1.O();
    for(int i = 0; i < buf0.size(); i ++) buf0[i] = p0maxNext<num_t>(off.col(i) - pm->col(i) / num_t(int(1) << (_P_BIT_ * 2)) );
    for(int i = 0; i < buf1.size(); i ++) buf1[i] = p0maxNext<num_t>(pm->col(i) / num_t(int(1) << (_P_BIT_ * 2)) );
    for(int i = 0; i < buf0.size(); i ++) M += v[i] * buf0[i] * buf1[i] * (* pm)(pm->rows() - 1, i);
    SimpleVector<SimpleVector<num_t> > wwb, wwf;
    wwb.entity = delta<SimpleVector<num_t> >(delta<SimpleVector<num_t> >(fb->next(v).entity));
    wwf.entity = delta<SimpleVector<num_t> >(delta<SimpleVector<num_t> >(ff->next(v).entity));
    pair<SimpleVector<SimpleVector<num_t> >, num_t> wb(normalizeS<num_t>(wwb));
    pair<SimpleVector<SimpleVector<num_t> >, num_t> wf(normalizeS<num_t>(wwf));
    (*s0) += unOffsetHalf<num_t>(pGuarantee<num_t, 0>(offsetHalf<num_t>(wb.first), string(""))) * wb.second;
    (*s0) += unOffsetHalf<num_t>(pGuarantee<num_t, 0>(offsetHalf<num_t>(wf.first), string(""))) * wf.second;
    for(int i = 1; i < pm->rows(); i ++) pm->row(i - 1) = pm->row(i);
    pm->row(pm->rows() - 1) = ((*s1) += (*s0) / num_t(int(2)));
  }
  {
    SimpleVector<SimpleVector<num_t> > wwb, wwf;
    wwb.entity = delta<SimpleVector<num_t> >(delta<SimpleVector<num_t> >(fb->next(   *bb ).entity));
    wwf.entity = delta<SimpleVector<num_t> >(delta<SimpleVector<num_t> >(ff->next(- (*bb)).entity));
    pair<SimpleVector<SimpleVector<num_t> >, num_t> wb(normalizeS<num_t>(wwb));
    pair<SimpleVector<SimpleVector<num_t> >, num_t> wf(normalizeS<num_t>(wwf));
    (*s0) += unOffsetHalf<num_t>(pGuarantee<num_t, 0>(offsetHalf<num_t>(wb.first), string(""))) * wb.second;
    (*s0) += unOffsetHalf<num_t>(pGuarantee<num_t, 0>(offsetHalf<num_t>(wf.first), string(""))) * wf.second;
    for(int i = 1; i < pm->rows(); i ++) pm->row(i - 1) = pm->row(i);
    pm->row(pm->rows() - 1) = ((*s1) += (*s0) / num_t(int(2)));
  }
  return M;
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
  int m('r');
  const num_t sq2(sqrt(num_t(int(2)) ));
  const num_t bmqpi(binMargin<num_t>(num_t().quatpi()));
  for(int i = 0; i < _P_MLEN_; i ++)
    printf("%d:", pnextcacher<num_t>(i + 1, 1).size());
  printf("mem usage temporal efficiency nil: %d, %d, %d, %d\n", sq2.m, sq2.e, bmqpi.m, bmqpi.e);
  printf("mode? (n for number | r for prng | k for keyboard [a-z] | d for pdata.h)\n");
  int ctr(0);
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
      if(num_t(int(0)) < next(nx)) ctr ++;
      break;
    } case 'r': {
      if(num_t(int(0)) <
        next(num_t(random() & 0x1fff) / num_t(0x1fff) - num_t(int(1)) /
          num_t(int(2)) ) ) ctr ++;
      break;
    } case 'd': {
      if(sizeof(pdata) / sizeof(int) <= lc) goto bbbreak;
      if(num_t(int(0)) <
        next(num_t(pdata[lc]) / num_t(999) - num_t(int(1)) /
          num_t(int(2)) ) ) ctr ++;
      break;
    } case 'k': {
      const int x(efi_cons_getc(0));
      if('a' <= x && x <= 'z') {
        if(num_t(int(0)) <
          next(num_t(x - 'a') / num_t('z' - 'a') - num_t(int(1)) /
            num_t(int(2)) ) ) ctr ++;
      } else lc --;
      break;
    } }
    const int per10000(num_t(ctr) / num_t(max(int(lc - 10), int(1))) * num_t(int(10000)));
    printf("%c%d: %d%c%d\r\n\0", m, lc, per10000 / 100, '.', per10000 % 100);
  }
 bbbreak:
  return EFI_SUCCESS;
}

