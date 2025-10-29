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

num_t next(const num_t& x, idFeeder<SimpleVector<num_t> >& fp, idFeeder<SimpleVector<num_t> >& fm, SimpleVector<num_t>& bb, SimpleVector<num_t>& bp, SimpleVector<num_t>& bm, SimpleVector<num_t>& bs, SimpleVector<num_t>& s0, SimpleVector<num_t>& s1, SimpleMatrix<num_t>& buf) {
  if(! bb.size()) {
    bb.entity.resize(_VDIM_, num_t(int(0)));
    bp.entity.resize(_VDIM_, num_t(int(0)));
    bm.entity.resize(_VDIM_, num_t(int(0)));
    bs.entity.resize(_VDIM_ * 2, num_t(int(0)));
    s0.entity.resize(_VDIM_ * 2, num_t(int(0)));
    s1.entity.resize(_VDIM_ * 2, num_t(int(0)));
    fp.next(  bb);
    fm.next(- bb);
  }
  SimpleVector<num_t> v(_VDIM_);
  for(int i = 0; i < v.size(); i ++)
    v[i] = random() & 1 ? x : - x;
  bb = v * num_t(int(2)) - bb;
  if(! fp.full) {
    fp.next(  v);
    fp.next(  bb);
    fm.next(  v);
    fm.next(- bb);
    return num_t(int(0));
  }
  fp.next(v);
  fm.next(v);
  printf("4");
  vector<SimpleVector<num_t> > ddfp(delta<SimpleVector<num_t> >(delta<SimpleVector<num_t> >(fp.res.entity) ));
  vector<SimpleVector<num_t> > ddfm(delta<SimpleVector<num_t> >(delta<SimpleVector<num_t> >(fm.res.entity) ));
  for(int i = 0; i < ddfp.size(); i ++) ddfp[i] = offsetHalf<num_t>(ddfp[i] /= num_t(int(4)) );
  for(int i = 0; i < ddfm.size(); i ++) ddfm[i] = offsetHalf<num_t>(ddfm[i] /= num_t(int(4)) );
  SimpleVector<SimpleVector<num_t> > rdfp;
  SimpleVector<SimpleVector<num_t> > rdfm;
  rdfp.entity = move(ddfp);
  rdfm.entity = move(ddfm);
  SimpleVector<num_t> p0(unOffsetHalf<num_t>(pGuarantee<num_t, 0>(rdfp, string("") )) );
  printf("\b3");
  SimpleVector<num_t> m0(unOffsetHalf<num_t>(pGuarantee<num_t, 0>(rdfm, string("") )) );
  fp.next(  bb);
  fm.next(- bb);
  printf("\b2");
  ddfp = delta<SimpleVector<num_t> >(delta<SimpleVector<num_t> >(fp.res.entity) );
  ddfm = delta<SimpleVector<num_t> >(delta<SimpleVector<num_t> >(fm.res.entity) );
  for(int i = 0; i < ddfp.size(); i ++) ddfp[i] = offsetHalf<num_t>(ddfp[i] /= num_t(int(4)) );
  for(int i = 0; i < ddfm.size(); i ++) ddfm[i] = offsetHalf<num_t>(ddfm[i] /= num_t(int(4)) );
  rdfp.entity = move(ddfp);
  rdfm.entity = move(ddfm);
  SimpleVector<num_t> p1(unOffsetHalf<num_t>(pGuarantee<num_t, 0>(rdfp, string("") )) );
  printf("\b1");
  SimpleVector<num_t> m1(unOffsetHalf<num_t>(pGuarantee<num_t, 0>(rdfm, string("") )) );
  if(buf.rows() != 3 || buf.cols() != p0.size() * 2) return num_t(int(0));
  for(int i = 1; i < buf.rows(); i ++) buf.row(i - 1) = buf.row(i);
  buf.row(buf.rows() - 1).setVector(0, v * num_t(int(2)) - (bp + bm));
  buf.row(buf.rows() - 1).setVector(v.size(), bp + bm);
  s1 += (s0 += buf.row(buf.rows() - 1) - bs);
  num_t rr(int(0));
  for(int i = 0; i < s1.size() / 2; i ++)
    rr += (s1[i] + s1[i + s1.size() / 2]) * s1[i + s1.size() / 2];
  for(int i = 0; i < bs.size(); i ++)
    bs[i] = pnextcacher<num_t>(buf.rows(), 1).dot(buf.col(i));
  for(int i = 1; i < buf.rows(); i ++) buf.row(i - 1) = buf.row(i);
  buf.row(buf.rows() - 1).setVector(0, v * num_t(int(2)) - (p0 + m0));
  buf.row(buf.rows() - 1).setVector(v.size(), p0 + m0);
  s1 += (s0 += buf.row(buf.rows() - 1) - bs);
  for(int i = 0; i < bs.size(); i ++)
    bs[i] = pnextcacher<num_t>(buf.rows(), 1).dot(buf.col(i));
  bp = p1;
  bm = m1;
  printf("\b");
  return rr;
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
  idFeeder<SimpleVector<num_t> > fp(_P_MLEN_);
  idFeeder<SimpleVector<num_t> > fm(_P_MLEN_);
  SimpleVector<num_t> bb, bp, bm, bs, s0, s1;
  SimpleMatrix<num_t> mbuf(3, _VDIM_ * 2);
  mbuf.O();
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
      if(num_t(int(0)) <
        next(nx, fp, fm, bb, bp, bm, bs, s0, s1, mbuf)) pctr ++;
      break;
    } case 'r': {
      if(num_t(int(0)) <
        next(num_t(random() & 0x1fff) / num_t(0x1fff) - num_t(int(1)) /
          num_t(int(2)), fp, fm, bb, bp, bm, bs, s0, s1, mbuf) )
            pctr ++;
      break;
    } case 'd': {
      if(sizeof(pdata) / sizeof(int) <= lc) goto bbbreak;
      if(num_t(int(0)) <
        next(num_t(pdata[lc]) / num_t(999) - num_t(int(1)) /
          num_t(int(2)), fp, fm, bb, bp, bm, bs, s0, s1, mbuf) )
            pctr ++;
      break;
    } case 'k': {
      const int x(efi_cons_getc(0));
      if('a' <= x && x <= 'z') {
        if(num_t(int(0)) <
          next(num_t(x - 'a') / num_t('z' - 'a') - num_t(int(1)) /
            num_t(int(2)), fp, fm, bb, bp, bm, bs, s0, s1, mbuf) )
              pctr ++;
      } else lc --;
      break;
    } }
    int per10000(num_t(pctr) / num_t(0 < lc ? lc : 1) * num_t(int(10000)));
    printf("%d: %d%c%d\r\n\0", lc, (per10000 / 100 ) % 100, '.', per10000 % 100);
  }
 bbbreak:
  return EFI_SUCCESS;
}

