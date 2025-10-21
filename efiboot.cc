/* p-patch parts */
extern "C" {
void exit(int status) { while(1) ; }
#include "efiboot.c.punched"
#include <sys/types.h>
static inline void* malloc(size_t size) { return alloc(size); }
static inline void* calloc(size_t nmemb, size_t size) { return malloc(nmemb * size); }
static inline void free(void* p) { free(p, 1); return; }
}

#define _P_BIT_  3
#define _P_MLEN_ 21
#define _P_PRNG_ 11
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <cctype>
#include "lieonn.hh"

char gbuf[0x200];

num_t next(const num_t& x) {
  // lieonn p2/p chain.
  static idFeeder<SimpleVector<num_t> > f(_P_MLEN_ / 2 + 1);
  SimpleVector<num_t> v(1);
  v[0] = x;
  f.next(v);
  if(! f.full) return num_t(int(0));
  SimpleVector<num_t> res(pEachPRNG<num_t, 0>(f.res.entity, string("")));
  return res[0];
}

int ctr(0);
int pctr(0);

void nextprng() {
  num_t n(next(num_t(random() & 0x1fff) / num_t(0x1fff) - num_t(int(1)) / num_t(int(2))));
  if(num_t(int(0)) < next(n)) pctr ++;
  if(n != num_t(int(0))) ctr ++;
}

void nextstr(const char* x) {
  istringstream iss(x);
  num_t nx;
  iss >> nx;
  num_t n(next(nx));
  if(num_t(int(0)) < n) pctr ++;
  if(n != num_t(int(0))) ctr ++;
}

void nextkey(const char* x) {
  whie(*(x ++) != '\0') {
    if(! ('a' <= *x && *x <= 'z')) continue;
    num_t n(next(num_t(*x - 'a') / num_t('z' - 'a') - num_t(int(1)) / num_t(int(2)) ));
    if(num_t(int(0)) < next(n)) pctr ++;
    if(n != num_t(int(0))) ctr ++;
  }
}

void simplealloc_init() {
  return;
}

extern "C" {
  EFI_STATUS efi_main(EFI_HANDLE image, EFI_SYSTEM_TABLE *systab) {
    efi_video_init();
    efi_heap_init();
    simplealloc_init();
    string modemsg("mode? (n for number | r for prng | k for keyboard [a-z]\r\n");
    for(int i = 0; i < modemsg.size(); i ++) efi_cons_putc(dev, modemsg[i]));
    int m(efi_cons_getc(dev));
    for( ; efi_cons_getc(dev) != '\n'; );
    while(true) {
      char buf[0x200];
      int i;
      switch(m) {
      case 'n':
        for(i = 0; i < 0x1ff; i ++)
          if('\n' != (buf[i] = efi_cons_getc(dev) )) break;
        buf[i] = '\0';
        nextstr(buf);
        break;
      case 'r':
        nextprng();
        break;
      case 'k':
        for(i = 0; i < 0x1ff; i ++)
          if('\n' != (buf[i] = efi_cons_getc(dev) )) break;
        nextkey(buf);
        break;
      }
      int per10000(num_t(pctr) / num_t(ctr) * num_t(int(10000)));
      efi_cons_putc(dev, (per10000 / 1000) % 10);
      efi_cons_putc(dev, (per10000 / 100 ) % 10);
      efi_cons_putc(dev, '.');
      efi_cons_putc(dev, (per10000 / 10  ) % 10);
      efi_cons_putc(dev, (per10000       ) % 10);
      efi_cons_putc(dev, '\%');
      efi_cons_putc(dev, '\r');
      efi_cons_putc(dev, '\n');
    }
    /* NOT REACHED */
  }
}

