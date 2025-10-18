/* p-patch parts */
#include "lieonn.hh"

char gbuf[0x200];

num_t next(const num_t& x) {
  // lieonn p2/p chain.
}

char* nextstr(const char* x) {
  istringstream iss(x);
  num_t nx;
  iss >> nx;
  num_t res(next(nx));
  stringstream ss(res);
  int i;
  for(i = 0; i < gbuf.size() && i < 0x1ff; i ++) gbuf[i] = ss.str[i];
  gbuf[i] = '\0';
  return gbuf;
}

void simplealloc_init() {
  return;
}

extern "C" {
  EFI_STATUS efi_main(EFI_HANDLE image, EFI_SYSTEM_TABLE *systab) {
    efi_video_init();
    efi_heap_init();
    simplealloc_init();
    while(true) {
      char buf[0x200];
      int i;
      for(i = 0; i < 0x1ff; i ++) if(! (buf[i] = efi_cons_getc(dev) )) break;
      buf[i] = '\0';
      char* pn(nextstr(buf));
      for( ; pn; pn ++) efi_cons_putc(dev, pn);
    }
    /* NOT REACHED */
  }
}

