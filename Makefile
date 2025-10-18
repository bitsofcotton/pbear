CXX=	clang++
#CXX=	eg++
#CXX=	c++

# Thanks to: /sys/arch/amd64/stand/efiboot/Makefile and so on.
# compiler flags.
#CXXFLAGS+=	-O0 -mtune=generic -gfull
#CXXFLAGS+=	-Ofast -mtune=native -gfull
#CXXFLAGS+=	-O3 -mtune=native -g3
CXXFLAGS+=	-Oz -mtune=native -gfull
#CXXFLAGS+=	-O2 -mtune=native -gfull
#CXXFLAGS+=	-O0 -mtune=native -gfull
#CXXFLAGS+=	-O2 -g3
#CXXFLAGS+=      -D_LIBCPP_ENABLE_ASSERTIONS
CXXFLAGS+=	-std=c++11
#CXXFLAGS+=	-std=gnu++98
LDFLAGS+=	-nostdlib -T -Bsymbolic -shared --pack-dyn-relocs=none
CXXFLAGS+=	-ffreestanding -fshort-wchar -fPIC -mno-red-zone
CXXFLAGS+=	-nostdinc -fno-builtin -Wno-pointer-sign -fno-stack-protector
CXXFLAGS+=	-fpack-struct
AFLAGS+=	-pipe -fPIC
OBJFMT=		efi-app-x86_64
LDSCRIT=	ldscript.amd64
LINKADDR=0x40120
LOADADDR=0x40000
HEAP_LIMIT=0xA0000
BOOTREL=0x60000
BOOTMAGIC=0xc001d00d

# lieonn.hh compile options
#CXXFLAGS+=	-D_P_BIT_=1
#CXXFLAGS+=	-D_P_MLEN_=0
#CXXFLAGS+=	-D_P_PRNG_=11
#CXXFLAGS+=	-D_ARCFOUR_

CLEANFILES= *.o pbear.efi

clean:
	@rm -rf ${CLEANFILES}

all:	pbear.efi

pbear.efi:	start_amd64.o run_i386.o mdrandom.o efi.o

