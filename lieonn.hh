/* You can use one of the both BSD 3-Clause License or GNU Lesser General Public License 3.0 for this source. */
/*
BSD 3-Clause License

Copyright (c) 2013-2025, bitsofcotton (kazunobu watatsu)
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

// N.B. there exists jemalloc/mimalloc or so.
//      if some of the performance regression occures, we can use them as
//      a vast performance increase.
#if !defined(_SIMPLELIN_)

// --- N.B. start approximate Lie algebra on F_2^k. ---
// N.B. start ifloat
// Double int to new int class.
template <typename T, int bits> class DUInt {
public:
  inline DUInt(const int& src = 0) {
    const int abssrc(src < 0 ? - src : src);
    e[0]   = T(abssrc);
    e[1]  ^= e[1];
    if(abssrc != src)
      *this = - *this;
  }
  inline DUInt(const T& src) {
    const T abssrc(src < T(int(0)) ? - src : src);
    e[0]   = abssrc;
    e[1]  ^= e[1];
    if(abssrc != src)
      *this = - *this;
  }
  inline DUInt(const DUInt<T,bits>& src) { *this = src; }
  inline DUInt(const DUInt<DUInt<T,bits>,bits*2>& src) { *this = src; }
  inline DUInt(DUInt<T,bits>&& src) { *this = src; }
  inline ~DUInt() { ; }
  inline DUInt<T,bits>& operator ++ () {
    ++ e[0];
    if(!e[0]) ++ e[1];
    return *this;
  }
  inline DUInt<T,bits>  operator ++ (int32_t) {
    const DUInt<T,bits> work(*this);
    ++ *this;
    return work;
  }
  inline DUInt<T,bits>& operator -- () {
    if(!e[0]) -- e[1];
    -- e[0];
    return *this;
  }
  inline DUInt<T,bits>  operator -- (int32_t) {
    const DUInt<T,bits> work(*this);
    -- *this;
    return work;
  }
  inline DUInt<T,bits>  operator -  () const {
    DUInt<T,bits> work(~ *this);
    return ++ work;
  }
  inline DUInt<T,bits>  operator +  (const DUInt<T,bits>& src) const {
    DUInt<T,bits> work(*this);
    return work += src;
  }
  inline DUInt<T,bits>& operator += (const DUInt<T,bits>& src) {
    // N.B. assembler can boost dramatically this code. but not here.
    const T e0(max(e[0], src.e[0]));
    e[0] += src.e[0];
    if(e[0] < e0)
      e[1] ++;
    e[1] += src.e[1];
    return *this;
  }
  inline DUInt<T,bits>  operator -  (const DUInt<T,bits>& src) const {
    DUInt<T,bits> work(*this);
    return work -= src;
  }
  inline DUInt<T,bits>& operator -= (const DUInt<T,bits>& src) {
    return *this += - src;
  }
  inline DUInt<T,bits>  operator *  (const DUInt<T,bits>& src) const {
    DUInt<T,bits> result;
    result ^= result;
    for(int i = 0; i < 2 * bits; i ++)
      if(int(src >> i) & 1)
        result += *this << i;
    // N.B.
    //   If we work with multiply with table and summing up with simple window,
    //   and the parallel condition, we can reduce better:
    //     with bit pair [x1, x2, ..., xn], [y1, y2, ..., yn],
    //     make table [[x1*y1, ..., x1*yn],...,[xn*y1, ... xn*yn]],
    //     then counter orthogonal sum-up with parallel.
    //     we get [z1, ... zn] == [x1*y1, ..., sum_i+j=k(x_i*y_j), ..., xn*yn],
    //     then sum-up with certain bit adder and fixing one by one:
    //     r1 := x1*y1, s1 := ((x1*y1) >> 1) + z2, r2 := s2 & 1, ... and so on.
    return result;
  }
  inline DUInt<T,bits>& operator *= (const DUInt<T,bits>& src) {
    return *this = *this * src;
  }
  inline DUInt<T,bits>  operator /  (const DUInt<T,bits>& src) const {
    DUInt<T,bits> work(*this);
    return work /= src;
  }
  inline DUInt<T,bits>& operator /= (const DUInt<T,bits>& src) {
    const DUInt<T,bits> one(int(1));
    if(! src) return *this;
    if(! *this)
      return *this;
    DUInt<T,bits> cache(*this);
    *this ^= *this;
    for(int i = 2 * bits - 1; 0 <= i; i --)
      if((cache >> i) >= src) {
        *this |= one << i;
        cache -= src << i;
      }
    return *this;
    // N.B. if we works with newton's method, better speed will be gained.
  }
  inline DUInt<T,bits>  operator %  (const DUInt<T,bits>& src) const {
    return *this - ((*this / src) * src);
  }
  inline DUInt<T,bits>& operator %= (const DUInt<T,bits>& src) {
    return *this = *this % src;
  }
  inline DUInt<T,bits>  operator << ( const int& b)            const {
    DUInt<T,bits> work(*this);
    return work <<= b;
  }
  inline DUInt<T,bits>& operator <<= (const int& b) {
    if(! b)                return *this;
    else if(b < 0)         return *this >>= (- b);
    else if(b >= bits * 2) return *this ^= *this;
    else if(b > bits) {
      e[1]  = e[0] << (b - bits);
      e[0] ^= e[0];
    } else if(b == bits) {
      e[1]  = e[0];
      e[0] ^= e[0];
    } else {
      e[1] <<= b;
      e[1]  |= e[0] >> (bits - b);
      e[0] <<= b;
    }
    return *this;
  }
  inline DUInt<T,bits>  operator >> ( const int& b)            const {
    DUInt<T,bits> work(*this);
    return work >>= b;
  }
  inline DUInt<T,bits>& operator >>= (const int& b) {
    if(! b)                return *this;
    else if(b < 0)         return *this <<= (- b);
    else if(b >= bits * 2) return *this ^= *this;
    else if(b > bits) {
      e[0]  = e[1] >> (b - bits);
      e[1] ^= e[1];
    } else if(b == bits) {
      e[0]  = e[1];
      e[1] ^= e[1];
    } else {
      e[0] >>= b;
      e[0]  |= e[1] << (bits - b);
      e[1] >>= b;
    }
    return *this;
  }
  inline DUInt<T,bits>  operator &  (const DUInt<T,bits>& src) const {
    DUInt<T,bits> work(*this);
    return work &= src;
  }
  inline DUInt<T,bits>& operator &= (const DUInt<T,bits>& src) {
    e[0] &= src.e[0]; e[1] &= src.e[1];
    return *this;
  }
  inline DUInt<T,bits>  operator |  (const DUInt<T,bits>& src) const {
    DUInt<T,bits> work(*this);
    return work |= src;
  }
  inline DUInt<T,bits>& operator |= (const DUInt<T,bits>& src) {
    e[0] |= src.e[0]; e[1] |= src.e[1];
    return *this;
  }
  inline DUInt<T,bits>  operator ^  (const DUInt<T,bits>& src) const {
    DUInt<T,bits> work(*this);
    return work ^= src;
  }
  inline DUInt<T,bits>& operator ^= (const DUInt<T,bits>& src) {
    e[0] ^= src.e[0]; e[1] ^= src.e[1];
    return *this;
  }
  inline DUInt<T,bits>  operator ~  ()                         const {
    DUInt<T,bits> work;
    work.e[0] = ~ e[0]; work.e[1] = ~ e[1];
    return work;
  }
  inline DUInt<T,bits>& operator =  (const DUInt<T,bits>& src) {
    e[0] = src.e[0]; e[1] = src.e[1];
    return *this;
  }
  inline DUInt<T,bits>& operator =  (const DUInt<DUInt<T,bits>,bits*2>& src) {
    return *this = src.e[0];
  }
  inline DUInt<T,bits>& operator =  (DUInt<T,bits>&& src) {
    e[0] = move(src.e[0]); e[1] = move(src.e[1]);
    return *this;
  }
  inline bool           operator <  (const DUInt<T,bits>& src) const {
    if(e[1]) return e[1] != src.e[1] ? e[1] < src.e[1] : e[0] < src.e[0];
    return bool(src.e[1]) || e[0] < src.e[0];
  }
  inline bool           operator <= (const DUInt<T,bits>& src) const {
    return *this < src || *this == src;
  }
  inline bool           operator >  (const DUInt<T,bits>& src) const {
    return ! (*this <= src);
  }
  inline bool           operator >= (const DUInt<T,bits>& src) const {
    return ! (*this < src);
  }
  inline bool           operator == (const DUInt<T,bits>& src) const {
    return ! (*this != src);
  }
  inline bool           operator != (const DUInt<T,bits>& src) const {
    return (*this ^ src).operator bool();
  }
  inline bool           operator && (const DUInt<T,bits>& src) const {
    return this->operator bool() && src.operator bool();
  }
  inline bool           operator || (const DUInt<T,bits>& src) const {
    return this->operator bool() || src.operator bool();
  }
  inline bool           operator !    () const {
    return ! this->operator bool();
  }
  inline                operator bool () const {
    return e[0] || e[1];
  }
  inline                operator int  () const {
    return int(e[0]);
  }
  inline                operator T    () const {
    return e[0];
  }

  T e[2];
};

// add sign.
template <typename T, int bits> class Signed : public T {
public:
  inline Signed() { ; }
  inline Signed(const int& src) {
    T tsrc(src);
    *this = reinterpret_cast<const Signed<T,bits>&>(tsrc);
  }
  inline Signed(const T& src) {
    *this = reinterpret_cast<const Signed<T,bits>&>(src);
  }
  inline Signed(const Signed<T,bits>& src) {
    *this = src;
  }
  inline bool operator != (const Signed<T,bits>& src) const {
    return dynamic_cast<const T&>(*this) != dynamic_cast<const T&>(src);
  }
  inline bool operator <  (const Signed<T,bits>& src) const {
    const int mthis(int(*this >> (bits - 1)));
    const int msrc( int(src   >> (bits - 1)));
    if(mthis ^ msrc) return mthis;
    if(mthis)
      return - dynamic_cast<const T&>(src) < - dynamic_cast<const T&>(*this);
    return dynamic_cast<const T&>(*this) < dynamic_cast<const T&>(src);
  }
  inline bool operator <= (const Signed<T,bits>& src) const {
    return ! (*this > src);
  }
  inline bool operator >  (const Signed<T,bits>& src) const {
    return ! (*this < src) && *this != src;
  }
  inline bool operator >= (const Signed<T,bits>& src) const {
    return ! (*this < src);
  }
  inline      operator double () const {
    if(*this < Signed<T,bits>(T(int(0))) ) {
      Signed<T,bits> mthis(- *this);
      T work(dynamic_cast<const T&>(mthis));
      return - double(work);
    }
    return double(dynamic_cast<const T&>(*this));
  }
};

// integer to integer float part.
template <typename T, typename W, int bits, typename U> class SimpleFloat {
public:
  inline SimpleFloat() {
    s |= (1 << NaN) | (1 << INF);
  }
  template <typename V> inline SimpleFloat(const V& src) {
    const V vzero(int(0));
    s ^= s;
    m  = T(int(src < vzero ? - src : src));
    e ^= e;
    s |= safeAdd(e, normalize(m));
    if(src < vzero)
      s |= 1 << SIGN;
    ensureFlag();
  }
  inline SimpleFloat(const SimpleFloat<T,W,bits,U>& src) { *this = src; }
  inline SimpleFloat(SimpleFloat<T,W,bits,U>&& src) { *this = src; }
  inline ~SimpleFloat() { ; }
  inline SimpleFloat<T,W,bits,U>  operator -  () const {
    SimpleFloat<T,W,bits,U> work(*this);
    work.s ^= 1 << SIGN;
    return work;
  }
  inline SimpleFloat<T,W,bits,U>  operator +  (const SimpleFloat<T,W,bits,U>& src) const {
    SimpleFloat<T,W,bits,U> work(*this);
    return work += src;
  }
         SimpleFloat<T,W,bits,U>& operator += (const SimpleFloat<T,W,bits,U>& src) {
    if((s |= src.s & (1 << NaN)) & (1 << NaN)) return *this;
    if(s & (1 << INF)) {
      if((src.s & (1 << INF)) && (s ^ src.s) & (1 << SIGN)) s |= 1 << NaN;
      return *this;
    }
    if(src.s & (1 << INF)) return *this = src;
    if(! m) return *this = src;
    if(! src.m) return *this;
    if(! ((s ^ src.s) & (1 << SIGN))) {
      if(e >= src.e) {
        m >>= 1;
        s |= safeAdd(e, 1);
        U se(e);
        if(! safeAdd(se, - src.e) && se < U(bits))
          m += se ? src.m >> int(se) : src.m;
      } else
        return *this = src + *this;
    } else {
      if(e > src.e) {
        U se(e);
        if(! safeAdd(se, - src.e) && se < U(bits))
          m -= se ? src.m >> int(se) : src.m;
      } else if(e == src.e) {
        if(m >= src.m)
          m -= src.m;
        else
          return *this = src + *this;
      } else
        return *this = src + *this;
    }
    s |= safeAdd(e, normalize(m));
    return ensureFlag();
  }
  inline SimpleFloat<T,W,bits,U>  operator -  (const SimpleFloat<T,W,bits,U>& src) const {
    SimpleFloat<T,W,bits,U> work(*this);
    return work -= src;
  }
  inline SimpleFloat<T,W,bits,U>& operator -= (const SimpleFloat<T,W,bits,U>& src) {
    s ^= 1 << SIGN;
    *this += src;
    s ^= 1 << SIGN;
    return *this;
  }
  inline SimpleFloat<T,W,bits,U>  operator *  (const SimpleFloat<T,W,bits,U>& src) const {
    SimpleFloat<T,W,bits,U> work(*this);
    return work *= src;
  }
         SimpleFloat<T,W,bits,U>& operator *= (const SimpleFloat<T,W,bits,U>& src) {
    s ^= src.s & (1 << SIGN);
    if((s |= src.s & (1 << NaN)) & (1 << NaN)) return *this;
    if((! m) || (! src.m)) {
      s |= 1 << DWRK;
      return ensureFlag();
    }
    if((s |= src.s & (1 << INF)) & (1 << INF)) return *this;
    W mm(W(m) * W(src.m));
    s |= safeAdd(e, src.e);
    s |= safeAdd(e, normalize(mm));
    s |= safeAdd(e, U(bits));
    m  = T(mm >> bits);
    return ensureFlag();
  }
  inline SimpleFloat<T,W,bits,U>  operator /  (const SimpleFloat<T,W,bits,U>& src) const {
    SimpleFloat<T,W,bits,U> work(*this);
    return work /= src;
  }
         SimpleFloat<T,W,bits,U>& operator /= (const SimpleFloat<T,W,bits,U>& src) {
    s ^= src.s & (1 << SIGN);
    if((s |= src.s & (1 << NaN)) & (1 << NaN)) return *this;
    if(! (s & (1 << INF)) && (src.s & (1 << INF))) {
      s |= 1 << DWRK;
      return ensureFlag();
    }
    if(s & (1 << INF)) {
      if(src.s & (1 << INF)) s |= 1 << NaN;
      return *this;
    }
    if(! src.m) {
      s |= 1 << NaN;
      return *this;
    }
    if(! m) return *this;
    W mm((W(m) << bits) / W(src.m));
    s |= safeAdd(e, - src.e);
    s |= safeAdd(e, normalize(mm));
    m  = T(mm >> bits);
    return ensureFlag();
  }
  inline SimpleFloat<T,W,bits,U>  operator %  (const SimpleFloat<T,W,bits,U>& src) const {
    return *this - (*this / src).absfloor() * src;
  }
  inline SimpleFloat<T,W,bits,U>& operator %= (const SimpleFloat<T,W,bits,U>& src) {
    return *this = *this % src;
  }
  inline SimpleFloat<T,W,bits,U>  operator <<  (const U& b) const {
    SimpleFloat<T,W,bits,U> work(*this);
    return work <<= b;
  }
  inline SimpleFloat<T,W,bits,U>& operator <<= (const U& b) {
    if(s & ((1 << INF) | (1 << NaN))) return *this;
    s |= safeAdd(e, b);
    return ensureFlag();
  }
  inline SimpleFloat<T,W,bits,U>  operator >>  (const U& b) const {
    SimpleFloat<T,W,bits,U> work(*this);
    return work >>= b;
  }
  inline SimpleFloat<T,W,bits,U>& operator >>= (const U& b) {
    if(s & ((1 << INF) | (1 << NaN))) return *this;
    s |= safeAdd(e, - b);
    return ensureFlag();
  }
  inline SimpleFloat<T,W,bits,U>& operator =  (const SimpleFloat<T,W,bits,U>& src) {
    s = src.s;
    e = src.e;
    m = src.m;
    return *this;
  }
  inline SimpleFloat<T,W,bits,U>& operator =  (SimpleFloat<T,W,bits,U>&& src) {
    s = move(src.s);
    e = move(src.e);
    m = move(src.m);
    return *this;
  }
  inline bool             operator == (const SimpleFloat<T,W,bits,U>& src) const {
    return ! (*this != src);
  }
  inline bool             operator != (const SimpleFloat<T,W,bits,U>& src) const {
    return (((s | src.s) & ((1 << INF) | (1 << NaN))) ||
             (s != src.s || e != src.e || m != src.m)) &&
           ! (! m && ! src.m);
  }
  inline bool             operator <  (const SimpleFloat<T,W,bits,U>& src) const {
    if((s | src.s) & (1 << NaN)) return false;
    const unsigned char s_is_minus(s & (1 << SIGN));
    if(s_is_minus ^ (src.s & (1 << SIGN))) return s_is_minus;
    if(s & (1 << INF)) {
      if(src.s & (1 << INF)) return false;
      return s_is_minus;
    }
    if(src.s & (1 << INF)) return ! s_is_minus;
    if(m && src.m) {
      if(e < src.e) return ! s_is_minus;
      if(e == src.e) return s_is_minus ? src.m < m : m < src.m;
      return s_is_minus;
    }
    return !m ? (! src.m ? m != src.m : ! (src.s & (1 << SIGN))) : s_is_minus;
  }
  inline bool             operator <= (const SimpleFloat<T,W,bits,U>& src) const {
    return *this < src || *this == src;
  }
  inline bool             operator >  (const SimpleFloat<T,W,bits,U>& src) const {
    return ! (*this <= src);
  }
  inline bool             operator >= (const SimpleFloat<T,W,bits,U>& src) const {
    return ! (*this < src);
  }
  inline bool             operator !  () const {
    return ! m && isfinite(*this);
  }
  inline                  operator bool () const {
    return ! (!*this);
  }
  inline                  operator int  () const {
    return int(this->operator T());
  }
  inline                  operator T    () const {
    SimpleFloat<T,W,bits,U> deci(*this);
    if(deci.s & (1 << INF)) return T(int(0));
    if(deci.s & (1 << NaN)) return T(int(0));
    if(! deci.m) return T(int(0));
    if(U(bits) <= deci.e || (uzero() < deci.e && (deci.m << int(deci.e)) >> int(deci.e) != deci.m))
      return T(int(0));
    if(deci.e <= - U(bits)) return T(int(0));
    if(     deci.e < uzero()) deci.m >>= - int(deci.e);
    else if(uzero() < deci.e) deci.m <<=   int(deci.e);
    return s & (1 << SIGN) ? - deci.m : deci.m;
  }
  inline SimpleFloat<T,W,bits,U>  absfloor() const {
    if(uzero() <= e) return *this;
    if(e <= - U(bits)) return zero();
    SimpleFloat<T,W,bits,U> deci(*this);
    deci.m >>= - int(deci.e);
    deci.m <<= - int(deci.e);
    return deci;
  }
  inline SimpleFloat<T,W,bits,U>  absceil() const {
    const SimpleFloat<T,W,bits,U> fl(this->absfloor());
    if(*this - fl) {
      SimpleFloat<T,W,bits,U> pmone(one());
      pmone.s |= s & (1 << SIGN);
      return fl + pmone;
    }
    return fl;
  }
  inline SimpleFloat<T,W,bits,U>  abs()  const {
    SimpleFloat<T,W,bits,U> work(*this);
    work.s &= ~ (1 << SIGN);
    return work;
  }
         SimpleFloat<T,W,bits,U>  log()  const;
         SimpleFloat<T,W,bits,U>  exp()  const;
         SimpleFloat<T,W,bits,U>  sin()  const;
         SimpleFloat<T,W,bits,U>  cos()  const;
         SimpleFloat<T,W,bits,U>  atan() const;
  inline SimpleFloat<T,W,bits,U>  sqrt() const;
  
  uint64_t s;
  typedef enum {
    INF = 0,
    NaN = 1,
    SIGN = 2,
    DWRK = 3
  } state_t;
  T m;
  U e;
  const U uzero() const {
    return U(0);
  }
  const SimpleFloat<T,W,bits,U> zero()   const {
    return SimpleFloat<T,W,bits,U>(T(int(0)));
  }
  const SimpleFloat<T,W,bits,U> one()    const {
    return SimpleFloat<T,W,bits,U>(T(int(1)));
  }
  const SimpleFloat<T,W,bits,U> two()    const {
    return SimpleFloat<T,W,bits,U>(one() << U(1));
  }
  const SimpleFloat<T,W,bits,U> pi()     const {
    return SimpleFloat<T,W,bits,U>(quatpi() << U(2));
  }
  const SimpleFloat<T,W,bits,U> halfpi() const {
    return SimpleFloat<T,W,bits,U>(quatpi() << U(1));
  }
  const SimpleFloat<T,W,bits,U> quatpi() const;
  const SimpleFloat<T,W,bits,U> twopi()  const {
    return SimpleFloat<T,W,bits,U>(quatpi() << U(3));
  }
  const SimpleFloat<T,W,bits,U> sqrt2()  const {
    return SimpleFloat<T,W,bits,U>((one() << U(1)).sqrt());
  }
private:
  template <typename V> inline U normalize(V& src) const {
    V   bt(int(0));
    bt = ~ bt;
    int doff(sizeof(V) * 4);
    int b(0);
    for( ; bool(bt) && bool(doff); doff /= 2) {
      bt   = (~ V(int(0))) >> int(sizeof(V) * 8 - doff);
      bt <<= b + doff;
      if(bool(src & bt)) b += doff;
    }
    const U shift(sizeof(V) * 8 - b - 1);
    if(shift) src <<= shift;
    return - U(shift);
  }
  inline SimpleFloat<T,W,bits,U>& ensureFlag() {
    if(s & (1 << INF))
      s &= ~ (1 << DWRK);
    else if(! m || (s & (1 << DWRK))) {
      e ^= e;
      m ^= m;
      s &= ~ ((1 << DWRK) | (1 << INF));
    }
    return * this;
  }
  inline uint64_t safeAdd(U& dst, const U& src) {
    const U dst0(dst);
    dst += src;
    if((dst0 > uzero() && src > uzero() && dst <= uzero()) ||
       (dst0 < uzero() && src < uzero() && dst >= uzero()))
      return 1 << (dst0 < uzero() ? DWRK : INF);
    return 0;
  }
  inline int64_t residue2() const {
    if(uzero() < e || U(bits) <= - e) return 0;
    if(! e) return int64_t(int(m) & 1);
    return int64_t(int(m >> - int(e)) & 1);
  }

  // XXX: these are NOT threadsafe on first call.
  const vector<SimpleFloat<T,W,bits,U> >& exparray()    const;
  const vector<SimpleFloat<T,W,bits,U> >& invexparray() const;
};

SimpleFloat<uint32_t, uint64_t, _FLOAT_BITS_, int32_t >* sf_qpi;
template <typename T, typename W, int bits, typename U> const SimpleFloat<T,W,bits,U> SimpleFloat<T,W,bits,U>::quatpi() const {
  if(! sf_qpi) {
    sf_qpi = SimpleAllocator<SimpleFloat<T,W,bits,U> >().allocate(1);
    ::new ((void*)sf_qpi) SimpleFloat<T,W,bits,U>();
    * sf_qpi = one().atan();
  }
  return * sf_qpi;
}

template <typename T, typename W, int bits, typename U> SimpleFloat<T,W,bits,U> SimpleFloat<T,W,bits,U>::log() const {
  const SimpleFloat<T,W,bits,U> einv(one() / one().exp());
  const SimpleFloat<T,W,bits,U> one_einv(one() + einv);
  if((s & (1 << SIGN)) && m)
    return *this;
  if(s & ((1 << INF) | (1 << NaN)))
    return *this;
  if(! m) {
    SimpleFloat<T,W,bits,U> work(*this);
    work.s |= (1 << INF) | (1 << SIGN);
    return work;
  }
  if(einv <= *this && *this <= one_einv) {
    // ln(x) = (x - 1) - (x - 1)^2/2 + (x-1)^3/3- ...
    const SimpleFloat<T,W,bits,U> dx(*this - one());
          SimpleFloat<T,W,bits,U> x(dx);
          SimpleFloat<T,W,bits,U> before(one());
          SimpleFloat<T,W,bits,U> res(zero());
    for(int t = 1; (res - before).m; t ++, x *= dx) {
      const SimpleFloat<T,W,bits,U> abst(x / SimpleFloat<T,W,bits,U>(t));
      before = res;
      res   += (t % 2 ? abst : - abst);
    }
    return res;
  }
  const vector<SimpleFloat<T,W,bits,U> >& ea(exparray());
  const vector<SimpleFloat<T,W,bits,U> >& iea(invexparray());
  SimpleFloat<T,W,bits,U>  result(zero());
  SimpleFloat<T,W,bits,U>  work(*this);
  if(one_einv < work) {
    for(int i = min(ea.size(), iea.size()) - 1; 0 < i; i --)
      if(ea[i] <= work) {
        result += one() << U(i - 1);
        work   *= iea[i];
      }
    if(! (work <= one_einv)) {
      result += one();
      work   *= iea[1];
    }
  } else {
    for(int i = min(ea.size(), iea.size()) - 1; 0 < i; i --)
      if(work <= iea[i]) {
        result -= one() << U(i - 1);
        work   *= ea[i];
      }
    if(! (einv <= work)) {
      result -= one();
      work   *= ea[1];
    }
  }
  return result += work.log();
}

template <typename T, typename W, int bits, typename U> SimpleFloat<T,W,bits,U> SimpleFloat<T,W,bits,U>::exp() const {
  if(s & ((1 << INF) | (1 << NaN))) {
    if(! (s & (1 << NaN)) && (s & (1 << SIGN)))
      return zero();
    return *this;
  }
  if(this->abs() <= one()) {
    // exp(x) = 1 + x/1! + x^2/2! + ...
    SimpleFloat<T,W,bits,U> denom(one());
    SimpleFloat<T,W,bits,U> x(*this);
    SimpleFloat<T,W,bits,U> before(zero());
    SimpleFloat<T,W,bits,U> res(one());
    for(int t = 1; (res - before).m; t ++, x *= *this) {
      before = res;
      denom *= SimpleFloat<T,W,bits,U>(t);
      res   += x / denom;
    }
    return res;
  }
  const vector<SimpleFloat<T,W,bits,U> >& en(exparray());
  const vector<SimpleFloat<T,W,bits,U> >& ien(invexparray());
  SimpleFloat<T,W,bits,U> work(this->abs());
  SimpleFloat<T,W,bits,U> result(one());
  for(int i = 1; 0 <= i && i < min(en.size(), ien.size()) && work.absfloor(); i ++, work >>= U(1))
    if(work.residue2())
      result *= s & (1 << SIGN) ? ien[i] : en[i];
  if(work.absfloor()) {
    work.s |= 1 << INF;
    return work;
  }
  const SimpleFloat<T,W,bits,U> residue(*this - this->absfloor());
  return result *= residue.exp();
}

template <typename T, typename W, int bits, typename U> SimpleFloat<T,W,bits,U> SimpleFloat<T,W,bits,U>::sin() const {
  if(s & ((1 << INF) | (1 << NaN))) {
    SimpleFloat<T,W,bits,U> res(*this);
    res.s |= 1 << NaN;
    return res;
  }
  if(- one() <= *this && *this <= one()) {
    // sin(x) = x - x^3/3! + x^5/5! - ...
    const SimpleFloat<T,W,bits,U> sqx(*this * *this);
          SimpleFloat<T,W,bits,U> denom(one());
          SimpleFloat<T,W,bits,U> x(sqx * *this);
          SimpleFloat<T,W,bits,U> before(zero());
          SimpleFloat<T,W,bits,U> res(*this);
    for(int t = 1; (res - before).m; t ++, x *= sqx) {
      SimpleFloat<T,W,bits,U> tt(t);
      tt   <<= U(1);
      before = res;
      denom *= - tt * (tt + one());
      res   += x / denom;
    }
    return res;
  }
  if(- halfpi() <= *this && *this <= halfpi())
    return ((*this - quatpi()).cos() + (*this - quatpi()).sin()) / sqrt2();
  if(this->abs() == pi())
    return zero();
  if(- pi() <= *this && *this <= pi())
    return (halfpi() - *this).cos();
  if(- twopi() <= *this && *this <= twopi())
    return - (*this + pi()).sin();
  return (*this % twopi()).sin();
}

template <typename T, typename W, int bits, typename U> SimpleFloat<T,W,bits,U> SimpleFloat<T,W,bits,U>::cos() const {
  if(s & ((1 << INF) | (1 << NaN))) {
    SimpleFloat<T,W,bits,U> res(*this);
    res.s |= 1 << NaN;
    return res;
  }
  if(- one() <= *this && *this <= one()) {
    // cos(x) = 1 - x^2/2! + x^4/4! - ...
    const SimpleFloat<T,W,bits,U> sqx(*this * *this);
          SimpleFloat<T,W,bits,U> denom(one());
          SimpleFloat<T,W,bits,U> x(sqx);
          SimpleFloat<T,W,bits,U> before(zero());
          SimpleFloat<T,W,bits,U> res(one());
    for(int t = 1; (res - before).m; t ++, x *= sqx) {
      SimpleFloat<T,W,bits,U> tt(t);
      tt   <<= U(1);
      before = res;
      denom *= - tt * (tt - one());
      res   += x / denom;
    }
    return res;
  }
  if(this->abs() == halfpi())
    return zero();
  if(- halfpi() <= *this && *this <= halfpi())
    return ((*this - quatpi()).cos() - (*this - quatpi()).sin()) / sqrt2();
  if(this->abs() == pi())
    return - one();
  if(- pi() <= *this && *this <= pi())
    return (halfpi() - *this).sin();
  if(- twopi() <= *this && *this <= twopi())
    return - (*this + pi()).cos();
  return (*this % twopi()).cos();
}

template <typename T, typename W, int bits, typename U> SimpleFloat<T,W,bits,U> SimpleFloat<T,W,bits,U>::atan() const {
  if(s & ((1 << INF) | (1 << NaN))) {
    if(! (s & (1 << NaN)))
      return s & (1 << SIGN) ? - halfpi() : halfpi();
    return *this;
  }
  if(s & (1 << SIGN))
    return - (- *this).atan();
  const SimpleFloat<T,W,bits,U> half(one() >> U(1));
  const SimpleFloat<T,W,bits,U> four(one() << U(2));
  const SimpleFloat<T,W,bits,U> five((one() << U(2)) + one());
  if(- half <= *this && *this <= half) {
    // arctan(x) = x - x^3/3 + x^5/5 - ...
    const SimpleFloat<T,W,bits,U> sqx(*this * *this);
          SimpleFloat<T,W,bits,U> x(sqx * *this);
          SimpleFloat<T,W,bits,U> before(zero());
          SimpleFloat<T,W,bits,U> res(*this);
    for(int t = 1; (res - before).m; t ++, x *= sqx) {
      const SimpleFloat<T,W,bits,U> abst(x / ((SimpleFloat<T,W,bits,U>(t) << U(1)) + one()));
      before = res;
      res   += (t % 2 ? - abst : abst);
    }
    return res;
  }
  // N.B.
  //  atan(u) + atan(v) = atan((u + v) / (1 - uv)) mod pi, uv != 1.
  //    in u = 0.5, v = x - 0.5 case,
  //  atan(x / (1 - x / 2 + 1 / 4)) = atan(.5) + atan(x - .5) =
  //  atan(x / (1.25 - .5 * x)) 
  //  y := x / (1.25 - .5 * x) then,
  //  (1.25 - .5 * x) * y = x,
  //  (5 - 2x) * y = 4 x
  //  x = 5y / (4 + 2y),
  //     y - x = ((4 + 2y) * y - 5y) / (4 + 2y)
  //           = y * (2y - 1) / (4 + 2y)
  //     so 0 <= y and 0 < y case, this makes decreasing function.
  //       (v = x - .5 and 0 <= 2y - 1)
  if(- two() <= *this && *this <= two()) {
    const SimpleFloat<T,W,bits,U> atanhalf(half.atan());
    const SimpleFloat<T,W,bits,U>  v(five * *this / (four + (*this << U(1))) - half);
    return atanhalf + v.atan();
  }
  // N.B.
  //    in u = v case,
  //  2 atan(u) = atan(2 * u / (1 - u * u))
  //  2 atan(u) = atan(2 * u / (1 - u) / (1 + u))
  //            = atan(2 / (1 / u - u))
  //    in 2Y := 1 / u - u case,
  //            = atan(1 / Y),
  //  u^2 + 2Yu - 1 == 0, u = - Y \pm sqrt(Y^2 + 1)
  const SimpleFloat<T,W,bits,U> Y(one() / (*this));
  const SimpleFloat<T,W,bits,U> u((Y * Y + one()).sqrt() - Y);
  return u.atan() << U(1);
}

// XXX:
vector<SimpleFloat<uint32_t, uint64_t, _FLOAT_BITS_, int32_t > >* sf_ea;
vector<SimpleFloat<uint32_t, uint64_t, _FLOAT_BITS_, int32_t > >* sf_iea;
template <typename T, typename W, int bits, typename U> const vector<SimpleFloat<T,W,bits,U> >& SimpleFloat<T,W,bits,U>::exparray() const {
  if(! sf_ea) {
    sf_ea = SimpleAllocator<vector<SimpleFloat<T,W,bits,U> > >().allocate(1);
    ::new ((void*)sf_ea) vector<SimpleFloat<T,W,bits,U> >();
  }
  vector<SimpleFloat<T,W,bits,U> >& ebuf(* sf_ea);
  if(ebuf.size())
    return ebuf;
  ebuf.emplace_back(one());
  ebuf.emplace_back(ebuf[0].exp());
  for(int i = 1; 0 <= i; i ++) {
    const SimpleFloat<T,W,bits,U> en(ebuf[i] * ebuf[i]);
    if(en && isfinite(en))
      ebuf.emplace_back(en);
    else
      break;
  }
  return ebuf;
}

template <typename T, typename W, int bits, typename U> const vector<SimpleFloat<T,W,bits,U> >& SimpleFloat<T,W,bits,U>::invexparray() const {
  if(! sf_iea) {
    sf_iea = SimpleAllocator<vector<SimpleFloat<T,W,bits,U> > >().allocate(1);
    ::new ((void*)sf_iea) vector<SimpleFloat<T,W,bits,U> >();
  }
  vector<SimpleFloat<T,W,bits,U> >& iebuf(* sf_iea);
  if(iebuf.size())
    return iebuf;
  const vector<SimpleFloat<T,W,bits,U> >& ea(exparray());
  for(int i = 0; 0 <= i && i < ea.size(); i ++) {
    const SimpleFloat<T,W,bits,U> ien(one() / ea[i]);
    if(ien && isfinite(ien))
      iebuf.emplace_back(ien);
    else
      break;
  }
  return iebuf;
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U> SimpleFloat<T,W,bits,U>::sqrt() const {
  if(s & ((1 << INF) | (1 << NaN))) {
    SimpleFloat<T,W,bits,U> res(*this);
    if(s & (1 << SIGN)) res.s |= 1 << NaN;
    return res;
  }
  SimpleFloat<T,W,bits,U> res((this->log() >> U(1)).exp());
  // get better accuracy (is this enough?, double accuracy on one loop.)
  // newton's method: 0 == f'(x_n) (x_{n+1} - x_n) + f(x_n)
  //            x_{n+1} := x_n - f(x_n)/f'(x_n).
  //         where f(x) := x_n * x_n - *this
  if(! res)
    return res;
  return (res + *this / res) >> U(1);
}

// N.B. only to enname, better languages can omit these.
template <typename T, typename W, int bits, typename U> static inline bool isinf(const SimpleFloat<T,W,bits,U>& src) {
  return src.s & (1 << src.INF);
}

template <typename T, typename W, int bits, typename U> static inline bool isnan(const SimpleFloat<T,W,bits,U>& src) {
  return src.s & (1 << src.NaN);
}

template <typename T, typename W, int bits, typename U> static inline bool isfinite(const SimpleFloat<T,W,bits,U>& src) {
  return ! (src.s & ((1 << src.INF) | (1 << src.NaN)));
}

template <typename T, typename W, int bits, typename U> static inline SimpleFloat<T,W,bits,U> absfloor(const SimpleFloat<T,W,bits,U>& src) {
  return src.absfloor();
}

template <typename T, typename W, int bits, typename U> static inline SimpleFloat<T,W,bits,U> absceil(const SimpleFloat<T,W,bits,U>& src) {
  return src.absceil();
}

template <typename T, typename W, int bits, typename U> static inline SimpleFloat<T,W,bits,U> abs(const SimpleFloat<T,W,bits,U>& src) {
  return src.abs();
}

template <typename T, typename W, int bits, typename U> static inline SimpleFloat<T,W,bits,U> sqrt(const SimpleFloat<T,W,bits,U>& src) {
  return src.sqrt();
}

template <typename T, typename W, int bits, typename U> static inline SimpleFloat<T,W,bits,U> exp(const SimpleFloat<T,W,bits,U>& src) {
  return src.exp();
}

template <typename T, typename W, int bits, typename U> static inline SimpleFloat<T,W,bits,U> log(const SimpleFloat<T,W,bits,U>& src) {
  return src.log();
}

template <typename T, typename W, int bits, typename U> static inline SimpleFloat<T,W,bits,U> sin(const SimpleFloat<T,W,bits,U>& src) {
  return src.sin();
}

template <typename T, typename W, int bits, typename U> static inline SimpleFloat<T,W,bits,U> cos(const SimpleFloat<T,W,bits,U>& src) {
  return src.cos();
}

template <typename T, typename W, int bits, typename U> static inline SimpleFloat<T,W,bits,U> tan(const SimpleFloat<T,W,bits,U>& src) {
  return src.sin() / src.cos();
}

template <typename T, typename W, int bits, typename U> static inline SimpleFloat<T,W,bits,U> atan(const SimpleFloat<T,W,bits,U>& src) {
  return src.atan();
}

template <typename T, typename W, int bits, typename U> static inline SimpleFloat<T,W,bits,U> atan2(const SimpleFloat<T,W,bits,U>& y, const SimpleFloat<T,W,bits,U>& x) {
  SimpleFloat<T,W,bits,U> atan0(y.halfpi());
  if(! x && ! y)
    return x / y;
  else if(isfinite(x)) {
    if(! isfinite(y) )
      goto ensure;
    if(! x)
      goto ensure;
    const SimpleFloat<T,W,bits,U> yoverx((y / x).abs());
    if(! isfinite(yoverx) )
      goto ensure;
    const SimpleFloat<T,W,bits,U> atan00(yoverx.atan());
    if(! isfinite(atan00) )
      goto ensure;
    atan0 = atan00;
    goto ensure;
  } else if(isfinite(y)) {
    atan0 = x.zero();
    goto ensure;
  }
  return y;
 ensure:
  if(y.s & (1 << y.SIGN)) {
    if(x.s & (1 << x.SIGN))
      atan0 = - (x.pi() - atan0);
    else
      atan0 = - atan0;
  } else if(x.s & (1 << x.SIGN))
    atan0  = x.pi() - atan0;
  return atan0;
}

template <typename T, typename W, int bits, typename U> static inline SimpleFloat<T,W,bits,U> pow(const SimpleFloat<T,W,bits,U>& src, const SimpleFloat<T,W,bits,U>& dst) {
  if(! dst) {
    if(! src)
      return T(int(0)) / T(int(0));
    return dst.one();
  }
  return exp(log(src) * dst);
}

// N.B. start class complex part:
template <typename T> class Complex {
public:
  inline Complex() { ; }
  inline Complex(const T& real, const T& imag = T(int(0))) {
    _real = real; _imag = imag;
  }
  inline Complex(const Complex<T>& s) { *this = s; }
  inline Complex(Complex<T>&& s) { *this = s; }
  inline Complex(T&& real) {
    const T zero(0);
    _real = move(real);
    _imag = zero;
  }
  inline Complex(T&& real, T&& imag) {
    _real = move(real);
    _imag = move(imag);
  }
  inline ~Complex() { ; }
  inline Complex<T>  operator ~  ()                    const {
    return Complex<T>(  _real, - _imag);
  }
  inline Complex<T>  operator -  ()                    const {
    return Complex<T>(- _real, - _imag);
  }
  inline Complex<T>  operator +  (const Complex<T>& s) const {
    Complex<T> result(*this);
    return result += s;
  }
  inline Complex<T>& operator += (const Complex<T>& s) {
    _real += s._real;
    _imag += s._imag;
    return *this;
  }
  inline Complex<T>  operator -  (const Complex<T>& s) const {
    Complex<T> result(*this);
    return result -= s;
  }
  inline Complex<T>& operator -= (const Complex<T>& s) {
    _real -= s._real;
    _imag -= s._imag;
    return *this;
  }
  inline Complex<T>  operator *  (const T& s)          const {
    Complex<T> result(*this);
    return result *= s;
  }
  inline Complex<T>& operator *= (const T& s) {
    _real *= s;
    _imag *= s;
    return *this;
  }
  inline Complex<T>  operator *  (const Complex<T>& s) const {
    return Complex<T>(_real * s._real - _imag * s._imag,
                      _real * s._imag + _imag * s._real);
  }
  inline Complex<T>& operator *= (const Complex<T>& s) {
    return (*this) = (*this) * s;
  }
  inline Complex<T>  operator /  (const T& s)          const {
    Complex<T> result(*this);
    return result /= s;
  }
  inline Complex<T>& operator /= (const T& s) {
    _real /= s;
    _imag /= s;
    return *this;
  }
  inline Complex<T>  operator /  (const Complex<T>& s) const {
    return (*this * (~ s)) / (s._real * s._real + s._imag * s._imag);
  }
  inline Complex<T>& operator /= (const Complex<T>& s) {
    return *this = *this / s;
  }
  inline bool        operator == (const Complex<T>& s) const {
    return !(*this != s);
  }
  inline bool        operator != (const Complex<T>& s) const {
    return (_real != s._real) || (_imag != s._imag);
  }
  inline bool        operator !  ()                    const {
    return !_real && !_imag;
  }
  inline Complex<T>  operator &  (const Complex<T>& s) const {
    Complex<T> result(*this);
    return result &= s;
  }
  inline Complex<T>& operator &= (const Complex<T>& s) {
    _real &= s._real;
    _imag &= s._imag;
    return *this;
  }
  inline Complex<T>  operator |  (const Complex<T>& s) const {
    Complex<T> result(*this);
    return result |= s;
  }
  inline Complex<T>& operator |= (const Complex<T>& s) {
    _real |= s._real;
    _imag |= s._imag;
    return *this;
  }
  inline Complex<T>  operator ^  (const Complex<T>& s) const {
    Complex<T> result(*this);
    return result ^= s;
  }
  inline Complex<T>& operator ^= (const Complex<T>& s) {
    _real ^= s._real;
    _imag ^= s._imag;
    return *this;
  }
  inline bool        operator && (const Complex<T>& s) const {
    return *this && s;
  }
  inline bool        operator || (const Complex<T>& s) const {
    return *this || s;
  }
  inline Complex<T>& operator =  (const Complex<T>& s) {
    _real = s._real;
    _imag = s._imag;
    return *this;
  }
  inline Complex<T>& operator =  (Complex<T>&& s) {
    _real = move(s._real);
    _imag = move(s._imag);
    return *this;
  }
  inline T&          operator [] (const size_t& i) {
    if(i) return _imag;
    return _real;
  }
  inline             operator bool () const {
    return ! (! *this);
  }
  inline             operator T    () const {
    return this->_real;
  }
  const Complex<T> i() const {
    return Complex<T>(Complex<T>(T(int(0)), T(int(1))));
  }
  inline T  abs() const {
    return sqrt(_real * _real + _imag * _imag);
  }
  inline T  arg() const {
    return atan2(_imag, _real);
  }
  inline T& real() {
    return _real;
  }
  inline T& imag() {
    return _imag;
  }
  inline const T& real() const {
    return _real;
  }
  inline const T& imag() const {
    return _imag;
  }
  T _real;
  T _imag;
};

template <typename T> static inline T abs(const Complex<T>& s) {
  return s.abs();
}

template <typename T> static inline T arg(const Complex<T>& s) {
  return s.arg();
}

template <typename T> static inline const T& real(const Complex<T>& s) {
  return s.real();
}

template <typename T> static inline const T& imag(const Complex<T>& s) {
  return s.imag();
}

template <typename T> static inline Complex<T> exp(const Complex<T>& s) {
  return Complex<T>(exp(s.real())) * Complex<T>(cos(s.imag()), sin(s.imag()));
}

template <typename T> static inline Complex<T> log(const Complex<T>& s) {
  // N.B. main branch
  return Complex<T>(log(abs(s)), arg(s));
}

template <typename T> static inline Complex<T> pow(const Complex<T>& s, const Complex<T>& p) {
  if(abs(s) == T(int(0))) {
    return T(int(0));
  }
  return exp(log(s) * p);
}

template <typename T> static inline Complex<T> sqrt(const Complex<T>& s) {
  return exp(log(s) * Complex<T>(T(int(1)) / T(int(2))));
}

template <typename T> static inline Complex<T> csin(const Complex<T>& s) {
  return (exp(Complex<T>(T(int(0)), s)) - exp(Complex<T>(T(int(0)), - s))) / Complex<T>(T(int(0)), T(int(2)));
}

template <typename T> static inline Complex<T> ccos(const Complex<T>& s) {
  return (exp(Complex<T>(T(int(0)), s)) + exp(Complex<T>(T(int(0)), - s))) / T(int(2));
}

template <typename T> static inline Complex<T> ctan(const Complex<T>& s) {
  return csin(s) / ccos(s);
}

template <typename T> static inline Complex<T> ccsc(const Complex<T>& s) {
  return Complex<T>(T(int(1))) / csin(s);
}

template <typename T> static inline Complex<T> csec(const Complex<T>& s) {
  return Complex<T>(T(int(1))) / ccos(s);
}

template <typename T> static inline T ccot(const T& s) {
  return Complex<T>(T(int(1))) / ctan(s);
}

#if _FLOAT_BITS_ == 32
    typedef uint32_t myuint;
    typedef int32_t  myint;
    typedef SimpleFloat<myuint, uint64_t, 32, myint> myfloat;
//    typedef SimpleFloat<myuint, DUInt<uint32_t, 32>, 32, myint> myfloat;
#elif _FLOAT_BITS_ == 64
    typedef uint64_t myuint;
    typedef int64_t  myint;
    typedef SimpleFloat<myuint, unsigned __int128, 64, myint> myfloat;
#elif _FLOAT_BITS_ == 128
    typedef DUInt<uint64_t, 64> uint128_t;
    typedef Signed<uint128_t, 128> int128_t;
    typedef uint128_t myuint;
    typedef int128_t  myint;
    typedef SimpleFloat<myuint, DUInt<myuint, 128>, 128, myint> myfloat;
#else
#  error cannot handle float
#endif

template <typename T> using complexC = Complex<T>;
#define complex(T) complexC<T>
#define complexctor(T) complex(T)

// N.B. start simplelin.
#if defined(_SIMPLEALLOC_)
#if defined(_OPENMP)
#error SimpleAllocator not supported
#endif
unsigned long long last;
unsigned long long lastptr;
unsigned long long sam_upper;
unsigned long long *v_alloc;
int *in_use;
// N.B. around 20% endurance, so making here into binary tree can reduce some.
template <typename T> class SimpleAllocator {
public:
  typedef T* pointer;
  typedef const T* const_pointer;
  typedef T value_type;
  inline SimpleAllocator() { }
  template <typename U> inline SimpleAllocator(const SimpleAllocator<U>&) { }
  inline ~SimpleAllocator() { }
  inline T* allocate(size_t n) {
    n *= sizeof(T);
    n  = (n + _SIMPLEALLOC_ - 1) / _SIMPLEALLOC_ * _SIMPLEALLOC_;
    if(M_ALLOC <= lastptr) { printf("pool full.\n"); for(;;) ; }
    v_alloc[lastptr] = n;
    in_use[lastptr ++] = 1;
    last += n;
    if(! (last < sam_upper)) { printf("memory full\n"); for(;;) ;}
    return reinterpret_cast<T*>(last - n);
  }
  inline void deallocate(T* p, size_t n) {
    size_t pp(reinterpret_cast<size_t>(p));
    size_t work(last);
    bool   flag(false);
    int    i;
    for(i = lastptr - 1; 0 <= i; i --) {
      if((work -= v_alloc[i]) == pp) break;
      flag = flag || in_use[i];
    }
    if(work != pp) { printf("free bug\n"); for(;;) ; }
    in_use[i] ^= in_use[i];
    if(! flag) {
      last = pp;
      lastptr = i;
    }
  }
  inline void destroy(T* p) { p->~T(); }
};
#endif

template <typename T> class SimpleVector {
public:
  inline SimpleVector() { ; }
  inline SimpleVector(const int& size) {
    this->entity.resize(size);
  }
  inline SimpleVector(const SimpleVector<T>& other) { *this = other; }
  inline SimpleVector(SimpleVector<T>&& other) { *this = other; }
  inline ~SimpleVector() { ; }
  inline       SimpleVector<T>  operator -  () const {
    SimpleVector<T> res(entity.size());
#if defined(_OPENMP)
#pragma omp simd
#endif
    for(int i = 0; i < entity.size(); i ++)
      res.entity[i] = - entity[i];
    return res;
  }
  inline       SimpleVector<T>  operator +  (const SimpleVector<T>& other) const {
    SimpleVector<T> res(*this);
    return res += other;
  }
  inline       SimpleVector<T>& operator += (const SimpleVector<T>& other) {
#if defined(_OPENMP)
#pragma omp simd
#endif
    for(int i = 0; i < entity.size(); i ++)
      entity[i] += other.entity[i];
    return *this;
  }
  inline       SimpleVector<T>  operator -  (const SimpleVector<T>& other) const {
    SimpleVector<T> res(*this);
    return res -= other;
  }
  inline       SimpleVector<T>& operator -= (const SimpleVector<T>& other) {
    return *this += - other;
  }
  inline       SimpleVector<T>  operator *  (const T& other) const {
    SimpleVector<T> res(*this);
    return res *= other;
  }
  inline       SimpleVector<T>& operator *= (const T& other) {
#if defined(_OPENMP)
#pragma omp simd
#endif
    for(int i = 0; i < entity.size(); i ++)
      entity[i] *= other;
    return *this;
  }
  inline       SimpleVector<T>  operator /  (const T& other) const {
    SimpleVector<T> res(*this);
    return res /= other;
  }
  inline       SimpleVector<T>& operator /= (const T& other) {
#if defined(_OPENMP)
#pragma omp simd
#endif
    for(int i = 0; i < entity.size(); i ++)
      entity[i] /= other;
    return *this;
  }
  inline       SimpleVector<T>& operator =  (const SimpleVector<T>& other) {
    entity = other.entity;
    return *this;
  }
  inline       SimpleVector<T>& operator =  (SimpleVector<T>&& other) {
    entity = move(other.entity);
    return *this;
  }
  inline       bool             operator == (const SimpleVector<T>& other) const {
    return ! (*this != other);
  }
  inline       bool             operator != (const SimpleVector<T>& other) const {
    if(entity.size() != other.entity.size()) return true;
    for(int i = 0; i < entity.size(); i ++)
      if(entity[i] != other.entity[i]) return true;
    return false;
  }
  template <typename U> inline T dot(const SimpleVector<U>& other) const {
    SimpleVector<T> work(other.size());
#if defined(_OPENMP)
#pragma omp simd
#endif
    for(int i = 0; i < entity.size(); i ++)
      work[i] = entity[i] * other.entity[i];
    T res(work[0]);
    for(int i = 1; i < entity.size(); i ++)
      res += work[i];
    return res;
  }
  inline       T&               operator [] (const int& idx) {
    return entity[idx];
  }
  inline const T&               operator [] (const int& idx) const {
    return entity[idx];
  }
  template <typename U> inline SimpleVector<U> real() const {
    SimpleVector<U> result(entity.size());
#if defined(_OPENMP)
#pragma omp simd
#endif
    for(int i = 0; i < entity.size(); i ++)
      result.entity[i] = U(entity[i].real());
    return result;
  }
  template <typename U> inline SimpleVector<U> imag() const {
    SimpleVector<U> result(entity.size());
#if defined(_OPENMP)
#pragma omp simd
#endif
    for(int i = 0; i < entity.size(); i ++)
      result.entity[i] = U(entity[i].imag());
    return result;
  }
  template <typename U> inline SimpleVector<U> cast() const {
    SimpleVector<U> result(entity.size());
#if defined(_OPENMP)
#pragma omp simd
#endif
    for(int i = 0; i < entity.size(); i ++)
      result.entity[i] = U(entity[i]);
    return result;
  }
  inline const int size() const {
    return entity.size();
  }
  inline       void resize(const int& size) {
    entity.resize(size);
    return;
  }
  inline       SimpleVector<T>  subVector(const int& i, const int& s) const {
    SimpleVector<T> res(s);
#if defined(_OPENMP)
#pragma omp simd
#endif
    for(int ii = i; ii < i + s; ii ++)
      res[ii - i] = (*this)[ii];
    return res;
  }
  inline       SimpleVector<T>& setVector(const int& i, const SimpleVector<T>& d) {
#if defined(_OPENMP)
#pragma omp simd
#endif
    for(int ii = i; ii < i + d.size(); ii ++)
      (*this)[ii] = d[ii - i];
    return *this;
  }
  inline       SimpleVector<T>& O(const T& r = T(int(0))) {
    return I(r);
  }
  inline       SimpleVector<T>& I(const T& r = T(int(1))) {
#if defined(_OPENMP)
#pragma omp simd
#endif
    for(int i = 0; i < size(); i ++)
      (*this)[i] = r;
    return *this;
  }
  inline       SimpleVector<T>& ek(const int& i, const T& r = T(int(1))) {
    const T zero(0);
#if defined(_OPENMP)
#pragma omp simd
#endif
    for(int ii = 0; ii < size(); ii ++)
      (*this)[ii] = ii == i ? r : zero;
    return *this;
  }
  inline       SimpleVector<T>  reverse() {
    SimpleVector<T> res(entity.size());
#if defined(_OPENMP)
#pragma omp simd
#endif
    for(int i = 0; i < res.size(); i ++)
      res[i] = entity[entity.size() - 1 - i];
    return res;
  }
#if defined(_SIMPLEALLOC_)
  vector<T, SimpleAllocator<T> > entity;
#else
  vector<T> entity;
#endif
};

template <typename T> class SimpleMatrix {
public:
  inline SimpleMatrix() { ecols = 0; }
  inline SimpleMatrix(const int& rows, const int& cols) {
    entity.resize(rows);
    for(int i = 0; i < entity.size(); i ++)
      entity[i].resize(cols);
    ecols = cols;
  }
  inline SimpleMatrix(const SimpleMatrix<T>& other) { *this = other; }
  inline SimpleMatrix(SimpleMatrix<T>&& other) { *this = other; }
  inline ~SimpleMatrix() { ; }
  inline       SimpleMatrix<T>  operator -  () const {
    SimpleMatrix<T> res(entity.size(), ecols);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < entity.size(); i ++)
      res.entity[i] = - entity[i];
    return res;
  }
  inline       SimpleMatrix<T>  operator +  (const SimpleMatrix<T>& other) const {
    SimpleMatrix<T> res(*this);
    return res += other;
  }
  inline       SimpleMatrix<T>& operator += (const SimpleMatrix<T>& other) {
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < entity.size(); i ++)
      entity[i] += other.entity[i];
    return *this;
  }
  inline       SimpleMatrix<T>  operator -  (const SimpleMatrix<T>& other) const {
    SimpleMatrix<T> res(*this);
    return res -= other;
  }
  inline       SimpleMatrix<T>& operator -= (const SimpleMatrix<T>& other) {
    return *this += - other;
  }
  inline       SimpleMatrix<T>  operator *  (const T& other) const {
    SimpleMatrix<T> res(*this);
    return res *= other;
  }
  inline       SimpleMatrix<T>& operator *= (const T& other) {
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < entity.size(); i ++)
      entity[i] *= other;
    return *this;
  }
  inline       SimpleMatrix<T>  operator *  (const SimpleMatrix<T>& other) const {
    SimpleMatrix<T> derived(other.transpose());
    SimpleMatrix<T> res(entity.size(), other.ecols);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < entity.size(); i ++) {
            SimpleVector<T>& resi(res.entity[i]);
      const SimpleVector<T>& ei(entity[i]);
      for(int j = 0; j < other.ecols; j ++)
        resi[j] = ei.dot(derived.entity[j]);
    }
    return res;
  }
  inline       SimpleMatrix<T>& operator *= (const SimpleMatrix<T>& other) {
    return *this = *this * other;
  }
  inline       SimpleVector<T>  operator *  (const SimpleVector<T>& other) const {
    SimpleVector<T> res(entity.size());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < entity.size(); i ++)
      res[i] = entity[i].dot(other);
    return res;
  }
  inline       SimpleMatrix<T>  operator /  (const T& other) const {
    SimpleMatrix<T> res(*this);
    return res /= other;
  }
  inline       SimpleMatrix<T>& operator /= (const T& other) {
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < entity.size(); i ++)
      entity[i] /= other;
    return *this;
  }
  inline       SimpleMatrix<T>& operator =  (const SimpleMatrix<T>& other) {
    ecols  = other.ecols;
    entity = other.entity;
    return *this;
  }
  inline       SimpleMatrix<T>& operator =  (SimpleMatrix<T>&& other) {
    ecols  = move(other.ecols);
    entity = move(other.entity);
    return *this;
  }
  inline       bool             operator == (const SimpleMatrix<T>& other) const {
    return ! (*this != other);
  }
  inline       bool             operator != (const SimpleMatrix<T>& other) const {
    if(entity.size() != other.entity.size() || ecols != other.ecols)
      return true;
    for(int i = 0; i < entity.size(); i ++)
      if(entity[i] != other.entity[i]) return true;
    return false;
  }
  inline       T&               operator () (const int& y, const int& x) {
    return entity[y][x];
  }
  inline const T&               operator () (const int& y, const int& x) const {
    return entity[y][x];
  }
  inline       SimpleVector<T>& row(const int& y) {
    return entity[y];
  }
  inline const SimpleVector<T>& row(const int& y) const {
    return entity[y];
  }
  inline const SimpleVector<T>  col(const int& x) const {
    SimpleVector<T> res(entity.size());
    for(int i = 0; i < entity.size(); i ++)
      res[i] = entity[i][x];
    return res;
  }
  inline       void             setCol(const int& x, const SimpleVector<T>& other) {
    for(int i = 0; i < entity.size(); i ++)
      entity[i][x] = other[i];
    return;
  }
  // N.B. transpose : exhaust of the resource, so Eigen library handles better.
  inline       SimpleMatrix<T>  transpose() const {
    SimpleMatrix<T> res(ecols, entity.size());
    for(int i = 0; i < ecols; i ++) {
      SimpleVector<T>& resi(res.entity[i]);
      for(int j = 0; j < entity.size(); j ++)
        resi[j] = entity[j][i];
    }
    return res;
  }
  inline       SimpleMatrix<T>  subMatrix(const int& y, const int& x, const int& h, const int& w) const {
    SimpleMatrix<T> res(h, w);
    for(int i = y; i < y + h; i ++)
      for(int j = x; j < x + w; j ++)
        res(i - y, j - x) = (*this)(i, j);
    return res;
  }
  inline       SimpleMatrix<T>& setMatrix(const int& y, const int& x, const SimpleMatrix<T>& d) {
    for(int i = y; i < y + d.rows(); i ++)
      for(int j = x; j < x + d.cols(); j ++)
        (*this)(i, j) = d(i - y, j - x);
    return *this;
  }
  inline       SimpleMatrix<T>& O(const T& r = T(int(0))) {
    for(int i = 0; i < rows(); i ++)
      for(int j = 0; j < cols(); j ++)
        (*this)(i, j) = r;
    return *this;
  }
  inline       SimpleMatrix<T>& I(const T& r = T(int(1))) {
    const T zero(0);
    for(int i = 0; i < rows(); i ++)
      for(int j = 0; j < cols(); j ++)
        (*this)(i, j) = (i == j ? r : zero);
    return *this;
  }
  inline       T                determinant(const bool& nonzero = false) const;
  inline       SimpleMatrix<T>  inverse() const {
    // XXX: extremely slow implementation.
    SimpleMatrix<T> result(entity.size(), ecols);
    result.I();
    for(int i = 0; i < result.cols(); i ++)
      result.setCol(i, solve(result.col(i)));
    return result;
  }
  inline       SimpleVector<T>  solve(SimpleVector<T> other) const;
  inline       SimpleVector<T>  solveN(SimpleVector<T> other) const;
  inline       SimpleVector<T>  projectionPt(const SimpleVector<T>& other) const;
  inline       SimpleMatrix<T>& fillP(const vector<int>& idx);
  inline       SimpleMatrix<T>  QR() const;
  inline       SimpleMatrix<T>  SVDleft1d() const;
  inline       pair<SimpleMatrix<T>, SimpleMatrix<T> > SVD1d() const;
  inline       SimpleMatrix<T> SVD() const;
  inline       pair<pair<SimpleMatrix<T>, SimpleMatrix<T> >, SimpleMatrix<T> > SVD(const SimpleMatrix<T>& src) const;
  inline       SimpleVector<T>  zeroFix(const SimpleMatrix<T>& A, vector<pair<T, int> > fidx);
  inline       SimpleVector<T>  inner(const SimpleVector<T>& bl, const SimpleVector<T>& bu) const;
  template <typename U> inline SimpleMatrix<U> real() const {
    SimpleMatrix<U> res(entity.size(), ecols);
    for(int i = 0; i < entity.size(); i ++)
      for(int j = 0; j < ecols; j ++)
        res(i, j) = U(entity[i][j].real());
    return res;
  }
  template <typename U> inline SimpleMatrix<U> imag() const {
    SimpleMatrix<U> res(entity.size(), ecols);
    for(int i = 0; i < entity.size(); i ++)
      for(int j = 0; j < ecols; j ++)
        res(i, j) = U(entity[i][j].imag());
    return res;
  }
  template <typename U> inline SimpleMatrix<U> cast() const {
    SimpleMatrix<U> res(entity.size(), ecols);
    for(int i = 0; i < entity.size(); i ++)
      for(int j = 0; j < ecols; j ++)
        res(i, j) = U(entity[i][j]);
    return res;
  }
  inline const int rows() const {
    return entity.size();
  }
  inline const int cols() const {
    return ecols;
  }
  inline       void resize(const int& rows, const int& cols) {
    ecols = cols;
    entity.resize(rows);
    for(int i = 0; i < entity.size(); i ++)
      entity[i].resize(ecols);
    return;
  }
  myfloat      epsilon() const {
#if defined(_PERSISTENT_) && ! defined(_FLOAT_BITS_)
    // N.B. conservative.
    const myfloat eps(sqrt(myfloat(int(1)) >> myint((sizeof(size_t) * 16) - 1)));
    // static const myfloat eps(myfloat(int(1)) >> myint((sizeof(size_t) * 16) - 1));
#elif defined(_FLOAT_BITS_)
    // N.B. conservative.
    const myfloat eps(sqrt(myfloat(int(1)) >> myint(_FLOAT_BITS_ - 1)));
    // static const myfloat eps(myfloat(int(1)) >> myint(_FLOAT_BITS_ - 1));
#else
    // N.B. conservative.
    const myfloat eps(sqrt(std::numeric_limits<myfloat>::epsilon()));
    // static const myfloat eps(std::numeric_limits<myfloat>::epsilon());
#endif
    return eps;
  }
  // this isn't better idea for faster calculations.
#if defined(_SIMPLEALLOC_)
  vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > > entity;
#else
  vector<SimpleVector<T> > entity;
#endif
  int ecols;
};

template <typename T> inline T SimpleMatrix<T>::determinant(const bool& nonzero) const {
  T det(1);
  SimpleMatrix<T> work(*this);
  for(int i = 0; i < entity.size(); i ++) {
    int xchg = i;
    for(int j = i + 1; j < entity.size(); j ++)
      if(abs(work.entity[j][i]) > abs(work.entity[xchg][i]))
        xchg = j;
    swap(work.entity[i], work.entity[xchg]);
    const SimpleVector<T>& ei(work.entity[i]);
    const T& eii(ei[i]);
    if(! nonzero || ! i || pow(abs(det), T(int(1)) / T(int(i))) * epsilon() <= abs(eii))
      det *= eii;
    if(ei.dot(ei) * epsilon() < eii * eii) {
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
      for(int j = i + 1; j < entity.size(); j ++) {
        const T ratio(work.entity[j][i] / eii);
        work.entity[j] -= ei * ratio;
      }
    }
  }
  return det;
}

template <typename T> inline SimpleVector<T> SimpleMatrix<T>::solve(SimpleVector<T> other) const {
  if(! (0 <= entity.size() && 0 <= ecols && entity.size() == ecols && entity.size() == other.size()) ) return SimpleVector<T>();
  SimpleMatrix<T> work(*this);
  for(int i = 0; i < entity.size(); i ++) {
    int xchg = i;
    for(int j = i + 1; j < entity.size(); j ++)
      if(abs(work.entity[j][i]) > abs(work.entity[xchg][i]))
        xchg = j;
    swap(work.entity[i], work.entity[xchg]);
    swap(other[i], other[xchg]);
    const SimpleVector<T>& ei(work.entity[i]);
    const T& eii(ei[i]);
    if(ei.dot(ei) * epsilon() < eii * eii) {
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
      for(int j = i + 1; j < entity.size(); j ++) {
        const T ratio(work.entity[j][i] / eii);
        work.entity[j] -= ei       * ratio;
        other[j]       -= other[i] * ratio;
      }
    }
  }
  for(int i = entity.size() - 1; 0 <= i; i --) {
    if(work.entity[i][i] == T(int(0))) continue;
    const T buf(other[i] / work.entity[i][i]);
    if(!isfinite(buf) || isnan(buf)) {
      continue;
    }
    other[i] = buf;
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int j = i - 1; 0 <= j; j --)
      other[j] -= other[i] * work.entity[j][i];
  }
  return other;
}

template <typename T> inline SimpleVector<T> SimpleMatrix<T>::solveN(SimpleVector<T> other) const {
  if(! (0 <= entity.size() && 0 <= ecols && entity.size() == ecols && entity.size() == other.size()) ) return SimpleVector<T>();
  SimpleVector<T> res;
  return res;
}

template <typename T> inline SimpleVector<T> SimpleMatrix<T>::projectionPt(const SimpleVector<T>& other) const {
  // also needs class or this->transpose() * (*this) == I assertion is needed.
  SimpleMatrix<T> work(entity.size(), ecols);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < work.rows(); i ++)
    work.row(i) = entity[i] * entity[i].dot(other);
  SimpleVector<T> res(ecols);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < other.size(); i ++) {
    res[i] = T(int(0));
    for(int j = 0; j < entity.size(); j ++)
      res[i] += work(j, i);
  }
  return res;
}

template <typename T> inline SimpleMatrix<T>& SimpleMatrix<T>::fillP(const vector<int>& idx) {
  int ii(0);
  for(int j = 0; j < cols() && ii < idx.size(); j ++) {
    SimpleVector<T> ek(cols());
    ek.ek(j);
    ek -= this->projectionPt(ek);
    const T n2(ek.dot(ek));
    if(n2 <= epsilon()) continue;
    this->row(idx[ii ++]) = ek / sqrt(n2);
  }
  return *this;
}

template <typename T> inline SimpleMatrix<T> SimpleMatrix<T>::QR() const {
  const T norm2(norm2M(*this));
  if(! isfinite(norm2)) return *this;
  SimpleMatrix<T> Q(min(this->rows(), this->cols()), this->rows());
  Q.O();
  vector<int> residue;
  residue.reserve(Q.rows());
  for(int i = 0; i < Q.rows(); i ++) {
    SimpleVector<T> work(this->col(i));
    work -= Q.projectionPt(work);
    const T n2(work.dot(work));
    if(n2 <= norm2 * epsilon())
      residue.emplace_back(i);
    else
      Q.row(i) = (work /= sqrt(n2));
  }
  return Q.fillP(residue);
}

template <typename T> inline SimpleMatrix<T> SimpleMatrix<T>::SVDleft1d() const {
  // N.B. A = QR, (S - lambda I)x = 0 <=> R^t Q^t U = R^-1 Q^t U Lambda
  //        <=> R^t Q^t U Lambda' = R^-1 Q^t U Lambda'^(- 1)
  //      A := R^t, B := Q^t U, C := Lambda'
  //          (A + A^-t)*B*(C + C^-1) - A^-t B C - A B C^-1 = 2 A B C
  //        <=> (A + A^-t) * B * (C + C^-1) = (2I + 2A^-tA^-1) * ABC
  //        <=> B = A^-1 (2I + 2A^-(t+1))^-1 (A + A^-t) * B *
  //                (C + C^-1) * C^-1
  // N.B. since S is symmetric, singular value on SS^t = QRR^tQ^t is
  //      same square root as singular value on R.
  // N.B. this is INCOMPLETE, so geometrical non separable ones cannot be
  //      separated by this SVD function but not affects to the most of our
  //      repositories.
  const SimpleMatrix<T> S(*this * transpose());
  const SimpleMatrix<T> Qt(S.QR());
  const SimpleMatrix<T> A((Qt * S).transpose());
  const SimpleMatrix<T> A1t(A * A.transpose());
        SimpleMatrix<T> Left(A.inverse() * (SimpleMatrix<T>(A.rows(), A.cols()).I(T(int(2))) + A1t.inverse() * T(int(2))).inverse() * (A + A.transpose().inverse()));
        SimpleMatrix<T> Right(SimpleMatrix<T>(Left.rows(), Left.cols()).O());
  for(int i = 0; i < Right.rows(); i ++)
    Right(i, i) = A(i, i) + T(int(1));
  Left  /= sqrt(norm2M(Left));
  Right /= sqrt(norm2M(Right));
  // N.B. now we have B = Left * B * Right.
  const int p(absceil(sqrt(- log(epsilon()))));
  for(int i = 0; i < p; i ++) {
    Left  *= Left;
    Right *= Right;
  }
  return (Left * Qt * /* U * */ Right).QR() * Qt;
}

template <typename T> inline pair<SimpleMatrix<T>, SimpleMatrix<T> > SimpleMatrix<T>::SVD1d() const {
  if(this->rows() < this->cols()) {
    SimpleMatrix<T> R(this->transpose().SVDleft1d().transpose());
    return make_pair(((*this) * R).QR(), move(R));
  }
  SimpleMatrix<T> L(this->SVDleft1d());
  return make_pair(move(L), (L * (*this)).transpose().QR().transpose());
}

// XXX: O(n^4) over all, we need O(n^3) methods they make SVD1d as SVDnd.
template <typename T> inline SimpleMatrix<T> SimpleMatrix<T>::SVD() const {
  SimpleMatrix<T> sym((*this) * this->transpose());
  SimpleMatrix<T> res(sym);
  res.I();
  for(int i = 0; i <= sym.rows() + 1; i ++) {
    SimpleMatrix<T> svd(sym.SVD1d());
    sym = svd * sym * sym.transpose().SVD1d().transpose();
    res = svd * res;
  }
  return res;
}

template <typename T> inline pair<pair<SimpleMatrix<T>, SimpleMatrix<T> >, SimpleMatrix<T> > SimpleMatrix<T>::SVD(const SimpleMatrix<T>& src) const {
  const T norm2(max(norm2M(*this), norm2M(src)));
  if(! isfinite(norm2)) return *this;
  // refered from : https://en.wikipedia.org/wiki/Generalized_singular_value_decomposition .
  SimpleMatrix<T> C(this->rows() + src.rows(), this->cols());
  C.setMatrix(0, 0, *this);
  C.setMatrix(this->rows(), 0, src);
  const SimpleMatrix<T> P(C.SVD());
  const SimpleMatrix<T> Qt(C.transpose().SVD().transpose());
  const SimpleMatrix<T> D(P * C * Qt.transpose());
  SimpleMatrix<T> P1(this->rows(), this->cols());
  SimpleMatrix<T> P2(src.rows(), this->cols());
  for(int i = 0; i < P1.rows(); i ++)
    P1.row(i) = P.col(i);
  for(int i = 0; i < P2.rows(); i ++)
    P2.row(i) = P.col(i + P1.rows());
  SimpleMatrix<T> U1(P1.SVD());
  SimpleMatrix<T> Wt(P1.transpose().SVD().transpose());
  SimpleMatrix<T> U2(Wt * P2.transpose());
  vector<int> fill;
  fill.reserve(U2.rows());
  for(int i = 0; i < U2.rows(); i ++) {
    const T n2(U2.row(i).dot(U2.row(i)));
    if(n2 <= epsilon())
      fill.emplace_back(i);
    else
      U2.row(i) /= sqrt(n2);
  }
  return make_pair(make_pair(move(U1), move(U2.fillP(fill))), (Wt * D).transpose().QR() * Qt);
}

template <typename T> inline SimpleVector<T> SimpleMatrix<T>::zeroFix(const SimpleMatrix<T>& A, vector<pair<T, int> > fidx) {
  // N.B. we now have |[A -bb] [x t]| <= 1 condition.
  // N.B. there's no difference |[A - bb] [x t]|^2 <= 1 condition in this.
  //      but not with mixed condition.
  const SimpleMatrix<T> R((*this) * A);
  SimpleVector<T> one(this->cols());
  one.I(T(int(1)));
  for(int i = 0; i < fidx.size(); i ++) {
    one[fidx[i].second] = - fidx[i].first;
    fidx[i].first = - T(int(1));
  }
  // N.B. we now have: Q [R [x t] ] <= {0, 1}^m cond.
  //      so this isn't normalized by super spherical ones.
  //      but it's normalized by elliptic like ones.
  SimpleMatrix<T> Pb(*this);
  const SimpleVector<T> on(projectionPt(one));
  fidx.reserve(fidx.size() + this->cols());
  for(int i = 0; i < this->cols(); i ++)
    fidx.emplace_back(make_pair(abs(on[i]), i));
  vector<bool> fixed;
  fixed.resize(this->cols(), false);
  // sort by: |<Q^t(1), q_k>|, we subject to minimize each, to do this,
  //   maximize minimum q_k orthogonality.
  for(int i = 0; i < this->rows() - 1; i ++) {
    int idx(- 1);
    T   M(int(0));
    for(int j = 0; j < fidx.size(); j ++) if(fidx[j].first < T(int(0))) idx = j;
    if(idx < 0)
      for(int j = 0; j < fidx.size(); j ++)
        if(fidx[j].first != T(int(0)) && (idx < 0 || fidx[idx].first < M)) M = fidx[idx = j].first;
    const int& iidx(fidx[idx].second);
    const SimpleVector<T>  orth(this->col(iidx));
    const T n2(orth.dot(orth));
    if(T(int(0)) < fidx[idx].first &&
       fidx[idx].first < sqrt(one.dot(one)) * epsilon()) {
      break;
    }
    // N.B. rank(*this) on call is max rank normally, should not be singular.
    //      however, masp calls with rank isn't max cases.
    if(n2 <= epsilon())
      continue;
    fixed[fidx[idx].second] = true;
    Pb = *this;
    // N.B. O(mn) can be written into O(lg m + lg n) in many core cond.
    for(int j = 0; j < this->cols(); j ++)
      this->setCol(j, this->col(j) - orth * this->col(j).dot(orth) / n2);
    if(T(int(0)) < fidx[idx].first) {
      if(one.size() != fidx.size()) fidx.resize(one.size());
      const SimpleVector<T> on(projectionPt(one));
      for(int j = 0; j < this->cols(); j ++) if(fixed[j]) { fidx[j].first = T(int(0)); fidx[j].second = j; } else { fidx[j].first = abs(on[j]); fidx[j].second = j; }
    } else fidx[idx].first = T(int(1));
    i ++;
  }
  {
    const SimpleVector<T> on(projectionPt(one));
    if(sqrt(on.dot(on)) < sqrt(one.dot(one)) * epsilon()) *this = move(Pb);
  }
  // N.B. now we have fix indices to be P R [x 1] * t == 0.
  return R.solve((*this) * one);
}

template <typename T> inline SimpleVector<T> SimpleMatrix<T>::inner(const SimpleVector<T>& bl, const SimpleVector<T>& bu) const {
  // |(2 / bu) A x - 1 - bl / bu| <= |1 - bl / bu|
  // <=> with (-A, -bu, -bl), |bl| <= |bu|, |(2 / bu) A x - 2| <= 2(1 - bl / bu)
  SimpleVector<T> bU(bu);
  SimpleVector<T> bL(bl);
  SimpleMatrix<T> A(*this);
  vector<pair<T, int> > fidx;
  for(int i = 0; i < bU.size(); i ++) {
    if(abs(bu[i]) < abs(bl[i])) {
      bU[i] = - bl[i];
      bL[i] = - bu[i];
      A.row(i) = - this->row(i);
    } else if(bu[i] == bl[i])
      fidx.emplace_back(make_pair(- T(int(bu[i] == T(0) ? 0 : 1)), i));
    A.row(i) /= (T(2) * bU[i] - bL[i]) / T(2);
  }
  // N.B. in zeroFix, we get linear Invariant s.t. |Ax| <= 1 possible enough.
        SimpleVector<T> res(A.QR().zeroFix(A, fidx));
  const SimpleVector<T> z(*this * res * T(int(4)));
        T    t(int(1));
  for(int i = 0; i < z.size(); i ++)
    if(bu[i] * z[i] < T(int(0))) // N.B.: infeasible.
      continue;
    else if(z[i] != T(int(0)))
      t = bl[i] * z[i] < T(int(0)) ? min(t, bu[i] / z[i]) :
                  min(t, min(bu[i] / z[i], bl[i] / z[i]));
  return res *= t;
}

template <typename T> static inline T norm2M(const SimpleMatrix<T>& m) {
  T norm2(m.row(0).dot(m.row(0)));
  for(int i = 1; i < m.rows(); i ++)
    norm2 = max(norm2, m.row(i).dot(m.row(i)));
  return norm2;
}

template <typename T> static inline SimpleMatrix<T> log(const SimpleMatrix<T>& m) {
  const int cut(- log(SimpleMatrix<T>().epsilon()) / log(T(int(2))) * T(int(2)) );
  SimpleMatrix<T> res(m.rows(), m.cols());
  const T c(sqrt(norm2M(m)) * T(2));
  const SimpleMatrix<T> residue(SimpleMatrix<T>(m.rows(), m.cols()).I() - m / c);
        SimpleMatrix<T> buf(residue);
  res.I(log(c));
  for(int i = 1; 0 < i && i < cut; i ++) {
    res -= buf / T(i);
    buf *= residue;
  }
  return res;
}

template <typename T> static inline SimpleMatrix<T> logSym(const SimpleMatrix<T>& x, const SimpleMatrix<T>& b) {
  // N.B. Ux Lx Uxt == X := B^A == Ub Ua exp(La) Lb_k Uat Ubt.
  const SimpleMatrix<T> Ub(b.SVD());
  const SimpleMatrix<T> Ubt(b.transpose().SVD());
  const SimpleMatrix<T> Lb(Ub * b * Ubt.transpose());
        SimpleMatrix<T> Ux(x.SVD());
        SimpleMatrix<T> Uxt(x.transpose().SVD());
  const SimpleMatrix<T> Lx(Ux * x * Uxt.transpose());
  // N.B. Lx == [[La log(Lb_k)]]
  SimpleVector<T> Lawork(Ux.rows());
  Lawork.O();
  for(int i = 0; i < Lawork.size(); i ++)
    // XXX stub:
    Lawork[i] = log(Lx(i, i) / Lb(i / (x.rows() / b.rows())));
  // N.B. Lawork.subVector... == another subVector in Ua diag(La) Uat condition.
  // XXX: might be a wrong method.
  SimpleMatrix<T> UUb(x.rows(), x.cols());
  SimpleMatrix<T> UUbt(UUb);
  UUb.O();
  UUbt.O();
  for(int i = 0; i < Ub.rows(); i ++)
    for(int j = 0; j < Ub.cols(); j ++) {
      UUb.setMatrix( i * (x.rows() / Ub.rows()),  j * (x.cols() / Ub.cols()),
        SimpleMatrix<T>(Ub.rows(), Ub.cols()).I(Ub(i, j)) );
      UUbt.setMatrix(i * (x.rows() / Ubt.rows()), j * (x.cols() / Ubt.cols()),
        SimpleMatrix<T>(Ubt.rows(), Ubt.cols()).I(Ubt(i, j)) );
    }
  Ux  = Ub.inverse()  * Ux;
  Uxt = Ubt.inverse() * Uxt;
  // N.B. Ux == Ua sqrt(diag(Lawork)) and same on right side.
  return x;
}

template <typename T> static inline SimpleMatrix<T> exp01(const SimpleMatrix<T>& m) {
  SimpleMatrix<T> res(m.rows(), m.cols());
  const int cut(- log(SimpleMatrix<T>().epsilon()) / log(T(int(2))) * T(int(2)) );
  SimpleMatrix<T> buf(m);
  res.I();
  for(int i = 1; 0 < i && i < cut; i ++) {
    res += buf;
    buf *= m / T(i + 1);
  }
  return res;
}

template <typename T> static inline SimpleMatrix<T> pow(const SimpleMatrix<T>& m, const T& p) {
  return exp(log(m) * p);
}

template <typename T> static inline SimpleMatrix<complex(T) > dft(const int& size0) {
  const int size(abs(size0));
  if(! size) {
    const SimpleMatrix<complex(T) > m0;
    return m0;
  }
  SimpleMatrix<complex(T) > edft( size, size);
  SimpleMatrix<complex(T) > eidft(size, size);
    const T Pi(T(4) * atan2(T(1), T(1)));
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < edft.rows(); i ++)
      for(int j = 0; j < edft.cols(); j ++) {
        const T theta(- T(int(2)) * Pi * T(i) * T(j) / T(edft.rows()));
        const T c(cos(theta));
        const T s(sin(theta));
        edft( i, j) = complexctor(T)(c,   s);
        eidft(i, j) = complexctor(T)(c, - s) / complexctor(T)(T(size));
      }
  return size0 < 0 ? eidft : edft;
}

// N.B. integrate(diff) isn't get original but is reasonable on IDFT*DFT meaning.
template <typename T> static inline SimpleMatrix<T> diff(const int& size0) {
  const int size(abs(size0));
  if(! size) {
    const SimpleMatrix<T> m0;
    return m0;
  }
  SimpleMatrix<T> dd;
  SimpleMatrix<T> ii;
    // N.B. if we return recursive each size diff,
    //      taylor series should be broken.
    SimpleMatrix<complex(T) > DD(dft<T>(size));
    SimpleMatrix<complex(T) > II(dft<T>(size));
    const T  Pi(T(4) * atan2(T(1), T(1)));
    // N.B. we should start this loop with i == 1 on integrate(diff) or inverse.
    //      we also should start with i == 0 on taylor series.
    //      we select latter one.
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < DD.rows(); i ++)
      DD.row(i) *= complexctor(T)(T(int(0)), - T(int(2)) * Pi * T(i) / T(DD.rows()));
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 1; i < II.rows(); i ++)
      II.row(i) /= complexctor(T)(T(int(0)), - T(int(2)) * Pi * T(i) / T(DD.rows()));
    // N.B. if we apply DD onto 1 / (1 / f(x)) graph, it's reverse order.
    //      if we average them, it's the only 0 vector.
    // N.B. there exists also completely correct differential matrix,
    //      but it's also be only 0 vector.
    //      (because it's a imaginary part on originals.)
    // N.B. in continuous function, we don't divide dd by &pi;.
    //      (d/dx sum exp(2 Pi i x theta / N) exp(- 2 Pi i y theta / N) f(y)).
    // N.B. from some numerical test, sign on DD, II is reverse side.
    // N.B. In discrete function, we had be chosen ||dd|| := 1, and not now.
    //      (because the (left or right) differential itself cannot be
    //       larger than |x_{k+1}-x_k| < ||x|| in discrete,
    //       sum_0^1 - 2 Pi i (theta/n)^2/2 -&gt; Pi)
    // N.B. if we make plain differential with no error on cosine curve,
    //      it causes constant 0 vector.
    // N.B. if we don't take this real operator, we cannot get better accuracy
    //      result on taylor series.
    dd =   (dft<T>(- size) * DD).template real<T>();
    ii = - (dft<T>(- size) * II).template real<T>();
  return size0 < 0 ? ii : dd;
}

template <typename T> static inline SimpleVector<complex(T) > taylorc(const int& size, const T& step, const T& stepw) {
  const int step00(max(int(0), min(size - 1, int(absfloor(step)))));
  const T   residue0(step - T(step00));
  const int step0(step00 == size - 1 || abs(residue0) <= T(int(1)) / T(int(2)) ? step00 : step00 + 1);
  const T   residue(step - T(step0));
  if(residue == T(int(0))) return SimpleVector<complex(T) >(size).ek(step0);
  const T   residuem(residue - (step - stepw));
  // N.B. following code is equivalent to exp each dft.
  //      this improves both accuracy and speed.
  // N.B. We don't need to matter which sign dft/idft uses till the sign
  //      we multiply is bonded to the transformation.
  const T Pi(T(4) * atan2(T(1), T(1) ));
  SimpleVector<complex(T) > res(dft<T>(- size).row(step0));
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < res.size(); i ++)
    res[i] *= step != stepw ? 
      (i ? exp(complexctor(T)(T(int(0)), - T(int(2)) * Pi * T(i) * residue / T(res.size()) ))
          / complexctor(T)(T(int(0)), - T(int(2)) * Pi * T(i) / T(res.size()) )
      - exp(complexctor(T)(T(int(0)), - T(int(2)) * Pi * T(i) * residuem / T(res.size()) ))
          / complexctor(T)(T(int(0)), - T(int(2)) * Pi * T(i) / T(res.size()) ) :
        complexctor(T)(T(int(0))) )
      : exp(complexctor(T)(T(int(0)), - T(int(2)) * Pi * T(i) * residue / T(res.size()) ));
  return dft<T>(size).transpose() * res;
}

template <typename T> static inline SimpleVector<T> taylor(const int& size, const T& step, const T& stepw) {
  return taylorc<T>(size, step, stepw).template real<T>();
}

template <typename T> static inline SimpleVector<T> taylor(const int& size, const T& step) {
  return taylor<T>(size, step, step);
}

// N.B. we only need cosine value on invariant, so normalize them into S^n.
template <typename T> static inline SimpleVector<T> linearInvariant(const SimpleMatrix<T>& in) {
  vector<pair<T, int> > sute;
  SimpleVector<T> res(in.QR().zeroFix(in, sute));
  const T nres(res.dot(res));
  return ! isfinite(nres) || nres <= T(int(0)) ? res.O() :
    res /= sqrt(nres);
}

// --- N.B. start small only to enname functions ---
// N.B. functions between R and [0,1], ]0,1[.
myfloat* bm_sqe;
myfloat* bm_denom;
template <typename T> static inline T binMargin(const T& in) {
  if(! bm_sqe) {
    bm_sqe   = SimpleAllocator<T>().allocate(1);
    bm_denom = SimpleAllocator<T>().allocate(1);
    ::new ((void*)bm_sqe) T();
    ::new ((void*)bm_denom) T();
    *bm_sqe    = sqrt(SimpleMatrix<T>().epsilon());
    *bm_denom = T(int(1)) + sqrt(*bm_sqe);
  }
  // N.B. better 0 handling, {0, 1} vanished before.
  T res(in + *bm_sqe);
  // N.B. CPU float glitch.
  res /= *bm_denom;
  return res;
}

template <typename T> static inline SimpleVector<T> binMargin(const SimpleVector<T>& in) {
  SimpleVector<T> res(in);
  for(int i = 0; i < res.size(); i ++) res[i] = binMargin<T>(res[i]);
  return res;
}

template <typename T> static inline T offsetHalf(const T& in, const T& o = T(int(1)) ) {
  return (in + o) / (T(int(1)) + o);
}

template <typename T> static inline SimpleVector<T> offsetHalf(const SimpleVector<T>& in, const T& o = T(int(1)) ) {
  SimpleVector<T> res(in);
  for(int i = 0; i < res.size(); i ++) res[i] = offsetHalf<T>(res[i], o);
  return res;
}

#if defined(_SIMPLEALLOC_)
template <typename T> static inline vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > > offsetHalf(const vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > >& in, const T& o = T(int(1)) ) {
  vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > > res(in);
#else
template <typename T> static inline vector<SimpleVector<T> > offsetHalf(const vector<SimpleVector<T> >& in, const T& o = T(int(1)) ) {
  vector<SimpleVector<T> > res(in);
#endif
  for(int i = 0; i < res.size(); i ++) res[i] = offsetHalf<T>(res[i], o);
  return res;
}

template <typename T> static inline SimpleVector<SimpleVector<T> > offsetHalf(const SimpleVector<SimpleVector<T> >& in, const T& o = T(int(1)) ) {
  SimpleVector<SimpleVector<T> > res;
  res.entity = offsetHalf<T>(in.entity, o);
  return res;
}

template <typename T> static inline T unOffsetHalf(const T& in, const T& o = T(int(1)) ) {
  return in * (T(int(1)) + o) - o;
}

template <typename T> static inline SimpleVector<T> unOffsetHalf(const SimpleVector<T>& in, const T& o = T(int(1)) ) {
  SimpleVector<T> res(in);
  for(int i = 0; i < res.size(); i ++) res[i] = unOffsetHalf<T>(res[i], o);
  return res;
}

#if defined(_SIMPLEALLOC_)
template <typename T> static inline vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > > unOffsetHalf(const vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > >& in, const T& o = T(int(1)) ) {
  vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > > res(in);
#else
template <typename T> static inline vector<SimpleVector<T> > unOffsetHalf(const vector<SimpleVector<T> >& in, const T& o = T(int(1)) ) {
  vector<SimpleVector<T> > res(in);
#endif
  for(int i = 0; i < res.size(); i ++) res[i] = unOffsetHalf<T>(res[i], o);
  return res;
}

template <typename T> static inline SimpleVector<SimpleVector<T> > unOffsetHalf(const SimpleVector<SimpleVector<T> >& in, const T& o = T(int(1)) ) {
  SimpleVector<SimpleVector<T> > res;
  res.entity = unOffsetHalf<T>(in.entity, o);
  return res;
}

template <typename T> static inline T clipBin(const T& in) {
  const T zero(int(0));
  const T one(int(1));
  return max(zero, min(one, in));
}

template <typename T> static inline SimpleVector<T> clipBin(const SimpleVector<T>& in) {
  SimpleVector<T> res(in);
  for(int i = 0; i < res.size(); i ++) res[i] = clipBin<T>(res[i]);
  return res;
}

template <typename T> static inline vector<SimpleVector<T> > clipBin(const vector<SimpleVector<T> >& in) {
  vector<SimpleVector<T> > res(in);
  for(int i = 0; i < res.size(); i ++) res[i] = clipBin<T>(res[i]);
  return res;
}

template <typename T> static inline T cutBin(const T& in) {
  const T zero(int(0));
  const T one(int(1));
  T res(in - absfloor(in));
  return res <= zero ? res += one : res;
}

template <typename T> static inline T R2bin(const T& in) {
  return offsetHalf<T>(atan(- in) / atan(T(int(1))) / T(int(2)));
}

template <typename T> static inline SimpleVector<T> R2bin(const SimpleVector<T>& in) {
  SimpleVector<T> res(in);
  for(int i = 0; i < res.size(); i ++) res[i] = R2bin<T>(res[i]);
  return res;
}

template <typename T> static inline T bin2R(const T& in) {
  return - tan(max(- T(int(1)) + sqrt(SimpleMatrix<T>().epsilon()),
               min(  T(int(1)) - sqrt(SimpleMatrix<T>().epsilon()),
                 unOffsetHalf<T>(in) )) * atan(T(int(1))) * T(int(2)) );
}

template <typename T> static inline SimpleVector<T> bin2R(const SimpleVector<T>& in) {
  SimpleVector<T> res(in);
  for(int i = 0; i < res.size(); i ++) res[i] = bin2R<T>(res[i]);
  return res;
}

template <typename T> T sgn(const T& x) {
  const T zero(0);
  const T one(1);
  const T mone(- T(int(1)));
  return x != zero ? (zero < x ? one : mone) : zero;
}

template <typename T> static inline T expscale(const T& x) {
  return sgn<T>(x) * (exp(abs(x)) - T(int(1))) /
    (exp(T(int(1))) - T(int(1)));
}

template <typename T> static inline SimpleVector<T> expscale(const SimpleVector<T>& x) {
  SimpleVector<T> res(x);
  for(int i = 0; i < res.size(); i ++) res[i] = expscale<T>(res[i]);
  return res;
}

#if defined(_SIMPLEALLOC_)
template <typename T> static inline vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > > expscale(const vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > >& x) {
  vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > > res(x);
#else
template <typename T> static inline vector<SimpleVector<T> > expscale(const vector<SimpleVector<T> >& x) {
  vector<SimpleVector<T> > res(x);
#endif
  for(int i = 0; i < res.size(); i ++) res[i] = expscale<T>(res[i]);
  return res;
}

template <typename T> static inline T logscale(const T& x) {
  return sgn<T>(x) * log(abs(x) + T(int(1))) / log(T(int(2)));
}

template <typename T> static inline SimpleVector<T> logscale(const SimpleVector<T>& x) {
  SimpleVector<T> res(x);
  for(int i = 0; i < res.size(); i ++) res[i] = logscale<T>(res[i]);
  return res;
}

#if defined(_SIMPLEALLOC_)
template <typename T> static inline vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > > logscale(const vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > >& x) {
  vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > > res(x);
#else
template <typename T> static inline vector<SimpleVector<T> > logscale(const vector<SimpleVector<T> >& x) {
  vector<SimpleVector<T> > res(x);
#endif
  for(int i = 0; i < res.size(); i ++) res[i] = logscale<T>(res[i]);
  return res;
}

#if defined(_SIMPLEALLOC_)
template <typename X> static inline vector<X, SimpleAllocator<X> > skipX(const vector<X, SimpleAllocator<X> >& in, const int& step = 1) {
  vector<X, SimpleAllocator<X> > res;
#else
template <typename X> static inline vector<X> skipX(const vector<X>& in, const int& step = 1) {
  vector<X> res;
#endif
  res.resize((in.size() + step - 1) / step);
  for(int i = (in.size() - 1) % step, ii = 0; i < in.size();
          i += step, ii ++) res[ii] = in[i];
  return res;
}

#if defined(_SIMPLEALLOC_)
template <typename X> static inline vector<X, SimpleAllocator<X> > delta(const vector<X, SimpleAllocator<X> >& in0) {
  vector<X, SimpleAllocator<X> > in(in0);
#else
template <typename X> static inline vector<X> delta(const vector<X>& in0) {
  vector<X> in(in0);
#endif
  for(int i = 1; i < in.size(); i ++) in[i] = in0[i] - in0[i - 1];
  return in;
}

template <typename T> static inline T getImgPt(const T& y, const T& h) {
  T yy(y % (2 * h));
  if(yy < 0)
    yy = - yy;
  if(yy >= h)
    yy = h - (yy - h);
  return yy % h;
}

// N.B. start raw prediction operations.
// N.B. please refer bitsofcotton/randtools.
template <typename T> static inline pair<SimpleVector<T>, T> makeProgramInvariant(const SimpleVector<T>& in, const T& index = - T(int(1)) ) {
  SimpleVector<T> res(in.size() + (T(int(0)) <= index ? 2 : 1));
  res.setVector(0, in);
  res[in.size()] = T(int(1));
  if(T(int(0)) <= index)
    res[in.size() + 1] = T(index);
  T ratio(0);
  for(int i = 0; i < res.size(); i ++)
    ratio += log(res[i] = binMargin<T>(res[i]));
  // N.B. x_1 ... x_n == 1.
  // <=> x_1 / (x_1 ... x_n)^(1/n) ... == 1.
  ratio = isfinite(ratio) ? exp(- ratio / T(res.size())) : T(int(1));
  return make_pair(res *= ratio, ratio);
}

template <typename T> static inline T revertProgramInvariant(const pair<T, T>& in) {
  return cutBin<T>(in.second == T(0) ?
    sgn<T>(in.second) / SimpleMatrix<T>().epsilon() : in.first / in.second);
}

template <typename T> static inline SimpleVector<T> revertProgramInvariant(const pair<SimpleVector<T>, T>& in) {
  SimpleVector<T> res(in.first);
  for(int i = 0; i < in.first.size(); i ++)
    res[i] = revertProgramInvariant<T>(make_pair(res[i], in.second));
  return res;
}

template <typename T, bool nonlinear> static inline T revertByProgramInvariant(SimpleVector<T> work, const SimpleVector<T>& invariant) {
  const T one(int(1));
  const T two(int(2));
  const int idx(work.size() - 1);
  if(invariant[idx] == T(int(0))) return work[idx - 1];
  work[idx] = T(int(0));
  if(nonlinear) {
    // N.B. t == revertProgramInvariant<T>(t0 + t, work2.second)).
    //      t == (t0 + t) / s, s == geometric average on work
    //      t^((n+1)/n) * s0 - t - t0 = 0.
    // newton's method with:
    //       f'(t) == (n+1)/n t^(1/n) - 1.
    //       f (t) == t^((n+1)/n) * s0 - t - t0
    // => t_{k+1} == t_k - (t^((n+1)/n)*s0-t-t0)/((n+1)/n*t^(1/n)-1)
          pair<SimpleVector<T>, T> vdp(makeProgramInvariant<T>(work, T(int(1))));
    const T nvdp(sqrt(vdp.first.dot(vdp.first)));
    vdp.first  /= nvdp;
    vdp.second /= nvdp;
    const T one(int(1));
    const int loop(T(int(2)) * sqrt(- log(SimpleMatrix<T>().epsilon()) /
      log(T(int(2))) ));
          complex(T) t(T(int(0)));
    const complex(T) t0(invariant.dot(vdp.first));
          complex(T) s0(T(int(0)));
    const T n(vdp.first.size());
    for(int i = 0; i < vdp.first.size(); i ++) if(i != idx)
      s0 += log(complexctor(T)(binMargin<T>(vdp.first[i])));
    s0  = exp(s0 /= complexctor(T)(T(vdp.first.size() - 1)));
    for(int i = 0; i <= loop; i ++)
      t -= (pow(t, complexctor(T)((n + one) / n)) * s0 - t - t0) /
        (complexctor(T)((n + one) / n) * pow(t, complexctor(T)(one / n)) -
          complexctor(T)(one) );
    work[idx] = (t.real() /= vdp.second);
  } else work[idx] = - invariant.dot(work) / invariant[idx];
  bool shown(false);
  work[idx] = cutBin<T>(work[idx]);
  return isfinite(work[idx]) ? work[idx] : T(int(0));
}

// N.B. F_2 case #f fixation on R^indim into R^result fixation.
//      if we object F_3 or more accuracy, we need to increase the result.
// N.B. there's also analogy to this as {0,1}^m*n operator orthogonalization
//      invariant residue maximum dimensions. (sqrt(3!) <~ (#{0,1})^2).
// N.B. also the binary 2 operand opereator is described in R^4 vector
//      as a invariant.
// N.B. this is for output is binary case especially sign bit on the
//      information amount on any p-adics but in [0, 1[.
static inline int ind2vd(const int& indim) {
  const int y(indim / 2 * (indim / 2 - 1));
        int varlen(4);
  for( ; 0 <= varlen && varlen < sizeof(int) * 8 &&
    (varlen - 1) * (varlen - 1) << varlen <= y; varlen ++) ;
  return max(-- varlen, int(4));
}

template <typename T> static inline SimpleVector<T> minsq(const int& size) {
  const T xsum(size * (size - 1) / 2);
  const T xdot(size * (size - 1) * (2 * size - 1) / 6);
  const T denom(xdot * T(size) - xsum * xsum);
  SimpleVector<T> s(size);
  for(int i = 0; i < s.size(); i ++)
    s[i] = (T(i) * T(size) - xsum) / denom;
  return s;
}

vector<vector<SimpleVector<myfloat> > > *pncr_cp;
template <typename T> const SimpleVector<T>& pnextcacher(const int& size, const int& step) {
  if(! pncr_cp) {
    pncr_cp = SimpleAllocator<vector<vector<SimpleVector<T> > > >().allocate(1);
    ::new ((void*)pncr_cp) vector<vector<SimpleVector<T> > >();
  }
  vector<vector<SimpleVector<T> > >& cp(*pncr_cp);
  if(cp.size() <= size)
    cp.resize(size + 1, vector<SimpleVector<T> >());
  if(cp[size].size() <= step)
    cp[size].resize(step + 1, SimpleVector<T>());
  if(cp[size][step].size()) return cp[size][step];
  return cp[size][step] = (dft<T>(- size) * (dft<T>(size * 2).subMatrix(0, 0, size, size * 2) * taylorc<T>(size * 2, T(step < 0 ? step * 2 : (size + step) * 2 - 1), T(step < 0 ? step * 2 + 2 : (size + step) * 2 - 3)) )).template real<T>();
}

// N.B. f in C1 case, F(z,theta) := complex(f(z+~z),f(z-~z)*tan(theta)) in C1,
//      z in C. For each theta, exists F: holomorphic function at some axis
//      that real axis is same as f. so |theta-pi/4|<=epsilon case,
//      we have a area, with guzmer's inequation and edit integrate path.
//
//      f have laurent series and as a upper bound of coefficinets, we can cut
//      of them with some error. if there's no essential singular point around
//      {x | |x-a|<1}, the correct result might be gained.
//
//      also we use weak differential taylor series on them.
//
// N.B. we can calculate real F as
//   exp(Sum log(z-|z|cis(pi/4+t_k)))
//     Sum((f(z+~z)+i*f(z-~z)*tan(pi/4+t_k))/(z-|z|cis(pi/4+t_k)))
//   =: f(z+~z)*g(z)+f(z-~z)*h(z)
//   so replace on (1+i)t, F(x):=(f(x)*real(g(x))-f(x)*imag(h(x))) + imaginary
//   F is holomorphic on some small range (f in C0) and
//   real(F)=f(x)*some G(x,f). so we can apply this condition with cauchy's
//   integrate theorem on ja.wikipedi.org link:
//   doi:10.1090/S0002-9947-1900-1500519-7 C1 condition to C0 condition with
//   some small area.
// N.B. also this is left differential one and periodical right differential one
//   average.
// N.B. if there exists f/(linear transform) =: g is reasonable, g have a
//   structure x+=Ax. with A.cols==n or n-markov.
// N.B. also with 2^x:=[1,x0,...,xn,x0&x1,...,x(n-1)&xn,...,x0&...&xn] form,
//   the operation '&' and '~' can be described as each taylor series also
//   2^y:=A*(2^x) A in R^(N*N), 2^x in {0,1}^N.
//   this concludes the recursive structure as:
//     Sum A_k*cosh(a_k x)+B_k*sinh(b_k x) because A^n calculation.
template <typename T> static inline T p0next(const SimpleVector<T>& in) {
  return pnextcacher<T>(in.size(), 1).dot(in);
}

template <typename T, T (*f)(const SimpleVector<T>&)> static inline T invNext(const SimpleVector<T>& in) {
  const T zero(int(0));
  const T one(int(1));
  SimpleVector<T> ff(in);
  for(int i = 0; i < in.size(); i ++) if(in[i] == zero) return in[in.size() - 1];
  else ff[i] = one / in[i];
  const T pn(f(ff));
  if(pn == zero) return in[in.size() - 1];
  return one / pn;
}

// N.B. some of the essential point hack.
myfloat* npole_M;
template <typename T, T (*f)(const SimpleVector<T>&)> static inline T northPoleNext(const SimpleVector<T>& in) {
  const T zero(int(0));
  const T one(int(1));
  if(! npole_M) {
    npole_M = SimpleAllocator<T>().allocate(1);
    ::new ((void*)npole_M) T();
    * npole_M = atan(one / sqrt(SimpleMatrix<T>().epsilon()));
  }
  const T& M(* npole_M);
  SimpleVector<T> ff(in);
  for(int i = 0; i < in.size(); i ++)
    if(! isfinite(in[i]) || in[i] == zero) return in[in.size() - 1];
    else {
      ff[i] = atan(in[i]);
      // N.B. we avoid right hand side, it's harmless.
      // ff[i] = atan(one / ff[i]);
      // assert(- M < ff[i] && ff[i] < M);
    }
  T work(f(ff));
  // if(! isfinite(work) || work == zero) return in[in.size() - 1];
  if(! isfinite(work)) return in[in.size() - 1];
  // work = tan(max(- M, min(M, one / tan(max(- M, min(M, work))))));
  work = tan(max(- M, min(M, work)));
  if(isfinite(work)) return work;
  return in[in.size() - 1];
}

// N.B. we can add some of the conditions on x_next := integrate^x f(x_now)
//   form with x'_next := integrate^x (f_x'now - alpha) + beta transforms.
template <typename T, bool avg, T (*f)(const SimpleVector<T>&)> static inline T sumCNext(const SimpleVector<T>& in) {
  SimpleVector<T> ff(in);
  for(int i = 1; i < ff.size(); i ++)
    ff[i] += ff[i - 1];
  if(! avg) return f(ff) - ff[ff.size() - 1];
  const T A(ff[ff.size() - 1] / T(ff.size()));
  for(int i = 0; i < ff.size(); i ++)
    ff[i] = in[i] - A;
  return f(ff) + A;
}

// N.B. Sum(d_k)/Sum(d_(k-1)) - 1 with i-axis plotted Sum f'/f goes to near
//      log(f), once goes log(f) + i pi/2, the series can be arg(z) depend one.
template <typename T, T (*f)(const SimpleVector<T>&)> static inline T logCNext(const SimpleVector<T>& in) {
  const T zero(int(0));
  const T one(int(1));
  SimpleVector<T> ff(in);
  if(ff[0] == zero) return in[in.size() - 1];
  for(int i = 1; i < ff.size(); i ++)
    if((ff[i] += ff[i - 1]) == zero) return in[in.size() - 1];
  SimpleVector<T> gg(ff.size() - 1);
  gg.O();
  for(int i = 1; i < ff.size(); i ++)
    if(! isfinite(gg[i - 1] = ff[i] / ff[i - 1] - one)) return in[in.size() - 1];
  return f(gg) * ff[ff.size() - 1];
}

template <typename T> static inline T p0max0next(const SimpleVector<T>& in) {
  // N.B. on existing taylor series in surface.
  return (sumCNext<T, true, p0next<T> >(in) +
    invNext<T, sumCNext<T, true, p0next<T> > >(in)) / T(int(2));
}

template <typename T> static inline T p0maxNext(const SimpleVector<T>& in) {
  // N.B. we only handle Riemann measurable and R(finite)-valued functions.
  //      so worse structures are handled by p01next or p012next.
  // N.B. o-minimal
  //      (https://ja.wikipedia.org/wiki/%E5%AE%9F%E9%96%89%E4%BD%93
  //       (2022/03/19)) continuous structure causes dim K == 1,2,4,8 real
  //      closed field if they're semi-ordered one.
  // N.B. we treat periodical part as non aligned complex arg part.
  // N.B. we make the prediction on (delta) summation also take average as
  //      origin of input.
  return in.size() <= 6 ?
    sumCNext<T, true, sumCNext<T, false, northPoleNext<T,
      p0max0next<T> > > >(in) :
    sumCNext<T, true, sumCNext<T, false, logCNext<T, logCNext<T,
      northPoleNext<T, p0max0next<T> > > > > >(in);
  // N.B. we need only once P0DFT in general because associative condition
  //      is necessary for input ordering even we work with sedenion.
  //      also this eliminates one dimension per each of complex-formed input
  //      on f as a continuous thing but this needs huge memory to run also
  //      if the original predictor is linear, they're only transparent id.
  //      transformation.
  // N.B. on any R to R into reasonable C^omega.
  // N.B. either there's plenty of a space to extend this with
  //      uparrow, downarrow operations they causes the result in H\C if second
  //      operand is in C\R.
  // N.B. however we don't need this in normal condition because if the
  //      prediction itself is linear, they doesn't attach the result.
  // return sumCNext<T, true, sumCNext<T, false, logCNext<T, logCNext<T,
  //   P0DFT<T, p0max0next<T> > > > > >(in);
}

// Get invariant structure that
// [0,1[-register computer with deterministic calculation.
// cf. bitsofcotton/randtools extract main part:
//   Xor_k And_m Xor_n x_{k,m,n}*x_n == any operation on {0,1}^(dim x)
//   because of pattern matching also (a xor b) xor (a and b) == a or b.
//   so they're Sum_k det diag (X_k x) == det diag (Y x) in first digit.
//   this is done by counter diagonal method and LDLt:
//  integrate X_0 and X_1 : det diag x + det diag X' x (max rank is always
//   artificially created.) == det diag x + det diag (LDL^t x),
//   det diag L^-t x' + det diag LD x', in the x' =: [x'', 1, x''_reverse]
//   condition, one dimension down, repeat them causes det diag Y x.
//  also we can do analytical calculus on them causes det diag (Yx) ==
//   d/d(x_1) S det diag Yx d(x_1) == d/d(x_1) ... S ... <y,x> d(x_1) ... .
//   with repeat, we get <y,x>(x_1...)^m in first digit we get.
//  also with negated gate, x_1...\bar(x_1)... is constant (usually
//   x_k := {1,1/2}) so <y,x> in first digit part is what we need.
//  this context can be applied to p-adics so usually {1,(p-1)/p,...,1/p}
//   elements.
//  however, this concludes #f count up collision, so ind2vd makes better
//   dimension we need when it's observed and fixed.
//  so the condition might came from external R^3n to R^4n matrices.
//  also this is the analogy to toeplitz matrix inversion with singular one.
// N.B. if the function has internal states variable to be projected into
//      series, they're looked as <a,x>+<b,y>==<a,x>==0, y is internal states.
//      so this causes A*x==B*y, so increasing dimension causes ok result.
//      however, we're in invariant condition (de)?compression destroys,
//      so {x,y} in R^varlen is upper bound of variables however they causes
//      some matrix timing attacks.
// N.B. there's also trivial invariant : if((forall k, x_k==a_0k) or ...)
//      return 1; program. this is also in the condition but the dimension
//      easily vanished. so when we met them we use:
//      ||Ax-1*x'||<epsilon condition with increased varlen.
template <typename T> vector<T> p01nextM(const SimpleVector<T>& in) {
  const T zero(0);
  const T one(1);
  const T two(2);
  const T nin(sqrt(in.dot(in)));
  const int varlen(ind2vd(in.size()));
  if(! isfinite(nin) || nin == zero) {
    vector<T> res;
    res.resize(max(int(1), int(in.size()) - varlen + 1), zero);
    return res;
  }
  SimpleMatrix<T> invariants(max(int(1), int(in.size()) - varlen + 1),
    varlen + 2);
  invariants.O();
  for(int i0 = varlen * 2 + 1; i0 < invariants.rows(); i0 ++) {
    SimpleMatrix<T> toeplitz(i0, invariants.cols());
    for(int i = 0; i < toeplitz.rows(); i ++) {
      SimpleVector<T> work(in.subVector(i, varlen));
      work[work.size() - 1] = in[i + varlen - 1];
      toeplitz.row(i) = makeProgramInvariant<T>(work,
        T(i + 1) / T(toeplitz.rows() + 1) ).first;
    }
    // N.B. this untangles input stream into invariant but the accuracy
    //      we make the hypothesis:
    //      ||invariant made stream||_sup <~ ||f||_sup / varlen!.
    //      this is because it's toeplitz made stream.
    const int ii0(i0 - (varlen * 2 + 1));
    const SimpleVector<T> lt(linearInvariant<T>(toeplitz));
    invariants.row(ii0) = lt;
  }
  vector<T> res;
  res.reserve(invariants.rows());
  for(int j = 1; j <= invariants.rows(); j++) {
    SimpleVector<T> invariant(invariants.cols());
    if(invariants.rows() <= 1) invariant = move(invariants.row(0));
    else {
      // N.B. the QR decomposition mixes as average on 2 or latter index.
      //      so the continuity concerns 1st index also orthogonalizing
      //      causes one dimension theta continuity concerns.
      //      in raw meaning, taking the invariant causes ||p_0||<<epsilon
      //      however we should take ||p_0'||/||p_0|| for this transformation
      //      is valid or not. in weak differential meaning, it's sliding and
      //      add some continuity by higher frequency as a average.
      //      so in very roughly this includes any of the cases except
      //      for the condition invariant-next's linearInvariant fixes
      //      the whole invariant by orthogonalization case, but this
      //      condition is rarely satisfied by PRNG blended streams.
      invariant.O();
      for(int i = 0; i < invariants.cols(); i ++)
        invariant[i] = p0maxNext<T>(invariants.col(i).subVector(0, j));
    }
    SimpleVector<T> work(varlen);
    for(int i = 1; i < work.size(); i ++)
      work[i - 1] = in[i - work.size() + in.size()];
    work[work.size() - 1] = T(int(0));
    res.emplace_back(revertByProgramInvariant<T, true>(work, invariant));
  }
  // efi_cons_putc(0, '.');
  return res;
}

template <typename T> static inline T p01next(const SimpleVector<T>& in) {
  vector<T> res(p01nextM<T>(in));
  return res[res.size() - 1];
}

// N.B. class-capsules for serial stream.
template <typename T> class idFeeder {
public:
  inline idFeeder(const int& size = 0) {
    res.resize(size);
    res.O();
    full = size ? 0 : 2;
    t = 0;
  }
  inline ~idFeeder() { ; }
  inline const SimpleVector<T>& next(const T& in) {
    if(full == 2) {
      res.entity.emplace_back(in);
      ++ t;
      return res;
    }
    if(t < res.size())
      res[t] = in;
    else {
      for(int i = 1; i < res.size(); i ++)
        res[i - 1] = res[i];
      res[res.size() - 1] = in;
    }
    if(res.size() <= ++ t) full = 1;
    return res;
  }
  SimpleVector<T> res;
  char full;
  int  t;
};
 
template <typename T> static inline pair<SimpleVector<SimpleVector<T> >, T> normalizeS(const SimpleVector<SimpleVector<T> >& in) {
  pair<SimpleVector<SimpleVector<T> >, T> res;
  res.second = T(int(0));
  for(int i = 0; i < in.size(); i ++) for(int j = 0; j < in[i].size(); j ++)
    res.second = max(res.second, abs(in[i][j]));
  res.first = in;
  if(T(int(0)) < res.second)
    for(int i = 0; i < in.size(); i ++)
      for(int j = 0; j < in[i].size(); j ++)
        res.first[i][j] /= res.second;
  return res;
}

template <typename T> static inline vector<vector<SimpleMatrix<T> > > normalize(const vector<vector<SimpleMatrix<T> > >& data, const T& upper = T(1)) {
  T MM(0), mm(0);
  bool fixed(false);
  for(int kk = 0; kk < data.size(); kk ++)
    for(int k = 0; k < data[kk].size(); k ++)
      for(int i = 0; i < data[kk][k].rows(); i ++)
        for(int j = 0; j < data[kk][k].cols(); j ++)
          if(! fixed || (isfinite(data[kk][k](i, j)) &&
               ! isinf(data[kk][k](i, j)) && ! isnan(data[kk][k](i, j)))) {
            if(! fixed)
              MM = mm = data[kk][k](i, j);
            else {
              MM = max(MM, data[kk][k](i, j));
              mm = min(mm, data[kk][k](i, j));
            }
            fixed = true;
          }
  vector<vector<SimpleMatrix<T> > > result(data);
  if(MM == mm) {
    for(int kk = 0; kk < data.size(); kk ++)
      for(int k = 0; k < data[kk].size(); k ++)
        for(int i = 0; i < data[kk][k].rows(); i ++)
          for(int j = 0; j < data[kk][k].cols(); j ++)
            result[kk][k](i, j) = T(int(1)) / T(int(2));
    return result;
  } else if(! fixed)
    return result;
  for(int kk = 0; kk < data.size(); kk ++)
    for(int k = 0; k < data[kk].size(); k ++)
      for(int i = 0; i < data[kk][k].rows(); i ++)
        for(int j = 0; j < data[kk][k].cols(); j ++) {
          if(isfinite(result[kk][k](i, j)) && ! isinf(data[kk][k](i, j)) && ! isnan(result[kk][k](i, j)))
            result[kk][k](i, j) -= mm;
          else
            result[kk][k](i, j)  = T(0);
          result[kk][k](i, j) *= upper / (MM - mm);
        }
  return result;
}

template <typename T> static inline vector<SimpleMatrix<T> > normalize(vector<SimpleMatrix<T> >& data, const T& upper = T(1)) {
  vector<vector<SimpleMatrix<T> > > work;
  work.emplace_back(move(data));
  vector<SimpleMatrix<T> > res(normalize<T>(work, upper)[0]);
  data = move(work[0]);
  return res;
}

template <typename T> static inline vector<SimpleMatrix<T> > normalize(const vector<SimpleMatrix<T> >& in, const T& upper = T(1)) {
  vector<SimpleMatrix<T> > d(in);
  return normalize<T>(d, upper);
}

template <typename T> static inline SimpleMatrix<T> normalize(SimpleMatrix<T>& data, const T& upper = T(1)) {
  vector<SimpleMatrix<T> > work;
  work.emplace_back(move(data));
  SimpleMatrix<T> res(normalize<T>(work, upper)[0]);
  data = move(work[0]);
  return res;
}

template <typename T> static inline SimpleMatrix<T> normalize(const SimpleMatrix<T>& in, const T& upper = T(1)) {
  SimpleMatrix<T> d(in);
  return normalize<T>(d, upper);
}

#if defined(_SIMPLEALLOC_)
template <typename T> static inline vector<vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > >, SimpleAllocator<vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > > > > normalize(const vector<vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > >, SimpleAllocator<vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > > > >& in, const T& upper = T(1)) {
#else
template <typename T> static inline vector<vector<SimpleVector<T> > > normalize(const vector<vector<SimpleVector<T> > >& in, const T& upper = T(1)) {
#endif
  vector<vector<SimpleMatrix<T> > > w;
  w.resize(in.size());
  for(int i = 0; i < in.size(); i ++) {
    w[i].resize(in[i].size(), SimpleMatrix<T>(1, in[i][0].size()).O());
    for(int j = 0; j < in[i].size(); j ++)
      w[i][j].row(0) = in[i][j];
  }
  vector<vector<SimpleMatrix<T> > > res(normalize<T>(w, upper));
  w.resize(0);
  vector<vector<SimpleVector<T> > > v;
  v.resize(res.size());
  for(int i = 0; i < res.size(); i ++) {
    v[i].resize(res[i].size(), SimpleVector<T>(res[i][0].cols()).O());
    for(int j = 0; j < v[i].size(); j ++)
      v[i][j] = res[i][j].row(0);
  }
  return v;
}

#if defined(_SIMPLEALLOC_)
template <typename T> static inline vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > > normalize(const vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > >& in, const T& upper = T(1)) {
#else
template <typename T> static inline vector<SimpleVector<T> > normalize(const vector<SimpleVector<T> >& in, const T& upper = T(1)) {
#endif
  SimpleMatrix<T> w;
  w.resize(in.size(), in[0].size());
  w.entity = in;
  return normalize<T>(w, upper).entity;
}

template <typename T> static inline SimpleVector<T> normalize(const SimpleVector<T>& in, const T& upper = T(1)) {
  SimpleMatrix<T> w;
  w.resize(1, in.size());
  w.row(0) = in;
  return normalize<T>(w, upper).row(0);
}

template <typename T> static inline SimpleMatrix<T> rotate(const SimpleMatrix<T>& d, const T& theta) {
  const T c(cos(theta));
  const T s(sin(theta));
  const int h0(abs(int(c * T(d.rows()) - s * T(d.cols()))));
  const int h1(h0 + abs(int(s * T(d.cols()))) * 2);
  const int w0(abs(int(s * T(d.rows()) + c * T(d.cols()))));
  const int w1(w0 + abs(int(s * T(d.rows()))) * 2);
  SimpleMatrix<T> res(h0 < d.rows() ? h1 : h0,
                      w0 < d.cols() ? w1 : w0);
  const T offy(h0 < d.rows() ? abs(int(s * T(d.cols()))) : 0);
  const T offx(w0 < d.cols() ? abs(int(s * T(d.rows()))) : 0);
  res.O();
  const int diag(absceil(sqrt(T(res.rows() * res.rows() +
                             res.cols() * res.cols()) )) );
  for(int j = - diag; j < diag; j ++)
    for(int k = - diag; k < diag; k ++) {
      const int yy(c * T(j) - s * T(k) + offy);
      const int xx(s * T(j) + c * T(k) + offx);
      if(0 <= yy && yy < res.rows() &&
         0 <= xx && xx < res.cols()) {
        const int dyy(getImgPt<int>(j, d.rows()));
        const int dxx(getImgPt<int>(k, d.cols()));
        {
          res(yy, xx) = res(min(yy + 1, int(res.rows()) - 1), xx) =
            res(yy, min(xx + 1, int(res.cols()) - 1)) =
            res(min(yy + 1, int(res.rows()) - 1),
                min(xx + 1, int(res.cols()) - 1)) =
              d(dyy, dxx);
        }
      }
    }
  return res;
}

template <typename T> static inline SimpleMatrix<T> center(const SimpleMatrix<T>& dr, const SimpleMatrix<T>& d) {
  SimpleMatrix<T> res(d.rows(), d.cols());
  for(int i = 0; i < res.rows(); i ++)
    for(int j = 0; j < res.cols(); j ++)
      res(i, j) = dr(max(int(0), min(i + (dr.rows() - d.rows()) / 2, dr.rows() - 1)),
                     max(int(0), min(j + (dr.cols() - d.cols()) / 2, dr.cols() - 1)));
  return res;
}

template <typename T, bool useful> static inline SimpleVector<T> bitsG(const SimpleVector<T>& d, const int& b) {
  SimpleVector<T> res;
  if(b < 0) {
    res.resize(d.size() / abs(b));
    res.O();
    for(int i = 0; i < res.size(); i ++)
      for(int j = 0; j < abs(b); j ++)
        res[i] += (useful ? offsetHalf<T>(sgn<T>(unOffsetHalf<T>(
          d[i * abs(b) + j]))) : d[i * abs(b) + j])
            * pow(T(int(2)), - T(j + 1));
  } else {
    res.resize(d.size() * b);
    res.O();
    for(int i = 0; i < d.size(); i ++)
      for(int j = 0; j < b; j ++) {
        T shift(d[i] << j);
        shift -= absfloor(shift);
        res[i * b + j] = useful ? shift : T(int(shift * T(int(2)) )) / T(int(2));
      }
  }
  return res;
}

#if defined(_SIMPLEALLOC_)
template <typename T, bool useful> static inline vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > > bitsG(const vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > >& d, const int& b) {
  vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > > res(d);
#else
template <typename T, bool useful> static inline vector<SimpleVector<T> > bitsG(const vector<SimpleVector<T> >& d, const int& b) {
  vector<SimpleVector<T> > res(d);
#endif
  for(int i = 0; i < res.size(); i ++) res[i] = bitsG<T, useful>(res[i], b);
  return res;
}

// N.B. start ddpmopt
// N.B. predict with discrete pseudo Riemann-Stieljes condition.
#if defined(_SIMPLEALLOC_)
template <typename T, int nprogress> vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > > pRS(const vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > >& intran0, const string& strloop) {
#else
template <typename T, int nprogress> vector<SimpleVector<T> > pRS(const vector<SimpleVector<T> >& intran0, const string& strloop) {
#endif
  if(intran0[0].size() < 12) {
#if defined(_SIMPLEALLOC_)
    if(! intran0[0].size()) return vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > >();
    vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > > res;
#else
    if(! intran0[0].size()) return vector<SimpleVector<T> >();
    vector<SimpleVector<T> > res;
#endif
    res.resize(intran0[0].size(), SimpleVector<T>(intran0.size()).O());
    for(int j = 0; j < res.size(); j ++) {
      for(int i = 0; i < res[j].size(); i ++)
        res[j][i] = p0maxNext<T>(intran0[i].subVector(0, j + 1));
    }
    return res;
  }
  SimpleVector<T> seconds(intran0[0].size());
  seconds.O();
#if defined(_SIMPLEALLOC_)
  vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > > intran(intran0);
#else
  vector<SimpleVector<T> > intran(intran0);
#endif
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < intran0[0].size(); i ++) {
    idFeeder<T> work0(intran0.size());
    for(int j = 0; j < intran0.size(); j ++) work0.next(intran[j][i]);
    pair<SimpleVector<T>, T> work(makeProgramInvariant<T>(work0.res));
    for(int j = 0; j < intran.size(); j ++) intran[j][i] = move(work.first[j]);
    seconds[i] = move(work.second);
  }
  T M(int(0));
  for(int i = 0; i < intran.size(); i ++)
    for(int j = 0; j < intran[i].size(); j ++)
      M = max(abs(intran[i][j]), M);
  for(int i = 0; i < intran.size(); i ++) intran[i] /= M;
  seconds /= M;
#if defined(_SIMPLEALLOC_)
  vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > > p;
#else
  vector<SimpleVector<T> > p;
#endif
  
  // N.B. p01next calls p0maxNext implicitly, this needs to be cached single
  //      threaded process on first call.
  vector<T> p0(p01nextM<T>(intran[0]));
  p.resize(p0.size(), SimpleVector<T>(intran.size()).O());
  for(int i = 0; i < p0.size(); i ++) p[i][0] = p0[i];
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int j = 1; j < p[0].size(); j ++) {
    vector<T> pp(p01nextM<T>(intran[j]));
    for(int i = 0; i < pp.size(); i ++) p[i][j] = pp[i];
  }
  const T nseconds(sqrt(seconds.dot(seconds)) );
  const vector<T> peconds(p01nextM<T>(seconds / nseconds));
  for(int i = 0; i < p.size(); i ++)
    p[i] = unOffsetHalf<T>(revertProgramInvariant<T>(make_pair(
      makeProgramInvariant<T>(normalize<T>(p[i])).first, peconds[i] * nseconds)
        ).subVector(0, intran.size()) );
  return p;
}

// N.B. we add some Lebesgue part by cutting input horizontal.
#if defined(_SIMPLEALLOC_)
template <typename T, int nprogress> vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > > pLebesgue(const vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > >& in, const int& horizontal, const string& strloop) {
#else
template <typename T, int nprogress> vector<SimpleVector<T> > pLebesgue(const vector<SimpleVector<T> >& in, const int& horizontal, const string& strloop) {
#endif
  if(in.size() <= horizontal * horizontal)
#if defined(_SIMPLEALLOC_)
    return vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > >();
#else
    return vector<SimpleVector<T> >();
#endif
#if defined(_SIMPLEALLOC_)
  vector<vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > >, SimpleAllocator<vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > > > > reform;
#else
  vector<vector<SimpleVector<T> > > reform;
#endif
  reform.resize(horizontal);
  for(int i = 0; i < reform.size(); i ++) {
    reform[i].resize(in[0].size());
    for(int j = 0; j < reform[i].size(); j ++)
      reform[i][j].entity.reserve(in.size() - horizontal * horizontal + 1);
  }
  for(int i = 0; i <= in.size() - horizontal * horizontal; i ++) {
    vector<vector<vector<T> > > les;
    les.resize(in[0].size());
    for(int j = 0; j < les.size(); j ++) {
      les[j].resize(horizontal);
      for(int k = 0; k < les[j].size(); k ++)
        les[j][k].reserve(horizontal * horizontal);
    }
    for(int j = i; j < i + horizontal * horizontal; j ++)
      for(int k = 0; k < les.size(); k ++)
        les[k][int(binMargin<T>(in[j][k]) * T(horizontal) )].emplace_back(in[j][k]);
    for(int k = 0; k < in[0].size(); k ++) {
      int Mtot(0);
      for(int j = 0; j < horizontal; j ++)
        Mtot = max(Mtot, int(les[k][j].size()));
      for(int j = 0; j < horizontal; j ++) {
        T sum(int(0));
        for(int n = 0; n < les[k][j].size(); n ++) sum += les[k][j][n];
        reform[j][k].entity.emplace_back(binMargin<T>(sum / T(Mtot) *
          T(horizontal) / T(j + 1) ) );
      }
    }
  }
#if defined(_SIMPLEALLOC_)
  vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > > res;
#else
  vector<SimpleVector<T> > res;
#endif
  res.resize(in.size(), SimpleVector<T>(in[0].size()).O());
  for(int i = 0; i < reform.size(); i ++) {
    T n2(int(0));
    for(int j = 0; j < reform[i].size(); j ++)
      n2 += reform[i][j].dot(reform[i][j]);
    if(n2 == T(int(0))) continue;
    // N.B. try 3 of the context from single input stream possible enough
    //      then take arithmetic average.
#if defined(_SIMPLEALLOC_)
    vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > > p0;
    vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > > p1;
    vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > > p2;
#else
    vector<SimpleVector<T> > p0;
    vector<SimpleVector<T> > p1;
    vector<SimpleVector<T> > p2;
#endif
    p0 = pRS<T, nprogress>(reform[i], string(""));
    p1 = logscale<T>(pRS<T, nprogress>(expscale<T>(reform[i]), string("") ));
    p1 = expscale<T>(pRS<T, nprogress>(logscale<T>(reform[i]), string("") ));
    for(int i = 0; i < p0.size(); i ++) {
      p0[i] += p1[i];
      p0[i] += p2[i];
      p0[i] *= T(i + 1) / T(horizontal) / T(int(3));
    }
    if(p0.size() < res.size()) res.resize(p0.size());
    for(int i = 0; i < res.size(); i ++) res[i] += p0[i];
    // efi_cons_putc(0, 'L');
  }
  for(int i = 0; i < res.size(); i ++) res[i] /= T(int(reform.size()));
  return res;
}

// N.B. add some sectional measurement part.
#if defined(_SIMPLEALLOC_)
template <typename T, int nprogress> static inline vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > > pSectional(const vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > >& in, const string& strloop) {
#else
template <typename T, int nprogress> static inline vector<SimpleVector<T> > pSectional(const vector<SimpleVector<T> >& in, const string& strloop) {
#endif
  const int n(absfloor(sqrt(max(T(int(0)), log(T(in.size())) / log(T(int(2))) )) ));
  const int range(max(int(2), n));
  const int sectional(range * range);
#if defined(_SIMPLEALLOC_)
  if(! in.size()) return vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > >();
#else
  if(! in.size()) return vector<SimpleVector<T> >();
#endif
  if(in.size() <= sectional * 2)
    return unOffsetHalf<T>(pLebesgue<T, nprogress>(in, range, strloop));
#if defined(_SIMPLEALLOC_)
  vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > > res(pLebesgue<T, nprogress>(in, range, strloop) );
#else
  vector<SimpleVector<T> > res(pLebesgue<T, nprogress>(in, range, strloop) );
#endif
  for(int j = 0; j < res.size(); j ++) {
    res[j] *= T(sectional);
    for(int i = 1; i < sectional; i ++)
      res[j] -= in[in.size() - i + j - res.size() + 1];
  }
  return unOffsetHalf<T>(res);
}

// N.B. the result somehow not offsetted so we offset to 0.
#if defined(_SIMPLEALLOC_)
template <typename T, int nprogress> vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > > pPolish(const vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > >& in, const string& strloop) {
  vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > > resp;
  vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > > resm;
#else
template <typename T, int nprogress> vector<SimpleVector<T> > pPolish(const vector<SimpleVector<T> >& in, const string& strloop) {
  vector<SimpleVector<T> > resp;
  vector<SimpleVector<T> > resm;
#endif
  resp = pSectional<T, nprogress>(in,  string("+") + strloop);
  {
#if defined(_SIMPLEALLOC_)
    vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > > inm(in);
#else
    vector<SimpleVector<T> > inm(in);
#endif
    for(int i = 0; i < inm.size(); i ++)
      inm[i] = offsetHalf<T>(- unOffsetHalf<T>(inm[i]));
    resm = pSectional<T, nprogress>(inm, string("-") + strloop);
    for(int i = 0; i < resm.size(); i ++) resm[i] = - resm[i];
  }
  for(int i = 0; i < resp.size(); i ++) {
    resp[i] += resm[i];
    resp[i] /= T(int(2));
  }
  // efi_cons_putc(0, 'r');
  return resp;
}

#if ! defined(_P_BIT_)
// N.B. we don't get 100% result on each prediction, so upper bit broken case,
//      lower bits says nothing.
#define _P_BIT_ 3
#endif

// N.B. p01next output meaning only binary we make hypothesis.
#if defined(_SIMPLEALLOC_)
template <typename T, int nprogress> static inline vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > > pGuaranteeM(const SimpleVector<SimpleVector<T> >& in, const string& strloop) {
  vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > > res(pPolish<T, nprogress>(
#else
template <typename T, int nprogress> static inline vector<SimpleVector<T> > pGuaranteeM(const SimpleVector<SimpleVector<T> >& in, const string& strloop) {
  vector<SimpleVector<T> > res(pPolish<T, nprogress>(
#endif
    bitsG<T, true>(in.entity, abs(_P_BIT_)), strloop) );
  for(int i = 0; i < res.size(); i ++)
    res[i] = bitsG<T, true>(offsetHalf<T>(res[i]), - abs(_P_BIT_) );
  efi_cons_putc(0, 'g');
  return res;
}

template <typename T, int nprogress> static inline SimpleVector<T> pGuarantee(const SimpleVector<SimpleVector<T> >& in, const string& strloop) {
  vector<SimpleVector<T> > res(pGuaranteeM<T, nprogress>(in, strloop));
  return res[res.size() - 1];
}

#if ! defined(_P_MLEN_)
// N.B. avoid exhaust of calculation.
#define _P_MLEN_ 21
#endif

// N.B. append pseudo-measureable condition into original input
//      stream but the predictor isn't depend pseudo-things.
//      also add whole context length markov feeding.
#if defined(_SIMPLEALLOC_)
template <typename T, int nprogress> SimpleVector<T> pAppendMeasure(const vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > >& in, const string& strloop) {
#else
template <typename T, int nprogress> SimpleVector<T> pAppendMeasure(const vector<SimpleVector<T> >& in, const string& strloop) {
#endif
  const int realin(_P_MLEN_ ? min(int(in.size()), int(_P_MLEN_)) : int(in.size()) );
#if defined(_SIMPLEALLOC_)
  vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > > pp;
  vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > > pm;
  vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > > p;
  vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > > q;
  vector<SimpleVector<T>, SimpleAllocator<SimpleVector<T> > > r;
#else
  vector<SimpleVector<T> > pp;
  vector<SimpleVector<T> > pm;
  vector<SimpleVector<T> > p;
  vector<SimpleVector<T> > q;
  vector<SimpleVector<T> > r;
#endif
  {
    SimpleVector<SimpleVector<T> > workp;
    workp.entity.reserve(realin * 2 + 1);
    SimpleVector<T> b(in[0].size());
    b.O();
    for(int i = 0; i < realin; i ++) {
      SimpleVector<T> uo(unOffsetHalf<T>(in[i - realin + in.size()]));
      workp.entity.emplace_back(b);
      workp.entity.emplace_back(uo);
      b = uo * T(int(2)) - b;
    }
    workp.entity.emplace_back(b);
    workp.entity = delta<SimpleVector<T> >(delta<SimpleVector<T> >(workp.entity));
    pair<SimpleVector<SimpleVector<T> >, T> wp(normalizeS<T>(workp));
    pp = unOffsetHalf<T>(pGuaranteeM<T, nprogress>(offsetHalf<T>(
      wp.first), string("-)") + strloop));
    for(int i = 0; i < pp.size(); i ++) pp[i] *= wp.second;
  }
  {
    SimpleVector<SimpleVector<T> > workm;
    workm.entity.reserve(realin * 2 + 1);
    SimpleVector<T> b(in[0].size());
    b.O();
    for(int i = 0; i < realin; i ++) {
      SimpleVector<T> uo(unOffsetHalf<T>(in[i - realin + in.size()]));
      workm.entity.emplace_back(- b);
      workm.entity.emplace_back(uo);
      b = uo * T(int(2)) - b;
    }
    workm.entity.emplace_back(- b);
    workm.entity = delta<SimpleVector<T> >(delta<SimpleVector<T> >(workm.entity));
    pair<SimpleVector<SimpleVector<T> >, T> wm(normalizeS<T>(workm));
    pm = unOffsetHalf<T>(pGuaranteeM<T, nprogress>(offsetHalf<T>(
      wm.first), string("-)") + strloop));
    for(int i = 0; i < pm.size(); i ++) pm[i] *= wm.second;
  }
  p.reserve(pp.size());
  q.reserve(pp.size());
  for(int i = 0; i < pp.size(); i ++) {
    p.emplace_back((pp[i] + pm[i]) / T(int(2)));
    q.emplace_back((i ^ pp.size()) & 1 ?
      SimpleVector<T>(p[i].size()).O() :
      unOffsetHalf<T>(in[(i - int(pp.size())) / 2 + in.size()]) );
  }
  r.reserve(p.size() - 2);
  for(int i = 3; i <= p.size(); i ++) {
    r.emplace_back(SimpleVector<T>(p[i - 1]).O());
    for(int j = 0; j < p[i - 1].size(); j ++) {
      idFeeder<T> buf0(3);
      idFeeder<T> buf1(3);
      for(int k = 0; k < 3; k ++) buf0.next(q[k - 3 + i][j] - p[k - 3 + i][j]);
      for(int k = 0; k < 3; k ++) buf1.next(p[k - 3 + i][j]);
      r[i - 3][j] = p0maxNext<T>(buf0.res) + p0maxNext<T>(buf1.res);
    }
  }
  for(int i = 1; i < r.size(); i ++) r[i] += r[i - 1];
  for(int i = 1; i < r.size(); i ++) {
    r[0] += r[i];
/*
    if(((i ^ r.size()) & 1) && i < r.size() - 1) {
      int sum(int(0));
      for(int j = 0; j < r[0].size(); j ++)
        sum += int(sgn<T>(r[0][j] *
          unOffsetHalf<T>(in[(i - int(r.size())) / 2 + in.size()][j]) ));
      int stat(T(sum) / T(r[0].size()) * T(int(10000)) );
      printf("%d%c%d%c\n\0", stat / 100, '.', stat % 100, '\%');
    }
*/
  }
  return r[0];
}

// N.B. each pixel each bit prediction with PRNG blended stream.
//      this breaks whole image context as each pixel going to 1 a.s. vs
//      jammer result causes 0 a.s. on our machine, don't know why (should be
//      1 a.s. on whole pixel context for each).
#if defined(_P_PRNG_)
template <typename T, int nprogress> static inline SimpleVector<T> pEachPRNG(const vector<SimpleVector<T> >& in0, const string& strloop) {
  for(int i = 0; i < in0.size(); i ++)
    pnextcacher<T>(i + 1, 1);
  const vector<SimpleVector<T> > in(bitsG<T, false>(in0, abs(_P_BIT_)) );
  SimpleVector<T> out;
  out.resize(in[0].size());
  out.O();
  for(int i = 0; i < in[0].size(); i ++) {
    printf("%d, %d\n", last, lastptr);
    vector<SimpleVector<T> > work;
    work.resize(in.size());
    for(int j = 0; j < in.size(); j ++) {
      work[j].resize(_P_PRNG_);
      for(int k = 0; k < work[j].size(); k ++) work[j][k] = offsetHalf<T>(
#if defined(_ARCFOUR_)
        arc4random() & 1 ?
#else
        random() & 1 ?
#endif
        - unOffsetHalf<T>(in[j][i]) : unOffsetHalf<T>(in[j][i]) );
    }
    SimpleVector<T> w(pAppendMeasure<T, 0>(work, string(" ") +
      to_string(i) + string("/") + to_string(in[0].size()) + strloop) );
    for(int j = 0; j < w.size(); j ++) out[i] += w[j];
    int sign(0);
    for(int j = 0; j < work[0].size(); j ++) sign += 
#if defined(_ARCFOUR_)
        arc4random() & 1;
#else
        random() & 1;
#endif
    if(sign < work[0].size() / 2)
      out[i] = offsetHalf<T>(- unOffsetHalf<T>(out[i]));
  }
  return bitsG<T, true>(normalize<T>(out), - abs(_P_BIT_) );
}
#endif

// N.B. we make the first hypothesis as the stream is calculatable by *single*
//      function as n-markov also which measureability-appendant stream can
//      seep out the original stream from information amount reason.
// N.B. layers:
//       | function           | layer# | [wsp1] | data amount* | time*(***)   |
//       +-----------------------------------------------------+--------------+
//       | pEachPRNG                   | -1  | w | _P_PRNG_    | _P_PRNG_
//       | pAppendMeasure              | 0   | w | in          | in
//       | pGuarantee                  | 1   | w | _P_BITS_    | _P_BITS_
//       | pPolish                     | *   | w | 2           | 2
//       | pSectional                  | 2   | w |             |
//       | pLebesgue                   | 3   | w | range^2     | range
//       | divide by program invariant | 4+  | s | +unit       | +O(GL)
//       | burn invariant by p0next    | 5++ | s | +unit       | +O(GL+L^3)
//       | makeProgramInvariant        | 6+  | p |             | +O(GL)
//       | linearInvariant             | -   | - | -           |
//       |   - QR decomposition        | 8+* | s | > 4!        | O(GL)
//       |   - orthogonalization       | 10+*| p | > 4!        | O(4!)
//       |   - solve                   | 12* | p | +> (4 * 4)  | +O(4^3)
//       | T::operator *,/             | 13+ | 1 |             |
//       | T::operator +,-             | 14+ | 1 |             |
//       | T::bit operation            | 15+ | 1 |             |
// *(++) | sumCNext                    | +0  | s |             |
//       | sumCNext                    | +1  | s |             |
//       | logCNext                    | +2  | s |             |
//       | logCNext                    | +3  | s |             |
//       | northPoleNext               | +4  | s |             |
//       | invNext                     | +5  | s |             |
//       | sumCNext                    | +6  | s |             |
//       | pnext                       | +7  | s | +once(dft)  |
//       | integrate-diff in taylorc   | +8  | p | +once(dft)  |
//       | exp to shift   in taylorc   | +9  | p | +once(dft)  |
//       | dft                         | +10 | p | +once(dft)  |
//       | exp-log complex operation   | +11 | 1 | +once(taylor) |
//       | T::operator *,/             | +12 | 1 |             |
//       | T::operator +,-             | +13 | 1 |             |
//       | T::bit operation            | +14 | 1 |             |
// (***) time order ratio, L for input stream length, G for input vector size,
//       stand from arithmatic operators. ind2varlen isn't considered.
// N.B. we need O(L^2*G) calculation time whole.

// N.B. rewrote 2025/08/30, last update 2025/10/17:
// N.B. generally speaking, the raw input stream with high entropy predictor
//      dislikes to predict in general meanings because {ok,ng,invariant}
//      each 1/3 condition. however, there's at least 2 hole to the condition.
//      they're (i) attach input stream as we can treat well (ii) learn large
//      set of data stream and predict with internal of explicitly described
//      numerical stream form. we select (i) condition also worked well for now.
//      this is from the stream and predictor better tangled condition with
//      good temperature. however, the tangle condition can slide to run away
//      from some of the reason we don't know why but the tanglement on the
//      calculation space structure - numerical series abstract structure
//      slip. so the predictor condition is for now, not the ever or never.
// N.B. the codes eliminated from this header the reason why:
// (i)  difference, summation concerns they have the condition s.t.
//      the noise of the prediction stream can spread entire of the result
//      so we eliminated: PdeltaOnce, Ppersistent, Pprogression, (P0DFT) and so.
//      also they're mostly equivalent to input timing skip concerns.
//      this is because A_0 ... A_k B x_0 matrix described form.
//      also patternized xor-filter have the timing related conditions.
//      however, skip conditions are to avoid jammer things, so we eliminated
//      from core predictors.
//      this include input stream non-linear scaling (even in x-axis).
//      also distant step predictions aren't fight with skip concerns.
// (ii) non linear function transformations we target is only exp/log scale.
//      this is because d^e/dx^e == dx condition and f^-1(f(x)) == x condition.
//      cf. (arctan(logscale))-n times chain causes y=x into sigmoid-like graph.
// (iii)predict twice or more by one predictor often causes clear edge but the
//      gulf things. this also includes {ok,ng,invariant}'s invariant condition
//      retry.
// (iv) ad-hoc layer implementations also inspired by numerical test isn't
//      useful for generic predictor because it's only ad-hoc to specific
//      numerical series. also after burner is.
// (v)  we don't need LoEM unstable case implementation on input stream attached
//      case because it's verbose.
// N.B. something XXX result descripton
// (00) there might exist non Lebesgue measureable condition discrete stream.
//      this is: there's no unique function on the range but AFTER all the
//      data is treated (observed), this condition never satisfied.
//      so this is the which is the latter chase. either if original stream
//      is something attached, this condition can be avoided well.
// (01) the prediction fail is come from first continuity hypothesis
//      satisfied or not. AFTER the whole stream is given context,
//      we can avoid such of the conditions with certain error.
// (02) (de)?compression concerns can jam out on N calculation matters.
//      we cannot avoid this other than verifying after the phenomenon
//      also having a verifiability of low of excluded middle based on
//      our calculation based on our conscious uniqueness.
// (03) might have once coded as obs. concerns. when we implement binary
//      they means we select one of the #f causes the jammer can jam out
//      our invariant condition. so if there's universal invariant,
//      once jammer targets us, they slip to non universal ones.
// (04) so the universal invariant condition needs to be hide from attacker
//      the binary tree or method itself to continue their effects really
//      grip on them. otherwise, we should use such a invariant from
//      the things we really trust from bottom of our hearts but this needs
//      a priori description on the stream however there exists the jammer
//      for any of the predictor, the description seems unfavorable.
// N.B. tips around jammer
// (-1) any of the predictor they have a jammer to them.
// (00) after of all, the dynamic jammer can be avoided if the predictor entropy
//      exceeds jammers one, so some of the first short range, the predictor
//      exceeds the jammer somehow.
// (01) however, the predictor entropy can be counted by program binary size
//      in some of the layer, so graphics predictor seems to have the quantity
//      so. either, once algorithm is coded as exist, they have upper entropy
//      size n bit-input, n bit-output, n^3 bit as a optimization result.
//      instead of the fixation of code optimization, we use optimization result
//      to get orthogonal to input stream condition or pivot to get high
//      frequency result.
// (02) the predictor vs. jammer made stream concludes the saturated input.
//      also the condition is the which side bore first chase.
//      if the saturated result we get, we should separate something on input.
// (03) there can be 0 invariant chain, so they can be caused by move average
//      they caused return to average works very well.
//      this is because <x,a> == 0, <[x,x+],[a,0]+[0,a]> == 0 chain in rough.
// (04) after some conversation with gemini around 2025/07, the jammers
//      they have internal optimization calculation to output to saturate
//      the stream intent condition or so also have the internal made intentions
//      to jam out the stream.
// N.B. there's plenty of the room to implement the predictor which is
//      saturating F_2^4 #f, the bra, ket condition indirect access.
//      either Riemann-measureable conditions' smaller than |dx| counting
//      causes whole function reference.
// N.B. things not implemented nor tested but abandoned.
// (01) brute force change state/output functions on (de)?compressed stream.
//      they are equivalent to p01next, p012next partially also we cannot
//      test because of their size on the memory.
//      the brute force condition isn't mean directly the things finding
//      pseudo-patternized ones because of (de)compression condition they
//      breaks LoEM causes some resonance on the stream.
// (02) predictor which shirking many much of the continuous functions:
//        this is with taking multiplication invariant on f,
//        S f(x) dx = S det(J((1,g0,...)/(1,x0,...)) dx0 ...
//        retaking their addition invariant as det(...) == 0, the given function
//        g0 ... should fit them also they describes much of continuities.
//        this can flatten N when our N is something infected.
//        also this is the analogy {1,x,x^2,...} on p0next meaning.
//      so the ongoing proceding machine learnings finds such a orthogonal
//      egg functions from input streams without the hypotehsis on PDE
//      structures. so we drop to implement them.
// (03) untangle by DFT or Wavelet triple. this is because R^R untangle
//      one by one causes Wavelet(Wavelet(Fourier+Discrete)+Discrete)+Discrete
//      causes only a combination ordinal, we need Discrete part separation
//      other than dft/mWavelet in fact. however, we drop this implementation.
#define _SIMPLELIN_
#endif

