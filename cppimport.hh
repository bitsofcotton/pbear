#if !defined(_SIMPLE_VECTOR_)

#include <new>

template <typename T> T&& move(T& x) { return static_cast<T&&>(x); }

template <typename T> class SimpleAllocator;

template <typename T, class _Allocator = SimpleAllocator<T> > class vector {
private:
  inline void initmem() { sz ^= sz; sz_alloc ^= sz_alloc; entity = reinterpret_cast<T*>(size_t(0)); }
public:
  inline vector<T, _Allocator>() { initmem(); }
  inline vector<T, _Allocator>(const int s) { initmem(); resize(s); }
  inline vector<T, _Allocator>(const vector<T, _Allocator>& x) { initmem(); *this = x; }
  inline vector<T, _Allocator>(vector<T, _Allocator>&& x) { initmem(); *this = x; }
  inline ~vector<T, _Allocator>() { if(sz_alloc || reinterpret_cast<size_t>(entity) ) alloc.deallocate(entity, sz_alloc); initmem(); }
  inline void reserve(const size_t& nsz) {
    if(nsz < sz) return;
    if(! entity || sz_alloc != nsz) {
      T* nentity(alloc.allocate(nsz));
      for(int i = 0; i < nsz; i ++) ::new ((void*)(&nentity[i])) T();
      // for(char* ptr(reinterpret_cast<char*>(nentity)); ptr < reinterpret_cast<char*>(&nentity[nsz + 1]); ++ ptr) *ptr = 0;
      if(entity) {
        for(int i = 0; i < sz; i ++) nentity[i] = entity[i];
        alloc.deallocate(entity, sz_alloc);
      }
      entity = nentity;
    }
    sz_alloc = nsz;
  }
  inline void resize(const size_t& nsz, const T& obj = T()) { 
    if(sz_alloc < nsz) reserve(nsz);
    if(sz < nsz) for(int i = sz; i < nsz; i ++) entity[i] = obj;
    else for(int i = nsz; i < sz; i ++) entity[i] = obj;
    sz = nsz;
  }
  inline void emplace_back(const T& obj) { T nobj(obj); emplace_back(move(nobj)); }
  inline void emplace_back(T&& obj) {
    if(! sz_alloc) reserve(20);
    if(sz_alloc <= sz + 1) reserve(sz_alloc * 2);
    entity[sz ++] = obj;
  }
  inline vector<T, _Allocator>& operator = (const vector<T, _Allocator>& x) {
    if(size() != x.size()) resize(x.size());
    for(int i = 0; i < x.size(); i ++) entity[i] = x.entity[i];
    return *this;
  }
  inline vector<T, _Allocator>& operator = (vector<T, _Allocator>&& x) {
    if(sz_alloc || reinterpret_cast<size_t>(entity)) alloc.deallocate(entity, sz_alloc);
    sz       = x.sz;
    sz_alloc = x.sz_alloc;
    entity   = x.entity;
    alloc    = x.alloc;
    x.sz    ^= x.sz;
    x.sz_alloc ^= x.sz_alloc;
    x.entity    = reinterpret_cast<T*>(size_t(0));
    return *this;
  }
  inline const T& operator [] (const size_t& idx) const { return entity[idx]; }
  inline T& operator [] (const size_t& idx) { return entity[idx]; }
  inline size_t size() const { return sz; }
  size_t sz;
  size_t sz_alloc;
  T*  entity;
  _Allocator alloc;
};

template <typename T, typename U> class pair {
public:
  T first;
  U second;
};

template <typename T, typename U> static inline pair<T, U> make_pair(const T& f, const U& s) {
  pair<T, U> res;
  res.first = f;
  res.second = s;
  return res;
}

template <typename T> static inline T abs(const T& x) { return x < T(int(0)) ? - x : x; }

template <typename T> static inline void swap(T& x, T& y) {
  const T xx(y);
  x = y; y = xx;
  return;
}

// stub:
class string {
public:
  inline string() { ; }
  inline string(const char* x) { ; }
  inline string operator + (const string& x) { return string(); }
};

// stub:
static inline string to_string(const int& x) { return string(); }
static inline string to_string(const size_t& x) { return string(); }

// N.B. thanks to musl-1.2.3/src/prng/rand.c
inline int random() {
  static unsigned long long seed(1234);
  seed = 6364136223846793005ULL * seed + 1;
  return int(seed >> 33);
}

#define _SIMPLE_VECTOR_
#endif

