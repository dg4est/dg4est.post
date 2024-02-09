/**
 * \file    memory.hxx
 * \author  akirby
 */

#ifndef MEMORY_HXX
#define MEMORY_HXX

/* system header files */
#include <memory>

/* library header files */
#include <occa.hpp>

namespace dg {

using properties_t = occa::json;
using device_t = occa::device;
using kernel_t = occa::kernel;
using stream_t = occa::stream;

template<typename T>
class memory {
  template <typename U> friend class memory;

  private:
    using size_t = std::size_t;
    using ptrdiff_t = std::ptrdiff_t;

    std::shared_ptr<T[]> shrdPtr;
    size_t lngth;
    size_t offset;
    bool allocated=false;

  public:
    memory() :
        shrdPtr{nullptr},
        lngth{0},
        offset{0},
        allocated{false} {}

    memory(const size_t lngth_) :
        shrdPtr(new T[lngth_]),
        lngth{lngth_},
        offset{0},
        allocated{true} {}

    memory(const size_t lngth_,
           const T val) :
        shrdPtr(new T[lngth_]),
        lngth{lngth_},
        offset{0},
        allocated{true} {
        for (size_t i=0;i<lngth;++i) {
          shrdPtr[i] = val;
        }
    }

    /* conversion constructor */
    template<typename U>
    memory(const memory<U> &m):
        shrdPtr{std::reinterpret_pointer_cast<T[]>(m.shrdPtr)},
        lngth{m.lngth*sizeof(T)/sizeof(U)},
        offset{m.offset*sizeof(T)/sizeof(U)},
        allocated{true} {
    }

    memory(const memory<T> &m)=default;
    memory& operator = (const memory<T> &m)=default;
   ~memory()=default;

    void malloc(const size_t lngth_) {
        *this = memory<T>(lngth_);
    }

    void malloc(const size_t lngth_, const T val) {
        *this = memory<T>(lngth_, val);
    }

    int newmalloc(const size_t lngth_) {
        if(!allocated) {*this = memory<T>(lngth_); return 1;}
        return 0;
    }

    int newmalloc(const size_t lngth_, const T val) {
        if(!allocated) {*this = memory<T>(lngth_, val); return 1;}
        return 0;
    }

    void calloc(const size_t lngth_) {
        *this = memory<T>(lngth_, T{0});
    }

    void realloc(const size_t lngth_) {
        memory<T> m(lngth_);
        const ptrdiff_t cnt = std::min(lngth, lngth_);
        m.copyFrom(*this, cnt);
        *this = m;
    }

    memory& swap(memory<T> &m) {
        std::swap(shrdPtr, m.shrdPtr);
        std::swap(lngth, m.lngth);
        std::swap(offset, m.offset);
        return *this;
    }

    void fill(const T val) {
        for (size_t i=0;i<lngth;++i) {
          shrdPtr[i] = val;
        }
    }

    T* ptr() {
        if (shrdPtr != nullptr) {
            return shrdPtr.get()+offset;
        } else {
            return nullptr;
        }
    }
    const T* ptr() const {
        if (shrdPtr != nullptr) {
            return shrdPtr.get()+offset;
        } else {
            return nullptr;
        }
    }

    T* begin() {return ptr();}
    T* end() {return ptr() + length();}

    size_t length() const {
        return lngth;
    }

    size_t size() const {
        return lngth*sizeof(T);
    }

    size_t use_count() const {
        return shrdPtr.use_count();
    }

    T& operator[](const ptrdiff_t idx) const {
        return shrdPtr[idx+offset];
    }

    bool operator == (const memory<T> &other) const {
        return (shrdPtr==other.shrdPtr && offset==other.offset);
    }

    bool operator != (const memory<T> &other) const {
        return (shrdPtr!=other.shrdPtr || offset!=other.offset);
    }

    memory<T> operator + (const ptrdiff_t offset_) const {
        return slice(offset_);
    }

    memory<T>& operator += (const ptrdiff_t offset_) {
        *this = slice(offset_);
        return *this;
    }

    memory<T> slice(const ptrdiff_t offset_,
                    const ptrdiff_t count = -1) const {
        memory<T> m(*this);
        m.offset = offset + offset_;
        m.lngth = (count==-1)
                ? (lngth - offset_)
                : count;
        return m;
    }

    /* copy from raw pointer */
    void copyFrom(const T* src,
                  const ptrdiff_t count = -1,
                  const ptrdiff_t offset_ = 0) {

        const ptrdiff_t cnt = (count==-1) ? lngth : count;

        std::copy(src,
                  src+cnt,
                  ptr()+offset_);
    }

    /* copy from memory */
    void copyFrom(const memory<T> src,
                  const ptrdiff_t count = -1,
                  const ptrdiff_t offset_ = 0) {
        const ptrdiff_t cnt = (count==-1) ? lngth : count;

        std::copy(src.ptr(),
                  src.ptr()+cnt,
                  ptr()+offset_);
    }

    /* copy to raw pointer */
    void copyTo(T *dest,
                const ptrdiff_t count = -1,
                const ptrdiff_t offset_ = 0) const {
        const ptrdiff_t cnt = (count==-1) ? lngth : count;

        std::copy(ptr()+offset_,
                  ptr()+offset_+cnt,
                  dest);
    }

    /* copy to memory */
    void copyTo(memory<T> dest,
                const ptrdiff_t count = -1,
                const ptrdiff_t offset_ = 0) const {
        const ptrdiff_t cnt = (count==-1) ? lngth : count;

        std::copy(ptr()+offset_,
                  ptr()+offset_+cnt,
                  dest.ptr());
    }

    memory<T> clone() const {
        memory<T> m(lngth);
        m.copyFrom(*this);
        return m;
    }

    void free() {
        shrdPtr = nullptr;
        lngth=0;
        offset=0;
        allocated=false;
    }
};

/* extern declare common instantiations for faster compilation */
template class memory<int>;
template class memory<float>;
template class memory<double>;

/* dg::deviceMemory is a wrapper around occa::memory */
template<typename T>
class deviceMemory: public occa::memory {
  public:
    deviceMemory() = default;
    deviceMemory(const deviceMemory<T> &m)=default;
    deviceMemory(occa::memory m) :
        occa::memory(m)
    {
        if (isInitialized()) {
            if (occa::dtype::get<T>() == occa::dtype::none) {
                occa::memory::setDtype(occa::dtype::byte);
            } else {
                occa::memory::setDtype(occa::dtype::get<T>());
            }
        }
    }

    /* conversion constructor */
    template<typename U>
    deviceMemory(const deviceMemory<U> &m) :
        occa::memory(m)
    {
        if (isInitialized()) {
            if (occa::dtype::get<T>() == occa::dtype::none) {
                occa::memory::setDtype(occa::dtype::byte);
            } else {
                occa::memory::setDtype(occa::dtype::get<T>());
            }
        }
    }

    deviceMemory<T>& operator = (const deviceMemory<T> &m)=default;
   ~deviceMemory()=default;

    T* ptr() {
        return static_cast<T*>(occa::memory::ptr());
    }
    const T* ptr() const {
        return static_cast<const T*>(occa::memory::ptr());
    }

    size_t length() const {
        return size()/sizeof(T);
    }

    T& operator[](const ptrdiff_t idx) {
        return ptr()[idx];
    }

    deviceMemory<T> operator + (const ptrdiff_t offset) const {
        if (isInitialized()) {
            return deviceMemory<T>(occa::memory::operator+(offset));
        } else {
            return deviceMemory<T>();
        }
    }

    deviceMemory<T>& operator += (const ptrdiff_t offset) {
        *this = deviceMemory<T>(occa::memory::slice(offset));
        return *this;
    }

    /* copy from dg::memory */
    void copyFrom(const dg::memory<T> src,
                  const ptrdiff_t count = -1,
                  const ptrdiff_t offset = 0,
                  const properties_t &props = properties_t()) {
        const ptrdiff_t cnt = (count==-1) ? length() : count;
        if(cnt==0) return;

        occa::memory::copyFrom(src.ptr(),
                               cnt*sizeof(T),
                               offset*sizeof(T),
                               props);
    }

    void copyFrom(const dg::memory<T> src,
                  const properties_t &props) {
        if(length()==0) return;

        occa::memory::copyFrom(src.ptr(),
                               length()*sizeof(T),
                               0,
                               props);
    }

    /* copy from dg::deviceMemory */
    void copyFrom(const deviceMemory<T> src,
                  const ptrdiff_t count = -1,
                  const ptrdiff_t offset = 0,
                  const properties_t &props = properties_t()) {
        const ptrdiff_t cnt = (count==-1) ? length() : count;
        if(cnt==0) return;

        occa::memory::copyFrom(src,
                               cnt*sizeof(T),
                               offset*sizeof(T),
                               0,
                               props);
    }

    void copyFrom(const deviceMemory<T> src,
                  const properties_t &props) {
        if(length()==0) return;

        occa::memory::copyFrom(src,
                               length()*sizeof(T),
                               0,
                               0,
                               props);
    }

    /* copy to dg::memory */
    void copyTo(dg::memory<T> dest,
                const ptrdiff_t count = -1,
                const ptrdiff_t offset = 0,
                const properties_t &props = properties_t()) const {
        const ptrdiff_t cnt = (count==-1) ? length() : count;
        if(cnt==0) return;

        occa::memory::copyTo(dest.ptr(),
                             cnt*sizeof(T),
                             offset*sizeof(T),
                             props);
    }

    void copyTo(dg::memory<T> dest,
                const properties_t &props) const {
        if(length()==0) return;

        occa::memory::copyTo(dest.ptr(),
                             length()*sizeof(T),
                             0,
                             props);
    }

    /* copy to dg::deviceMemory */
    void copyTo(deviceMemory<T> dest,
                const ptrdiff_t count = -1,
                const ptrdiff_t offset = 0,
                const properties_t &props = properties_t()) const {
        const ptrdiff_t cnt = (count==-1) ? length() : count;
        if(cnt==0) return;

        occa::memory::copyTo(dest,
                             cnt*sizeof(T),
                             0,
                             offset*sizeof(T),
                             props);
    }

    void copyTo(deviceMemory<T> dest,
                const properties_t &props) const {
        if(length()==0) return;

        occa::memory::copyTo(dest,
                             length()*sizeof(T),
                             0,
                             0,
                             props);
    }
};

/* declare common instantiations for faster compilation */
template class deviceMemory<int>;
template class deviceMemory<float>;
template class deviceMemory<double>;

} //namespace dg

#endif /* MEMORY_HXX */
