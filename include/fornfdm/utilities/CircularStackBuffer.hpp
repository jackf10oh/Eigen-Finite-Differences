// CircularStackBuffer.hpp
//
// Circular buffer of a fixed size on the stack 
//
// JAF 5/2/2026 

#ifndef FORNFDM_UTILS_CIRCULARSTACKBUFFER_H
#define FORNFDM_UTILS_CIRCULARSTACKBUFFER_H 

#include<algorithm>
#include<array>
#include<cstdint>
#include<cstddef>
#include<iterator>
#include<stdexcept>
#include<type_traits>
#include<utility>
#include<initializer_list>
#include<memory> // std::addressof

namespace fornfdm{
namespace utils{

template<class T, std::size_t N>
class CircularStackBuffer
{
  static_assert(N>0, "size 0 CircularStackBuffer is UB."); 
  private:
    // Implemenations -------------------------------------------------- 
    template<bool isConst>
    class IteratorImpl
    { 
      // Friends -----------
      template<bool> friend class IteratorImpl;

      public:
        // Type Defs --------------- 
        using iterator_concept = std::random_access_iterator_tag;
        using iterator_category = std::random_access_iterator_tag; 
        using value_type = T;  
        using reference = typename std::conditional<isConst, const T&, T&>::type;
        using pointer = typename std::conditional<isConst, const T*, T*>::type; 
        using difference_type = std::ptrdiff_t; 

      private:
        // Member Data ------------ 
        pointer m_ptr; 
        difference_type m_offset; 

      public:
        // Constructor ----------------
        IteratorImpl(pointer ptr, difference_type offset)
          : m_ptr(ptr), m_offset(offset)
        {}

        IteratorImpl(const IteratorImpl& other)=default;

        // const from non const
        template<bool B = isConst, typename = std::enable_if_t<B>>
        IteratorImpl(const IteratorImpl<false>& other)
          : m_ptr(other.m_ptr), m_offset(other.m_offset)
        {}

        // Member Functions ----------- 
        IteratorImpl& operator++(){++m_offset; return *this; }
        IteratorImpl operator++(int){auto tmp = *this; ++m_offset; return tmp; }
        IteratorImpl& operator+=(difference_type n){m_offset+=n; return *this; }
        IteratorImpl operator+(difference_type n) const { return IteratorImpl(m_ptr, m_offset + n); }
        friend IteratorImpl operator+(difference_type n, IteratorImpl it){ return IteratorImpl(it.m_ptr, it.m_offset + n);}

        IteratorImpl& operator--(){--m_offset; return *this; }
        IteratorImpl operator--(int){auto tmp = *this; --m_offset; return tmp; }
        IteratorImpl& operator-=(difference_type n){m_offset-=n; return *this; }
        IteratorImpl operator-(difference_type n) const { return IteratorImpl(m_ptr, m_offset - n); }

        reference operator[](difference_type n) const 
        {
          auto idx = (m_offset+n) % static_cast<difference_type>(N); 
          if(idx<0) idx += static_cast<difference_type>(N); 
          return m_ptr[idx];  
        }

        reference operator*() const 
        { 
          auto idx = m_offset % static_cast<difference_type>(N); 
          if(idx<0) idx += static_cast<difference_type>(N); 
          return m_ptr[idx]; 
        }

        pointer operator->() const
        {
          return std::addressof(operator*());
        }

        difference_type operator-(const IteratorImpl& other) const { return m_offset - other.m_offset; }

        bool operator<(const IteratorImpl& other) const { return m_offset<other.m_offset; }
        bool operator<=(const IteratorImpl& other) const { return m_offset<=other.m_offset; }
        bool operator>(const IteratorImpl& other) const { return m_offset>other.m_offset; }
        bool operator>=(const IteratorImpl& other) const { return m_offset>=other.m_offset; }
        bool operator==(const IteratorImpl& other) const { return m_offset==other.m_offset; }
        bool operator!=(const IteratorImpl& other) const { return m_offset!=other.m_offset; }
    };

  public:
    // Type Defs -------------------------------------------
      using value_type = T; 
      using reference = T&; 
      using const_reference = const T&; 
      using difference_type = std::ptrdiff_t; 
      using size_type = std::size_t; 
      using pointer = T*; 
      using const_pointer = const T*; 
      using const_iterator = IteratorImpl<true>; 
      using iterator = IteratorImpl<false>;
      using reverse_iterator = std::reverse_iterator<iterator>;  
      using const_reverse_iterator = std::reverse_iterator<const_iterator>;

  private:  
    // Member Data ------------------------------------------
    size_type m_head;
    std::array<T,N> m_data;  

  public:
    // Constructors + Destructor =============================
    
    // default 
    CircularStackBuffer()
      : m_head(0), m_data()
    {}; 

    // from variadic args 
    template<typename... Args>
    CircularStackBuffer(Args&&... args)
      : m_data{ std::forward<Args>(args)... }, m_head(sizeof...(args))
    {
      static_assert(sizeof...(args) <= N, "must construct from <= N elements");
    }

    // copy 
    CircularStackBuffer(const CircularStackBuffer& other)=default; 

    // move
    CircularStackBuffer(CircularStackBuffer&& other)=default;

    // TODO copy/move from difference size, convertible U to T

    // destructor 
    ~CircularStackBuffer()=default; 

    // Member Functions -------------------------------------- 
    reference operator[](size_type idx) noexcept
    {
        size_type i = m_head + idx;
        if (i >= N) i -= N;
        return m_data[i];
    }
    const_reference operator[](size_type idx) const noexcept
    {
        size_type i = m_head + idx;
        if (i >= N) i -= N;
        return m_data[i];
    }
    reference at(size_type idx)
    { 
      if (idx >= N) throw std::out_of_range("CircularStackBuffer::at: index out of range");
      size_type i = m_head + idx;
      if (i >= N) i -= N;
      return m_data[i];
    }
    const_reference at(size_type idx) const
    { 
      if (idx >= N) throw std::out_of_range("CircularStackBuffer::at: index out of range");
      size_type i = m_head + idx;
      if (i >= N) i -= N;
      return m_data[i];
    }
    reference front() noexcept { return m_data[m_head]; }
    const_reference front() const noexcept { return m_data[m_head]; }
    reference back() noexcept { return m_data[(m_head==0) ? N-1 : m_head-1]; }
    const_reference back() const noexcept { return m_data[(m_head==0) ? N-1 : m_head-1]; }
    constexpr size_type size() const noexcept { return N; } 
    constexpr size_type max_size() const noexcept { return N; } 
    constexpr bool empty() const noexcept {return false; /*N!=0 asserted earlier */}    
    iterator begin() noexcept { return iterator(m_data.data(), static_cast<difference_type>(m_head)); }
    iterator end() noexcept { return iterator(m_data.data(), static_cast<difference_type>(m_head+N)); }
    const_iterator begin() const noexcept { return cbegin(); }
    const_iterator end() const noexcept { return cend(); }
    const_iterator cbegin() const noexcept { return const_iterator(m_data.data(), static_cast<difference_type>(m_head)); }
    const_iterator cend() const noexcept { return const_iterator(m_data.data(), static_cast<difference_type>(m_head+N)); }
    reverse_iterator rbegin() noexcept { return std::make_reverse_iterator(end());}      
    const_reverse_iterator rbegin() const noexcept { return crbegin();}      
    const_reverse_iterator crbegin() const noexcept { return std::make_reverse_iterator(cend());}      
    reverse_iterator rend() noexcept { return std::make_reverse_iterator(begin());}  
    const_reverse_iterator rend() const noexcept { return crend(); }  
    const_reverse_iterator crend() const noexcept { return std::make_reverse_iterator(cbegin());}  
    void swap(CircularStackBuffer& other) noexcept(noexcept(std::swap(std::declval<T&>(), std::declval<T&>()))) { std::swap(m_head, other.m_head); m_data.swap(other.m_data);}
    pointer data() noexcept { return m_data.data(); }
    const_pointer data() const noexcept { return m_data.data(); }
    void fill(const T& val){ m_data.fill(val); m_head=0; }
    void clear(){ m_data.fill(T()); m_head=0; }

    // shifts everything so that m_head is back at 0. 
    size_type offset() const noexcept { return m_head; }
    void linearize()
    {
      std::rotate(m_data.begin(), std::next(m_data.begin(), m_head), m_data.end()); 
      m_head = 0; 
    }  
    void rotate(difference_type r=1) noexcept
    {
      if constexpr(N <= 1) return;
      difference_type n = static_cast<difference_type>(N); 
      difference_type mod = ((static_cast<difference_type>(m_head)+r) % n);
      m_head = (mod >= 0) ? mod : (mod + n); 
    }
    void push_back(const T& val) noexcept(noexcept(std::declval<T&>() = std::declval<const T&>()))
    {
      if constexpr(N==1){
        back() = val; 
      }
      else{
        rotate(1); 
        back() = val; 
      }
    }
    void push_back(T&& val) noexcept(noexcept(std::declval<T&>() = std::declval<T&&>()))
    {
      if constexpr(N==1){
        back() = std::move(val); 
      }
      else{
        rotate(1); 
        back() = std::move(val); 
      }
    }
    template<class... Args>
    void emplace_back(Args&&... args)
    { 
      if constexpr(N==1){
        back() = T(std::forward<Args>(args)...);
      }
      else{
        rotate(1); 
        back() = T(std::forward<Args>(args)...); 
      }
    }

    // Operators ----------------------- 
    // lval assignment
    CircularStackBuffer& operator=(const CircularStackBuffer& other)=default;

    // rval assignment
    CircularStackBuffer& operator=(CircularStackBuffer&& other)=default;

    // TODO copy/move from difference size, convertible U to T
}; 

} // end namespace utils 
} // end namespace fornfdm 

#endif // CircularStackBuffer.hpp 
