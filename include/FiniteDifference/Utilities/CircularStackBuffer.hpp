// CircularStackBuffer.hpp
//
// Circular buffer of a fixed size on the stack 
//
// JAF 5/2/2026 

#ifndef FDM_UTILS_CIRCULARSTACKBUFFER_H
#define FDM_UTILS_CIRCULARSTACKBUFFER_H 

namespace fdm{
namespace utils{

template<class T, std::size_t N>
class CircularStackBuffer
{
  private:
    // Implemenations -------------------------------------------------- 
    template<bool isConst>
    class IteratorImpl
    { 
      public:
        // Type Defs --------------- 
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

        // Member Functions ----------- 
        IteratorImpl& operator++(){++m_offset; return *this; }
        IteratorImpl operator++(int){auto tmp = *this; ++m_offset; return tmp; }
        IteratorImpl& operator+=(difference_type n){m_offset+=n; return *this; }
        IteratorImpl operator+(difference_type n) const { return IteratorImpl(m_ptr, m_offset + n); }

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

  private:  
    // Member Data ------------------------------------------
    std::size_t m_head;
    std::array<T,N> m_data;  

  public:
    // Constructors + Destructor =============================
    
    // default 
    CircularStackBuffer()=default; 

    // from variadic args 
    template<class... Args>
    CircularStackBuffer(Args&&... args)
      : m_head(sizeof...(args) % N), m_data{std::forward<Args>(args)...}
    {} 

    // copy 
    CircularStackBuffer(const CircularStackBuffer& other)=default; 

    // destructor 
    ~CircularStackBuffer()=default; 

    // Member Functions -------------------------------------- 
    reference operator[](size_type idx){ return m_data[(m_head+idx)%N]; }
    const_reference operator[](size_type idx) const { return m_data[(m_head+idx)%N]; }
    reference at(size_type idx){ return m_data.at((m_head+idx)%N); }
    const_reference at(size_type idx) const { return m_data.at((m_head+idx)%N); }
    reference front(){ return m_data[m_head]; }
    const_reference front() const { return m_data[m_head]; }
    reference back(){ return m_data[(m_head==0) ? N-1 : m_head-1]; }
    const_reference back() const { return m_data[(m_head==0) ? N-1 : m_head-1]; }
    constexpr size_type size() const noexcept { return N; } 
    constexpr size_type max_size() const noexcept { return N; } 
    constexpr bool empty() const noexcept {return N==0;}    
    iterator begin(){ return iterator(m_data.data(), m_head); }
    iterator end(){ return iterator(m_data.data(), m_head+N); }
    const_iterator begin() const { return cbegin(); }
    const_iterator end() const { return cend(); }
    const_iterator cbegin() const { return const_iterator(m_data.data(), m_head); }
    const_iterator cend() const { return const_iterator(m_data.data(), m_head+N); }
    reverse_iterator rbegin(){ return std::make_reverse_iterator(end());}      
    reverse_iterator rend(){ return std::make_reverse_iterator(begin());}  
    
    // shifts everything so that m_head is back at 0. 
    void linearize()
    {
      std::rotate(m_data.begin(), std::next(m_data.begin(), m_head), m_data.end()); 
      m_head = 0; 
    }  
    void rotate(difference_type r=1)
    {
      if constexpr(N <= 1) return; 
      auto n = static_cast<difference_type>(N);
      m_head = static_cast<size_type>((static_cast<difference_type>(m_head)+r%n+n) % n); 
    }
    void push_back(const T& val){
      if constexpr(N==0){
        return; 
      }
      else if constexpr(N==1){
        back().assign(val); 
      }
      else{
        rotate(1); 
        back().assign(val); 
      }
    }
    void push_back(T&& val){
      if constexpr(N==0){
        return; 
      }
      else if constexpr(N==1){
        back() = std::move(val); 
      }
      else{
        rotate(1); 
        back() = std::move(val); 
      }
    }
    template<class... Args>
    void emplace_back(Args&&... args){ 
      if constexpr(N==0){
        return; 
      }
      else if constexpr(N==1){
        back() = T(std::forward<Args>(args...))
      }
      else{
        rotate(1); 
        back() = T(std::forward<Args>(args...)); 
      }
    }
}; 

} // end namespace utils 
} // end namespace fdm 

#endif // CircularStackBuffer.hpp 
