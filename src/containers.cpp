#include <iostream>
#include <vector>
#include <stdexcept>
#include <ranges>
#include <span>


// https://www.modernescpp.com/index.php/expression-templates

template<typename T, typename VEC= std::vector<T> >
class MAT {

    VEC d;

    public:

    MAT(const std::size_t n) : d(n) {}

    MAT(const std::size_t n, const double init) : d(n, init) {}

    MAT(const VEC& other) : d(other) {}

    friend std::ostream& operator<<(std::ostream& os, const MAT     & obj) {
        for(int i = 0; i< obj.d.size(); i++) {
                os << obj.d[i] << std::endl;
        }
        return os;
    }

    template <typename T2, typename R2>
    MAT& operator=(const MAT<T2, R2>& other) {
        for(int i = 0; i < d.size(); i++) {
            d[i] = other[i];
        }

      return *this;
    }

    std::size_t size() const{
        return d.size();
    }

    T operator[](const std::size_t i) const{
        return d[i];
    }

    T& operator[](const std::size_t i) {
        return d[i];
    }

    const VEC& data() const { 
        return d; 
    }

    VEC& data() { 
        return d; 
    }
};


template<typename T, typename OP1, typename OP2>
class MAT_PLUS {
    const OP1& op1;
    const OP2& op2;

    public:

    MAT_PLUS(const OP1& a, const OP2& b): op1(a), op2(b){}

    T operator[](const std::size_t i) const{
        return op1[i] + op2[i];
    }

    std::size_t size() const {
        return op1.size();
    }
};


template<typename T, typename OP1, typename OP2>
class MAT_TIME {
    const OP1& op1;
    const OP2& op2;

    public:

    MAT_TIME(const OP1& a, const OP2& b): op1(a), op2(b){}

    T operator[](const std::size_t i) const{
        return op1[i] * op2[i];
    }

    std::size_t size() const {
        return op1.size();
    }
};


template<typename T, typename R1, typename R2>
MAT<T, MAT_PLUS<T, R1, R2> > operator+(
    const MAT<T, R1>&a, const MAT<T, R2>&b) {
        return MAT<T, MAT_PLUS<T, R1, R2> >(MAT_PLUS<T, R1, R2 >(a.data(), b.data()));
}


template<typename T, typename R1, typename R2>
MAT<T, MAT_TIME<T, R1, R2> > operator*(
    const MAT<T, R1>&a, const MAT<T, R2>&b) {
        return MAT<T, MAT_TIME<T, R1, R2> >(MAT_TIME<T, R1, R2 >(a.data(), b.data()));
}


int main() {

  MAT<double> x(10,5.4);
  MAT<double> y(10,10.3);

  MAT<double> result(10);
  
  result= x+x + y*y;
  
  std::cout << result << std::endl;
} 