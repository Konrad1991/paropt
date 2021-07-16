//#include "Rcpp.h"
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



// Matrix class
class M {
    public:
        std::vector<std::vector<double> > data;
        std::span<std::vector<double> > sp;


    M(int a, int b) {
        data.resize(a);

        for(int i = 0; i < a; i++) {
            data[i].resize(b);
        }

        sp = data;
    }

    friend std::ostream& operator<<(std::ostream& os, const M& obj) {
        for(int i = 0; i< obj.data.size(); i++) {
            for(int j = 0; j < obj.data[i].size(); j++) {
                os << obj.data[i][j] << ' ';
            }
            os << ' ' << std::endl;
        }
        return os;
    }

    void operator=(double inp) {
        for(int i = 0; i < data.size(); i++) {
            for(int j = 0; j < data[i].size(); j++) {
                data[i][j] = inp;
            }
        }
    }

    void operator=(M inp) {

        if(inp.data.size() != this -> data.size()) {
            throw std::invalid_argument("size incompatible");
        } 

        for(int i = 0; i < inp.data.size(); i++) {
            if(inp.data[i].size() != this -> data[i].size()) {
                throw std::invalid_argument("size incompatible");
            }
        }


        for(int i = 0; i < data.size(); i++) {
            for(int j = 0; j < data[i].size(); j++) {
                this -> data[i][j] = inp.data[i][j];
            }
        }

    }

    
    M operator+(M inp) {
        if(inp.data.size() != this -> data.size()) {
            throw std::invalid_argument("size incompatible");
        } 

        for(int i = 0; i < inp.data.size(); i++) {
            if(inp.data[i].size() != this -> data[i].size()) {
                throw std::invalid_argument("size incompatible");
            }
        }

        for(int i = 0; i < data.size(); i++) {
            for(int j = 0; j < data[i].size(); j++) {
                this -> sp[i][j] += inp.sp[i][j];
            }
        }

        return *this;
    }
    

    M operator+(double inp) {

        for(int i = 0; i < data.size(); i++) {
            for(int j = 0; j < data[i].size(); j++) {
                this -> data[i][j] += inp;
            }
        }

        return *this;
    }


    M operator-(M inp) {
        if(inp.data.size() != this -> data.size()) {
            throw std::invalid_argument("size incompatible");
        } 

        for(int i = 0; i < inp.data.size(); i++) {
            if(inp.data[i].size() != this -> data[i].size()) {
                throw std::invalid_argument("size incompatible");
            }
        }

        for(int i = 0; i < data.size(); i++) {
            for(int j = 0; j < data[i].size(); j++) {
                this -> data[i][j] -= inp.data[i][j];
            }
        }

        return *this;
    }

    M operator-(double inp) {

        for(int i = 0; i < data.size(); i++) {
            for(int j = 0; j < data[i].size(); j++) {
                this -> data[i][j] -= inp;
            }
        }

        return *this;
    }


    M operator*(M inp) {
        if(inp.data.size() != this -> data.size()) {
            throw std::invalid_argument("size incompatible");
        } 

        for(int i = 0; i < inp.data.size(); i++) {
            if(inp.data[i].size() != this -> data[i].size()) {
                throw std::invalid_argument("size incompatible");
            }
        }

        for(int i = 0; i < data.size(); i++) {
            for(int j = 0; j < data[i].size(); j++) {
                this -> data[i][j] *= inp.data[i][j];
            }
        }

        return *this;
    }

    M operator*(double inp) {

        for(int i = 0; i < data.size(); i++) {
            for(int j = 0; j < data[i].size(); j++) {
                this -> data[i][j] *= inp;
            }
        }

        return *this;
    }


    M operator/(M inp) {
        if(inp.data.size() != this -> data.size()) {
            throw std::invalid_argument("size incompatible");
        } 

        for(int i = 0; i < inp.data.size(); i++) {
            if(inp.data[i].size() != this -> data[i].size()) {
                throw std::invalid_argument("size incompatible");
            }
        }

        for(int i = 0; i < data.size(); i++) {
            for(int j = 0; j < data[i].size(); j++) {
                this -> data[i][j] /= inp.data[i][j];
            }
        }

        return *this;
    }

    M operator/(double inp) {

        for(int i = 0; i < data.size(); i++) {
            for(int j = 0; j < data[i].size(); j++) {
                this -> data[i][j] /= inp;
            }
        }

        return *this;
    }


    double operator()(int a, int b) {

        if((a < 0) | (b < 0)) {
            throw std::invalid_argument("size incompatible");
        }

        if(a > data.size()) {
            throw std::invalid_argument("size incompatible");
        }

        if(b > data[a].size()) {
            throw std::invalid_argument("size incompatible");
        }

        return this -> data[a][b];
    }



};


int main() {

    M m1(40,50);

    m1 = 4.5;

    M m2(40,50);

    m2 = 6.;

    //m1 = m1 - m2 - m2;


  MAT<double> x(10,5.4);
  MAT<double> y(10,10.3);

  MAT<double> result(10);
  
  result= x+x + y*y;
  
  std::cout << result << std::endl;
} 