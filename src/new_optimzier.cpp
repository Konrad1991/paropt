//#include "Rcpp.h"
#include <iostream>
#include <vector>
#include <stdexcept>
#include <ranges>
#include <span>


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
                this -> data[i][j] += inp.data[i][j];
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

    M m1(4,5);

    m1 = 4.5;

    M m2(4,5);

    m2 = 6.;

    m1 = m1 + m2/3. - .5;

    std::cout << m1 << std::endl;


    std::cout << m1(2,3) << std::endl;


    std::vector<double> vec{1,2,3};

    auto v = std::span{vec};

} 