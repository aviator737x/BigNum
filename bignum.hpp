#ifndef BIGNUM_BIGNUM_H
#define BIGNUM_BIGNUM_H

#include <string>
#include <ostream>
#include <istream>
#include <vector>
#include <limits>
#include <sstream>

namespace mp {
    struct bignum {
        bignum() {
            digits_.push_back(0);
        };
        bignum(uint32_t num) {
            digits_.push_back(num);
        };
        explicit bignum(const std::string& num) {
            bignum result;
            for (size_t i = 0; i < num.size(); ++i) {
                result *= 10;
                bignum temp(static_cast<uint32_t >(num[i] - '0'));
                result += temp;
            }
            swap(result, *this);
        };
        
        bignum(const bignum& num) = default;
        bignum &operator=(bignum rhs) {
            swap(*this, rhs);
            return *this;
        }
        
        explicit operator bool() const {
            return digits_.size() > 1 || digits_[0];
        }
        
        explicit operator uint32_t () const {
            return digits_[0];
        }
        
        std::string to_string() const {
            std::string str;
            bignum num = *this;
            
            str = std::to_string(num.devide(10));
            while (num) {
                str = std::to_string(num.devide(10)) + str;
            }
            return str;
        }
        
        bignum &operator+=(const bignum& num) {
            bool carry = false;
            uint64_t base = 4294967296;
            
            for (size_t i=0; i < std::max(digits_.size(), num.digits_.size()) || carry; ++i) {
                if (i == digits_.size())
                    digits_.push_back(0);
                uint64_t tmp = (uint64_t)digits_[i] + (uint64_t)(carry ? 1 : 0) + (uint64_t)(i < num.digits_.size() ? num.digits_[i] : 0);
                carry = (tmp >= base);
                if (carry) {
                    tmp -= base;
                }
                digits_[i] = static_cast<uint32_t >(tmp);
            }
            
            return *this;
        }
        bignum &operator*=(const bignum& num) {
            uint64_t base = 4294967296;
            
            std::vector<uint32_t> c(digits_.size() + num.digits_.size(), 0);
            for (size_t i = 0; i < digits_.size(); ++i) {
                uint64_t carry = 0;
                for (size_t j = 0; j < num.digits_.size() || carry; ++j) {
                    uint64_t cur = c[i+j] + (uint64_t)digits_[i] * (uint64_t)(j < num.digits_.size() ? num.digits_[j] : 0) + carry;
                    c[i+j] = (uint32_t)(cur % base);
                    carry = (uint32_t)(cur / base);
                }
            }
            while (c.size() > 1 && c.back() == 0)
                c.pop_back();
            
            digits_ = c;
            return *this;
        }
        
    private:
        uint64_t devide(uint64_t b) {
            uint64_t base = 4294967296;
            uint64_t carry = 0;
            for (int64_t i = digits_.size() - 1; i >= 0; --i) {
                uint64_t cur = (uint64_t)digits_[i] + carry * base;
                digits_[i] = uint32_t (cur / b);
                carry = uint32_t (cur % b);
            }
            while (digits_.size() > 1 && digits_.back() == 0)
                digits_.pop_back();
            return carry;
        }
        
        friend std::ostream& operator<<(std::ostream& os, const bignum& num);
        friend std::istream& operator>>(std::istream& is, bignum& num);
        friend void swap(bignum& lhs, bignum& rhs);
        
        std::vector<uint32_t> digits_;
    };
    
    bignum operator+(const bignum& lhs, const bignum& rhs) {
        bignum res = lhs;
        res+=rhs;
        return res;
    }
    bignum operator*(const bignum& lhs, const bignum& rhs) {
        bignum res = lhs;
        res*=rhs;
        return res;
    }
    
    void swap(bignum& lhs, bignum& rhs) {
        std::swap(lhs.digits_, rhs.digits_);
    }
    
    std::ostream& operator<<(std::ostream& os, const bignum& num) {
        os << num.to_string();
        return os;
    }
    std::istream& operator>>(std::istream& is, bignum& num) {
        std::string str;
        is >> str;
        num = bignum(str);
        return is;
    }
    
    
    
    struct polinomial {
        explicit polinomial(const std::string & polinom) {
            coefficients.resize(1000);
            uint32_t coeff = 0;
            std::vector<int> digits;
            for (int i = 0; i < polinom.size(); ++i) {
                if (polinom[i] >= '0' && polinom[i] <= '9')
                {
                    digits.push_back(polinom[i] - '0');
                }
                else if (polinom[i] == '+') {
                    continue;
                }
                else {
                    i += 3;
                    std::vector<int> degree;
                    
                    while (i < polinom.size() && polinom[i] >= '0' && polinom[i] <= '9')
                    {
                        degree.push_back(polinom[i] - '0');
                        i++;
                    }
                    int deg = 0;
                    int base = 1;
                    for (int j = 0; j < degree.size(); ++j) {
                        deg += degree[degree.size() - 1 - j] * base;
                        base *= 10;
                    }
                    uint32_t numb = 0;
                    base = 1;
                    for (int j = 0; j < digits.size(); ++j) {
                        numb += digits[digits.size() - 1 - j] * base;
                        base *= 10;
                    }
                    if (deg + 1 >= coefficients.size()) {
                        coefficients.resize(coefficients.size() * 2);
                    }
                    
                    if (deg > degr)
                        degr = deg;
                    coefficients[deg] = numb;
                    digits.clear();
                    
                }
            }
        }
        
        uint32_t at(size_t i) {
            if (i > degr)
                degr = i;
            return coefficients[i];
        }
        
        uint32_t at(size_t i) const {
            if (i > degr)
                return 0;
            return coefficients[i];
        }
        
        template <typename T>
        T operator() (T x) {
            std::vector<T> b(degr+1);
            b[degr] = coefficients[degr-1];
            for (int i = degr - 1; i >= 0; --i) {
                b[i] = coefficients[i] + b[i+1] * x;
            }
            if (degr > 0) {
                b[0] = coefficients[0] + b[1] * x;
            }
            return b[0];
        }
        
        
    private:
        int degr = 0;
        std::vector<uint32_t> coefficients;
    };
}


#endif //BIGNUM_BIGNUM_H


