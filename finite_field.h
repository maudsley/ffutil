#ifndef __FINITE_FIELD_H__
#define __FINITE_FIELD_H__

#include <exception>
#include <vector>

class ffelement_exception : public std::logic_error {
public:
  ffelement_exception(const std::string& what) : std::logic_error(what) {
  }
};

class ffelement {
public:

  ffelement() {
  }

  ffelement(const std::vector<uint8_t>& p) : p_(p) {
    v_ = std::vector<uint8_t>(1, 0); /* zero */
  }

  ffelement(const std::vector<uint8_t>& v, const std::vector<uint8_t>& p) : v_(v), p_(p) {
  }

  ffelement(const uint8_t v0, const std::vector<uint8_t>& p) : p_(p) {
    v_.push_back(v0);
  }

  ffelement(const uint8_t v0, const uint8_t v1, const std::vector<uint8_t>& p) : p_(p) {
    v_.push_back(v0);
    v_.push_back(v1);
  }

  ffelement(const uint8_t v0, const uint8_t v1, const uint8_t v2, const std::vector<uint8_t>& p) : p_(p) {
    v_.push_back(v0);
    v_.push_back(v1);
    v_.push_back(v2);
  }

  ffelement(const ffelement& rhs) {
    *this = rhs;
  }

  ffelement& operator =(const ffelement& rhs) {
    v_ = rhs.v_;
    p_ = rhs.p_;
    return *this;
  }

  size_t degree() const {
    size_t d = 0;
    for (size_t i = 0; i < v_.size(); ++i) {
      for (size_t j = 0; j < 8; ++j) {
        if ((v_[i] >> j) & 1) {
          d = i * 8 + j;
        }
      }
    }
    return d;
  }

  bool operator ==(const ffelement& rhs) const {
    if (!vectors_match(p_, rhs.p_)) {
      return false;
    }
    return vectors_match(v_, rhs.v_);
  }

  bool operator !=(const ffelement& rhs) const {
    if (!vectors_match(p_, rhs.p_)) {
      return false;
    }
    return !vectors_match(v_, rhs.v_);
  }

  bool operator <(const ffelement& rhs) const {
    if (!vectors_match(p_, rhs.p_)) {
      return false;
    }
    size_t degree_lhs = degree();
    size_t degree_rhs = rhs.degree();
    if (degree_lhs < degree_rhs) {
      return true;
    } else if (degree_lhs > degree_rhs) {
      return false;
    } else { /* same degree */
      for (size_t i = 0; i < v_.size(); ++i) {
        if (v_[i] < rhs.v_[i]) {
          return true;
        }
      }
      return false;
    }
  }

  bool operator >(const ffelement& rhs) const {
    if (!vectors_match(p_, rhs.p_)) {
      return false;
    }
    size_t degree_lhs = degree();
    size_t degree_rhs = rhs.degree();
    if (degree_lhs > degree_rhs) {
      return true;
    } else if (degree_lhs < degree_rhs) {
      return false;
    } else { /* same degree */
      for (size_t i = 0; i < v_.size(); ++i) {
        if (v_[i] > rhs.v_[i]) {
          return true;
        }
      }
      return false;
    }
  }

  bool operator <=(const ffelement& rhs) const {
    if (!vectors_match(p_, rhs.p_)) {
      return false;
    }
    return !(*this > rhs);
  }

  bool operator >=(const ffelement& rhs) const {
    if (!vectors_match(p_, rhs.p_)) {
      return false;
    }
    return !(*this < rhs);
  }

  ffelement operator +(const ffelement& rhs) const {
    if (!vectors_match(p_, rhs.p_)) {
      throw ffelement_exception("Polynomials do not match.");
    }
    ffelement r;
    if (v_.size() > rhs.v_.size()) {
      r = *this;
      for (size_t i = 0; i < rhs.v_.size(); ++i) {
        r.v_[i] ^= rhs.v_[i];
      }
    } else {
      r = rhs;
      for (size_t i = 0; i < rhs.v_.size(); ++i) {
        r.v_[i] ^= v_[i];
      }
    }
    return r;
  }

  ffelement& operator +=(const ffelement& rhs) {
    *this = *this + rhs;
    return *this;
  }

private:

  static ffelement full_divide(const ffelement& lhs, const ffelement& rhs, ffelement* remainder) {
    if (!vectors_match(lhs.p_, rhs.p_)) {
      throw ffelement_exception("Polynomials do not match.");
    }
    if (lhs.degree() < rhs.degree()) {
      if (remainder) {
        *remainder = lhs;
      }
      return ffelement(lhs.p_); /* zero */
    }
    
  }

  static ffelement mul_no_reduction(const ffelement& lhs, const ffelement& rhs) {
    if (!vectors_match(lhs.p_, rhs.p_)) {
      throw ffelement_exception("Polynomials do not match.");
    }
    ffelement r(lhs.p_);
    ffelement k = rhs;
    for (size_t i = 0; i < lhs.v_.size(); ++i) {
      for (size_t j = 0; j < 8; ++j) {
        if ((rhs.v_[i] >> j) & 1) {
          r += k;
        }
        k = mul_by_x2_no_reduction(k);
      }
    }
    return r;
  }

  static ffelement mul_by_x2_no_reduction(const ffelement& arg) {
    ffelement r = arg;
    int carry = 0;
    for (size_t i = 0; i < r.v_.size(); ++i) {
      int next_carry = r.v_[i] >> 7;
      r.v_[i] <<= 1;
      r.v_[i] |= carry;
      carry = next_carry;
    }
    if (carry) {
      r.v_.push_back(1);
    }
    return r;
  }

  static bool vectors_match(const std::vector<uint8_t>& lhs, const std::vector<uint8_t>& rhs) {
    if (lhs.size() != rhs.size()) {
      return false;
    }
    for (size_t i = 0; i < lhs.size(); ++i) {
      if (lhs[i] != rhs[i]) {
        return false;
      }
    }
    return true;
  }

  std::vector<uint8_t> v_; /* element */
  std::vector<uint8_t> p_; /* ideal */

};

#endif /* __FINITE_FIELD_H__ */
