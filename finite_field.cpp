#include "finite_field.h"

ffelement::ffelement() {
}

ffelement::ffelement(const std::vector<uint8_t>& p) : p_(p) {
  v_ = std::vector<uint8_t>(1, 0); /* zero */
}

ffelement::ffelement(const std::vector<uint8_t>& v, const std::vector<uint8_t>& p) : v_(v), p_(p) {
}

ffelement::ffelement(const uint8_t v0, const std::vector<uint8_t>& p) : p_(p) {
  v_.push_back(v0);
}

ffelement::ffelement(const uint8_t v0, const uint8_t v1, const std::vector<uint8_t>& p) : p_(p) {
  v_.push_back(v0);
  v_.push_back(v1);
}

ffelement::ffelement(const uint8_t v0, const uint8_t v1, const uint8_t v2, const std::vector<uint8_t>& p) : p_(p) {
  v_.push_back(v0);
  v_.push_back(v1);
  v_.push_back(v2);
}

ffelement::ffelement(const ffelement& rhs) {
  *this = rhs;
}

ffelement& ffelement::operator =(const ffelement& rhs) {
  v_ = rhs.v_;
  p_ = rhs.p_;
  return *this;
}

bool ffelement::is_zero() const {
  for (size_t i = 0; i < v_.size(); ++i) {
    if (v_[i]) {
      return false;
    }
  }
  return true;
}

bool ffelement::is_one() const {
  return v_.size() >= 1 && v_[0] == 1;
}

size_t ffelement::degree() const {
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

std::string ffelement::to_string() const {
  if (is_zero()) {
    return "{0}";
  }
  std::stringstream stream;
  stream << "{";
  size_t term_count = 0;
  size_t d = degree() + 1;
  for (size_t i = 0; i < d; ++i) {
    size_t j = d - i - 1;
    size_t byte = j / 8;
    size_t bit = j % 8;
    size_t value = (v_[byte] >> bit) & 1;
    if (value) {
      stream << (i != 0 ? " + " : "");
      switch (j) {
      case 0:
        stream << "1";
        break;
      case 1:
        stream << "x";
        break;
      default:
        stream << "x^" << j;
        break;
      }
      term_count++;
    }
  }
  if (!term_count) {
    stream << "0";
  }
  stream << "}";
  return stream.str();
}

bool ffelement::operator ==(const ffelement& rhs) const {
  if (!vectors_match(p_, rhs.p_)) {
    return false;
  }
  return vectors_match(v_, rhs.v_);
}

bool ffelement::operator !=(const ffelement& rhs) const {
  if (!vectors_match(p_, rhs.p_)) {
    return false;
  }
  return !vectors_match(v_, rhs.v_);
}

bool ffelement::operator <(const ffelement& rhs) const {
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

bool ffelement::operator >(const ffelement& rhs) const {
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

bool ffelement::operator <=(const ffelement& rhs) const {
  if (!vectors_match(p_, rhs.p_)) {
    return false;
  }
  return !(*this > rhs);
}

bool ffelement::operator >=(const ffelement& rhs) const {
  if (!vectors_match(p_, rhs.p_)) {
    return false;
  }
  return !(*this < rhs);
}

ffelement ffelement::operator +(const ffelement& rhs) const {
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
    for (size_t i = 0; i < v_.size(); ++i) {
      r.v_[i] ^= v_[i];
    }
  }
  return r;
}

ffelement& ffelement::operator +=(const ffelement& rhs) {
  *this = *this + rhs;
  return *this;
}

ffelement ffelement::operator *(const ffelement& rhs) const {
  if (!vectors_match(p_, rhs.p_)) {
    throw ffelement_exception("Polynomials do not match.");
  }
  ffelement r;
  ffelement n = mul_no_reduction(*this, rhs);
  ffelement poly(p_, p_);
  full_divide(n, poly, &r);
  //std::cout << "Multiplying: " << to_string() << std::endl;
  //std::cout << " by: " << rhs.to_string() << std::endl;
  //std::cout << " got: " << r.to_string() << std::endl;
  return r;
}

ffelement& ffelement::operator *=(const ffelement& rhs) {
  *this = *this * rhs;
  return *this;
}

ffelement ffelement::inverse() const {
  ffelement inv;
  ffelement poly(p_, p_);
  ffelement r = gcd(*this, poly, &inv, 0);
  if (!r.is_one()) {
    throw ffelement_exception("Element has no inverse");
  }
  return inv;
}

ffelement ffelement::gcd(const ffelement& a, const ffelement& b, ffelement* p, ffelement* q) {
  /* sort */
  bool swap = false;
  ffelement wa, wb;
  if (a < b) {
    wa = a;
    wb = b;
  } else {
    wa = b;
    wb = a;
    swap = true;
  }

  /* special case where either inputs are zero */
  if (wa.is_zero()) {
    if (p) {
      *p = ffelement(a.p_);
    }
    if (q) {
      *q = ffelement(a.p_);
    }
    return ffelement(a.p_); /* zero */
  }

  /* init bezout coefficients */
  ffelement wp[2], wq[2];
  wp[0] = wq[1] = monomial(0, a.p_); /* one */
  wp[1] = wq[0] = ffelement(a.p_); /* zero */

  //std::cout << "Computing GCD of " << std::endl;
  //std::cout << " " << wa.to_string() << std::endl;
  //std::cout << " " << wb.to_string() << std::endl;

  /* extended Euclidean algorithm */
  ffelement r;
  do {
    ffelement d = full_divide(wa, wb, &r);
    ffelement u = wp[0] + d * wp[1];
    ffelement v = wq[0] + d * wq[1];
    wa = wb;
    wb = r;
    wp[0] = wp[1];
    wq[0] = wq[1];
    wp[1] = u;
    wq[1] = v;
    //std::cout << "Remainder is: " << r.to_string() << std::endl;
  } while (!r.is_zero());

  if (swap) {
    ffelement t = wp[0];
    wp[0] = wq[0];
    wq[0] = t;
  }

  if (p) {
    *p = wp[0];
  }

  if (q) {
    *q = wq[0];
  }

  //std::cout << "GCD is: " << wa.to_string() << std::endl;
  //std::cout << "P is: " << wp[0].to_string() << std::endl;
  //std::cout << "Q is: " << wq[0].to_string() << std::endl;

  return wa;
}

ffelement ffelement::full_divide(const ffelement& lhs, const ffelement& rhs, ffelement* remainder) {
  if (!vectors_match(lhs.p_, rhs.p_)) {
    throw ffelement_exception("Polynomials do not match.");
  }
  if (lhs.degree() < rhs.degree()) {
    if (remainder) {
      *remainder = lhs;
    }
    return ffelement(lhs.p_); /* zero */
  }
  //std::cout << "Dividing " << lhs.to_string() << std::endl;
  //std::cout << " by " << rhs.to_string() << std::endl;
  ffelement q(lhs.p_);
  ffelement w = lhs;
  while (w.degree() >= rhs.degree()) {
    //std::cout << "Working is: " << w.to_string() << std::endl;
    size_t delta = w.degree() - rhs.degree();
    ffelement m = monomial(delta, lhs.p_);
    //std::cout << "Monomial is: " << m.to_string() << std::endl;
    q += m;
    //std::cout << "Quotient is: " << q.to_string() << std::endl;
    ffelement s = mul_no_reduction(rhs, m);
    //std::cout << "Will eliminate: " << s.to_string() << std::endl;
    w += s;
    if (w.is_zero()) {
      break;
    }
    if (q.is_zero()) {
      break;
    }
  }
  //std::cout << "Final quotient is: " << q.to_string() << std::endl;
  //std::cout << "Final remainder is: " << w.to_string() << std::endl;
  if (remainder) {
    *remainder = w;
  }
  return q;
}

ffelement ffelement::mul_no_reduction(const ffelement& lhs, const ffelement& rhs) {
  if (!vectors_match(lhs.p_, rhs.p_)) {
    throw ffelement_exception("Polynomials do not match.");
  }
  ffelement r(lhs.p_);
  ffelement k = rhs;
  for (size_t i = 0; i < lhs.v_.size(); ++i) {
    for (size_t j = 0; j < 8; ++j) {
      if ((lhs.v_[i] >> j) & 1) {
        r += k;
      }
      k = mul_by_x_no_reduction(k);
    }
  }
  return r;
}

ffelement ffelement::mul_by_x_no_reduction(const ffelement& arg) {
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

ffelement ffelement::monomial(const size_t degree, const std::vector<uint8_t>& p) {
  ffelement r(p);
  if (!degree) {
    r.v_[0] = 1;
    return r;
  } else {
    size_t bytes = degree / 8;
    size_t index = degree % 8;
    r.v_ = std::vector<uint8_t>(bytes + 1, 0);
    r.v_[bytes] |= 1 << index;
    return r;
  }
}

bool ffelement::vectors_match(const std::vector<uint8_t>& lhs, const std::vector<uint8_t>& rhs) {
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
