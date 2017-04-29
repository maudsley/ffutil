#ifndef __FINITE_FIELD_H__
#define __FINITE_FIELD_H__

#include <stdexcept>
#include <vector>
#include <iostream>
#include <sstream>

class ffelement_exception : public std::logic_error {
public:
  ffelement_exception(const std::string& what) : std::logic_error(what) {
  }
};

class ffelement {
public:

  ffelement();
  explicit ffelement(const std::vector<uint8_t>& p);
  ffelement(const std::vector<uint8_t>& v, const std::vector<uint8_t>& p);
  ffelement(const uint8_t v0, const std::vector<uint8_t>& p);
  ffelement(const uint8_t v0, const uint8_t v1, const std::vector<uint8_t>& p);
  ffelement(const uint8_t v0, const uint8_t v1, const uint8_t v2, const std::vector<uint8_t>& p);
  ffelement(const ffelement& rhs);
  ffelement& operator =(const ffelement& rhs);

  bool is_zero() const;
  bool is_one() const;
  size_t degree() const;
  std::string to_string() const;

  bool operator ==(const ffelement& rhs) const;
  bool operator !=(const ffelement& rhs) const;
  bool operator <(const ffelement& rhs) const;
  bool operator >(const ffelement& rhs) const;
  bool operator <=(const ffelement& rhs) const;
  bool operator >=(const ffelement& rhs) const;

  ffelement operator +(const ffelement& rhs) const;
  ffelement& operator +=(const ffelement& rhs);

  ffelement operator *(const ffelement& rhs) const;
  ffelement& operator *=(const ffelement& rhs);

  ffelement inverse() const;

  static ffelement gcd(const ffelement& a, const ffelement& b, ffelement* p, ffelement* q);
  static ffelement full_divide(const ffelement& lhs, const ffelement& rhs, ffelement* remainder);
  static ffelement mul_no_reduction(const ffelement& lhs, const ffelement& rhs);
  static ffelement mul_by_x_no_reduction(const ffelement& arg);
  static ffelement monomial(const size_t degree, const std::vector<uint8_t>& p);

  static bool vectors_match(const std::vector<uint8_t>& lhs, const std::vector<uint8_t>& rhs);

  std::vector<uint8_t> v_; /* element */
  std::vector<uint8_t> p_; /* polynomial (ideal) */

};

#endif /* __FINITE_FIELD_H__ */
