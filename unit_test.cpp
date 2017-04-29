#include <gtest/gtest.h>
#include <iostream>
#include "finite_field.h"

uint8_t polynomials[][3] = {
  {0x07, 0x00, 0x00}, /* x^2 + x + 1 */
  {0x0b, 0x00, 0x00}, /* x^3 + x + 1 */
  {0x83, 0x00, 0x00}, /* x^7 + x + 1 */
  {0x1b, 0x20, 0x00}, /* x^13 + x^4 + x^3 + x + 1 */
  {0x09, 0x00, 0x02} /* x^17 + x^3 + 1*/
};

std::vector<std::vector<uint8_t>> get_polynomials() {
  std::vector<std::vector<uint8_t>> p;
  for (size_t i = 0; i < sizeof(polynomials)/sizeof(*polynomials); ++i) {
    size_t p_size = sizeof(polynomials[i]) / sizeof(*polynomials[i]);
    p.push_back(std::vector<uint8_t>(polynomials[i], polynomials[i] + p_size));
  }
  return p;
}

TEST(ffutil, inverse_1) {
  std::vector<std::vector<uint8_t>> p = get_polynomials();
  ffelement a(123, p[p.size() - 1]);
  ffelement b = a * a.inverse();
  EXPECT_TRUE(b.is_one());
}

TEST(ffutil, inverse_many) {
  std::vector<std::vector<uint8_t>> p = get_polynomials();
  for (size_t i = 0; i < p.size(); ++i) {
    for (size_t a0 = 1; a0 < 256; ++a0) {
      for (size_t a1 = 0; a1 < 256; ++a1) {
        ffelement a(a0, p[i]);
        ffelement poly(p[i], p[i]);
        if (a < poly) {
          ffelement b = a * a.inverse();
          EXPECT_TRUE(b.is_one());
        }
      }
    }
  }
}

TEST(ffutil, divide_eq) {
  std::vector<std::vector<uint8_t>> p = get_polynomials();
  ffelement a(123, p[p.size() - 1]);
  ffelement b(123, p[p.size() - 1]);
  ffelement q, r;
  q = ffelement::full_divide(a, b, &r);
  EXPECT_EQ(1, q.v_[0]);
  EXPECT_EQ(0, r.v_[0]);
}

TEST(ffutil, divide_1) {
  std::vector<std::vector<uint8_t>> p = get_polynomials();
  ffelement a(123, 1, p[p.size() - 1]);
  ffelement b(3, p[p.size() - 1]);
  ffelement q, r;
  q = ffelement::full_divide(a, b, &r);
  ffelement t = ffelement::mul_no_reduction(q, b) + r;
  EXPECT_EQ(t, a);
}

TEST(ffutil, divide_2) {
  std::vector<std::vector<uint8_t>> p = get_polynomials();
  ffelement a(123, 1, 128, p[p.size() - 1]);
  ffelement b(7, p[p.size() - 1]);
  ffelement q, r;
  q = ffelement::full_divide(a, b, &r);
  ffelement t = ffelement::mul_no_reduction(q, b) + r;
  EXPECT_EQ(t, a);
}

TEST(ffutil, divide_many) {
  std::vector<std::vector<uint8_t>> p = get_polynomials();
  for (size_t a0 = 0; a0 < 256; ++a0) {
    for (size_t a1 = 0; a1 < 2; ++a1) {
      ffelement a(a0, a1, p[p.size() - 1]);
      for (size_t b0 = 0; b0 < 256; ++b0) {
        for (size_t b1 = 0; b1 < 2; ++b1) {
          ffelement b(b0, b1, p[p.size() - 1]);
          ffelement q, r;
          q = ffelement::full_divide(a, b, &r);
          ffelement t = ffelement::mul_no_reduction(q, b) + r;
          EXPECT_EQ(t, a);
        }
      }
    }
  }
}

TEST(ffutil, degree) {
  std::vector<std::vector<uint8_t>> p = get_polynomials();
  ffelement a(1, p[p.size() - 1]);
  EXPECT_EQ(0U, a.degree());
  ffelement b(3, p[p.size() - 1]);
  EXPECT_EQ(1U, b.degree());
  ffelement c(1, 1, p[p.size() - 1]);
  EXPECT_EQ(8U, c.degree());
}

TEST(ffutil, mul_by_x) {
  std::vector<std::vector<uint8_t>> p = get_polynomials();
  ffelement t = ffelement::mul_by_x_no_reduction(ffelement(1, p[p.size() - 1]));
  EXPECT_EQ(1U, t.v_.size());
  EXPECT_EQ(2, t.v_[0]);
  t = ffelement::mul_by_x_no_reduction(ffelement(63, p[p.size() - 1]));
  EXPECT_EQ(1U, t.v_.size());
  EXPECT_EQ(126, t.v_[0]);
  t = ffelement::mul_by_x_no_reduction(ffelement(0xc8, 0x01, p[p.size() - 1]));
  EXPECT_EQ(2U, t.v_.size());
  EXPECT_EQ(0x90, t.v_[0]);
  EXPECT_EQ(0x03, t.v_[1]);
  t = ffelement::mul_by_x_no_reduction(ffelement(128, p[p.size() - 1]));
  EXPECT_EQ(2U, t.v_.size());
  EXPECT_EQ(0, t.v_[0]);
  EXPECT_EQ(1, t.v_[1]);
}

TEST(ffutil, monomial) {
  std::vector<std::vector<uint8_t>> p = get_polynomials();
  ffelement a = ffelement::monomial(4, p[p.size() - 1]);
  EXPECT_EQ(1U, a.v_.size());
  EXPECT_EQ(16, a.v_[0]);
  ffelement b = ffelement::monomial(15, p[p.size() - 1]);
  EXPECT_EQ(2U, b.v_.size());
  EXPECT_EQ(128, b.v_[1]);
  ffelement c = ffelement::monomial(0, p[p.size() - 1]);
  EXPECT_EQ(1U, c.v_.size());
  EXPECT_EQ(1, c.v_[0]);
}

TEST(ffutil, compare_ineq_1) {
  std::vector<std::vector<uint8_t>> p = get_polynomials();
  ffelement a(0x00, p[0]);
  ffelement b(0x01, p[0]);
  EXPECT_LT(a, b);
  EXPECT_GT(b, a);
}

TEST(ffutil, compare_ineq_2) {
  std::vector<std::vector<uint8_t>> p = get_polynomials();
  ffelement a(0xff, 0xf7, p[p.size() - 1]);
  ffelement b(0xff, 0xff, p[p.size() - 1]);
  EXPECT_LT(a, b);
  EXPECT_GT(b, a);
}

TEST(ffutil, compare_ineq_3) {
  std::vector<std::vector<uint8_t>> p = get_polynomials();
  ffelement a(0x7f, 0xff, p[p.size() - 1]);
  ffelement b(0xff, 0xff, p[p.size() - 1]);
  EXPECT_LT(a, b);
  EXPECT_GT(b, a);
}

TEST(ffutil, compare_eq_1) {
  std::vector<std::vector<uint8_t>> p = get_polynomials();
  ffelement a(0x01, p[0]);
  ffelement b(0x01, p[0]);
  EXPECT_EQ(a, b);
  EXPECT_GE(a, b);
  EXPECT_LE(a, b);
}

TEST(ffutil, compare_eq_2) {
  std::vector<std::vector<uint8_t>> p = get_polynomials();
  ffelement a(0x01, p[1]);
  ffelement b(0x01, p[1]);
  EXPECT_EQ(a, b);
  EXPECT_GE(a, b);
  EXPECT_LE(a, b);
}

TEST(ffutil, compare_eq_3) {
  std::vector<std::vector<uint8_t>> p = get_polynomials();
  ffelement a(0xff, 0xff, p[p.size() - 1]);
  ffelement b(0xff, 0xff, p[p.size() - 1]);
  EXPECT_EQ(a, b);
  EXPECT_GE(a, b);
  EXPECT_LE(a, b);
}

TEST(ffutil, compare_ne_1) {
  std::vector<std::vector<uint8_t>> p = get_polynomials();
  ffelement a(0xff, 0xff, p[p.size() - 1]);
  ffelement b(0xf7, 0xff, p[p.size() - 1]);
  EXPECT_NE(a, b);
}

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
