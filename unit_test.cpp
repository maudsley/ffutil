#include <gtest/gtest.h>
#include <iostream>
#include "finite_field.h"

uint8_t polynomials[][3] = {
  {0x00, 0x00, 0x07}, /* x^2 + x + 1 */
  {0x00, 0x00, 0x0b}, /* x^3 + x + 1 */
  {0x00, 0x00, 0x83}, /* x^7 + x + 1 */
  {0x00, 0x20, 0x1b}, /* x^13 + x^4 + x^3 + x + 1 */
  {0x02, 0x00, 0x09} /* 1 + x^3 + x^17 */
};

std::vector<std::vector<uint8_t>> get_polynomials() {
  std::vector<std::vector<uint8_t>> p;
  for (size_t i = 0; i < sizeof(polynomials)/sizeof(*polynomials); ++i) {
    size_t p_size = sizeof(polynomials[i]) / sizeof(*polynomials[i]);
    p.push_back(std::vector<uint8_t>(polynomials[i], polynomials[i] + p_size));
  }
  return p;
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
  ffelement a(0xff, 0xf7, p[4]);
  ffelement b(0xff, 0xff, p[4]);
  EXPECT_LT(a, b);
  EXPECT_GT(b, a);
}

TEST(ffutil, compare_ineq_3) {
  std::vector<std::vector<uint8_t>> p = get_polynomials();
  ffelement a(0x7f, 0xff, p[4]);
  ffelement b(0xff, 0xff, p[4]);
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
  ffelement a(0xff, 0xff, p[4]);
  ffelement b(0xff, 0xff, p[4]);
  EXPECT_EQ(a, b);
  EXPECT_GE(a, b);
  EXPECT_LE(a, b);
}

TEST(ffutil, compare_ne_1) {
  std::vector<std::vector<uint8_t>> p = get_polynomials();
  ffelement a(0xff, 0xff, p[4]);
  ffelement b(0xf7, 0xff, p[4]);
  EXPECT_NE(a, b);
}

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
