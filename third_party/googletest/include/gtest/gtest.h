#pragma once

#include <cmath>
#include <exception>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace testing {

using TestFunction = void (*)();

struct TestSpec {
    std::string suite;
    std::string name;
    TestFunction function = nullptr;
};

class FatalTestFailure : public std::exception {
public:
    const char* what() const noexcept override { return "fatal test failure"; }
};

void RegisterTest(const char* suite, const char* name, TestFunction function);
void AddFailure(const char* file, int line, const std::string& message, bool fatal);
void InitGoogleTest(int* argc, char** argv);
int RunAllTests();

namespace internal {

class TestRegistrar {
public:
    TestRegistrar(const char* suite, const char* name, TestFunction function) {
        RegisterTest(suite, name, function);
    }
};

template <typename T>
std::string ToString(const T& value) {
    std::ostringstream oss;
    oss << value;
    return oss.str();
}

}  // namespace internal

}  // namespace testing

#define RUN_ALL_TESTS() ::testing::RunAllTests()

#define TEST(SuiteName, TestName)                                                                                   \
    static void SuiteName##_##TestName##_Impl();                                                                    \
    static ::testing::internal::TestRegistrar SuiteName##_##TestName##_registrar(#SuiteName, #TestName,            \
                                                                                   &SuiteName##_##TestName##_Impl); \
    static void SuiteName##_##TestName##_Impl()

#define GTEST_EXPECT_COND(cond, fatal)                                                      \
    do {                                                                                     \
        if (!(cond)) {                                                                       \
            std::ostringstream _gtest_oss;                                                   \
            _gtest_oss << "Expectation failed: " << #cond;                                 \
            ::testing::AddFailure(__FILE__, __LINE__, _gtest_oss.str(), (fatal));           \
        }                                                                                    \
    } while (0)

#define GTEST_EXPECT_BINARY(op, lhs, rhs, fatal)                                             \
    do {                                                                                      \
        const auto _gtest_lhs = (lhs);                                                       \
        const auto _gtest_rhs = (rhs);                                                       \
        if (!(_gtest_lhs op _gtest_rhs)) {                                                   \
            std::ostringstream _gtest_oss;                                                   \
            _gtest_oss << "Expectation failed: " << #lhs << " " << #op << " " << #rhs;      \
            ::testing::AddFailure(__FILE__, __LINE__, _gtest_oss.str(), (fatal));           \
        }                                                                                     \
    } while (0)

#define EXPECT_TRUE(cond) GTEST_EXPECT_COND((cond), false)
#define EXPECT_FALSE(cond) GTEST_EXPECT_COND(!(cond), false)
#define ASSERT_TRUE(cond) GTEST_EXPECT_COND((cond), true)
#define ASSERT_FALSE(cond) GTEST_EXPECT_COND(!(cond), true)

#define EXPECT_EQ(lhs, rhs) GTEST_EXPECT_BINARY(==, lhs, rhs, false)
#define EXPECT_NE(lhs, rhs) GTEST_EXPECT_BINARY(!=, lhs, rhs, false)
#define EXPECT_LT(lhs, rhs) GTEST_EXPECT_BINARY(<, lhs, rhs, false)
#define EXPECT_LE(lhs, rhs) GTEST_EXPECT_BINARY(<=, lhs, rhs, false)
#define EXPECT_GT(lhs, rhs) GTEST_EXPECT_BINARY(>, lhs, rhs, false)
#define EXPECT_GE(lhs, rhs) GTEST_EXPECT_BINARY(>=, lhs, rhs, false)

#define ASSERT_EQ(lhs, rhs) GTEST_EXPECT_BINARY(==, lhs, rhs, true)
#define ASSERT_NE(lhs, rhs) GTEST_EXPECT_BINARY(!=, lhs, rhs, true)

#define EXPECT_NEAR(lhs, rhs, abs_error)                                                           \
    do {                                                                                            \
        const auto _gtest_lhs = (lhs);                                                             \
        const auto _gtest_rhs = (rhs);                                                             \
        const auto _gtest_tol = (abs_error);                                                       \
        if (std::fabs((_gtest_lhs) - (_gtest_rhs)) > (_gtest_tol)) {                              \
            std::ostringstream _gtest_oss;                                                         \
            _gtest_oss << "Expectation failed: |" << #lhs << " - " << #rhs << "| <= "        \
                       << #abs_error << " (lhs=" << _gtest_lhs << ", rhs=" << _gtest_rhs       \
                       << ", tol=" << _gtest_tol << ")";                                        \
            ::testing::AddFailure(__FILE__, __LINE__, _gtest_oss.str(), false);                   \
        }                                                                                           \
    } while (0)
