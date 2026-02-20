#include "gtest/gtest.h"

#include <cstddef>

namespace testing {
namespace {

std::vector<TestSpec>& Registry() {
    static auto* registry = new std::vector<TestSpec>();
    return *registry;
}

int& CurrentTestFailureCount() {
    static int failureCount = 0;
    return failureCount;
}

int& FailedTestCount() {
    static int failedTestCount = 0;
    return failedTestCount;
}

}  // namespace

void RegisterTest(const char* suite, const char* name, TestFunction function) {
    Registry().push_back(TestSpec{suite, name, function});
}

void AddFailure(const char* file, int line, const std::string& message, bool fatal) {
    std::cerr << file << ":" << line << ": Failure\n" << message << "\n";
    ++CurrentTestFailureCount();
    if (fatal) {
        throw FatalTestFailure();
    }
}

void InitGoogleTest(int* argc, char** argv) {
    (void)argc;
    (void)argv;
}

int RunAllTests() {
    FailedTestCount() = 0;

    std::cout << "[==========] Running " << Registry().size() << " tests" << std::endl;
    for (const auto& test : Registry()) {
        CurrentTestFailureCount() = 0;
        std::cout << "[ RUN      ] " << test.suite << "." << test.name << std::endl;

        try {
            test.function();
        } catch (const FatalTestFailure&) {
            // Already reported.
        } catch (const std::exception& ex) {
            std::cerr << "Unhandled exception: " << ex.what() << std::endl;
            ++CurrentTestFailureCount();
        } catch (...) {
            std::cerr << "Unhandled unknown exception" << std::endl;
            ++CurrentTestFailureCount();
        }

        if (CurrentTestFailureCount() == 0) {
            std::cout << "[       OK ] " << test.suite << "." << test.name << std::endl;
        } else {
            std::cout << "[  FAILED  ] " << test.suite << "." << test.name << " ("
                      << CurrentTestFailureCount() << " failures)" << std::endl;
            ++FailedTestCount();
        }
    }

    std::cout << "[==========] " << Registry().size() << " tests ran" << std::endl;
    if (FailedTestCount() == 0) {
        std::cout << "[  PASSED  ] All tests passed" << std::endl;
    } else {
        std::cout << "[  FAILED  ] " << FailedTestCount() << " tests" << std::endl;
    }

    return FailedTestCount();
}

}  // namespace testing
