#include "ddlsolver.hpp"

int main() {
  Variable X{"X"};
  Predicate P0("P0", {X});
  Predicate P1("P1", {X});
  Predicate P2("P2", {X});
  RuleTbl KB{};
  KB.insert(1, {{P0}, P1});
  KB.insert(2, {{P1}, P2});
  Solver s1{{}, &KB, {}};
  
  // The above two should be logically equivalent
  auto ret = s1.backward_chain(P2);
  std::cout << "Theorems that need to hold for KB1\n";
  for (const auto &x : ret) {
    std::cout << "[";
    std::for_each(x.cbegin(), x.cend(),
                  [](const Formula &x) { std::cout << x.toString() << " "; });
    std::cout << "]\n";
  }

  // Make another knowledge base
  RuleTbl KB2{};
  KB2.insert(1, {{P1, P0}, P2});
  Solver s2{{}, &KB2, {}};

  std::cout << "Theorems that need to hold for KB2\n";
  ret = s2.backward_chain(P2);
  for (const auto &x : ret) {
    std::cout << "[";
    std::for_each(x.cbegin(), x.cend(),
                  [](const Formula &x) { std::cout << x.toString() << " "; });
    std::cout << "]\n";
  }

  // Example 2
  Variable A{"A"};
  Predicate A1("A1", {A});
  Predicate A2("A2", {A});
  Predicate A3("A3", {A});
  Predicate A4("A4", {A});
  RuleTbl KB3{};
  KB3.insert(1, {{A4}, A3});
  KB3.insert(2, {{A2}, A3});
  KB3.insert(3, {{A1}, A2});

  Solver s3{{}, &KB3, {}};
  ret = s3.backward_chain(A3);
  std::cout << "Theorems that need to hold for KB3\n";
  for (const auto &x : ret) {
    std::cout << "[";
    std::for_each(x.cbegin(), x.cend(),
                  [](const Formula &x) { std::cout << x.toString() << " "; });
    std::cout << "]\n";
  }

  RuleTbl KB4{};
  KB4.insert(1, {{A1, A4}, A3});
  Solver s4{{}, &KB4, {}};
  ret = s4.backward_chain(A3);
  std::cout << "Theorems that need to hold for KB4\n";
  for (const auto &x : ret) {
    std::cout << "[";
    std::for_each(x.cbegin(), x.cend(),
                  [](const Formula &x) { std::cout << x.toString() << " "; });
    std::cout << "]\n";
  }

  // Example 3
  RuleTbl KB5{};
  KB5.insert(1, {{A2}, A1});
  KB5.insert(2, {{A4}, PNot(A1)});
  Solver s5{{}, &KB5, {{2, 1}}};
  ret = s5.backward_chain(A1);
  std::cout << "Theorems that need to hold for KB5\n";
  for (const auto &x : ret) {
    std::cout << "[";
    std::for_each(x.cbegin(), x.cend(),
                  [](const Formula &x) { std::cout << x.toString() << " "; });
    std::cout << "]\n";
  }

  RuleTbl KB6{};
  KB6.insert(1, {{A2, PNot(A4)}, A1});
  Solver s6{{}, &KB6, {}};
  ret = s6.backward_chain(A1);
  std::cout << "Theorems that need to hold for KB6\n";
  for (const auto &x : ret) {
    std::cout << "[";
    std::for_each(x.cbegin(), x.cend(),
                  [](const Formula &x) { std::cout << x.toString() << " "; });
    std::cout << "]\n";
  }

  return 0;
}
