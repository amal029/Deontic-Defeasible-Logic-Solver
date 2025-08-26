#include "ddlsolver.hpp"
#include <algorithm>

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
  auto ret = s1.build_and_or_tree(P2);
  std::cout << "Theorems that need to hold for KB1\n";
  for (const auto &x : ret) {
    std::cout << "[";
    std::for_each(x.begin(), x.end(),
                  [](const Formula &x) { std::cout << x.toString() << " "; });
    std::cout << "]\n";
  }

  // // Make another knowledge base
  RuleTbl KB2{};
  KB2.insert(1, {{P1, P0}, P2});
  Solver s2{{}, &KB2, {}};

  std::cout << "Theorems that need to hold for KB2\n";
  ret = s2.build_and_or_tree(P2);
  for (const auto &x : ret) {
    std::cout << "[";
    std::for_each(x.begin(), x.end(),
                  [](const Formula &x) { std::cout << x.toString() << " "; });
    std::cout << "]\n";
  }

  return 0;
}
