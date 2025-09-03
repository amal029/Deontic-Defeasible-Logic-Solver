#include "ddlsolver.hpp"
#include <iostream>

int main() {
  // Variables for the predicates
  Variable x{"x"}, y{"y"}, z{"z"};
  // Predicates themseleves
  Predicate ip1{"p1", {x, y}};
  Predicate ip2{"p2", {x, z}};
  Predicate ip3{"p3", {x, y, z}};
  // Atoms
  Atom A{"A"}, B{"B"}, C{"C"}, D{"D"}, E{"E"};
  // The predicates in the fact
  Formula fpp1{{"p1", {A, B}}};
  Formula fpp11{{"p1", {C, D}}};
  Formula fpp2{{"p2", {C, E}}};
  // Now add the implication to the rule table
  RuleTbl KB{};
  KB.insert(1, {{ip1, ip2}, ip3});

  // The solver
  Solver s{{&fpp1, &fpp11, &fpp2}, &KB, {}};
  std::cout << s.toString() << "\n";

  // Performing backward inference with a given goal.
  if (s.check(Predicate{"p3", {C, D, E}}))
    std::cout << s.toString() << "\n";

  // Performing forward inference without a given goal.
  if (s.check())
    std::cout << s.toString() << "\n";

  return 0;
}
