#include "ddlsolver.hpp"

int main() {
  // First declare all the literals
  Atom A{"A"}, B{"B"}, C{"C"}, D{"D"}, L{"L"}, M{"M"}, T{"T"}, X{"X"};
  // Now the nots if any
  PNot ND{D}, NX{X};
  // An example of a deontic
  OBL ob = OBL(B), od = OBL(D);

  // Declare the formulas to be used
  Formula ra1{Formula{A}};
  Formula ra2 {Formula{ob}};
  Formula ra3{Formula{ob}}, ra4{Formula{C}};
  Formula ra5{L};
  Formula ra6{A};
  Formula ra7{M};
  Formula ra8{od};
  Formula ra9{T};
  Formula ra10{D};
  Formula ra11{ND};
  Formula ra12{L};
  
  // Now the implications (defeasible rules)
  Implication r1{{ra1}, Formula(ob)};
  Implication r2{{ra2}, Formula{C}};
  Implication r3{{ra3, ra4}, Formula{D}};
  Implication r4{{ra5}, Formula{M}};
  Implication r5{{ra6, ra7, ra8}, Formula{T}};
  Implication r6{{ra9}, Formula{ND}};
  Implication r7{{ra10}, Formula{NX}};
  Implication r8{{ra11}, Formula{X}};
  Implication r9{{ra12}, Formula{od}};

  // Now the rule table
  RuleTbl rules{};
  rules.insert(1, std::move(r1));
  rules.insert(2, std::move(r2));
  rules.insert(3, std::move(r3));
  rules.insert(4, std::move(r4));
  rules.insert(5, std::move(r5));
  rules.insert(6, std::move(r6));
  rules.insert(7, std::move(r7));
  rules.insert(8, std::move(r8));
  rules.insert(9, std::move(r9));

  // Now the facts
  Formula f1{A}, f2{L};
  Solver s{{&f1, &f2}, &rules, {{6, 1}}};
  std::cout << s.toString() << "\n";
  // Now we can check if the rules are satisfied
  if (s.check())
    // Get the model after checking
    std::cout << s.toString() << "\n";
  // Try to copy a formula
  Formula ra1c = ra1;
  std::cout << ra1c.toString() << " " << ra1.toString() << "\n";
  // Substitute the proposition with another proposition
  Formula ra2c{ra1c.substituteProp("Z")};
  std::cout << ra2c.toString() << " " << ra1c.toString() << " "
            << ra1.toString() << "\n";
  return 0;
}
