#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <variant>
#include <vector>

// The literal atom
class Literal {
public:
  Literal(std::string atom) : atom(atom) {}
  Literal(const Literal &) = delete;
  Literal(Literal &&) = default;
  ~Literal() {}
  std::string toString() const { return atom; }

private:
  std::string atom;
};

// The not for literals
class PNot {
public:
  PNot(Literal *lit) : l(lit) {}
  PNot(const PNot &) = delete;
  PNot(PNot &&) = default;
  ~PNot() {}
  std::string toString() const { return "(Not " + l->toString() + ")"; }
  Literal *getLiteral() const { return l; }

private:
  Literal *l;
};

// The obligation deontic
class OBL {
public:
  OBL(Literal *l) : pformula(l) {};
  OBL(PNot *n) : pformula(n) {};
  OBL(const OBL &) = delete;
  OBL(OBL &&) = default;
  ~OBL() {}
  std::string toString() const {
    if (const auto p = std::get_if<Literal *>(&pformula)) {
      return "(O " + (**p).toString() + ")";
    } else {
      return "(O " + std::get<PNot *>(pformula)->toString() + ")";
    }
  }

private:
  std::variant<Literal *, PNot *> pformula;
};

class DNot {
public:
  DNot(OBL *o) : o(o) {}
  DNot(const DNot &) = delete;
  DNot(DNot &&) = default;
  ~DNot() {}
  std::string toString() const { return "(Not " + o->toString() + ")"; }
  OBL *getOBL() const { return o; }

private:
  OBL *o;
};

class Formula {
public:
  Formula(Literal *l) : formula(l) {}
  Formula(PNot *l) : formula(l) {}
  Formula(OBL *l) : formula(l) {}
  Formula(DNot *l) : formula(l) {}
  Formula(const Formula &) = delete;
  Formula(Formula &&) = default;
  ~Formula() {}
  std::string toString() const {
    std::string ss;
    if (const auto l = std::get_if<Literal *>(&formula)) {
      ss = (**l).toString();
    } else if (const auto p = std::get_if<PNot *>(&formula)) {
      ss = (**p).toString();
    } else if (const auto o = std::get_if<OBL *>(&formula)) {
      ss = (**o).toString();
    } else if (const auto d = std::get_if<DNot *>(&formula)) {
      ss = (**d).toString();
    }
    return ss;
  }

  // Is this a complement.
  bool isNot() const {
    return std::holds_alternative<PNot *>(formula) ||
           std::holds_alternative<DNot *>(formula);
  }
  bool isLiteral() const { return std::holds_alternative<Literal *>(formula); }
  bool isOBL() const { return std::holds_alternative<OBL *>(formula); }
  std::string getComplement() const {
    if (this->isLiteral()) {
      return PNot(std::get<Literal *>(formula)).toString();
    } else if (this->isOBL()) {
      return DNot(std::get<OBL *>(formula)).toString();
    }
    return "";
  }
  std::string getComplementInner() const {
    if (std::holds_alternative<PNot *>(formula)) {
      return Formula{std::get<PNot *>(formula)->getLiteral()}.toString();
    } else if (std::holds_alternative<DNot *>(formula)) {
      return Formula{std::get<DNot *>(formula)->getOBL()}.toString();
    }
    return "";
  }

private:
  std::variant<Literal *, PNot *, OBL *, DNot *> formula;
};

// The hash for the Formula
// struct FormulaHash {
//   std::size_t operator()(const Formula *&s) const noexcept {
//     return std::hash<std::string>{}(s->toString());
//   }
// };

// // Equality of formulas
// struct FormulaEq {
//   bool operator()(const Formula *&lhs, const Formula *&rhs) const {
//     return lhs->toString() == rhs->toString();
//   }
// };

// This is the defeasible rule
using Antecedent = std::unordered_set<Formula *>;
class Implication {
public:
  Implication(Antecedent &&antecedents, Formula &&l)
      : antecedents(std::move(antecedents)), consequent(std::move(l)) {}
  ~Implication() {}
  Implication(const Implication &) = delete;
  Implication(Implication &&) = default;
  const Formula *getConsequent() const { return &consequent; }
  const Antecedent &getAntecedents() const { return antecedents; }
  const bool getDone() const { return done; }
  void setDone(bool val = true) { done = val; }
  std::string toString() const {
    std::string toret = "{";
    for (auto it = antecedents.cbegin(); it != antecedents.cend(); ++it)
      toret += (*it)->toString() + " ";
    toret += "} => " + consequent.toString();
    return toret;
  }

private:
  Antecedent antecedents;
  Formula consequent;
  bool done = false;
};

// The rules type (hashtable rule# -> Formula)
class RuleTbl {
public:
  RuleTbl() {}
  RuleTbl(const RuleTbl &) = delete;
  RuleTbl(RuleTbl &&) = delete;
  RuleTbl(std::unordered_map<uint64_t, Implication> &&rules)
      : rules(std::move(rules)) {}
  void insert(uint64_t rnum, Implication &&rule) {
    assert(rnum > 0);
    auto res = rules.find(rnum);
    if (res != rules.end()) {
      std::cout << "Rule already in the rule table\n";
      assert(false);
    }
    rules.insert({rnum, std::move(rule)});
  }
  size_t size() const { return rules.size(); }
  auto begin() { return rules.begin(); }
  auto end() { return rules.end(); }
  std::string toString() {
    std::string toret = "{\n";
    for (auto it = rules.cbegin(); it != rules.cend(); ++it) {
      toret += std::to_string(it->first) + " : " + it->second.toString() + "\n";
    }
    toret += "};";
    return toret;
  }
  ~RuleTbl() {}

private:
  std::unordered_map<uint64_t, Implication> rules;
};

using Conc = std::pair<const Formula *, uint64_t>;

// // The hash for the pair
// struct ConcHash {
//   std::size_t operator()(const Conc &s) const noexcept {
//     auto fhash = std::hash<std::string>{}(s.first->toString());
//     auto shash = std::hash<int64_t>{}(s.second);
//     return fhash ^ (shash << 1);
//   }
// };

// // Equality of formulas
// struct ConcEq {
//   bool operator()(const Conc &lhs, const Conc &rhs) const {
//     bool l = lhs.first->toString() == rhs.first->toString();
//     bool r = lhs.second == rhs.second;
//     return l && r;
//   }
// };

enum class HL { HIGHER = 0, LOWER = 1 };

class Solver {
public:
  Solver(std::vector<const Formula *> &&facts, RuleTbl *tbl,
         std::vector<std::pair<uint64_t, uint64_t>> &&precedence)
      : ruletbl(tbl), facts(std::move(facts)),
        precedence(std::move(precedence)) {}
  bool check() {
    bool ret = true;
    // Run the algorithm here
    // 1. Push the facts into the processing queue
    for (const auto &x : facts) {
      processing.push_back({x, 0});
    }
#ifdef DEBUG
    std::cout << this->toString() << "\n";
#endif
    // 2. Now loop until the queue is done
    while (!processing.empty()) {
      Conc &temp = processing.front();
      processing.pop_front(); // remove the rule from processing.

      // Check if this rule has a higher precedence rule in the
      // conclusions already
      auto iter = check_precedence(std::get<1>(temp), HL::HIGHER);
      if (iter != conclusions.begin()) {
#ifdef DEBUG
        std::cout << "Found a higher precedence rule in conclusions\n";
        std::cout << "Not adding " << temp.first->toString()
                  << " to conclusions \n";
#endif
        // We have a higher precedence rule in the conclusions, so don't
        // add this one.
        continue;
      }
      iter = check_precedence(std::get<1>(temp), HL::LOWER);
      if (iter != conclusions.begin()) {
#ifdef DEBUG
        std::cout << "Found a lower precedence conclusion " << "\n";
#endif
        // a. First remove everything from start of conclusions to this
        // iterator and then move them to the remove list
        for (auto it = conclusions.begin(); it != iter; ++it)
          remove.push_back(*it);
#ifdef DEBUG
        std::cout << this->toString() << "\n";
#endif
        // Now remove all these from the conclusions.
        conclusions.erase(conclusions.begin(), iter);
#ifdef DEBUG
        std::cout << "after erasing from conclusions: " << this->toString()
                  << "\n";
#endif
        // backward removal
        backward_removal();
      }
      // Check if there is a contradiction before adding this conclusion.
      if (check_contradiction(temp.first)) {
        return false;
      }
      // Finally we can add this rule to the conclusions and trigger the
      // required implications.
      conclusions.push_back(std::move(temp));
#ifdef DEBUG
      std::cout << this->toString() << "\n";
#endif
      // Now enable all the possible implications
      for (auto it = ruletbl->begin(); it != ruletbl->end(); ++it) {
        // Now check if all antecedents are satisfied.
        if (!it->second.getDone() && check_antecedents(std::get<1>(*it))) {
          std::get<1>(*it).setDone();
          processing.push_back(
              {std::get<1>(*it).getConsequent(), std::get<0>(*it)});
        }
      }
    }
    return ret;
  }
  bool check_contradiction(const Formula *v) {
    std::string fstring;
    if (v->isNot()) {
      fstring = v->getComplementInner();
    } else if (v->isLiteral() || v->isOBL()) {
      fstring = v->getComplement();
    }
    // Now check if the conclusion has this fstring
    auto iter = std::find_if(
        conclusions.begin(), conclusions.end(),
        [&fstring](const auto &x) { return x.first->toString() == fstring; });
    return iter != conclusions.end();
  }

  std::string toString() const {
    std::string toret = "model(\n";
    toret += "Rule table: \n" + ruletbl->toString() + "\n";
    toret += "Facts: \n";
    std::for_each(facts.cbegin(), facts.cend(),
                  [&toret](const Formula *x) { toret += x->toString() + " "; });
    toret += "\nConclusions: ";
    std::for_each(
        conclusions.cbegin(), conclusions.cend(), [&toret](const Conc &x) {
          toret += x.first->toString() + " " + std::to_string(x.second) + ", ";
        });
    toret += "\nProcessing: ";
    std::for_each(
        processing.cbegin(), processing.cend(), [&toret](const Conc &x) {
          toret += x.first->toString() + " " + std::to_string(x.second) + ", ";
        });
    toret += "\nRemove: ";
    std::for_each(remove.cbegin(), remove.cend(), [&toret](const Conc &x) {
      toret += x.first->toString() + " " + std::to_string(x.second) + ", ";
    });
    toret += "\nPrecendence: ";
    std::for_each(precedence.cbegin(), precedence.cend(),
                  [&toret](const auto &x) {
                    toret += std::to_string(x.first) + ">" +
                             std::to_string(x.second) + ", ";
                  });

    toret += "\n)";
    return toret;
  }
  ~Solver() {}

private:
  // Removing from conclusions and processing ad long as remove exists
  void backward_removal() {
    while (!remove.empty()) {
      const auto &[temp, rnum] = remove.front();
      const std::string stemp{temp->toString()};
      remove.pop_front();
      // Get all the Formulas in conclusions/processing that are
      // triggerd by temp.
      for (auto it = ruletbl->begin(); it != ruletbl->end(); ++it) {
        Implication &trule = std::get<1>(*it);
        if (trule.getDone()) {
          // Loop through the antecedents and check if this rule has
          // temp.
#ifdef DEBUG
          std::cout << "checking in antecedent formula: " << stemp << "\n";
          std::cout << "rule number being checked: "
                    << std::to_string(it->first) << "\n";
#endif
          auto fit = std::find_if(
              trule.getAntecedents().begin(), trule.getAntecedents().end(),
              [&stemp](const Formula *x) { return x->toString() == stemp; });

          if (fit != trule.getAntecedents().end()) {
            // We have the rule. Now get iterator from the conclusion
            // and processing to remove it from there.
            uint64_t ruleNum = std::get<0>(*it);
            std::string consequent =
                std::get<1>(*it).getConsequent()->toString();
#ifdef DEBUG
            std::cout << "Trying to remove the consequent:" << consequent
                      << "\n";
#endif
            // Get the iterator from the conclusions
            auto cit = std::partition(conclusions.begin(), conclusions.end(),
                                      [&ruleNum, &consequent](const auto &x) {
                                        return (std::get<1>(x) == ruleNum &&
                                                std::get<0>(x)->toString() ==
                                                    consequent);
                                      });
#ifdef DEBUG
            std::string toret = "\nConclusions: ";
            std::for_each(conclusions.cbegin(), conclusions.cend(),
                          [&toret](const Conc &x) {
                            toret += x.first->toString() + " " +
                                     std::to_string(x.second) + ", ";
                          });
            std::cout << toret << "\n";
#endif
            if (cit != conclusions.begin()) {
              for (auto ccit = conclusions.begin(); ccit != cit; ++ccit)
                remove.push_back(*ccit);
              conclusions.erase(conclusions.begin(), cit);
            } else {
              // Remove from the processing queue if not in the
              // conclusion
              auto cit = std::partition(processing.begin(), processing.end(),
                                        [&ruleNum, &consequent](const auto &x) {
                                          return (std::get<1>(x) == ruleNum &&
                                                  std::get<0>(x)->toString() ==
                                                      consequent);
                                        });
              if (cit != processing.begin()) {
                for (auto ccit = processing.begin(); ccit != cit; ++ccit)
                  remove.push_back(*ccit);
                processing.erase(processing.begin(), cit);
              }
            }
#ifdef DEBUG
            toret = "\n Remove: ";
            std::for_each(remove.cbegin(), remove.cend(),
                          [&toret](const Conc &x) {
                            toret += x.first->toString() + " " +
                                     std::to_string(x.second) + ", ";
                          });
            std::cout << toret << "\n";

            toret = "\nConclusions: ";
            std::for_each(conclusions.cbegin(), conclusions.cend(),
                          [&toret](const Conc &x) {
                            toret += x.first->toString() + " " +
                                     std::to_string(x.second) + ", ";
                          });
            std::cout << toret << "\n";
#endif
            // FIXME: Remove it from processing queue too -- later
          }
        }
      }
    }
  }

  // Checking antecedent
  bool check_antecedents(const Implication &rule) {
    std::vector<std::string> cons;
    for (const auto &x : conclusions) {
      cons.push_back(std::get<0>(x)->toString());
    }
    size_t counter = 0;
    for (const auto &x : rule.getAntecedents()) {
      if (std::find(cons.cbegin(), cons.cend(), x->toString()) != cons.end())
        counter++;
    }
    return (counter == rule.getAntecedents().size());
  }

  // Checking if the conclusion already has a higher order rule
  std::vector<Conc>::iterator check_precedence(uint64_t rNum, HL lh) {
#ifdef DEBUG
    std::cout << "check precedence: " << rNum << "\n";
#endif
    std::vector<uint64_t> hls;
    for (size_t i = 0; i < precedence.size(); ++i) {
      if (lh == HL::HIGHER && (std::get<1>(precedence[i]) == rNum)) {
        // Getting hls
        hls.push_back(std::get<0>(precedence[i]));
      } else if (lh == HL::LOWER && (std::get<0>(precedence[i]) == rNum)) {
        // Getting lowers
        hls.push_back(std::get<1>(precedence[i]));
      }
    }
#ifdef DEBUG
    for (auto it = hls.cbegin(); it != hls.cend(); ++it)
      std::cout << *it << "\n";
#endif
    // Now check if there are any!
    if (hls.empty()) {
      return conclusions.begin();
    } else {
      // Get the projected rule numbers from conclusions
      // Now we need to check if the hls are in the conclusions set.
      auto res = std::partition(conclusions.begin(), conclusions.end(),
                                [&hls](auto &x) {
                                  return std::find(hls.begin(), hls.end(),
                                                   std::get<1>(x)) != hls.end();
                                });
      return res; // we have found a higher one
    }
  }
  RuleTbl *ruletbl;
  std::vector<const Formula *> facts;
  std::vector<std::pair<uint64_t, uint64_t>> precedence;
  std::vector<Conc> conclusions;
  std::deque<Conc> processing;
  std::deque<Conc> remove;
};

int main() {
  // First declare all the literals
  Literal A{"A"}, B{"B"}, C{"C"}, D{"D"}, L{"L"}, M{"M"}, T{"T"}, X{"X"};
  // Now the nots if any
  PNot ND{&D}, NX{&X};
  // An example of a deontic
  OBL ob = OBL(&B), od = OBL(&D);

  // Declare the formulas to be used
  Formula ra1 = Formula{&A};
  Formula ra2 = Formula{&ob};
  Formula ra3 = Formula{&ob}, ra4 = Formula{&C};
  Formula ra5{&L};
  Formula ra6{&A};
  Formula ra7{&M};
  Formula ra8{&od};
  Formula ra9{&T};
  Formula ra10{&D};
  Formula ra11{&ND};
  Formula ra12{&L};
  
  // Now the implications (defeasible rules)
  Implication r1{{&ra1}, Formula(&ob)};
  Implication r2{{&ra2}, Formula{&C}};
  Implication r3{{&ra3, &ra4}, Formula{&D}};
  Implication r4{{&ra5}, Formula{&M}};
  Implication r5{{&ra6, &ra7, &ra8}, Formula{&T}};
  Implication r6{{&ra9}, Formula{&ND}};
  Implication r7{{&ra10}, Formula{&NX}};
  Implication r8{{&ra11}, Formula{&X}};
  Implication r9{{&ra12}, Formula{&od}};

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
  Formula f1{&A}, f2{&L};
  Solver s{{&f1, &f2}, &rules, {{6, 1}}};
  std::cout << s.toString() << "\n";
  // Now we can check if the rules are satisfied
  if (s.check())
    // Get the model after checking
    std::cout << s.toString() << "\n";
  return 0;
}
