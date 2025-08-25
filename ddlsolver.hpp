#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <iterator>
#include <numeric>
#include <string>
#include <sys/types.h>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <variant>
#include <vector>

// The literal atom
class Atom {
public:
  Atom(const char *atom) : atom(atom) {}
  Atom(const Atom &) = default;
  Atom(Atom &&) noexcept = default;
  Atom &operator=(const Atom &other) = default;
  Atom &operator=(Atom &&other) = default;
  ~Atom() {}
  std::string toString() const { return atom; }

  Atom substituteProp(const Atom &x, const char *y) const {
    if (x.toString() == this->atom) {
      return Atom(y);
    } else
      return *this;
  }

private:
  std::string atom;
};

// We need the variable class/type
class Variable {
public:
  Variable(const char *name) : name(name) {}
  ~Variable() {}
  std::string toString() const { return "Var(" + name + ")"; }
  const std::string &getName() const { return name; }

private:
  std::string name;
};

// Now need the predicate class
using PredicateType = std::variant<Variable, Atom>;
class Predicate {
public:
  Predicate(std::string &&name, std::initializer_list<PredicateType> &&args) noexcept
    : args(std::move(args)), name(std::move(name)) {}
  Predicate(std::string name, std::vector<PredicateType> &&args)
    : args(std::move(args)), name(std::move(name)) {}

  ~Predicate() {}
  std::string toString() const {
    std::string toret{};
    size_t counter = 0;
    for (auto it = args.cbegin(); it != args.cend(); ++it) {
      if (std::holds_alternative<Variable>(*it)) {
        toret += std::get<Variable>(*it).toString();
      } else {
        toret += std::get<Atom>(*it).toString();
      }
      if (counter < args.size() - 1)
        toret += ", ";
      counter += 1;
    }
    return name + "(" + toret + ")";
  }

  bool hasVariables() const {
    // accumulate the args inside this to check if there are any vars in
    // here!
    return std::accumulate(args.cbegin(), args.cend(), false,
                           [](const bool &y, const PredicateType &x) {
                             return y || std::holds_alternative<Variable>(x);
                           });
  }

  bool hasVariable(const Variable &tocheck) const {
    // accumulate the args inside this to check if there are any vars in
    // here!
    return std::accumulate(
			   args.cbegin(), args.cend(), false,
			   [&tocheck](const bool &y, const PredicateType &x) {
			     return ((std::holds_alternative<Variable>(x)) &&
				     (std::get<Variable>(x).getName() == tocheck.getName())) ||
			       y;
			   });
  }

  Predicate subsVartoAtom(const Variable &x, const char *y) const {
    std::vector<PredicateType> rargs{args};
    for (size_t counter = 0; counter < rargs.size(); ++counter) {
      if (std::holds_alternative<Variable>(rargs[counter]) &&
          std::get<Variable>(rargs[counter]).getName() == x.getName()) {
        rargs[counter] = Atom(y);
      }
    }
    Predicate toret{name, std::move(rargs)};
#ifdef DEBUG
    std::cout << "return predicate " << toret.toString() << "\n";
#endif
    return toret;
  }

  Predicate substituteProp(const Atom &x, const char *y) const {
    std::vector<PredicateType> rargs;
    rargs.reserve(args.size());
    std::transform(args.begin(), args.end(), rargs.begin(),
                   [&x, &y](const PredicateType &z) -> PredicateType {
                     if (std::holds_alternative<Atom>(z)) {
                       return std::get<Atom>(z).substituteProp(x, y);
                     } else
                       return z;
                   });
    return Predicate{name, std::move(rargs)};
  }

  const std::string &getName() const { return name; }
  size_t getArity() const { return args.size(); }

  bool getMatchingPredicate(const Predicate &x) const {
    return (x.getName() == name) && (x.getArity() == getArity());
  }
  const std::vector<PredicateType> &getArgs() const { return args; }

private:
  std::vector<PredicateType> args;
  std::string name;
};

using LiteralType = std::variant<Atom, Predicate>;
// The not for literals
class PNot {
public:
  PNot(Atom lit) : l(lit) {}
  PNot(Predicate lit) : l(lit) {}
  PNot(const PNot &) = default;
  PNot &operator=(const PNot &) noexcept = default;
  PNot(PNot &&) noexcept = default;
  ~PNot() {}
  std::string toString() const {
    if (std::holds_alternative<Atom>(l)) {
      return std::get<Atom>(l).toString();
    } else
      return std::get<Predicate>(l).toString();
  }
  const LiteralType &getLiteral() const { return l; }
  PNot substituteProp(const Atom &x, const char *y) const {
    if (std::holds_alternative<Atom>(l)) {
      return PNot(std::get<Atom>(l).substituteProp(x, y));
    } else
      return PNot(std::get<Predicate>(l).substituteProp(x, y));
  }

  bool hasVariable(const Variable &tocheck) const {
    return (std::holds_alternative<Predicate>(l))
               ? std::get<Predicate>(l).hasVariable(tocheck)
               : false;
  }

  bool hasVariables() const {
    return (std::holds_alternative<Predicate>(l))
               ? std::get<Predicate>(l).hasVariables()
               : false;
  }

  Predicate getPredicate() const {
    if (std::holds_alternative<Predicate>(l)) {
      return std::get<Predicate>(l);
    } else
      assert(false);
  }

  bool getMatchingPredicate(const Predicate &x) const {
    return (std::holds_alternative<Predicate>(l))
               ? std::get<Predicate>(l).getMatchingPredicate(x)
               : false;
  }

  PNot subsVartoAtom(const Variable &x, const char *y) const {
    return (std::holds_alternative<Atom>(l))
               ? *this
               : PNot(std::get<Predicate>(l).subsVartoAtom(x, y));
  }

private:
  LiteralType l;
};

// The obligation deontic
class OBL {
public:
  OBL(Atom l) : pformula(l) {};
  OBL(PNot n) : pformula(n) {};
  OBL(const OBL &) noexcept = default;
  OBL &operator=(const OBL &) = default;
  OBL(OBL &&) = default;
  ~OBL() {}
  std::string toString() const {
    if (const auto p = std::get_if<Atom>(&pformula)) {
      return "(O " + (*p).toString() + ")";
    } else {
      return "(O " + std::get<PNot>(pformula).toString() + ")";
    }
  }

  OBL substituteProp(const Atom &x, const char *&y) const {
    return (std::holds_alternative<Atom>(pformula))
               ? OBL(std::get<Atom>(pformula).substituteProp(x, y))
               : OBL(std::get<PNot>(pformula).substituteProp(x, y));
  }

  bool hasVariable(const Variable &tocheck) const {
    return (std::holds_alternative<PNot>(pformula))
               ? std::get<PNot>(pformula).hasVariable(tocheck)
               : false;
  }

  bool hasVariables() const {
    return (std::holds_alternative<PNot>(pformula))
               ? std::get<PNot>(pformula).hasVariables()
               : false;
  }

  Predicate getPredicate() const {
    if (std::holds_alternative<PNot>(pformula)) {
      return std::get<PNot>(pformula).getPredicate();
    } else
      assert(false);
  }

  bool getMatchingPredicate(const Predicate &x) const {
    return (std::holds_alternative<PNot>(pformula))
               ? std::get<PNot>(pformula).getMatchingPredicate(x)
               : false;
  }

  OBL subsVartoAtom(const Variable &x, const char *y) const {
    return (std::holds_alternative<Atom>(pformula))
               ? *this
               : OBL(std::get<PNot>(pformula).subsVartoAtom(x, y));

    // if (std::holds_alternative<PNot>(pformula)) {
    //   std::get<PNot>(pformula).subsVartoAtom(x, y);
    // }
  }

private:
  std::variant<Atom, PNot> pformula;
};

class DNot {
public:
  DNot(OBL o) : o(o) {}
  DNot(const DNot &) = default;
  DNot &operator=(const DNot &) noexcept = default;
  DNot(DNot &&) noexcept = default;
  ~DNot() {}
  std::string toString() const { return "(Not " + o.toString() + ")"; }
  const OBL &getOBL() const { return o; }

  DNot substituteProp(const Atom &x, const char *y) const {
    return o.substituteProp(x, y);
  }

  bool hasVariables() const { return o.hasVariables(); }
  bool hasVariable(const Variable &tocheck) const {
    return o.hasVariable(tocheck);
  }
  bool getMatchingPredicate(const Predicate &x) const {
    return o.getMatchingPredicate(x);
  }
  DNot subsVartoAtom(const Variable &x, const char *y) const {
    return o.subsVartoAtom(x, y);
  }

  Predicate getPredicate() const { return o.getPredicate(); }

private:
  OBL o;
};

class Formula {
public:
  Formula(Atom l) : formula(std::move(l)) {}
  Formula(Predicate l) : formula(std::move(l)) {}
  Formula(PNot l) : formula(std::move(l)) {}
  Formula(OBL l) : formula(std::move(l)) {}
  Formula(DNot l) : formula(std::move(l)) {}
  Formula(const Formula &) = default;
  Formula &operator=(const Formula &) noexcept = default;
  Formula(Formula &&) noexcept = default;
  ~Formula() {}
  std::string toString() const {
    std::string ss;
    if (const auto l = std::get_if<Atom>(&formula)) {
      ss = (*l).toString();
    } else if (const auto p = std::get_if<PNot>(&formula)) {
      ss = (*p).toString();
    } else if (const auto o = std::get_if<OBL>(&formula)) {
      ss = (*o).toString();
    } else if (const auto d = std::get_if<DNot>(&formula)) {
      ss = (*d).toString();
    } else if (const auto d = std::get_if<Predicate>(&formula)) {
      ss = (*d).toString();
    }
    return ss;
  }

  // Get the predicate version of this formula
  Predicate getPredicate() const {
    if (std::holds_alternative<Predicate>(formula))
      return std::get<Predicate>(formula);
    else if (std::holds_alternative<PNot>(formula)) {
      return std::get<PNot>(formula).getPredicate();
    } else if (std::holds_alternative<OBL>(formula)) {
      return std::get<OBL>(formula).getPredicate();
    } else if (std::holds_alternative<DNot>(formula)) {
      return std::get<DNot>(formula).getPredicate();
    } else {
      assert(false);
    }
  }

  Formula substituteProp(const Atom &x, const char *y) const {
    if (std::holds_alternative<Atom>(formula)) {
      return Formula(std::get<Atom>(formula).substituteProp(x, y));
    } else if (std::holds_alternative<PNot>(formula)) {
      return Formula(std::get<PNot>(formula).substituteProp(x, y));
    } else if (std::holds_alternative<OBL>(formula)) {
      return Formula(std::get<OBL>(formula).substituteProp(x, y));
    } else if (std::holds_alternative<DNot>(formula)) {
      return Formula(std::get<DNot>(formula).substituteProp(x, y));
    } else {
      // This is the predicate case
      return Formula{std::get<Predicate>(formula).substituteProp(x, y)};
    }
  }

  // Is this a complement.
  bool isNot() const {
    return std::holds_alternative<PNot>(formula) ||
           std::holds_alternative<DNot>(formula);
  }
  bool isAtom() const { return std::holds_alternative<Atom>(formula); }
  bool isPredicate() const {
    return std::holds_alternative<Predicate>(formula);
  }
  bool isOBL() const { return std::holds_alternative<OBL>(formula); }

  Formula getComplementFormula() const {
    if (this->isAtom()) {
      return PNot(std::get<Atom>(formula));
    } else if (this->isOBL()) {
      return DNot(std::get<OBL>(formula));
    } else if (this->isPredicate()) {
      return PNot(std::get<Predicate>(formula));
    }
    else
      assert(false);
  }

  std::string getComplement() const {
    if (this->isAtom()) {
      return PNot(std::get<Atom>(formula)).toString();
    } else if (this->isOBL()) {
      return DNot(std::get<OBL>(formula)).toString();
    } else if (this->isPredicate()) {
      return PNot(std::get<Predicate>(formula)).toString();
    }
    return "";
  }

  Formula getComplementInnerFormula() const {
    if (std::holds_alternative<PNot>(formula)) {
      auto x = std::get<PNot>(formula).getLiteral();
      if (std::holds_alternative<Atom>(x)) {
        return std::get<Atom>(x);
      } else
        return std::get<Predicate>(x);
    } else if (std::holds_alternative<DNot>(formula)) {
      return std::get<DNot>(formula).getOBL();
    } else
      assert(false);
  }

  std::string getComplementInner() const {
    if (std::holds_alternative<PNot>(formula)) {
      auto x = std::get<PNot>(formula).getLiteral();
      if (std::holds_alternative<Atom>(x)) {
        return std::get<Atom>(x).toString();
      } else
        return std::get<Predicate>(x).toString();
    } else if (std::holds_alternative<DNot>(formula)) {
      return std::get<DNot>(formula).getOBL().toString();
    }
    return "";
  }

  bool hasVariable(const Variable &tocheck) const {
    if (std::holds_alternative<Atom>(formula)) {
      return false;
    } else if (std::holds_alternative<PNot>(formula)) {
      return std::get<PNot>(formula).hasVariable(tocheck);
    } else if (std::holds_alternative<OBL>(formula)) {
      return std::get<OBL>(formula).hasVariable(tocheck);
    } else if (std::holds_alternative<DNot>(formula)) {
      return std::get<DNot>(formula).hasVariable(tocheck);
    } else {
      return std::get<Predicate>(formula).hasVariable(tocheck);
    }
  }

  bool hasVariables() const {
    if (std::holds_alternative<Atom>(formula)) {
      return false;
    } else if (std::holds_alternative<PNot>(formula)) {
      return std::get<PNot>(formula).hasVariables();
    } else if (std::holds_alternative<OBL>(formula)) {
      return std::get<OBL>(formula).hasVariables();
    } else if (std::holds_alternative<DNot>(formula)) {
      return std::get<DNot>(formula).hasVariables();
    } else {
      return std::get<Predicate>(formula).hasVariables();
    }
  }

  bool getMatchingPredicate(const Predicate &x) const {
    if (std::holds_alternative<Atom>(formula)) {
      return false;
    } else if (std::holds_alternative<PNot>(formula)) {
      return std::get<PNot>(formula).getMatchingPredicate(x);
    } else if (std::holds_alternative<OBL>(formula)) {
      return std::get<OBL>(formula).getMatchingPredicate(x);
    } else if (std::holds_alternative<DNot>(formula)) {
      return std::get<DNot>(formula).getMatchingPredicate(x);
    } else {
      return std::get<Predicate>(formula).getMatchingPredicate(x);
    }
  }

  Formula subsVartoAtom(const Variable &x, const char *y) const {
    if (std::holds_alternative<PNot>(formula)) {
      return std::get<PNot>(formula).subsVartoAtom(x, y);
    } else if (std::holds_alternative<OBL>(formula)) {
      return std::get<OBL>(formula).subsVartoAtom(x, y);
    } else if (std::holds_alternative<DNot>(formula)) {
      return std::get<DNot>(formula).subsVartoAtom(x, y);
    } else if (std::holds_alternative<Predicate>(formula)) {
      return std::get<Predicate>(formula).subsVartoAtom(x, y);
    } else {
      return *this;
    }
  }

private:
  std::variant<Atom, PNot, OBL, DNot, Predicate> formula;
};

// The Formula hash
struct FormulaHash {
  std::size_t operator()(const Formula &s) const noexcept {
    return std::hash<std::string>{}(s.toString());
  }
};

// The Formula equality
struct FormulaEq {
  bool operator()(const Formula &lhs, const Formula &rhs) const {
    return lhs.toString() == rhs.toString();
  }
};

// This is the defeasible rule
using Antecedent = std::unordered_set<Formula, FormulaHash, FormulaEq>;
class Implication {
public:
  Implication(Antecedent &&antecedents, Formula &&l)
      : antecedents(std::move(antecedents)), consequent(std::move(l)) {}
  ~Implication() {}
  Implication(const Implication &) = default;
  Implication &operator=(const Implication &) noexcept = delete;
  Implication &operator=(Implication &&) noexcept = default;
  Implication(Implication &&) noexcept = default;
  const Formula *getConsequent() const { return &consequent; }
  const Antecedent &getAntecedents() const { return antecedents; }
  bool getDone() const { return done; }
  void setDone(bool val = true) { done = val; }
  std::string toString() const {
    std::string toret = "{";
    for (auto it = antecedents.cbegin(); it != antecedents.cend(); ++it)
      toret += (*it).toString() + " ";
    toret += "} => " + consequent.toString();
    return toret;
  }

  bool hasVariables() const {
    return std::accumulate(antecedents.cbegin(), antecedents.cend(), false,
                           [](const bool &val, const Formula &x) {
                             return val || x.hasVariables();
                           });
  }

  bool hasVariable(const Variable &tocheck) const {
    return std::accumulate(antecedents.cbegin(), antecedents.cend(), false,
                           [&tocheck](const bool &val, const Formula &x) {
                             return val || x.hasVariable(tocheck);
                           });
  }

  void getPredicateWithVar(std::vector<Predicate> &toret) const {
    auto it = std::find_if(antecedents.begin(), antecedents.end(),
                           [](const Formula &x) { return x.hasVariables(); });
    while (it != antecedents.end()) {
      toret.push_back(it->getPredicate());
      it = std::find_if(std::next(it), antecedents.end(),
                        [](const Formula &x) { return x.hasVariables(); });
    }
#ifdef DEBUG
    std::cout << "Got the predicates with vars in implication\n";
#endif
  }

  Implication subsVartoAtom(const Variable &v, const char *atom) const {
    // Go through all the antecedents and replace the variable with the
    // atom
    Antecedent rantecedents;
    for (auto it = antecedents.begin(); it != antecedents.end(); ++it) {
      if (it->hasVariable(v)) {
        rantecedents.insert(it->subsVartoAtom(v, atom));
      } else
        rantecedents.insert(*it);
    }
    Formula rcons = consequent.hasVariable(v)
                        ? Formula{consequent.subsVartoAtom(v, atom)}
                        : consequent;
    Implication toret{std::move(rantecedents), std::move(rcons)};
#ifdef DEBUG
    std::cout << "conclusion implication: " << this->toString() << "\n";
#endif
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
  const Implication &find(uint64_t key) { return rules.find(key)->second; }
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

enum class HL { HIGHER = 0, LOWER = 1 };

struct VariableHash {
  std::size_t operator()(const Variable &s) const noexcept {
    return std::hash<std::string>{}(s.toString());
  }
};

// The Formula equality
struct VariableEq {
  bool operator()(const Variable &lhs, const Variable &rhs) const {
    return lhs.toString() == rhs.toString();
  }
};

using VarMap = std::unordered_map<Variable, Atom, VariableHash, VariableEq>;

// Backward chaining for getting the facts (propositions and predicates)
// that are needed to satisfy the given goal.
struct Node {
  std::vector<size_t> edges;
  const Formula *f;
  bool fact;
  bool tofree = false;
  ~Node() {
    if (tofree)
      delete f;
  }
};

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
#ifdef DEBUG
      std::cout << "temp : " << this->toString() << "\n";
#endif
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
        if (!it->second.getDone() && it->second.hasVariables()) {
          // XXX: We perform unification here
          // Implication out{it->second}; // copy constructor call here
#ifdef DEBUG
          std::cout << "entering unify\n";
#endif
          if (unify(it->second)) {
            std::get<1>(*it)
                .setDone(); // setting the original implication as done
#ifdef DEBUG
            std::cout << "Adding the consequent: "
                      << it->second.getConsequent()->toString() << "\n";
#endif
            processing.push_back(
                {it->second.getConsequent(), std::get<0>(*it)});
          }
        }
        // Now check if all antecedents are satisfied.
        else if (!it->second.getDone() && check_antecedents(std::get<1>(*it))) {
          std::get<1>(*it).setDone();
          processing.push_back(
              {std::get<1>(*it).getConsequent(), std::get<0>(*it)});
        }
      }
    }
    return ret;
  }

  std::vector<Predicate> getConclusionPredicate(const Predicate &in,
                                                std::vector<Predicate> &toret) {
    auto it = std::find_if(
        conclusions.cbegin(), conclusions.cend(),
        [&in](const Conc &x) { return (x.first)->getMatchingPredicate(in); });
    while (it != conclusions.cend()) {
      toret.push_back(it->first->getPredicate());
      it =
          std::find_if(std::next(it), conclusions.cend(), [&in](const Conc &x) {
            return (x.first)->getMatchingPredicate(in);
          });
    }
    return toret;
  }

  // This can possibly be NP-hard
  bool unify(Implication &out) {
    bool toret = false;

    // Go through every predicate (with a variable) one by one
    std::vector<Predicate> ps;
    out.getPredicateWithVar(ps);
#ifdef DEBUG
    std::cout << "Predicates with Var in implication: \n";
    for_each(ps.begin(), ps.end(),
             [](const Predicate &x) { std::cout << x.toString() << "\n"; });
#endif
    if (ps.empty())
      return toret; // No predicates with variables to unify

    // Here we need to find the substitution for the variable(s)

    // 1. Get the predicate with same name and arity from the
    // conclusions.
    std::vector<std::vector<Predicate>> concps;
    concps.reserve(ps.size());
    for (const auto &x : ps) {
      std::vector<Predicate> res;
      getConclusionPredicate(x, res);
      if (res.empty())
        break; // Some predicate has no known fact
      else
        concps.push_back(res);
    }
#ifdef DEBUG
    std::cout << "Predicates in conclusions: \n";
    for_each(concps.begin(), concps.end(), [](const auto &x) {
      std::cout << "[";
      for_each(x.begin(), x.end(),
               [](const Predicate &y) { std::cout << y.toString() << " "; });
      std::cout << "]\n";
    });
#endif

    if (concps.size() != ps.size())
      return toret;

    // 2. Make the cartresian product if ps.size() > 1
    std::vector<std::vector<Predicate>> cartresianconcps;
    if (ps.size() > 1) {
      for (const auto &x : concps[0]) {
        cartresianconcps.push_back({x});
      }
      cartresianconcps =
          std::accumulate(concps.cbegin() + 1, concps.cend(), cartresianconcps,
                          [](const std::vector<std::vector<Predicate>> &f,
                             const std::vector<Predicate> &s)
                              -> std::vector<std::vector<Predicate>> {
                            std::vector<std::vector<Predicate>> res;
                            for (const Predicate &y : s) {
                              for (const std::vector<Predicate> &x : f) {
                                std::vector<Predicate> res1{x};
                                res1.push_back(y);
                                res.push_back(std::move(res1));
                              }
                            }
                            return res;
                          });
    } else
      cartresianconcps = std::move(concps);
#ifdef DEBUG
    std::cout << "Predicates in cartresian: \n";
    for_each(cartresianconcps.begin(), cartresianconcps.end(),
             [](const auto &x) {
               std::cout << "[";
               for_each(x.begin(), x.end(), [](const Predicate &y) {
                 std::cout << y.toString() << " ";
               });
               std::cout << "]\n";
             });
#endif
    // 3. Now we have the cartresian product of the conclusion
    // predicates. Now we can start performing substitution.
    for (const auto &x : cartresianconcps) {
      Implication temp{out};         // copy ctor
      assert(x.size() == ps.size()); // This has to hold!
      for (size_t i = 0; i < ps.size(); ++i) {
        const Predicate &concp = x[i];
        const Predicate &implp = ps[i];
#ifdef DEBUG
        std::cout << "comparing: impl predicate: " << implp.toString()
                  << " conc predicate :" << concp.toString() << "\n";
#endif
        VarMap subs;
        if (!getSubstitutes(concp, implp, subs))
          break;
#ifdef DEBUG
        for (const auto &[k, v] : subs) {
          std::cout << "key: " << k.toString() << " value: " << v.toString()
                    << "\n";
        }
#endif
        // Now replace all the variables with Atoms in the
        // Implication.
        for (const auto &[key, value] : subs)
          if (temp.hasVariable(key))
            temp = temp.subsVartoAtom(key, value.toString().c_str());
      }
#ifdef DEBUG
      std::cout << "after replacement implication: " << temp.toString() << "\n";
#endif
      // Now check if antecendets are satisfied
      if (check_antecedents(temp)) {
#ifdef DEBUG
        std::cout << "Antecedents met! Triggering the consequent\n";
#endif
        toret = true;
        out = std::move(temp);
        break;
      }
    }
    return toret;
  }

  // XXX: For each of the impl predicate variables replace the
  // conc predicate Atom. If the impl predicate arg is Atom then
  // it should match the conc predicate arg Atom.
  bool getSubstitutes(const Predicate &concp, const Predicate &implp,
                      VarMap &subs) {
    bool toret = true;
    assert(concp.getName() == implp.getName());
    assert(concp.getArity() == implp.getArity());
    for (size_t i = 0; i < concp.getArity(); ++i) {
      auto &concp_arg = concp.getArgs()[i];
      auto &implp_arg = implp.getArgs()[i];
      assert(!std::holds_alternative<Variable>(concp_arg));
      if (std::holds_alternative<Variable>(implp_arg)) {
        subs.insert({std::get<Variable>(implp_arg), std::get<Atom>(concp_arg)});
      } else {
        // The atoms should match
        toret &= (std::get<Atom>(concp_arg).toString() ==
                  std::get<Atom>(implp_arg).toString());
      }
    }
    return toret;
  }

  bool check_contradiction(const Formula *v) {
    std::string fstring;
    if (v->isNot()) {
      fstring = v->getComplementInner();
    } else if (v->isAtom() || v->isOBL()) {
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
  // Removing from conclusions and processing as long as remove exists
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
              [&stemp](const Formula &x) { return x.toString() == stemp; });

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
      if (std::find(cons.cbegin(), cons.cend(), x.toString()) != cons.end())
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

  // Arena for allocating nodes
  std::vector<Node> Arena;

  // Get the leaves from a given node
  void get_leaves(Node &node, std::vector<Node *> &leaves) {
    if (node.edges.empty()) {
      leaves.push_back(&node);
    } else {
      for (size_t i = 0; i < node.edges.size(); ++i)
        get_leaves(Arena[node.edges[i]], leaves);
    }
  }

  std::vector<uint64_t> get_sat_rules(const Node &node) const {
    std::vector<uint64_t> rules{};
    for (const auto &[k, v] : *ruletbl)
      if (v.getConsequent()->toString() == node.f->toString())
        rules.push_back(k);
    return rules; // Hope this does copy elision.
  }

  void get_higher_precedence_rules(uint64_t rule,
                                   std::vector<uint64_t> &hrules) {
    for (auto &[hi, i] : precedence) {
      if (i == rule) {
        hrules.push_back(hi);
      }
    }
  }

  // Process the node
  void process_node(Node &node, std::vector<Node *> &&leaves) {
    // First get all the leaves for this node
    leaves.clear();           // reuse this vector
    get_leaves(node, leaves); // leaves in the leaves vector
    // Now get the rules that makes this node satisfied
    std::vector<uint64_t> rules = get_sat_rules(node);

    // If there is no rule that is needed to satisfy this node' formula then it
    // has to be a fact (theorem)
    if (rules.empty()) {
      node.fact = true;
      return;
    }

    // Now get the antecedents from each of these rules
    for (const uint64_t key : rules) {
      const Antecedent &ant = ruletbl->find(key).getAntecedents();
      for (const Formula &f : ant) {
        Arena.push_back({{}, &f, false});
      }
      // Now attach the nodes just pushed to the arena
      for (size_t i = Arena.size() - ant.size(); i < Arena.size() - 1; ++i) {
        Arena[i].edges.push_back(i + 1);
      }
      // Attach it to each of the leaves.
      for (auto node : leaves) {
        node->edges.push_back(Arena.size() - ant.size());
      }
      // Now check if any rule has higher precedence than the current
      // rule?
      std::vector<uint64_t> higher_rules;
      get_higher_precedence_rules(key, higher_rules);
      // Now get the not of the higher precedence rule' consequent.
      for (auto i : higher_rules) {
        const Formula *v = ruletbl->find(i).getConsequent();
        if (v->isNot()) {
          const Formula *cons = new Formula(v->getComplementInnerFormula());
          Arena.push_back({{}, cons, false, true});
        } else if (v->isAtom() || v->isOBL()) {
          const Formula *cons = new Formula(v->getComplementFormula());
          Arena.push_back({{}, cons, false, true});
        }
      }
      // Now connect the higher precedence nodes one to another
      for (size_t i = Arena.size() - higher_rules.size(); i < Arena.size() - 1;
           ++i) {
        Arena[i].edges.push_back(i + 1);
      }
      // Connect the last antecedent node to the first higher precedence node
      Arena[Arena.size() - higher_rules.size() - 1].edges.push_back(
          Arena.size() - higher_rules.size());
    }
    // Now we are done with this node, so we proceed with dfs
    for (uint64_t x : node.edges) {
      process_node(Arena[x], std::move(leaves));
    }
  }

  void _get_facts(size_t index, std::vector<size_t> &facts,
                  std::vector<std::vector<size_t>> &ffacts) {

    if (Arena[index].fact) {
      facts.push_back(index);
    }

    // We have reached the leaf node following this path, so copy the
    // fact into the final facts vector.
    if (Arena[index].edges.empty()) {
      // copy the facts into the ffacts
      ffacts.push_back(facts);
      return;
    }

    for (size_t index : Arena[index].edges) {
      _get_facts(index, facts, ffacts);
      // Remove this node from the facts too
      auto it = std::find(facts.begin(), facts.end(), index);
      if (it != facts.end()) {
        facts.erase(it);
      }
    }
  }

  void get_facts(std::vector<std::vector<size_t>> &facts) {
    std::vector<size_t> mfacts;
    // Start from the root node.
    _get_facts(0, mfacts, facts);
  }

  // Backward chaining (reasoning) for getting the facts
  void build_and_or_tree(const Formula &goal) {
    // First make the node for the goal
    Arena.push_back(Node{{}, &goal, false});
    // Now just traverse the tree
    process_node(Arena[0], {});
    // Now just get the set of facts needed to prove the final goal.
    std::vector<std::vector<size_t>> facts;
    get_facts(facts);
  }
};
