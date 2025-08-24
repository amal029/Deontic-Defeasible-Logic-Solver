  * This is a solver for Deontic Defeasible Logic programs.
	1. The solver is written in C++ (std C++-17).
	2. The solver can be compiled using: ``g++ -O3 -Werror -Wall
	-Wpedantic -std=c++17 main.cpp -o solver && ./solver``
	3. Or using: ``g++ -O3 -Werror -Wall -Wpedantic -std=c++17
	mainfol.cpp -o solver && ./solver``
    4. The example is written in the main function.
	
* What is Deontic logic 

    Deontic logic is a modal logic with the primary modal construct
    _Obligation_ (``O``). The semantics of obligation is given in terms
    of Kripke' many world' model $M = <W, R, L>$, where $W$ is the set
    of worlds, $R \subseteq W \times W$ is the acceptance relation and
    $L : (W) \rightarrow p, p \in Prop$, where $Prop$ are a set of
    propositions.
	
	Then, $M, w |= O(p), if\ \forall w', s.t.\  R(w, w'), M, w' |= p$.
	
	We can also define, permission $P(p) = \neg O(\neg p)$ and Forbidden
    as $F(p) = \neg O(p)$.
	
* What is Defeasible logic

  Consider that your knowledge base (KB) has two rules: $r1: A \implies
  B$ and $r2: C \implies \neg B$. Moreover, also consider that both
  propositions $A$ and $C$ are facts given as input to your inference
  engine. Then, _classical_ propositional logic inference engine will
  give an output of _unsatisfied_.
  
  Defeasible logic gives precedence to logical rules in the KB. If $r1 >
  r2$, then $r1$ _defeats_ rule $r2$ to produce an output that
  proposition $B$ holds.
  
* Deontic and Defeasible logic can be used to infer legal rules, e.g.,
  contracts, driving road rules, building construction rules, etc.
  
* We are developing a tool suite for natural language translation to
  Deontic Defeasible logic via AI for runtime verification and
  compliance. This tool forms the part of the tool suite.
  
* Supported fragments of logic:
  1. Defeasible logic.
  2. Deontic Obligation operator.
  3. Predicates and propositions.
  4. Negation.
  5. Conjunction in the antecedents.
  6. Consequent should be a literal or a predicate.
  
* Algorithm
  1. Uses unification for solving predicates quantified with $\forall$
     in each KB rule.
  2. Backtracking for solving defeasible logic.
  3. See two example files: ``main.cpp`` and ``mainfol.cpp`` for
     propositional deontic defeasible logic and first order deontic
     defeasible logic examples, respectively.
	 
* Research questions:
  1. What is an efficient algorithm to prove equivalence of two
     (possibly syntactically different) deontic defeasible logic
     formulae?
  2. What is an efficient runtime verification and enforcement strategy
     for deontic defeasible logic formulae on a given model?
  3. How can one capture ethical rules and cultural rules described in
     Natural language into __temporal deontic defeasible logic__?
