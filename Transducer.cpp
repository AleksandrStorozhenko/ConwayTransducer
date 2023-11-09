#include <algorithm>
#include <assert.h>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <vector>

using namespace std;

struct Transducer {

  using letter = int;
  using state = int;

  int inputLetters, outputLetters;

  set<state> startNodes, finalNodes;

  vector<vector<vector<pair<letter, state>>>> table;

  // Transducer Struct Constructors

  Transducer(int inputLetters, int outputLetters, int stateCount = 0)
      : inputLetters(inputLetters), outputLetters(outputLetters) {}

  bool operator==(const Transducer &other) const {
    return inputLetters == other.inputLetters &&
           outputLetters == other.outputLetters &&
           startNodes == other.startNodes && finalNodes == other.finalNodes &&
           table == other.table;
  };

  void addEdge(letter inpchar, letter outchar, state A, state B) {
    for (state k = table.size(); k <= max(A, B); k++) {
      table.push_back(vector<vector<pair<letter, state>>>(inputLetters + 1));
    }

    table[A][inpchar].push_back({outchar, B});
  }

  // compute the compositon of this *then* T1
  Transducer compose(const Transducer &other) {
    // Map T indexes pairs of states as the states of the composed transducer.
    map<pair<state, state>, state> T;

    Transducer trsdComp = Transducer(this->inputLetters, other.outputLetters);

    // Multisource DFS
    vector<pair<state, state>> q;
    set<pair<state, state>> visited;

    // Helper function returns the index of a pair state. If the pair state is
    // met for the first time, assigns a new index and enqueues it for further
    // exploration.
    auto index = [&](pair<state, state> node) {
      auto [it, insertion_did_happen] = T.try_emplace(node, T.size());
      auto idx = it->second;
      if (insertion_did_happen) {
        q.push_back(node);
        if (finalNodes.find(node.first) != finalNodes.end() &&
            other.finalNodes.find(node.second) != other.finalNodes.end())
          trsdComp.finalNodes.insert(idx);
      }
      return idx;
    };

    for (auto i : startNodes)
      for (auto j : other.startNodes)
        trsdComp.startNodes.insert(index(pair(i, j)));

    while (!q.empty()) {
      auto node = q.back();
      q.pop_back();

      for (letter i = 0; i <= inputLetters; i++) {
        for (auto edge1 : table[node.first][i]) {
          for (auto edge2 : other.table[node.second][edge1.first]) {
            pair<state, state> next = {edge1.second, edge2.second};

            trsdComp.addEdge(i, edge2.first, T.at(node), index(next));
          }

          if (edge1.first == 0) {
            pair<state, state> next = {edge1.second, node.second};
            trsdComp.addEdge(i, 0, T.at(node), index(next));
          }
        }
      }

      for (auto edge1 : other.table[node.second][0]) {
        pair<state, state> next = {node.first, edge1.second};
        trsdComp.addEdge(0, edge1.first, T.at(node), index(next));
      }
    }

    return trsdComp;
  }

  Transducer reverse() {
    Transducer rev = Transducer(inputLetters, outputLetters);
    rev.startNodes = finalNodes;
    rev.finalNodes = startNodes;

    for (state i = 0; i < table.size(); i++)
      for (letter c = 0; c <= inputLetters; c++)
        for (auto edge : table[i][c])
          rev.addEdge(c, edge.first, edge.second, i);

    return rev;
  }

  template <typename It>
  void run(It first, It last, state node, vector<letter> &out,
           vector<vector<letter>> &ret) {
    if (first == last) {
      if (finalNodes.find(node) != finalNodes.end()) {
        ret.push_back(out);
      }
    }

    for (auto e : table[node][0]) {
      if (e.first)
        out.push_back(e.first);
      run(first, last, e.second, out, ret);
      if (e.first)
        out.pop_back();
    }

    if (first < last) {
      for (auto e : table[node][*first]) {
        if (e.first)
          out.push_back(e.first);
        run(first + 1, last, e.second, out, ret);
        if (e.first)
          out.pop_back();
      }
    }
  }

  // return all possible output of this given the input. Assumes that the set of
  // all output is indeed finite and that the transducer does not have a cycle
  // of epsilon-transitions.
  template <typename It> vector<vector<letter>> run(It first, It last) {
    vector<letter> out;
    vector<vector<letter>> ret;

    for (state s : startNodes) {
      run(first, last, s, out, ret);
    }

    return ret;
  }

  template <typename It> bool match(It first, It last, state node) {
    if (first == last && finalNodes.find(node) != finalNodes.end())
      return true;

    for (auto e : table[node][0]) {
      if (match(first, last, e.second))
        return true;
    }

    if (first < last) {
      for (auto e : table[node][*first]) {
        if (match(first + 1, last, e.second))
          return true;
      }
    }

    return false;
  }

  // return true if the input is accepted by this.
  template <typename It> bool match(It first, It last) {
    for (auto s : startNodes) {
      if (match(first, last, s))
        return true;
    }
    return false;
  }

  // closure of S under epsilon transitions
  set<state> closure(const set<state> &S) {
    assert(!S.empty());
    vector<state> queue(S.begin(), S.end());
    set<state> visited;

    while (!queue.empty()) {
      state node = queue.back();
      queue.pop_back();
      if (visited.find(node) != visited.end())
        continue;

      visited.insert(node);

      // iterate over epsilon transitions
      for (auto e : table[node][0])
        queue.push_back(e.second);
    }

    return visited;
  }

  Transducer determinize() {
    Transducer B = Transducer(inputLetters, outputLetters);
    vector<state> queue;
    map<set<state>, state> T;
    vector<set<state>> states;

    // Helper function returns the index of a set of states. If the set of
    // states is met for the first time, assigns a new index and explores it
    // further.
    auto index = [&](const set<state> &U) {
      auto [it, insertion_did_happen] = T.try_emplace(U, states.size());
      int idx = it->second;

      if (insertion_did_happen) {
        auto cl = closure(U);
        queue.push_back(idx);

        // is it a final state?
        for (state s : cl) {
          if (finalNodes.find(s) != finalNodes.end()) {
            B.finalNodes.insert(idx);
            break;
          }
        }

        states.emplace_back(std::move(cl));
      }

      return idx;
    };

    B.startNodes.insert(index(startNodes));

    while (!queue.empty()) {
      state Sidx = queue.back();
      queue.pop_back();

      for (letter c = 1; c <= inputLetters; c++) {
        set<state> T;
        for (state node : states[Sidx])
          for (auto el : table[node][c])
            T.insert(el.second);

        if (!T.empty()) {
          B.addEdge(c, 0, Sidx, index(T));
        }
      }
    }

    return B;
  }

  // Minimization of Recognisers
  Transducer zip_alphabet() {
    int newAlph = (inputLetters + 1) * (outputLetters + 1) - 1;

    auto rec = Transducer(newAlph, 0);
    rec.startNodes = startNodes;
    rec.finalNodes = finalNodes;

    for (state i = 0; i < table.size(); i++)
      for (letter j = 0; j <= this->inputLetters; j++)
        for (auto e : table[i][j])
          rec.addEdge(e.first * (inputLetters + 1) + j, 0, i, e.second);

    return rec;
  }

  Transducer unzip_alphabet(int inputLetters, int outputLetters) {
    auto split = Transducer(inputLetters, outputLetters);
    split.startNodes = startNodes;
    split.finalNodes = finalNodes;

    for (state i = 0; i < table.size(); i++) {
      for (letter j = 0; j <= this->inputLetters; j++) {
        for (auto e : table[i][j]) {
          letter inpLetter = j % (inputLetters + 1);
          letter outLetter = j / (inputLetters + 1);
          split.addEdge(inpLetter, outLetter, i, e.second);
        }
      }
    }

    return split;
  }

  Transducer minimize() {

    if (outputLetters) {
      // merge the alphabets
      auto merge = zip_alphabet();
      auto minimized = merge.reverse().determinize().reverse().determinize();
      // split the alphabets
      auto split =
          minimized.unzip_alphabet(this->inputLetters, this->outputLetters);
      return split;
    } else {
      return reverse().determinize().reverse().determinize();
    }
  }

  // Automata completion - The implemented algorithm already produces a
  // complete dfa - The function is therefore not necessary for our
  // implementation

  Transducer complement() {

    // The final determinization produces a complete dfa
    auto complement = minimize();

    // Flip the final nodes vector;
    set<state> invertFinal;

    for (state i = 0; i < complement.table.size(); i++) {
      if (complement.finalNodes.find(i) == complement.finalNodes.end())
        invertFinal.insert(i);
    }

    complement.finalNodes = invertFinal;

    return complement;
  }

  Transducer filter() {

    Transducer convert = Transducer(inputLetters, inputLetters);
    convert.startNodes = startNodes;
    convert.finalNodes = finalNodes;

    for (state i = 0; i < table.size(); i++) {
      for (letter j = 0; j <= convert.inputLetters; j++) {
        for (auto edge : table[i][j]) {
          convert.addEdge(j, j, i, edge.second);
        }
      }
    }

    return convert;
  }

  Transducer recognizer() {

    Transducer convert = Transducer(inputLetters, 0);
    convert.startNodes = startNodes;
    convert.finalNodes = finalNodes;

    for (state i = 0; i < table.size(); i++) {
      for (letter j = 0; j <= convert.inputLetters; j++) {
        for (auto edge : table[i][j]) {
          convert.addEdge(j, 0, i, edge.second);
        }
      }
    }

    return convert;
  }

  Transducer transpose() {

    auto transp = Transducer(outputLetters, inputLetters);
    transp.startNodes = startNodes;
    transp.finalNodes = finalNodes;

    for (state i = 0; i < table.size(); i++)
      for (letter j = 0; j <= inputLetters; j++)
        for (auto edge : table[i][j])
          transp.addEdge(edge.first, j, i, edge.second);

    return transp;
  }
};

bool hamiltonian_path(const vector<vector<size_t>> &adj, size_t element,
                      size_t size, vector<size_t> &path, set<size_t> &visited) {
  if (visited.find(element) != visited.end())
    return false;

  visited.insert(element);
  path.push_back(element);
  size--;

  if (size == 0)
    return true;

  for (auto next : adj[element]) {
    if (hamiltonian_path(adj, next, size, path, visited))
      return true;
  }

  visited.erase(element);
  path.pop_back();
  size++;

  return false;
}

// Definition of the fundemental transducers (Transducers that are not derived
// from the fundemental composition and complementation operations on
// transducers)

Transducer multimark() {

  Transducer multimark = Transducer(5, 5);
  multimark.startNodes.insert(0);
  multimark.finalNodes.insert(0);

  for (int i = 1; i <= multimark.inputLetters; i++) {
    multimark.addEdge(i, i, 0, 0);
  }

  multimark.addEdge(0, 5, 0, 0);

  return multimark.minimize();
}

Transducer singlemark() {

  Transducer singlemark = Transducer(4, 5);
  singlemark.startNodes.insert(0);
  singlemark.finalNodes.insert(0);
  singlemark.finalNodes.insert(3);

  for (int i = 1; i <= singlemark.inputLetters; i++) {
    singlemark.addEdge(i, i, 0, 1);
    singlemark.addEdge(i, i, 1, 1);
    singlemark.addEdge(i, i, 2, 2);
    singlemark.addEdge(i, i, 2, 3);
  }

  singlemark.addEdge(0, 5, 1, 2);

  return singlemark.minimize();
}

Transducer scissors() {
  Transducer scissors = Transducer(5, 4);
  scissors.startNodes.insert(0);
  scissors.finalNodes.insert(2);

  for (int i = 1; i <= scissors.inputLetters; i++) {
    scissors.addEdge(i, 0, 0, 0);
    scissors.addEdge(i, 0, 2, 2);
    if (i != 5) {
      scissors.addEdge(i, i, 1, 1);
    }
  }

  scissors.addEdge(5, 0, 0, 1);
  scissors.addEdge(5, 0, 1, 2);

  return scissors.minimize();
}

Transducer audioactive(bool augmented = false) {
  Transducer audioactiveT = Transducer(4 + augmented, 4 + augmented);

  for (int i = 0; i <= 3; i++) {
    audioactiveT.startNodes.insert(i);
    audioactiveT.finalNodes.insert(i);
  }

  for (int c = 1; c <= 4; ++c) {

    audioactiveT.addEdge(c, 1, c - 1, 5 * (c - 1) + 1 + 3);
    audioactiveT.addEdge(0, c, 5 * (c - 1) + 1 + 3, 5 * (c - 1) + 1 + 4);

    audioactiveT.addEdge(c, 2, c - 1, 5 * (c - 1) + 1 + 5);
    audioactiveT.addEdge(c, c, 5 * (c - 1) + 1 + 5, 5 * (c - 1) + 1 + 4);

    audioactiveT.addEdge(c, 3, c - 1, 5 * (c - 1) + 1 + 6);
    audioactiveT.addEdge(c, c, 5 * (c - 1) + 1 + 6, 5 * (c - 1) + 1 + 7);
    audioactiveT.addEdge(c, 0, 5 * (c - 1) + 1 + 7, 5 * (c - 1) + 1 + 4);
  }

  for (int c = 1; c <= 4; ++c) {
    for (int i = 0; i <= 3; ++i) {
      if (i != c - 1) {
        audioactiveT.addEdge(0, 0, 5 * (c - 1) + 1 + 4, i);
      }
    }
  }

  if (augmented) {
    for (int i = 0; i <= 3; i++) {
      audioactiveT.addEdge(5, 5, i, i);
    }
  }

  return audioactiveT.minimize();
}

int main(int argc, const char *argv[]) {

  // CONSTRUCTION OF THE SPLITTING RECOGNIZER

  // We begin by defining a series of useful "example" transducers
  // These transducers will serve as the fundamental building blocks for the
  // more complex automata, i.e all others will derive from these via
  // composition.
  auto multimarkT = multimark();
  auto singlemarkT = singlemark();
  auto scissorsT = scissors();
  // Notably, we define we define the audioactive transducer, modelizing the
  // derivation process for day-one sequences.
  auto audioactiveT = audioactive();
  // Finally, we define the augmented audioactive transducer, additionally
  // verifying whether there's a . (5, in our internal representation) between
  // 2 equal characters.
  auto aaT = audioactive(true);

  // From these 5 relatively simple automata, we will be able to prove the
  // Cosmological Theorem.

  // We begin by proving the Splitting Theorem.

  cout << "Splitting Theorem Proof" << endl;

  auto splitR = aaT.recognizer().minimize();

  int count = 0;
  while (true) {
    cout << "Composition order n: " << setw(2) << count++
         << " | Number of States: " << setw(3) << splitR.table.size() << endl;

    auto next = aaT.compose(splitR).minimize();

    if (splitR == next) {
      cout << "Fixed point found\n" << endl;
      break;
    }

    splitR = next;
  }

  // PROOF OF THE COSMOLOGICAL THEOREM

  cout << "Cosmological Theorem Proof" << endl;

  // From the splitting recognizer, we can construct a recognizer for the set
  // of irreducible words: As outlined in the paper, we begin by computing the
  // complement for the recogniser of the reducible words over W,
  // and further compose it with a filter for the input language of the
  // audioactive transducer to reject the words not in W.
  auto irreducibleR = singlemarkT.compose(splitR).complement();
  irreducibleR = audioactiveT.filter().compose(irreducibleR);

  // Prior to proving the main result, let us construct the irreducible factor
  // extractor (Constant minimization allows to assure better compositional
  // runtime)

  auto atomT = multimarkT.compose(splitR.filter())
                   .minimize()
                   .compose(scissorsT)
                   .minimize()
                   .compose(irreducibleR.filter())
                   .minimize();

  // We can now give the proof of the Cosmological Theorem !

  auto T = audioactiveT.compose(atomT).minimize();

  // we want to work with recognizers, not generators, so we transpose
  // everything
  auto Ttranspose = T.transpose();
  auto elementR = Ttranspose.recognizer().minimize();

  count = 1;
  while (true) {

    cout << "Composition order n: " << setw(2) << count++
         << " | Number of States: " << setw(3) << elementR.table.size() << endl;

    auto next = Ttranspose.compose(elementR).minimize();

    if (next == elementR) {
      cout << "Fixed point found\n" << endl;
      break;
    }

    elementR = next;
  }

  // element generator
  auto elementG = elementR.transpose();

  // COMPILATION OF THE ELEMENT TABLE

  // atomT is now a generator with finite output language, consisting of the 92
  // common elements and the 2 transuranic elements. We enumerate the language
  // and numer the elements so that the n-th element appears in the derication
  // of the (n+1)-th element.

  auto makesplitT = singlemarkT.compose(splitR.filter()).minimize();

  vector<Transducer::letter> empty;
  vector<vector<Transducer::letter>> elements =
      elementG.run(empty.begin(), empty.end());
  vector<vector<size_t>> derivation(elements.size());

  // compute the derivation of each element in terms of the other elements
  for (size_t i = 0; i < elements.size(); i++) {
    // a vector containing a vector containing the unique derivation of
    // elements[i]
    vector<vector<Transducer::letter>> queue =
        audioactiveT.run(elements[i].begin(), elements[i].end());
    assert(queue.size() == 1);

    // splits the derivation into irreducible factors
    while (!queue.empty()) {
      auto w = queue.back();
      queue.pop_back();
      auto trysplit = makesplitT.run(w.begin(), w.end());

      if (!trysplit.empty()) {
        auto split = trysplit[0];
        auto mark = find(split.begin(), split.end(), 5);
        queue.push_back(vector(mark + 1, split.end()));
        queue.push_back(vector(split.begin(), mark));
      } else {
        auto itfact = find(elements.begin(), elements.end(), w);
        assert(itfact != elements.end());
        derivation[i].push_back(itfact - elements.begin());
      }
    }
  }

  vector<size_t> path, inv;
  size_t uranium =
      find(elements.begin(), elements.end(), vector<Transducer::letter>{3}) -
      elements.begin();

  {
    // Recover Conway's ordering of the elements
    set<size_t> visited;
    assert(hamiltonian_path(derivation, uranium, 92, path, visited));
    reverse(path.begin(), path.end());

    for (size_t i = 0; i < elements.size(); i++) {
      if (elements[i].back() == 4)
        path.push_back(i);
    }

    assert(path.size() == 94);

    for (size_t i = 0; i < path.size(); i++)
      inv.push_back(i);

    sort(inv.begin(), inv.end(),
         [&](size_t i, size_t j) { return path[i] < path[j]; });
  }

  vector<string> names{
      "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg",
      "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr",
      "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
      "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
      "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
      "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf",
      "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po",
      "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu"};

  cout << "Element table" << endl;
  for (size_t i = 0; i < elements.size(); i++) {
    string elstr;
    for (auto a : elements[path[i]])
      elstr.append(1, (char)('0' + a));

    cout << i + 1 << '\t' << names[i] << '\t' << std::left << setw(50) << elstr
         << " (→";
    for (auto d : derivation[path[i]]) {
      cout << ' ' << names[inv[d]];
    }
    cout << ")" << endl;
  }

  return 0;
}
