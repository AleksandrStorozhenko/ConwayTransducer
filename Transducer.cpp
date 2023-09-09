#include <algorithm>
#include <assert.h>
#include <chrono>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <stack>
#include <string>
#include <vector>

using namespace std;

struct Transducer {

  int inputLetters, outputLetters;

  set<int> startNodes, finalNodes;

  vector<vector<vector<pair<int, int>>>> table;

  Transducer(int inputLetters, int outputLetters)
      : inputLetters(inputLetters), outputLetters(outputLetters){

                                    };

  Transducer(int stateCount, int inputLetters, int outputLetters,
             set<int> &startNodes)
      : inputLetters(inputLetters), outputLetters(outputLetters),
        startNodes(startNodes) {
    table.resize(stateCount);
    for (int i = 0; i < stateCount; ++i) {
      table[i].resize(inputLetters + 1);
    }
  };

  Transducer(int stateCount, int inputLetters, int outputLetters,
             set<int> &startNodes, set<int> &finalNodes)
      : Transducer(stateCount, inputLetters, outputLetters, startNodes) {
    this->finalNodes = finalNodes;
  };

  void addEdge(int inpchar, int outchar, int A, int B) {

    // check that both A and B are present in the array;

    for (int k = table.size(); k <= max(A, B); k++) {
      table.push_back(vector<vector<pair<int, int>>>(inputLetters + 1));
    }

    // basically if it's the first one then
    // table[A].resize(inputLetters+1);
    table[A][inpchar].push_back({outchar, B});
  }

  Transducer compose(const Transducer &T1) {

    map<pair<int, int>, int> T;

    Transducer trsdComp = Transducer(this->inputLetters, T1.outputLetters);

    vector<pair<int, int>> q;
    set<pair<int, int>> visited;

    // returns the index of a pair state. If the pair state is met for the first
    // time, assign a new index and enqueue it for further exploration.
    auto index = [&](pair<int, int> state) {
      auto [it, insertion_did_happen] = T.try_emplace(state, T.size());
      auto idx = it->second;
      if (insertion_did_happen) {
        q.push_back(state);
        if (this->finalNodes.find(state.first) != this->finalNodes.end() &&
            T1.finalNodes.find(state.second) != T1.finalNodes.end())
          trsdComp.finalNodes.insert(idx);
      }
      return idx;
    };

    for (auto i : this->startNodes)
      for (auto j : T1.startNodes)
        trsdComp.startNodes.insert(index(pair(i, j)));

    while (!q.empty()) {
      auto state = q.back();
      q.pop_back();

      for (int i = 0; i <= this->inputLetters; i++) {
        for (auto edge1 : table[state.first][i]) {
          for (auto edge2 : T1.table[state.second][edge1.first]) {
            pair<int, int> next = {edge1.second, edge2.second};

            trsdComp.addEdge(i, edge2.first, T.at(state), index(next));
          }

          if (edge1.first == 0) {
            pair<int, int> next = {edge1.second, state.second};
            trsdComp.addEdge(i, 0, T.at(state), index(next));
          }
        }
      }

      for (auto edge1 : T1.table[state.second][0]) {
        pair<int, int> next = {state.first, edge1.second};
        trsdComp.addEdge(0, edge1.first, T.at(state), index(next));
      }
    }

    return trsdComp;
  }

  Transducer transpose() {

    // the number of states remains
    Transducer transpose =
        Transducer(table.size(), this->inputLetters, this->outputLetters,
                   this->finalNodes, this->startNodes);

    // iterate over the nodes
    for (int i = 0; i < (int)table.size(); i++) {
      for (int c = 0; c <= this->inputLetters; c++) {
        for (auto edge : table[i][c]) {
          transpose.addEdge(c, edge.first, edge.second, i);
        }
      }
    }

    return transpose;
  }

  void backtrack(string word, int node, string &out, set<string> &res) {

    if (word.size() == 0) {
      if (count((this->finalNodes).begin(), (this->finalNodes).end(), node)) {
        res.insert(out);
      }
      return;
    }

    for (auto edge : table[node][(int)(word[0] - 48)]) {

      if (((int)(word[0]) - 48 == edge.first) and (edge.first == 0) &&
          node == edge.second) {
        continue;
      }

      if (edge.first != 0) {
        out += to_string(edge.first);
      }

      backtrack(word.substr(1), edge.second, out, res);
      backtrack("0" + word.substr(1), edge.second, out, res);

      if (edge.first != 0) {
        out = out.substr(0, out.size() - 1);
      }
    }
  }

  set<string> traverse(string word) {

    set<string> res;
    string out = "";

    for (auto sn : this->startNodes) {

      this->backtrack(word, sn, out, res);
      this->backtrack("0" + word, sn, out, res);
    }

    return res;
  }

  // Element Retrieval Backtrack Function

  void dfs(int node, vector<int> path, string &word, set<vector<int>> &res_path,
           set<string> &res_word, set<int> &visited) {

    if (count((this->finalNodes).begin(), (this->finalNodes).end(), node)) {

      res_path.insert(path);
      res_word.insert(word);

      int total_edges = 0;

      for (int i = 0; i < table.size(); i++) {
        total_edges += table[i].size();
        if (total_edges)
          break;
      }

      if (!total_edges)
        return;
    }

    for (int i = 0; i <= this->inputLetters; i++) {
      for (auto edge : table[node][i]) {

        if (visited.find(edge.second) != visited.end())
          continue;

        path.push_back(edge.second);
        string s = to_string(edge.first);
        word += s;
        visited.insert(edge.second);

        dfs(edge.second, path, word, res_path, res_word, visited);

        visited.erase(edge.second);
        word = word.substr(0, word.size() - 1);
        path.pop_back();
      }
    }
  }

  set<string> getElements() {

    // Backtrack through the graph to get the elements (There are no cycles), so
    // we just find all the possible paths through the graph

    set<vector<int>> res_path;
    set<string> res_word;

    for (auto s : this->startNodes) {
      set<int> visited;
      vector<int> path{s};
      string w = "";
      dfs(s, path, w, res_path, res_word, visited);
    }

    return res_word;
  }

  void next(int node, int letter, set<int> &nextStates) {

    for (auto el : table[node][letter]) {
      nextStates.insert(el.second);
    }
  }

  set<int> nextSet(set<int> &S, int letter) {

    // S - set of states; T - set of next states;

    set<int> T;

    for (auto q : S) {
      this->next(q, letter, T);
    }

    return this->closure(T);
  }

  set<int> closure(set<int> closure) {

    if (closure.size() != 0) {
      stack<int> S;
      set<int> visited;

      for (auto el : closure) {
        S.push(el);
      }

      while (!S.empty()) {

        int node = S.top();
        S.pop();
        visited.insert(node);

        for (auto e : table[node][0]) {

          int out_node = e.second;

          if (visited.find(out_node) == visited.end()) {
            closure.insert(out_node);
            S.push(out_node);
          }
        }
      }
    }

    return closure;
  }

  void explore(map<set<int>, int> &T, set<int> &S, Transducer &B) {

    for (int c = 1; c <= B.inputLetters; c++) {

      set<int> U = this->nextSet(S, c);

      if (T.find(U) != T.end()) {

        if (T.size() > B.table.size()) {

          int temp = B.table.size();

          B.table.resize(T.size());

          for (int i = temp; i < T.size(); i++) {
            B.table[i].resize(B.inputLetters + 1);
          }
        }

        B.table[T.at(S)][c].push_back({0, T.at(U)});
      } else {

        T.insert({U, T.size()});

        if (T.size() > B.table.size()) {

          int temp = B.table.size();
          B.table.resize(T.size());

          for (int i = temp; i < T.size(); i++) {
            B.table[i].resize(B.inputLetters + 1);
          }
        }

        B.table[T.at(S)][c].push_back({0, T.at(U)});

        this->explore(T, U, B);
      }
    }
  }

  Transducer determinize() {

    Transducer B = Transducer(this->inputLetters, this->outputLetters);

    set<int> I = this->closure(this->startNodes);

    map<set<int>, int> T;
    T.insert({I, T.size()});

    this->explore(T, I, B);

    // Start & Final Nodes Specification

    if (B.table.size()) {

      B.startNodes.insert(T.at(I));

      for (const auto &el : T) {
        for (int i = 0; i < table.size(); i++) {
          if (count(el.first.begin(), el.first.end(), i) &&
              count((this->finalNodes).begin(), (this->finalNodes).end(), i)) {
            B.finalNodes.insert(el.second);
          }
        }
      }
    }

    return B;
  }

  // Minimization of Recognisers

  Transducer mergeAlph() {

    int newAlph = ((this->inputLetters) + 1) * ((this->outputLetters) + 1) - 1;

    auto rec = Transducer(table.size(), newAlph, 0, this->startNodes,
                          this->finalNodes);

    for (int i = 0; i < table.size(); i++) {
      for (int j = 0; j <= this->inputLetters; j++) {
        for (auto e : table[i][j]) {
          rec.table[i][(e.first * (this->inputLetters + 1) + j)].push_back(
              {0, e.second});
        }
      }
    }

    return rec;
  }

  Transducer splitAlph(int inputLetters, int outputLetters) {

    auto split = Transducer(table.size(), inputLetters, outputLetters,
                            this->startNodes, this->finalNodes);

    for (int i = 0; i < table.size(); i++) {
      for (int j = 0; j <= this->inputLetters; j++) {
        for (auto e : table[i][j]) {

          int inpLetter = j % (inputLetters + 1);

          int outLetter = j / (inputLetters + 1);

          split.table[i][inpLetter].push_back({outLetter, e.second});
        }
      }
    }

    return split;
  }

  Transducer minimize() {

    if (this->outputLetters != 0) {
      // merge the alphabets
      auto merge = mergeAlph();
      auto minimized =
          merge.transpose().determinize().transpose().determinize();
      // split the alphabets
      auto split = minimized.splitAlph(this->inputLetters, this->outputLetters);
      return split;
    } else {
      return transpose().determinize().transpose().determinize();
    }
  }

  // Automata completion - The implemented algorithm already produces a complete
  // dfa - The function is therefore not necessary for our implementation

  Transducer complement() {

    // The final determinization produces a complete dfa
    auto complement = this->minimize();

    // Flip the final nodes vector;
    set<int> invertFinal;

    for (int i = 0; i < (complement.table).size(); i++) {
      if (complement.finalNodes.find(i) == complement.finalNodes.end())
        invertFinal.insert(i);
    }

    complement.finalNodes = invertFinal;

    return complement;
  }

  // It is sufficient to implement the conversion to a filter
  Transducer RtF() {

    Transducer convert =
        Transducer(table.size(), this->inputLetters, this->inputLetters,
                   this->startNodes, this->finalNodes);

    for (int i = 0; i < table.size(); i++) {
      for (int j = 0; j <= convert.inputLetters; j++) {
        for (auto edge : table[i][j]) {
          convert.addEdge(j, j, i, edge.second);
        }
      }
    }

    return convert;
  }

  Transducer FtR() {

    Transducer convert = Transducer(table.size(), this->inputLetters, 0,
                                    this->startNodes, this->finalNodes);

    for (int i = 0; i < table.size(); i++) {
      for (int j = 0; j <= convert.inputLetters; j++) {
        for (auto edge : table[i][j]) {
          convert.addEdge(j, 0, i, edge.second);
        }
      }
    }

    return convert;
  }

  Transducer invert() {

    Transducer invert =
        Transducer(table.size(), this->outputLetters, this->inputLetters,
                   this->startNodes, this->finalNodes);

    for (int i = 0; i < table.size(); i++) {
      for (int j = 0; j <= this->inputLetters; j++) {
        for (auto edge : table[i][j]) {
          invert.addEdge(edge.first, j, i, edge.second);
        }
      }
    }

    return invert;
  }

  bool languageEquality(int u, int v, vector<int> &inverse1,
                        vector<int> &inverse2, Transducer &T1) {

    // u - node of the first graph;
    // v - node of the second graph;

    inverse1[u] = v;
    inverse2[v] = u;

    for (int i = 1; i <= this->inputLetters; i++) {

      int n1, n2;

      if ((table[u][i].size() == (T1.table)[u][i].size())) {
        if (table[u][i].size()) {
          n1 = table[u][i][0].second;
          n2 = T1.table[v][i][0].second;
        } else {
          continue;
        }
      } else {
        return false;
      }

      if (inverse1[n1] == -1 && inverse2[n2] == -1) {
        // if we simply return the language inequality directly then just
        // whichever letter is executed first will fully determine
        if (!languageEquality(n1, n2, inverse1, inverse2, T1))
          return false;
      } else if (inverse1[n1] != n2 || inverse2[n2] != n1) {
        return false;
      }
    }
    return true;
  }

  void latex_table(ostream& out) {
    const string alph = "Sabcdefghijklmnopqrstuvwxyz";
    assert(table.size() <= alph.size());

    out << "\\begin{tabular}{ll}\\toprule\n";

    out << "states & \\texttt{";
    for(int s = 0; s < table.size(); s++) {
      out << alph[s] << "\\,";
    }
    out << "}\\\\\n";

    out << "initial & \\texttt{";
    for(int s = 0; s < table.size(); s++) {
      out << (startNodes.find(s) != startNodes.end() ? '+' : ' ') << "\\,";
    }
    out << "}\\\\\n";

    out << "final & \\texttt{";
    for(int s = 0; s < table.size(); s++) {
      out << (finalNodes.find(s) != finalNodes.end() ? '+' : ' ') << "\\,";
    }
    out << "}\\\\\n";

    out << "\\midrule \n";

    for(int i = 1; i <= inputLetters; i++) {
      out << "input $";
      if (i < 4)
        out << char('0' + i);
      else if(i == 4)
        out << 'd';
      else if(i == 5)
        out << "\\splt";

      out << "$ & " << "\\texttt{";

      for(auto t : table) {
        out << alph[t[i].at(0).second] << "\\,";
      }

      out << "}\\\\\n";
    }

    out << "\\bottomrule \\end{tabular}" << endl;
  }
};

// Implementation of the transducers used in the proof.

Transducer multimark() {

  set<int> sn{0};
  set<int> fn{0};

  Transducer multimark = Transducer(1, 5, 5, sn, fn);

  for (int i = 1; i <= multimark.inputLetters; i++) {
    multimark.addEdge(i, i, 0, 0);
  }

  multimark.addEdge(0, 5, 0, 0);

  return multimark;
}

Transducer singlemark() {

  set<int> sn{0};
  set<int> fn{0, 3};

  Transducer singlemark = Transducer(4, 4, 5, sn, fn);

  for (int i = 1; i <= singlemark.inputLetters; i++) {
    singlemark.addEdge(i, i, 0, 1);
    singlemark.addEdge(i, i, 1, 1);
    singlemark.addEdge(i, i, 2, 2);
    singlemark.addEdge(i, i, 2, 3);
  }

  singlemark.addEdge(0, 5, 1, 2);

  return singlemark;
}

Transducer scissors() {

  set<int> sn{0};
  set<int> fn{2};

  Transducer scissors = Transducer(3, 5, 4, sn, fn);

  for (int i = 1; i <= scissors.inputLetters; i++) {
    scissors.addEdge(i, 0, 0, 0);
    scissors.addEdge(i, 0, 2, 2);
    if (i != 5) {
      scissors.addEdge(i, i, 1, 1);
    }
  }

  scissors.addEdge(5, 0, 0, 1);
  scissors.addEdge(5, 0, 1, 2);

  return scissors;
}

Transducer audioactiveT(bool augmented = false) {

  // it makes sense to place them first i.e
  set<int> sn{0, 1, 2, 3};

  Transducer audioactiveT =
      Transducer(24, 4 + augmented, 4 + augmented, sn, sn);

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

  return audioactiveT;
}

Transducer splitRec() {

  auto aat = audioactiveT(true);
  auto split = aat.FtR().minimize();

  for (int i = 0; i < 9; i++) {
    split = aat.compose(split).minimize();
  }

  return split.minimize();
}

Transducer irredFactorRec() {

  auto sm = singlemark().minimize();
  auto sr = splitRec().minimize();
  auto iwr = audioactiveT().RtF().compose(sm.compose(sr).complement().minimize()).minimize();

  sr.latex_table(cout);
  iwr.latex_table(cout);
  return iwr;
}

Transducer irredFactorDer() {

  auto mmt = multimark();
  auto sf = splitRec().RtF();
  auto sc = scissors();
  auto isf = irredFactorRec().RtF();

  // definition of irredFactorDer - There's still the problem of minimization of
  // derivators

  auto ifd =
      mmt.compose(sf).minimize().compose(sc).minimize().compose(isf).minimize();

  return ifd;
}

// Theorem Proofs

Transducer Theorem2() {

  cout << "Splitting Theorem Proof" << endl;

  auto aat = audioactiveT(true);
  auto split = aat.FtR().minimize();
  auto splitPrev = split;

  int count = 0;
  while (true) {
    cout << "Composition degree n: " << count
         << " | Number of States: " << split.table.size() << endl;

    count++;
    // check for language equality
    split = aat.compose(split).minimize();

    if (split.table.size() == splitPrev.table.size()) {
      vector<int> inverse1((split.table.size()), -1);
      vector<int> inverse2((split.table.size()), -1);
      if (split.languageEquality(0, 0, inverse1, inverse2, splitPrev)) {
        return split.minimize();
      }
    }
    splitPrev = split;
  }
}

void hamiltonPath(string element, vector<string> &path, set<string> &visited,
                  set<vector<string>> &paths,
                  map<string, vector<string>> &adj) {

  if (path.size() == 92) {
    paths.insert(path);
    return;
  }

  for (auto el_next : adj[element]) {

    if (visited.find(element) == visited.end()) {
      path.push_back(element);
      visited.insert(element);
      hamiltonPath(el_next, path, visited, paths, adj);
      path.pop_back();
      visited.erase(element);
    }
  }
}

set<string> CosmologicalTheorem() {

  cout << "Cosmological Theorem Proof" << endl;

  auto at = audioactiveT();

  auto factorizer = irredFactorDer();

  // Deal with the inverted edges

  auto T = at.compose(factorizer).minimize();

  // This is supposed to be the generator - transposed to a recogniser
  // Instead of directly inverting the T automata it might be easier to do

  // Test which one is faster
  auto T_inv = T.invert();

  auto Tn = T_inv.FtR().minimize();

  auto Tn_prev = Tn;

  int count = 1;
  while (true) {

    cout << "Composition degree n: " << count
         << " | Number of States: " << Tn.table.size() << endl;
    count++;

    // Can still make it more efficient by
    Tn = T_inv.compose(Tn).minimize();

    if (Tn_prev.table.size() == Tn.table.size()) {

      cout << "Potential Isomorphism" << endl;

      if (Tn.table.size() == Tn_prev.table.size()) {
        vector<int> inverse1((Tn.table.size()), -1);
        vector<int> inverse2((Tn.table.size()), -1);
        if (Tn.languageEquality(0, 0, inverse1, inverse2, Tn_prev)) {
          cout << "Verified Isomorphism" << endl;
          break;
        }
      }
    }

    Tn_prev = Tn;
  }

  Tn = Tn.invert();

  set<string> words = Tn.getElements();

  cout << "Number of Elements: " << words.size() << endl;

  map<string, vector<string>> adj;

  for (auto word : words) {

    auto deriv = factorizer.traverse(*at.traverse(word).begin());

    for (auto el : deriv) {
      if (adj.count(word)) {
        adj[word].push_back(el);
      } else {
        vector<string> elts;
        elts.push_back(el);
        adj.insert({word, elts});
      }
    }
  }

  // So now we have the elements and want to construct an adjacency list from
  // them;

  // We also want the mapping between elements and atomic labelings;

  vector<string> per_elt_names{
      "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg",
      "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr",
      "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
      "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
      "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
      "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf",
      "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po",
      "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U"};

  set<string> visited;
  set<vector<string>> paths;

  string start = "3";
  vector<string> path;
  hamiltonPath(start, path, visited, paths, adj);

  auto it = paths.begin();
  advance(it, 2);
  path = *it;
  reverse(path.begin(), path.end());

  cout << "Common Elements (Conway's ordering)" << endl;
  for (int i = 0; i < path.size(); i++) {
    cout << i + 1 << '\t' << per_elt_names[i] << '\t'
         << path[i] << string(50 - path[i].size(), ' ') << " (→";
    for (auto el : adj[path[i]]) {
      auto itr = find(path.begin(), path.end(), el);
      cout << ' ' << per_elt_names[distance(path.begin(), itr)];
    }
    cout << ")" << endl;
  }
  cout << endl;

  return words;
}

int main(int argc, const char *argv[]) {

  auto start = chrono::steady_clock::now();

  Theorem2();
  cout << endl;
  CosmologicalTheorem();

  auto end = chrono::steady_clock::now();

  auto diff = end - start;

  cout << chrono::duration<double, milli>(diff).count() << " ms" << endl;
}
