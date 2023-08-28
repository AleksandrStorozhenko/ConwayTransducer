#include <iostream>
#include <vector>
#include <stack>
#include <set>
#include <map>
#include <algorithm>
#include <string>
#include <assert.h>

using namespace std;

#define endl "\n";

struct Transducer{
    
    int inputLetters, outputLetters;
    
    // at this point we might as well just use sets
    vector<int> startNodes, finalNodes;
    
    vector<vector<vector<pair<int, int>>>> table;

    // Use initializer lists for the constructors + Also pass optional parameters to condense all the constructors into 1.
    Transducer(int inputLetters, int outputLetters){
        // 0 -> e; 1 - 1; 2 - 2; 3 - 3; 4 - n; 5 - .;
        this -> inputLetters = inputLetters;
        this -> outputLetters = outputLetters;
        
    };

    Transducer(int stateCount, int inputLetters, int outputLetters, vector<int> &startNodes){
        this -> inputLetters = inputLetters;
        this -> outputLetters = outputLetters;
        
        (this -> table).resize(stateCount);
        
        for(int i = 0; i < stateCount; ++i){
            (this -> table)[i].resize((this->inputLetters)+1);
        }
        
        this -> startNodes = startNodes;
        
        vector<int> finalNodes;
        
        this -> finalNodes = finalNodes;
    };
    
    Transducer(int stateCount, int inputLetters, int outputLetters, vector<int> &startNodes, vector<int> &finalNodes){
        
        // 0 -> e; 1 - 1; 2 - 2; 3 - 3; 4 - n; 5 - .;
        
        this -> inputLetters = inputLetters;
        this -> outputLetters = outputLetters;
        
        (this -> table).resize(stateCount);
        
        for(int i = 0; i < stateCount; ++i){
            
            (this -> table)[i].resize((this->inputLetters)+1);
            
        }
        
        this -> startNodes = startNodes;

        this -> finalNodes = finalNodes;

    };
    
    // Separate the process of adding a new node from adding an edge
    void addEdge(int inpchar, int outchar, int A, int B){
        
        // check that both A and B are present in the array
    
        if(max(A, B) > ((int)(this->table).size() - 1)){
            (this->table).resize(max(A,B)+1);
        }
        
        
        //TODO: This can be optimized
        for(int i = 0; i < (this->table).size(); i++){
            (this->table)[i].resize(inputLetters+1);
        }
        
        // basically if it's the first one then
        // (this->table)[A].resize(inputLetters+1);
        (this->table)[A][inpchar].push_back({outchar, B});
        
    }
    
    Transducer compose(Transducer &T1){
        
        //number of nodes in the 2 transducers
        auto X = (this -> table).size();
        auto Y = T1.table.size();
        
        vector<int> initialNodes;
        
        for(auto node1: this->startNodes){
            for(auto node2: T1.startNodes){
                auto pos = node2 * X + node1;
                initialNodes.push_back(pos);
            }
        }
        
        //Define the composed transducer - the size should be set to the product of the number of nodes XY;
        Transducer trsdComp = Transducer(X*Y, this->inputLetters, T1.outputLetters, initialNodes);
        
        //iterate over the states of transducer 1
        for(int A = 0; A < X; A++){
            
            // iterate over the letters of transducer 1
            for(int c1 = 0; c1 <= this -> inputLetters; c1++){
                
                for(auto edge1: (this -> table)[A][c1]){

                    int c2 = edge1.first;
                    int B = edge1.second;
                    
                    if(c2 == 0){
                        for(int I = 0; I < Y; I++){
//                            cout << "Add Case 1 Epsilon Edge:" << "("<< c1 << 0 << I*X + A << I*X + B << ")" << endl;
                            trsdComp.addEdge(c1, 0, I * X + A, I*X + B);
                        }
                    }

                    //iterate over the states of transducer 2
                    for(int C = 0; C < Y; C++){
                        for(auto edge2: T1.table[C][c2]){

                            int c3 = edge2.first;
                            int D = edge2.second;
                            
//                            cout << "Add Edge:" << "("<< c1 << c3 << C*X + A << D*X + B << ")" << endl;
                            
                            trsdComp.addEdge(c1, c3, C * X + A, D * X + B);

                        }
                    }
                }
            }
        }
        
        for(int C = 0; C < Y; C++){
            for(auto e: T1.table[C][0]){
                int outChar = e.first;
                int D = e.second;
                for(int A = 0; A < X; A++){
//                    cout << "Add Case 2 Epsilon Edge:" << "("<< 0 << outChar << C*X + A << D*X + A << ")" << endl;
                    trsdComp.addEdge(0, outChar, C * X + A, D * X + A);
                }
            }
        }
        
        for(int i = 0; i < X; i++){
            for(int j = 0; j < Y; j++){
                
                // the trsdComp.table[j * X + i] entry is not necessarily defined;
                
                if(count((this->finalNodes).begin(), (this->finalNodes).end(), i) &&
                   count(T1.finalNodes.begin(), T1.finalNodes.end(), j) &&
                   trsdComp.table[j * X + i].size() != 0){
                    trsdComp.finalNodes.push_back(j * X + i);
                }
                
            }
        }

        return trsdComp;
    }
    
    Transducer transpose(){
        
        // the number of states remains
        Transducer transpose = Transducer((this->table).size(), this->inputLetters, this->outputLetters, this->finalNodes, this->startNodes);
        
        //iterate over the nodes
        for(int i = 0; i < (int)(this -> table).size(); i++){
            for(int c = 0; c <= this -> inputLetters; c++){
                for(auto edge: (this->table)[i][c]){
                    transpose.addEdge(c, edge.first, edge.second, i);
                }
            }
        }
        
        return transpose;
    }
    
    // Traversal
    
    void backtrack(string word, int node, string &out, set<string> &res){
        
        if(word.size() == 0){
            if(count((this->finalNodes).begin(), (this->finalNodes).end(), node)){
                res.insert(out);
            }
            return;
        }
        
        for(auto edge: (this->table)[node][(int)(word[0]-48)]){
            
            if(((int)(word[0])-48 == edge.first) and (edge.first == 0) && node == edge.second){
                continue;
            }
            
            if(edge.first != 0){
                out += to_string(edge.first);
            }
        
            backtrack(word.substr(1), edge.second, out, res);
            backtrack("0" + word.substr(1), edge.second, out, res);
            
            if(edge.first != 0){
                out = out.substr(0, out.size()-1);
            }
        }
    }
    
    set<string> traverse(string word){
        
        set<string> res;
        string out = "";
        
        for(auto sn: this->startNodes){
            
            this->backtrack(word, sn, out, res);
            this->backtrack("0" + word, sn, out, res);
        }

        return res;
    }
    
    // Element Retrieval Backtrack Function
    
//    void dfs(int node, vector<int> path, string &word, set<vector<int>> &res_path, set<string> &res_word){
//
//        // This is a particular case, in general we need to be able to continue
//        if(count((this->finalNodes).begin(), (this->finalNodes).end(), node)){
//
//            res_path.insert(path);
//            res_word.insert(word);
//
//            // Keep exploring until you exhaust all the possible paths
//            int total_edges = 0;
//
//            for(int i = 0; i < (this->table).size(); i++){
//                total_edges += (this->table)[i].size();
//                if(total_edges) break;
//            }
//
//            if(!total_edges) return;
//        }
//
//        for((char_inp, char_out, next_node) in self.G[node]){
//            if next_node in visited:
//                continue
//
//            path.append(next_node)
//            word += char_out
//            visited.add(next_node)
//
//            dfs(next_node, path, word)
//
//            visited.remove(next_node)
//            word = word[:-1]
//            path.pop()
//        }
//    }
//
//    set<string> getElements(){
//
//        // Backtrack through the graph to get the elements (There are no cycles), so we just find all the possible paths through the graph
//
//        set<vector<int>> res_path;
//        set<string> res_word;
//
//        for(auto s: this -> startNodes){
//            set<int> visited;
//            vector<int> path {s};
//            dfs(s, path, "");
//        }
//
//        return res_word;
//    }
    
    // Determinization of Recognisers
    
    void next(int node, int letter, vector<int> &nextStates){
        
        for(auto el: (this->table)[node][letter]){
            nextStates.push_back(el.second);
        }
    
    }
    
    vector<int> nextSet(vector<int> &S, int letter){

        // S - set of states; T - set of next states;

        vector<int> T;
        
        for(auto q: S){
            this->next(q, letter, T);
        }
        
        return this->closure(T);
    }
    
    vector<int> closure(vector<int> closure){
        
        if(closure.size() != 0){
            stack<int> S;
            set<int> visited;
            
            for(auto el: closure){
                S.push(el);
            }
            
            while(!S.empty()){
                
                int node = S.top();
                S.pop();
                visited.insert(node);
                
                for(auto e: (this->table)[node][0]){
                    
                    int out_node = e.second;
                    
                    // if the element is not in visited
                    if (visited.find(out_node) == visited.end()){
                        closure.push_back(out_node);
                        S.push(out_node);
                    }
                }
            }
        }
        
        return closure;
    }
    
    void explore(map<vector<int>, int> &T, vector<int> &S, Transducer &B){
        
        for(int c = 1; c <= B.inputLetters; c++){
            
            vector<int> U = this->nextSet(S, c);
            
            // check if all the elements of U are unique.
            set<int> USet( U.begin(), U.end());
            U.assign(USet.begin(), USet.end());
            
            assert(USet.size() == U.size());
            
            sort(U.begin(), U.end());
            
            if(T.find(U) != T.end()){

                if(T.size() > B.table.size()){
                    
                    int temp = B.table.size();
                    
                    B.table.resize(T.size());
                
                    for(int i = temp; i < T.size(); i++){
                        B.table[i].resize(B.inputLetters + 1);
                    }
                }
                
                B.table[T.at(S)][c].push_back({0, T.at(U)});
            }else{
                
                T.insert({U, T.size()});
                
                if(T.size() > B.table.size()){
                    
                    int temp = B.table.size();
                    B.table.resize(T.size());
                    
                    
                    for(int i = temp; i < T.size(); i++){
                        B.table[i].resize(B.inputLetters + 1);
                    }
                }
                
                B.table[T.at(S)][c].push_back({0, T.at(U)});
            
                if(U.size() != 0){
                    this->explore(T, U, B);
                }
            }
        }
    }
    
    Transducer determinize(){
        
        Transducer B = Transducer(this->inputLetters, this->outputLetters);
        
        vector<int> I = this -> closure(this -> startNodes);
        
        sort(I.begin(), I.end());
        
        map<vector<int>, int> T;
        T.insert({I, T.size()});
        
        this->explore(T, I, B);
        
        // Start & Final Nodes Specification
        
        if(B.table.size()){
            
            B.startNodes.push_back(T.at(I));
            
            for (const auto &el : T) {
                for(int i = 0; i < (this->finalNodes).size(); i++){
                    if(count(el.first.begin(), el.first.end(), i) && count((this->finalNodes).begin(), (this->finalNodes).end(), i)){
                        B.finalNodes.push_back(el.second);
                    }
                }
            }
        }
        
        return B;
    }
    
    // Minimization of Recognisers
    
    Transducer mergeAlph(){

        int newAlph = ((this -> inputLetters) + 1) * ((this -> outputLetters) + 1) - 1;
        
        auto rec = Transducer((this -> table).size(), newAlph, 0, this -> startNodes, this -> finalNodes);
        
        for(int i = 0; i < (this -> table).size(); i++){
            for(int j = 0; j <= this -> inputLetters; j++){
                for(auto e: (this->table)[i][j]){
                    rec.table[i][(e.first * (this->inputLetters + 1) + j)].push_back({0, e.second});
                }
            }
        }

        return rec;
    }
    
    Transducer splitAlph(int inputLetters, int outputLetters){
        
        auto split = Transducer((this -> table).size(), inputLetters, outputLetters, this -> startNodes, this -> finalNodes);
        
        for(int i = 0; i < (this -> table).size(); i++){
            for(int j = 0; j <= this -> inputLetters; j++){
                for(auto e: (this->table)[i][j]){
                    
                    int inpLetter = j % (inputLetters + 1);
                    
                    int outLetter = j / (inputLetters + 1);
                    
                    split.table[i][inpLetter].push_back({outLetter, e.second});
                }
            }
        }
        
        return split;
    }
    
    Transducer minimize(){
        
        if(this -> outputLetters != 0){
            // merge the alphabets
            auto merge = this -> mergeAlph();
            auto rev1 = merge.transpose();
            auto det1 = rev1.determinize();
            auto rev2 = det1.transpose();
            auto det2 = rev2.determinize();
            // split the alphabets
            auto split = det2.splitAlph(this->inputLetters, this->outputLetters);
            return split;
        }
        else{
            auto rev1 = this -> transpose();
            auto det1 = rev1.determinize();
            auto rev2 = det1.transpose();
            auto det2 = rev2.determinize();
            return det2;
        }
    }
    
    // Automata completion - The implemented algorithm already produces a complete dfa - The function is therefore not necessary for our implementation
    
    Transducer complete(){
        
        // the + 1 adds the dead state
        auto completeDFA = Transducer((int)(this->table).size() + 1, this->inputLetters, this->outputLetters, this -> startNodes, this -> finalNodes);
        
        // check if there are any vectors of the form table[i][c] such that their size is 0 and then make a transition to the "dead-state S" from them;
        
        int dead_state = (int)(this->table).size();
        
        // add an extra state
        
        for(int i = 0; i < (this->table).size(); i++){
            for(int c = 0; c <= this -> inputLetters; c++){
                if(table[i][c].size() == 0){
                    completeDFA.table[i][c].push_back({0, dead_state});
                }
            }
        }
        
        completeDFA.finalNodes = this -> finalNodes;
        
        return completeDFA;
        
    }
    
    Transducer complement(){
        
        // The final determinization produces a complete dfa
        auto complement = this -> minimize();
        
        // Flip the final nodes vector;
        
        vector<int> invertFinal;
        
        for(int i = 0; i < (complement.table).size(); i++){
            if(count((complement.finalNodes).begin(), (complement.finalNodes).end(), i) == 0){
                invertFinal.push_back(i);
            }
        }
        
        complement.finalNodes = invertFinal;
        
        return complement;
    }
    
    // It is sufficient to implement the conversion to a filter
    Transducer RtF(){
        
        Transducer convert = Transducer((this -> table).size(), this->inputLetters, this->inputLetters, this->startNodes, this->finalNodes);
        
        for(int i = 0; i < (this->table).size(); i++){
            for(int j = 0; j <= convert.inputLetters; j++){
                for(auto edge: (this -> table)[i][j]){
                    convert.addEdge(j, j, i, edge.second);
                }
            }
        }
        
        return convert;
    }
    
    Transducer FtR(){
        
        Transducer convert = Transducer((this -> table).size(), this->inputLetters, 0, this->startNodes, this->finalNodes);
        
        for(int i = 0; i < (this->table).size(); i++){
            for(int j = 0; j <= convert.inputLetters; j++){
                for(auto edge: (this -> table)[i][j]){
                    convert.addEdge(j, 0, i, edge.second);
                }
            }
        }
        
        return convert;
    }
    
    Transducer invert(){
        
        Transducer invert = Transducer((this -> table).size(), this->outputLetters, this->inputLetters, this->startNodes, this->finalNodes);
        
        for(int i = 0; i < (this->table).size(); i++){
            for(int j = 0; j <= this -> inputLetters; j++){
                for(auto edge: (this -> table)[i][j]){
                    invert.addEdge(edge.first, j, i, edge.second);
                }
            }
        }
        
        return invert;
    }
    
    bool languageEquality(){
        
        //finding an isomorphism between finite automata;
        
        return true;
        
    }
    
    // Check if there's a cycle in the directed graph;
    bool cycleCheck(){
        
        return true;
    }
    
};

// Implementation of the transducers used in the proof.

Transducer multimark(){
    
    // intialize a vector of bool instead ...
    vector<int> sn {0};
    vector<int> fn {0};
    
    Transducer multimark = Transducer(1, 5, 5, sn, fn);
    
    for(int i = 1; i <= multimark.inputLetters; i++){
        multimark.addEdge(i, i, 0, 0);
    }
    
    multimark.addEdge(0, 5, 0, 0);
    
    return multimark;
}

Transducer singlemark(){
    
    vector<int> sn {0};
    vector<int> fn {0, 3};
    
    Transducer singlemark = Transducer(4, 4, 5, sn, fn);
    
    for(int i = 1; i <= singlemark.inputLetters; i++){
        singlemark.addEdge(i, i, 0, 1);
        singlemark.addEdge(i, i, 1, 1);
        singlemark.addEdge(i, i, 2, 2);
        singlemark.addEdge(i, i, 2, 3);
    }
    
    singlemark.addEdge(0, 5, 1, 2);
    
    return singlemark;
}

Transducer scissors(){
    
    vector<int> sn {0};
    vector<int> fn {2};
    
    Transducer scissors = Transducer(3, 5, 4, sn, fn);
    
    for(int i = 1; i <= scissors.inputLetters; i++){
        scissors.addEdge(i, 0, 0, 0);
        scissors.addEdge(i, 0, 2, 2);
        if(i != 5){
            scissors.addEdge(i, i, 1, 1);
        }
    }
    
    scissors.addEdge(5, 0, 0, 1);
    scissors.addEdge(5, 0, 1, 2);
    
    return scissors;
}

// Use the traversal algorithm to check if the automata is doing what it is supposed to .;..

Transducer audioactiveT(){
    
    // it makes sense to place them first i.e
    vector<int> sn {0, 1, 2, 3};
    
    Transducer audioactiveT = Transducer(24, 4, 4, sn, sn);
    
    for(int c = 1; c <= 4; ++c){
        
        audioactiveT.addEdge(c, 1, c - 1, 5 * (c-1) + 1 + 3);
        audioactiveT.addEdge(0, c, 5 * (c-1) + 1 + 3, 5 * (c-1) + 1 + 4);
        
        audioactiveT.addEdge(c, 2, c - 1, 5 * (c-1) + 1 + 5);
        audioactiveT.addEdge(c, c, 5 * (c-1) + 1 + 5, 5 * (c-1) + 1 + 4);
        
        audioactiveT.addEdge(c, 3, c - 1, 5 * (c-1) + 1 + 6);
        audioactiveT.addEdge(c, c, 5 * (c-1) + 1 + 6, 5 * (c-1) + 1 + 7);
        audioactiveT.addEdge(c, 0, 5 * (c-1) + 1 + 7, 5 * (c-1) + 1 + 4);
    }
    
    for(int c = 1; c <= 4; ++c){
        for(int i = 0; i <= 3; ++i){
            if(i != c-1){
                audioactiveT.addEdge(0, 0, 5 * (c-1) + 1 + 4, i);
            }
        }
    }
    
    return audioactiveT;
}

Transducer augmentedAudioactiveT(){
    
    Transducer aat = audioactiveT();
    
    aat.inputLetters = 5; aat.outputLetters = 5;
    
    for(int i = 0; i <= 3; i++){
        aat.addEdge(5, 5, i, i);
    }
    
    return aat;
}

Transducer splitRec(){

    auto aat = augmentedAudioactiveT();
    auto split = aat.FtR().minimize();
    
    for(int i = 0; i < 9; i++){
        split = aat.compose(split).minimize();
    }
    
    return split;
}

Transducer irredFactorRec(){

    auto sm = singlemark().minimize();
    auto sr = splitRec().minimize();
    auto iwr = sm.compose(sr).complement();
    
    return iwr;
}

Transducer irredFactorDer(){

    auto mmt = multimark();
    auto sf = splitRec().RtF();
    auto sc = scissors();
    auto isf = irredFactorRec().RtF();

    // definition of irredFactorDer - There's still the problem of minimization of derivators
    
    auto ifd = mmt.compose(sf).minimize().compose(sc).minimize().compose(isf).minimize();

    return ifd;
}

// Theorem Proofs

Transducer Theorem2(){
    
    auto aat = augmentedAudioactiveT();
    auto split = aat.FtR().minimize();
    
    for(int i = 0; i < 9; i++){
        split = aat.compose(split).minimize();
    }
    
    return split;

}

set<string> CosmologicalTheorem(){

    auto aat = augmentedAudioactiveT();

    auto factorizer = irredFactorDer();
    
    // Deal with the inverted edges
    
    auto T = aat.compose(factorizer).minimize();
    
    // This is supposed to be the generator - transposed to a recogniser
    auto Tn = T.invert().FtR().minimize();
    
    auto Tn_prev = Tn;
    
    for(int i = 0; i < 26; i++){
        
        cout << "Iteration: " << i << endl;
        
        // Test which one is faster
        auto T_inv = T.invert();
        Tn = T_inv.compose(Tn).minimize();
//        Tn = fact_inv.compose(aat_inv).minimize().compose(Tn).minimize();
        
        if(Tn_prev.table.size() == Tn.table.size()){
            
            // perform the equivalence check
            cout << "Potential Isomorphism" << endl;
            
        }
        
        Tn_prev = Tn;
    }
    
    // for the final step what do we do ?
    
    set<string> a {"placeholder for elements"};
    return a;
}

int main(int argc, const char * argv[]){
    
    auto start = chrono::steady_clock::now();
    
    auto sr = splitRec();
    
    // examples
    
    auto res = sr.traverse("11132513211");
    
    if(res.size()){
        cout << "Valid" << endl;
    }else{
        cout << "Invalid" << endl;
    }
    
    res = sr.traverse("11135213211");
    
    if(res.size()){
        cout << "Valid" << endl;
    }else{
        cout << "Invalid" << endl;
    }
    
    auto end = chrono::steady_clock::now();

    // Store the time difference between start and end
    auto diff = end - start;
    
    cout << chrono::duration <double, milli> (diff).count() << " ms" << endl;
    
}

