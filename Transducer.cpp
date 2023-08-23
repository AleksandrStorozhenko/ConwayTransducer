#include <iostream>
#include <vector>
#include <stack>
#include <set>
#include <map>
#include <algorithm>

using namespace std;

#define endl "\n";

struct Transducer{
    
    int inputLetters, outputLetters;
    
    vector<bool> startNodes, finalNodes;
    
    vector<vector<vector<pair<int, int>>>> table;

    Transducer(int inputLetters, int outputLetters){
        // 0 -> e; 1 - 1; 2 - 2; 3 - 3; 4 - n; 5 - .;
        this -> inputLetters = inputLetters;
        this -> outputLetters = outputLetters;
        
    };

    Transducer(int stateCount, int inputLetters, int outputLetters, vector<bool> &startNodes){
        this -> inputLetters = inputLetters;
        this -> outputLetters = outputLetters;
        
        (this -> table).resize(stateCount);
        
        for(int i = 0; i < stateCount; ++i){
            (this -> table)[i].resize((this->inputLetters)+1);
        }
        
        this -> startNodes = startNodes;
        
        vector<bool> finalNodes;
        finalNodes.resize(stateCount);
        
        this -> finalNodes = finalNodes;
    };
    
    Transducer(int stateCount, int inputLetters, int outputLetters, vector<bool> &startNodes, vector<bool> &finalNodes){
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
        
//        (this->table)[A].resize(inputLetters+1);
        (this->table)[A][inpchar].push_back({outchar, B});
        
    }
    
    Transducer compose(Transducer &T1){
        
        //number of nodes in the 2 transducers
        auto X = (this -> table).size();
        auto Y = T1.table.size();
        
        vector<bool> initialNodes;
        initialNodes.resize(X*Y);
        
        for(auto node1: this->startNodes){
            for(auto node2: T1.startNodes){
                auto pos = node2 * X + node1;
                initialNodes[pos] = true;
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
        
        //TODO: Add Final Nodes
        
        vector<bool> finalNodes;
        finalNodes.resize(X*Y);

        for(int i = 0; i < X; i++){
            for(int j = 0; j < Y; j++){
                // the trsdComp.table[j * X + i] entry is not necessarily defined;
                if((this->finalNodes)[i] && T1.finalNodes[j] && trsdComp.table[j * X + i].size() != 0){
                    finalNodes[j * X + i] = true;
                }
            }
        }
        
        //TODO: Condense the number of states
        
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
    
    // Determinization of Recognisers
    
    vector<bool> next(int node, int letter){
        
        vector<bool> nextStates;
        
        nextStates.resize((this->table).size());
        
        for(auto el: (this->table)[node][letter]){
            nextStates[el.second] = true;
        }
        
        return nextStates;
    }
    
    vector<bool> nextSet(vector<bool> &S, int letter){

        // S - set of states; T - set of next states;

        vector<bool> T;
        T.resize((this->table).size());
        
        for(int q = 0; q < S.size(); q++){
            
            if(S[q]){
                vector<bool> nextStates = (this->next(q, letter));
                
                for(int i = 0; i < nextStates.size(); i++){
                    if(T[i] || nextStates[i]){
                        T[i] = true;
                    }
                }
            }
        }
        
        return this->closure(T);
    }
    
    vector<bool> closure(vector<bool> &t){
        
        //closure should be a copy of t i.e passed by value, rather than by reference
        vector<bool> closure = t;
        
        // iterate over the values to find a true
        
        int size = 0;
        for(auto el: closure){
            if(el){
                size = 1;
                break;
            }
        }
        
        if(size == 0){
            return closure;
        }
        
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
                    closure[out_node] = true;
                    S.push(out_node);
                }
            }
        }
        
        return closure;
    }
    
    void explore(map<vector<bool>, int> &T, vector<bool> &S, Transducer &B){
        
        for(int c = 1; c <= B.inputLetters; c++){
            
            vector<bool> U = this->nextSet(S, c);
            
            // check if the key is already present in the T map
            
            bool empty = true;
            
            for(int i = 0; i < U.size(); i++){
                if(U[i]){
                    empty = false;
                    break;
                }
            }
            
            // Because we haven't yet set the size - edge case
            
            if(T.find(U) != T.end()){

                // resizing of the vector table
                if(T.size() > B.table.size()){
                    
                    B.table.resize(T.size());
                    
                    //TODO: This can be optimized
                    
                    for(int i = 0; i < B.table.size(); i++){
                        B.table[i].resize(B.inputLetters + 1);
                    }
                }
                
                B.table[T.at(S)][c].push_back({0, T.at(U)});
            }else{
                
                T.insert({U, T.size()});
                
                // resizing of the vector table
                if(T.size() > B.table.size()){
                    
                    B.table.resize(T.size());
                    
                    //TODO: This can be optimized
                    
                    for(int i = 0; i < B.table.size(); i++){
                        B.table[i].resize(B.inputLetters + 1);
                    }
                }
                
                B.table[T.at(S)][c].push_back({0, T.at(U)});
                
                // this should technically ensure that U = {true, true} is only explored once
                if(!empty){
                    this->explore(T, U, B);
                }
            }
        }
    }
    
    Transducer determinize(){
        
        Transducer B = Transducer(this->inputLetters, this->outputLetters);
        
        vector<bool> I = this -> closure(this -> startNodes);
        
        map<vector<bool>, int> T;
        T.insert({I, T.size()});
        
        this->explore(T, I, B);
        
        // Start & Final Nodes Specification
        
        vector<bool> startNodes, finalNodes;
        
        if(B.table.size()){
            
            startNodes.resize(B.table.size());
            finalNodes.resize(B.table.size());
            
            startNodes[T.at(I)] = true;
            B.startNodes = startNodes;
            
            for (const auto &el : T) {
                for(int i = 0; i < (this->finalNodes).size(); i++){
                    if(el.first[i] && (this->finalNodes)[i]){
                        finalNodes[el.second] = true;
                    }
                }
            }
            
            B.finalNodes = finalNodes;
        }
        
        return B;
    }
    
    // Traversal (Helps with testing)
    
    void traverse(){
        
        // A simple backtracking function with the purpose of finding all the possible words in the path
        
    }
    
    // Minimization of Recognisers
    
    Transducer minimize(){
    
        auto rev1 = this -> transpose();
        auto det1 = rev1.determinize();
        auto rev2 = det1.transpose();
        auto det2 = rev2.determinize();
        
        return det2;
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
        
        vector<bool> invertFinal;
        invertFinal.resize((this -> finalNodes).size());
        
        for(int i = 0; i < invertFinal.size(); i++){
            if(!(this -> finalNodes)[i]){
                invertFinal[i] = true;
            }
        }
        
        complement.finalNodes = invertFinal;
        
        return complement;
    }
    
    Transducer convert(int type_in, int type_out){
        
        Transducer convert = Transducer(5, 5);
        // The function requires the specification of the 2 transducer types.
        
        return convert;
    }
    
    // Language equality check
    
    bool languageEquality(){
        
        //finding an isomorphism between finite automata;
        
        return true;
        
    }
    
    void getElements(){
        
    }
    
};

// Implementation of the transducers used in the proof.

Transducer multimark(){
    
    vector<bool> sn {true};
    vector<bool> fn {true};
    
    Transducer multimark = Transducer(1, 4, 5, sn, fn);
    
    for(int i = 1; i <= multimark.inputLetters; i++){
        multimark.addEdge(i, i, 0, 0);
    }
    
    multimark.addEdge(0, 5, 0, 0);
    
    return multimark;
}

Transducer singlemark(){
    
    vector<bool> sn {true, false, false, false};
    vector<bool> fn {true, false, false, true};
    
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
    
    vector<bool> sn {true, false, false};
    vector<bool> fn {false, false, true};
    
    Transducer scissors = Transducer(3, 5, 4, sn, fn);
    
    for(int i = 1; i <= scissors.inputLetters; i++){
        scissors.addEdge(i, 0, 0, 0);
        scissors.addEdge(i, 0, 2, 2);
        if(i != 5){
            scissors.addEdge(i, 0, 1, 1);
        }
    }
    
    scissors.addEdge(5, 0, 0, 1);
    scissors.addEdge(5, 0, 1, 2);
    
    return scissors;
}

Transducer audioactiveT(){
    
    // it makes sense to place them first i.e
    vector<bool> sn;
    sn.resize(24);
    for(int i = 0; i <= 3; ++i){
        sn[i] = true;
    }
    
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
    
    auto aat = augmentedAudioactiveT().minimize();
    
    //Convert to a recogniser
    auto split = aat.convert(2, 0).minimize();
    
    for(int i = 0; i < 9; ++i){
        split = aat.compose(split).minimize();
    }
    
    return split;
}

Transducer irredFactorRec(){

    auto sm = singlemark();
    auto sr = splitRec();
    auto iwr = sm.compose(sr).minimize().complement();
    
    return iwr;
}

Transducer irredFactorDer(){
    
    auto mmt = multimark();
    // Convert to a filter
    auto sf = splitRec().convert("Filt");
    auto sc = scissors();
    auto isf = irredFactorRec().convert("Filt");
    
    // definition of irredFactorDer
    auto ifd = mmt.compose(sf).minimize().compose(sc).minimize().compose(isf).minimize();
    
    return ifd;
}

// Theorem Proofs

void CosmologicalTheorem(){
    
    auto aat = augmentedAudioactiveT().minimize();
    
    // Multimark transducer
    auto mmt = multimark();
    // Split filter
    auto sf = splitRec().convert("Filt");
    // Scissors
    auto sc = scissors();
    // Irreducible string filter
    auto isf = irredFactorRec().convert("Filt");

    auto factorizer = irredFactorDer().minimize();

    auto T = aat.compose(factorizer);

    auto Tn = T.convert("Gen").minimize();

    auto Tn_prev = Tn;
    
    int counter = 0;
    
    while(true){
        
        counter++;
        
        Tn = Tn.compose(aat).minimize().compose(mmt).minimize().compose(sf).minimize().compose(sc).minimize().compose(isf).minimize();
        
        // preliminary check if the number of nodes is the same
        
        if(Tn.table.size() == Tn_prev.table.size()){
            if(Tn.languageEquality(Tn_prev)){
                cout << "The language stabilizes after: " << counter << " iterations" << endl;
                break;
            }
        }
        
        Tn_prev = Tn;
    }
    
    auto[paths, words] = Tn.getElements();
    
    cout << "Number of elements: " << words.size() << endl;
    
    return words;
}

int main(int argc, const char * argv[]){
    auto start = chrono::steady_clock::now();
    
    vector<bool> sn {true, false, false};
    
    vector<bool> fn {false, false, true};
    
    Transducer T1 = Transducer(3, 2, 0, sn, fn);
    
    T1.addEdge(1, 0, 0, 1);
    T1.addEdge(2, 0, 0, 3);

    T1.addEdge(1, 0, 1, 0);
    T1.addEdge(2, 0, 1, 3);

    T1.addEdge(1, 0, 2, 1);
    T1.addEdge(2, 0, 2, 4);

    T1.addEdge(1, 0, 3, 5);
    T1.addEdge(2, 0, 3, 5);

    T1.addEdge(1, 0, 4, 3);
    T1.addEdge(2, 0, 4, 3);

    T1.addEdge(1, 0, 5, 5);
    T1.addEdge(2, 0, 5, 5);
    
    auto dfa = T1.determinize();
    
    auto end = chrono::steady_clock::now();

    // Store the time difference between start and end
    auto diff = end - start;
    
    cout << chrono::duration <double, milli> (diff).count() << " ms" << endl;
    
}
