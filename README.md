# Conway’s Cosmological Theorem and Automata Theory: A C++ Implementation

This repository provides the C++ implementation supporting the computational proofs detailed in the paper *"Conway’s Cosmological Theorem and Automata Theory"*, which has been submitted to the *American Mathematical Monthly* for consideration. The work introduces a novel approach to Conway’s theorem using automata theory, providing a systematic understanding that was not achieved by previous *ad hoc* methods. This implementation complements the theoretical developments by offering a rigorous computational demonstration of the key results.

## Overview

We implement a key automaton-based structure—the **Transducer**—which models finite-state automata essential to proving both Conway’s **Cosmological Theorem** and the **Splitting Theorem**. The automata approach modernizes the previously manual methods used to prove these theorems, providing anruntime efficient computational demonstration that verifies Conway's result of audioactive sequence stabilization after 24 iterations.

## Core Concepts

### 1. **Transducer Object**
Central to the implementation, the Transducer struct models finite-state machines. This data structure supports a range of operations necessary for the automata-theoretic proof:

- **Transition Management**: The `addEdge` method allows for the definition of transitions between states.
- **Transducer Operations**: Functions like `compose`, `determinize`, `reverse`, and `minimize` enable the verification of the theorem's results through composition and minimization of transducers.

### 2. **Fundamental Transducers**

We begin by defining a series of useful "example" transducers. These transducers will serve as the fundamental building blocks for the more complex automata, i.e all others will derive from these via composition.
- **`Multimark()`**: Inserts marks nondeterministically into sequences.
- **`Singlemark()`**: Inserts a single mark deterministically.
- **`Scissors()`**: Extracts subwords from sequences, essential for proving the splitting properties.
- **`Audioactive(bool augmented)`**: Simulates the behavior of the look-and-say sequence, following the transformations detailed in Conway’s theorem. Depending on the value of the augmented boolean argument, the function constructs either the standard audioactive transducer or its augmented version.

From these 5 relatively simple automata, we will be able to prove the **Cosmological Theorem**.

### 3. **Cosmological Theorem: Computational Proof**

The `main()` function implements the proof of the **Cosmological Theorem** through the following steps:
- **Splitting Theorem** *(runtime: 50 ms)*
- **Cosmological Theorem** *(runtime: 150 ms)*

## Running the Code

### Prerequisites

- A C++ compiler such as `g++` is required to compile the source code.

### Compilation and Execution

1. Clone the repository:
    ```bash
    git clone https://github.com/AleksandrStorozhenko/ConwayTransducer
    ```
2. Navigate to the directory:
    ```bash
    cd ConwayTransducer
    ```
3. Compile the code:
    ```bash
    g++ -std=c++17 -o conway Transducer.cpp
    ```
4. Run the program:
    ```bash
    ./conway
    ```

## Contribution

This codebase complements the theoretical work described in the paper, providing a clear demonstration of automata theory in action. Contributions to further optimize or extend the code are welcomed; feel free to open an issue or submit a pull request.

For questions or detailed discussions, contact us via the repository’s issue tracker.

# References

This repository is based on the theoretical framework presented in the paper "Conway’s Cosmological Theorem and Automata Theory". We encourage readers interested in the detailed proofs and theoretical background to consult the original paper for a deeper understanding.
