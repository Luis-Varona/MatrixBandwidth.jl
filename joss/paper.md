---
title: "MatrixBandwidth.jl: Fast algorithms for matrix bandwidth minimization and recognition"
tags:
  - matrix bandwidth
  - sparse matrices
  - optimization
  - scientific computing
  - Julia
authors:
  - name: Luis M. B. Varona
    orcid: 0009-0003-7784-5415
    affiliation: "1,2,3"
affiliations:
  - name: Department of Politics and International Relations, Mount Allison University
    index: 1
  - name: Department of Mathematics and Computer Science, Mount Allison University
    index: 2
  - name: Department of Economics, Mount Allison University
    index: 3
date: 27 October 2025
bibliography: paper.bib
---

# Summary

The *bandwidth* of an $n \times n$ matrix $A$ is the minimum non-negative integer $k \in
\{0, 1, \ldots, n - 1\}$ such that $A_{i,j} = 0$ whenever $\lvert i - j \rvert > k$. Reordering the
rows and columns of a matrix to reduce its bandwidth has many practical applications in engineering
and scientific computing: it can improve performance when solving linear systems, approximating
partial differential equations, optimizing circuit layout, and more [@Maf14]. There are two variants
of this problem: *minimization*, which involves finding a permutation matrix $P$ such that the
bandwidth of $PAP^\mathsf{T}$ is minimized, and *recognition*, which entails determining whether
there exists a permutation matrix $P$ such that the bandwidth of $PAP^\mathsf{T}$ is less than or
equal to some fixed non-negative integer (an optimal permutation that fully minimizes the bandwidth
of $A$ is not required). Accordingly,
[MatrixBandwidth.jl](https://github.com/Luis-Varona/MatrixBandwidth.jl) offers fast algorithms for
matrix bandwidth minimization and recognition. Julia's combination of easy syntax and high
performance, along with its rapidly growing ecosystem for scientific computing, made it the ideal
language of choice for this project.

## Example

Consider the following $60 \times 60$ sparse matrix with initial bandwidth $51$:

\begin{figure}[H]
  \centering
  \includegraphics[height=2in]{assets/A.png}
  \caption{Original $60 \times 60$ matrix with bandwidth $51$}
  \label{fig:A}
\end{figure}

MatrixBandwidth.jl can both recognize whether the minimum bandwidth of $A$ is less than or equal to
some fixed integer (\autoref{fig:A_rec}) and actually minimize the bandwidth of $A$
(\autoref{fig:A_min}):

\begin{figure}[H]
  \begin{minipage}[b]{.475\textwidth}
    \centering
    \includegraphics[height=1.5in]{assets/A_rec.png}
    \caption{The matrix with bandwidth recognized as $\le 6$ via the Del Corso--Manzini algorithm}
    \label{fig:A_rec}
  \end{minipage}\hfill
  \begin{minipage}[b]{.475\textwidth}
    \centering
    \includegraphics[height=1.5in]{assets/A_min.png}
    \caption{The matrix with bandwidth minimized to $5$ via the Gibbs--Poole--Stockmeyer algorithm}
    \label{fig:A_min}
  \end{minipage}
\end{figure}

Note that since Gibbs&ndash;Poole&ndash;Stockmeyer is a heuristic algorithm, $5$ may not be the
*true* minimum bandwidth of $A$, but it is likely close.

## Algorithms

The following matrix bandwidth reduction algorithms are currently available:

- Minimization
  - Exact
    - Caprara&ndash;Salazar-González [@CS05]
    - Del Corso&ndash;Manzini [@DM99]
    - Del Corso&ndash;Manzini with perimeter search [@DM99]
    - Saxe&ndash;Gurari&ndash;Sudborough [@Sax80; @GS84]
    - Brute-force search
  - Heuristic
    - Gibbs&ndash;Poole&ndash;Stockmeyer [@GPS76]
    - Cuthill&ndash;McKee [@CM69]
    - Reverse Cuthill&ndash;McKee [@CM69; @Geo71]
- Recognition
  - Caprara&ndash;Salazar-González [@CS05]
  - Del Corso&ndash;Manzini [@DM99]
  - Del Corso&ndash;Manzini with perimeter search [@DM99]
  - Saxe&ndash;Gurari&ndash;Sudborough [@Sax80; @GS84]
  - Brute-force search

Recognition algorithms determine whether any row-and-column permutation of a matrix induces
bandwidth less than or equal to some fixed integer. Exact minimization algorithms always guarantee
optimal orderings to minimize bandwidth, while heuristic minimization algorithms produce
near-optimal solutions more quickly. Metaheuristic minimization algorithms employ iterative search
frameworks to find better solutions than heuristic methods (albeit more slowly); no such algorithms
are already implemented, but several (e.g., simulated annealing) are currently under development.

Thus far, the Caprara&ndash;Salazar-González algorithms are the only ones implemented that require
integer linear programming; it is for these that the [JuMP.jl](https://github.com/jump-dev/JuMP.jl)
package [@LDD+23] is included as a dependency.

# Statement of need

Many matrix bandwidth reduction algorithms exist in the literature, but implementations in the
open-source ecosystem are scarce, with those that do exist primarily tackling older, less efficient
algorithms. The [Boost](https://www.boost.org/) libraries in C++ [@LLS+01], the
[NetworkX](https://networkx.org/) library in Python [@Net25], and the MATLAB standard library
[@MAT25] all only implement the aforementioned reverse Cuthill&ndash;McKee algorithm from 1971.
In Julia, the only other relevant packages identified by the author are
[BandedMatrices.jl](https://github.com/JuliaLinearAlgebra/BandedMatrices.jl) [@Jul16] and
[SymRCM.jl](https://github.com/PetrKryslUCSD/SymRCM.jl) [@Krys20], both of which also only implement
reverse Cuthill&ndash;McKee as their sole bandwidth reduction algorithm.

Furthermore, not enough attention is given to recognition algorithms or exact minimization
algorithms. Although more performant modern alternatives are often neglected, at least reverse
Cuthill&ndash;McKee is a widely implemented method of approximating a minimal bandwidth ordering (as
noted above). However, no such functionality for recognition or exact minimization is widely
available, requiring researchers with such needs to fully re-implement these algorithms themselves.

These two gaps in the ecosystem not only make it difficult for researchers to benchmark and compare
new proposed algorithms but also preclude the application of the most performant modern algorithms
in real-life industry settings. MatrixBandwidth.jl aims to bridge this gap by presenting a unified
interface for matrix bandwidth reduction algorithms in Julia.

# Research applications

The author either has used or is using MatrixBandwidth.jl to do the following:

- Develop a new polynomial-time algorithm for "bandwidth $\le k$" recognition efficient for both
  small and large $k$, and benchmarking it against other approaches [@Sax80; @GS84]
- Speed up $k$-coherence checks of quantum states in many cases by confirming that the density
  matrix's minimum bandwidth is greater than $k$ [@JMP25]
- Compute the spectral graph property of "$S$-bandwidth" [@JP25] via the
  [SDiagonalizability.jl](https://github.com/GraphQuantum/SDiagonalizability.jl) package [@VJP25],
  which depends critically on MatrixBandwidth.jl for bandwidth recognition
- Investigate the precise performance benefits of reducing the propagation graph's bandwidth when
  training a recurrent neural network, building on @BvMM+19

The first three use cases rely on the recognition and exact minimization functionality unique to
MatrixBandwidth.jl (indeed, they largely motivated the package's development). The last (ongoing)
research project *could* be facilitated by SymRCM.jl instead, but the author intends to use more
performant metaheuristic minimization algorithms currently under development when producing the
final computational results, as well as use recognition algorithms to minimize bandwidth to various
target levels when quantifying performance improvements.

# Limitations

Currently, MatrixBandwidth.jl's core functions generically accept any input of the type
`AbstractMatrix{<:Number}`, not behaving any differently when given sparsely stored matrices (e.g.,
from the [SparseArrays.jl](https://github.com/JuliaSparse/SparseArrays.jl) standard library
package). Capabilities for directly handling graph inputs (aiming to reduce the matrix bandwidth of
a graph's adjacency) are also not available. Given that bandwidth reduction is often applied to
sparse matrices and graphs, this will be addressed in future releases.

Moreover, many of the algorithms only apply to structurally symmetric matrices (i.e., those whose
nonzero pattern is symmetric). However, this is a limitation of the algorithms themselves, not the
package's implementation. Future releases with metaheuristic algorithms will include more methods
that accept structurally asymmetric inputs.

# Conflict of interests

The author declares no conflict of interest.

# Acknowledgements

I owe much to my research supervisors&mdash;Nathaniel Johnston, Sarah Plosker, and Craig
Brett&mdash;for supporting and guiding me in my work. I would also like to thank Liam Keliher,
Peter Leli&egrave;vre, and Marco Cognetta for useful discussions. Finally, credit for
MatrixBandwidth.jl's telepathic-cat-and-turtle logo goes to Rebekka Jonasson.

# References
