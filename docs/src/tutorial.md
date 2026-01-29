```@meta
CurrentModule = MatrixBandwidth
```

# Tutorial

MatrixBandwidth.jl offers unified interfaces for both bandwidth minimization and bandwidth recognition via the `minimize_bandwidth` and `has_bandwidth_k_ordering` functions, respectively—the algorithm itself is specified as an argument. Comprehensive documentation for these two methods, as well as their output structs, is available in the [public API documentation](public_api.md). We go over some basic examples below.

## Getting started

Throughout these examples, we shall work with a random sparse, structurally symmetric matrix of size $60 \times 60$ with initial bandwidth $51$:

```julia-repl
julia> using Random, SparseArrays

julia> Random.seed!(8675309);

julia> A = sprand(60, 60, 0.01); A = A + A' # Ensure structural symmetry
60×60 SparseMatrixCSC{Float64, Int64} with 93 stored entries:
⎡⠀⠀⠀⠀⠠⠀⠀⠀⠀⠀⠀⠀⠠⠀⢀⠀⠀⠒⠀⠀⠀⠀⠀⠀⡀⠨⠀⠀⠀⠀⎤
⎢⠀⠀⠀⠀⠅⠀⠀⠀⠀⠀⠐⠀⠠⡀⢀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠂⠠⠀⠀⠀⠀⎥
⎢⠀⠂⠁⠁⢀⠐⠀⠂⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠠⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠁⠀⡀⠀⠀⠀⠀⠀⠀⠈⠁⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠐⠀⠀⠀⠀⠀⠀⠀⠀⠐⠀⠀⠈⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠐⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠁⠀⠄⡀⠀⠀⎥
⎢⠀⠂⠀⠢⠀⠀⠀⠀⠀⠀⠀⠄⠀⠀⠀⠀⢀⠀⠁⠀⠀⠀⠐⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠐⠀⠐⠀⠀⠇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠐⠁⠀⠀⎥
⎢⢠⠀⠀⠀⠀⢀⠀⠠⢀⠀⠀⠀⠀⠐⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠂⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠐⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠂⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠊⠀⠀⢠⠀⠀⠀⠀⠀⠠⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠐⠀⠀⠀⠀⠀⠀⠀⠀⣀⠀⠀⠁⠀⠀⠀⠄⠀⎥
⎢⡀⡈⠈⡀⠀⠀⠆⠀⠀⠀⠁⠀⠀⠀⠀⠀⠈⠀⠀⠀⠀⠀⠁⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠡⠀⠀⠔⠀⠀⠀⠐⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡀⠀⠁⠀⠀⠀⠠⠄⠁⎦
```

## Bandwidth minimization

To minimize the bandwidth of this matrix, we can use the `minimize_bandwidth` function along with any of the available minimization algorithms (a full list can be found in the [home page's **Algorithms** section](index.md#Algorithms)). For example, to use the reverse Cuthill&ndash;McKee algorithm, we can run:

```julia-repl
julia> using MatrixBandwidth

julia> res_minimize = minimize_bandwidth(A, Minimization.ReverseCuthillMcKee())
Results of Bandwidth Minimization Algorithm
 * Algorithm: Reverse Cuthill–McKee
 * Approach: heuristic
 * Minimum Bandwidth: 9
 * Original Bandwidth: 51
 * Matrix Size: 60×60

julia> A[res_minimize.ordering, res_minimize.ordering]
60×60 SparseMatrixCSC{Float64, Int64} with 93 stored entries:
⎡⠠⡢⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎤
⎢⠀⠈⠠⡢⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠈⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠠⠂⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠠⡢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠠⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠢⡀⠈⠀⠂⠀⢠⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠠⠀⠀⠀⠀⠘⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀⣀⠀⠀⠀⠣⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠉⠢⣀⡸⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠪⡢⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢊⠐⢠⡀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠲⢄⡱⢀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠐⠪⢂⡄⠀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠪⡢⎦
```

If no algorithm is explicitly specified, `minimize_bandwidth` defaults to the Gibbs–Poole–Stockmeyer algorithm:

```julia-repl
julia> res_minimize_default = minimize_bandwidth(A)
Results of Bandwidth Minimization Algorithm
 * Algorithm: Gibbs–Poole–Stockmeyer
 * Approach: heuristic
 * Minimum Bandwidth: 5
 * Original Bandwidth: 51
 * Matrix Size: 60×60

julia> A[res_minimize_default.ordering, res_minimize_default.ordering]
60×60 SparseMatrixCSC{Float64, Int64} with 93 stored entries:
⎡⠠⡢⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎤
⎢⠀⠈⠠⡢⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠈⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠠⠂⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠠⡢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠠⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠢⡀⠈⠀⢢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠠⣀⣀⠘⢄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⠤⠃⠠⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠦⢄⠑⠤⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠃⠠⡢⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠪⢂⡄⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⢊⡰⡀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⣄⠙⢤⠀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠓⠠⡢⎦
```

(We default to Gibbs–Poole–Stockmeyer because it is one of the most accurate heuristic algorithms—note how in this case, it produced a lower-bandwidth ordering than reverse Cuthill–McKee. Of course, if true optimality is required, an exact algorithm should be used instead.)

## Bandwidth recognition

Similarly, we can use the `has_bandwidth_k_ordering` function to determine whether a given matrix has bandwidth at most some fixed integer $k$, via any of the available recognition algorithms (again, a full list can be found in the [home page's **Algorithms** section](index.md#Algorithms)). For example, to determine whether our example matrix has bandwidth *at most* 3 via the Del Corso–Manzini algorithm, we can run:

```julia-repl
julia> res_recognize = has_bandwidth_k_ordering(A, 3, Recognition.SaxeGurariSudborough())
Results of Bandwidth Recognition Algorithm
 * Algorithm: Saxe–Gurari–Sudborough
 * Bandwidth Threshold k: 3
 * Has Bandwidth ≤ k Ordering: true
 * Original Bandwidth: 51
 * Matrix Size: 60×60

julia> A[res_recognize.ordering, res_recognize.ordering]
60×60 SparseMatrixCSC{Float64, Int64} with 93 stored entries:
⎡⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎤
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠠⠂⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠘⠀⡠⢀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠐⠊⠀⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⠀⠀⠠⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠢⡄⠉⡢⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠪⢀⠔⢂⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠰⡊⠈⠂⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠀⠊⠀⠠⡀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠢⣀⡸⢀⡀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠰⡀⠈⢠⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠒⣀⡸⢂⡀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠰⡀⠈⠠⡀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠢⠊⡠⎦
```

(In this case, though, it turns out that 3 *is* the true minimum bandwidth of the matrix, as can be verified by running `minimize_bandwidth` with any exact algorithm.)

As with `minimize_bandwidth`, if no algorithm is explicitly specified, `has_bandwidth_k_ordering` defaults to the Caprara–Salazar-González algorithm:

```julia-repl
julia> res_recognize_default = has_bandwidth_k_ordering(A, 6)
Results of Bandwidth Recognition Algorithm
 * Algorithm: Caprara–Salazar-González
 * Bandwidth Threshold k: 6
 * Has Bandwidth ≤ k Ordering: true
 * Original Bandwidth: 51
 * Matrix Size: 60×60

julia> A[res_recognize_default.ordering, res_recognize_default.ordering]
60×60 SparseMatrixCSC{Float64, Int64} with 93 stored entries:
⎡⠀⠀⠐⠡⢄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎤
⎢⠔⡀⠀⠀⠀⠐⢄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠑⢀⠀⠀⠀⠀⠐⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠑⢀⠀⠠⠂⠈⠑⢄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠈⢆⠀⠀⢀⠀⠲⠄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠑⢠⡀⡀⠈⠀⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠁⠄⠀⠀⠀⠐⠐⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢐⠀⠄⠁⡉⠁⠄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠇⠈⡀⠈⠈⠑⠄⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠁⢆⠀⠀⠀⠈⣀⣀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠁⠂⢠⠀⠀⠀⠐⠄⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⢀⠀⠄⠁⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠁⠀⠀⠀⠀⠀⢐⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⢀⠀⠀⠠⠀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠂⠀⠀⎦
```

## Additional utilities

Complementing our various bandwidth minimization and recognition algorithms, MatrixBandwidth.jl exports several additional core functions, including (but not limited to) `bandwidth` and `profile` to compute the original bandwidth and profile of a matrix:

```julia-repl
julia> using MatrixBandwidth

julia> using Random, SparseArrays

julia> Random.seed!(1234);

julia> A = sprand(50, 50, 0.02)
50×50 SparseMatrixCSC{Float64, Int64} with 49 stored entries:
⎡⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠀⠀⠀⠀⎤
⎢⠀⠀⠀⠁⠀⠀⠀⠀⠀⠀⠀⠂⢀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠁⠀⠀⠀⡀⠀⡀⠀⠀⠀⠄⠀⠀⠀⠄⠀⎥
⎢⠀⠐⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠄⠂⠀⠀⠐⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠈⠀⠀⠀⠀⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⢀⠀⠀⠀⠂⠀⠀⠀⠁⎥
⎢⡀⡀⢀⠄⠀⠁⠄⢀⠀⠀⠀⠀⠀⢀⠀⠀⠠⠀⠀⠀⠀⣀⠀⠀⠀⎥
⎢⠀⢀⠀⠀⠀⠊⠀⠄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠄⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠐⠀⠀⠀⠀⠀⠀⠀⠀⠀⠐⠀⠀⠀⠀⠀⠀⎥
⎢⠈⠀⢀⠀⠀⠀⢀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠂⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠂⠀⠀⠀⠀⠀⠀⠀⠀⠨⠐⠀⠀⢀⠀⠀⠀⠀⠀⠀⠀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎦

julia> bandwidth(A) # Bandwidth prior to any reordering of rows and columns
38

julia> profile(A) # Profile prior to any reordering of rows and columns
703
```

(Closely related to bandwidth, the *column profile* of a matrix is the sum of the distances from each diagonal entry to the farthest nonzero entry in that column, whereas the *row profile* is the sum of the distances from each diagonal entry to the farthest nonzero entry in that row. `profile(A)` computes the column profile of `A` by default, but it can also be used to compute the row profile.)
