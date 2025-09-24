using MatrixBandwidth
using NPZ
using Random
using SparseArrays

DST_DIR = joinpath(@__DIR__, "npz")

Random.seed!(8675309);

(n, k) = (60, 6)
A = sprand(n, n, 0.01);
A = A + A'

res_rec = has_bandwidth_k_ordering(A, k, Recognition.DelCorsoManzini())
res_min = minimize_bandwidth(A)

ordering_rec = res_rec.ordering
ordering_min = res_min.ordering

A_rec = A[ordering_rec, ordering_rec]
A_min = A[ordering_min, ordering_min]

npzwrite(joinpath(DST_DIR, "A.npz"), A)
npzwrite(joinpath(DST_DIR, "A_min.npz"), A_min)
npzwrite(joinpath(DST_DIR, "A_rec.npz"), A_rec)
