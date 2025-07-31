# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

@setup_workload begin
    A = BitMatrix([
        0 0 1 1 1;
        0 0 0 0 0;
        1 0 0 0 0;
        1 0 0 0 1;
        1 0 0 1 0
    ])

    @compile_workload begin
        bandwidth(A)
        profile(A)
        bandwidth_lower_bound(A)

        for decider_type in ALGORITHMS[:Recognition]
            has_bandwidth_k_ordering(A, 2, decider_type())
        end

        for solvers in values(ALGORITHMS[:Minimization]), solver_type in solvers
            minimize_bandwidth(A, solver_type())
        end
    end
end
