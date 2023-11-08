% Gets lowerbound and upperbound indices to subset dataMatrixCell

%{
© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are.
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare.
derivative works, distribute copies to the public, perform publicly and display publicly, and to permit.
others to do so.

Author: Adam Watts (acwatts@lanl.gov)
%}
function function_Subsetter(self)
assert(numel(self.boundsArray) == 2, 'Requires 2 Bounds')
assert(numel(unique(self.boundsArray)) == 2, 'Lower and Upper Bounds Must be Different from each Other')

lb = min(self.boundsArray);
ub = max(self.boundsArray);


if lb < min(self.d_grid)
    lb_index = 1;
elseif lb > max(self.d_grid)
    error('Lower Bound is outside the function domain. The Lower Bound exceeds that largest value in the domain')
else
    [~, idx] = min(abs(self.d_grid - lb));
    min_data = self.d_grid(idx);
    min_data = min_data(min_data >= lb); % Must be within the bounds
    min_data = min(min_data);

    lb_index = find(min_data==self.d_grid);
end

if ub > max(self.d_grid)
    ub_index = numel(self.d_grid);
elseif ub < min(self.d_grid)
    error('Upper Bound is Outside the function domain. The upper Bound is smaller than the smallest value in the domain ')
else
    [~, idx] = min(abs(self.d_grid - ub));
    max_data = self.d_grid(idx);
    max_data = max_data(max_data <= ub); % Must be within the bounds
    max_index = max(max_data);
    ub_index = find(max_index==self.d_grid );

end


self.lb_index = lb_index;
self.ub_index = ub_index;
end