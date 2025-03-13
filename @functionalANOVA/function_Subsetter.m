% Gets lowerbound and upperbound indices to subset dataMatrixCell
function function_Subsetter(self)
% Cleaner Implementation

lb = min(self.boundsArray);
ub = max(self.boundsArray);


% Check that lb is not greater than the maximum value in self.d_grid
if lb > max(self.d_grid)
    error('The lower bound (lb=%g) is greater than the maximum value in self.d_grid (%g).', lb, max(self.d_grid));
end

% Check that ub is not less than the minimum value in self.d_grid
if ub < min(self.d_grid)
    error('The upper bound (ub=%g) is less than the minimum value in self.d_grid (%g).', ub, min(self.d_grid));
end

% Get the indices corresponding to the values within the bounds
subset_idx = find(self.d_grid >= lb & self.d_grid <= ub);

% Create the subset of d_grid
subset_d_grid = self.d_grid(subset_idx);

% Find the local indices for the smallest and largest values in the subset
[~, local_min_idx] = min(subset_d_grid);
[~, local_max_idx] = max(subset_d_grid);

% Map the local indices back to the indices in the original self.d_grid array
global_min_idx = subset_idx(local_min_idx);
global_max_idx = subset_idx(local_max_idx);

self.lb_index = global_min_idx;
self.ub_index = global_max_idx;

end