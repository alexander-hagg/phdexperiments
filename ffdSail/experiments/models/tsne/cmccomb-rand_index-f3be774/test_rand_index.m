% Two arbitrary partitions
p1 = [1 1 1 1 0 0];
p2 = [1 1 1 0 0 0];

% Compute the unadjusted rand index
ri = rand_index(p1, p2)

% Compute the adjusted rand index
ri = rand_index(p1, p2, 'adjusted')
