function [dX] = cweq(~,X, mu,a, u)
    % Calculate mean motion:
    n = sqrt(mu/a^3);
    
    % Define linear state space model:
    A = [  0   0   0    1    0    0;
           0   0   0    0    1    0;
           0   0   0    0    0    1;
         3*n^2 0   0    0   2*n   0;
           0   0   0  -2*n   0    0;
           0   0 -n^2   0    0    0];
    
    B = [0 0 0;
         0 0 0;
         0 0 0;
         1 0 0;
         0 1 0;
         0 0 1];
    
    dX = A*X + B*u;
end