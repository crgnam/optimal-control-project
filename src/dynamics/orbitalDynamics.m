function [dX] = orbitalDynamics(~, X, mu, u)
    r = X(1:3);
    v = X(4:6);
    
    a = -(mu*r)/norm(r)^3 + u;
    dX = [v; a];
end