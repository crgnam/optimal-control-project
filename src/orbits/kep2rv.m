function [r,v] = kep2rv(a,e,i,omega,Omega,theta,mu, deg)
    if nargin == 7
        deg = false;
    end
    
    if deg
        theta = deg2rad(theta);
    end
    
    % Define Orbit Parameters
    h = sqrt(mu.*a.*(1-e.*e));
    r = ((h.*h)./mu).*(1./(1+(e.*cos(theta))));
    
    % Calculate R,V in PQW coordinates
    R_pqw(:,1) = r.*cos(theta);
    R_pqw(:,2) = r.*sin(theta);
    R_pqw(:,3) = 0;
    V_pqw(:,1) = (mu./h).*-sin(theta);
    V_pqw(:,2) = (mu./h).* (e+cos(theta));
    V_pqw(:,3) = 0;

    % Calcualte the rotation matrix portion:
    A = ea2rotmat('313',Omega,i,omega, deg);
    r = A'*R_pqw';
    v = A'*V_pqw';
end



