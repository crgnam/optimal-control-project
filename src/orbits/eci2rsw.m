function [T, rRSW, vRSW] = eci2rsw(rECI, vECI)
%% Function Description:
% Convert a position in ECI to RSW.
%   - Chris Gnam, 2019
%
% Reference: Fundamentals of Astrodynamics and Applications - Vallado
%
% ===============
% INPUT VARIABLES
% ===============
%   rECI   - ECI Position vector of reference frame (km)                                           
%   vECI   - ECI Velocity vector of reference frame (km/s)                                          
%
% ================
% OUTPUT VARIABLES
% ================
%   T     - Transformation matrix between ECI and this RSW
%   rRSW  - Hill's relative position vector (km)
%   vRSW  - Hill's relative velocity vector (km/s)

%% ECI 2 RSW Conversion:
rvec = rECI/norm(rECI);
wvec = cross(rECI, vECI)/norm(cross(rECI, vECI));
svec = cross(wvec, rvec)/norm(cross(wvec, rvec));

T = [rvec'; svec'; wvec'];
 
rRSW = T*rECI;
vRSW = T*vECI;

end