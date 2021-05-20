function [rHill,vHill] = eci2hill(rTgt, vTgt, rChase, vChase)
%% Function Description:
% Converts position and velocity from ECI J2000 to Hill's frame (using
% target + chase)
%   - Chris Gnam, 2019
%
% Reference: Fundamentals of Astrodynamics and Applications - Vallado
%
% ===============
% INPUT VARIABLES
% ===============
%   rTgt   - ECI Position vector of first (km)                                         
%   vTgt   - ECI Velocity vector of first (km/s)                                         
%   rChase - ECI Position vector of second (km)
%   vChase - ECI Velocity vector of second (km/s)
%
% ================
% OUTPUT VARIABLES
% ================
%   rHill  - Hill's relative position vector (km)
%   vHill  - Hill's relative velocity vector (km/s)
%
% (Radial, In-Track, Cross-Track)

%% ECI 2 Hill frame Conversion:
rTgtMag   = sqrt(sum(rTgt.^2,1));
rChaseMag = sqrt(sum(rChase.^2,1));
vTgtMag   = sqrt(sum(vTgt.^2,1));

RSW = eci2rsw(rTgt, vTgt);

%Use RSW rotation matrix to convert rChase and vChase to RSW
r_Chase_RSW = RSW*rChase;
v_Chase_RSW = RSW*vChase;

%Find Rotation angles to go from target to interceptor
phi_chase    = asin(r_Chase_RSW(3)/rTgtMag); % OLD --> asin(r_Chase_RSW(3)/rChaseMag)
lambda_chase = atan2(r_Chase_RSW(2),r_Chase_RSW(1));
CPC = cos(phi_chase);     
SPC = sin(phi_chase);
SLC = sin(lambda_chase);  
CLC = cos(lambda_chase);

%Find Position component rotations
rHill = cat(1, rChaseMag - rTgtMag, ...
               lambda_chase*rTgtMag, ...
               phi_chase*rTgtMag);
           
%Find the rotation matrix RSW->SEZ of chaser
RSW_SEZ(1,1) = SPC*CLC;  RSW_SEZ(1,2) = SPC*SLC;   RSW_SEZ(1,3) = -CPC;
RSW_SEZ(2,1) = -SLC;     RSW_SEZ(2,2) = CLC;       RSW_SEZ(2,3) = 0;
RSW_SEZ(3,1) = CPC*CLC;  RSW_SEZ(3,2) = CPC.*SLC;  RSW_SEZ(3,3) = SPC;

%Find the velocity component of positions using the angular rates in SEZ frame
v_Chase_SEZ = RSW_SEZ*v_Chase_RSW;
vHill = [v_Chase_SEZ(3);
         rTgtMag*(v_Chase_SEZ(2)/(rChaseMag*CPC)-vTgtMag/rTgtMag);
         -rTgtMag*v_Chase_SEZ(1)/rChaseMag];
end