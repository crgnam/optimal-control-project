function [rotmat] = ea2rotmat(sequence, rot1,rot2,rot3, deg)
    if nargin == 5
        assert(islogical(deg),'Input option for degrees must be a logical type')
        if deg
            rot1 = deg2rad(rot1);
            rot2 = deg2rad(rot2);
            rot3 = deg2rad(rot3);
        end
    end
    
    % Parse input:
    if nargin == 3
        sequence = '321';
    end
    
    % Generate basic rotaiton matrices:
    T{1} = @(rot1) [1       0            0;
                     0   cos(rot1)   sin(rot1);
                     0  -sin(rot1)   cos(rot1)];
    T{2} = @(rot2) [cos(rot2)  0  -sin(rot2);
                         0     1       0;
                    sin(rot2)  0   cos(rot2)];
    T{3} = @(rot3) [ cos(rot3)   sin(rot3)  0;
                    -sin(rot3)   cos(rot3)  0;
                         0           0      1];

    % Calculate overall rotation matrix:
    rotmat = T{str2double(sequence(3))}(rot3)*...
             T{str2double(sequence(2))}(rot2)*...
             T{str2double(sequence(1))}(rot1);
end