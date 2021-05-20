function [val] = cost2(u,X0,mu,dt,dt_u,num_steps,a,target)
    u = reshape(u,3,[]);
    
    % Number of runs to complete a burn
    n = dt_u/dt;

    % Initial setup:
    X = zeros(6,num_steps);
    X(:,1) = X0;
    kk = 1;
    for ii = 1:num_steps
        if mod(ii,n) == 0
            u_in = u(:,kk);
            kk = kk+1;
        else
            u_in = [0;0;0];
        end
        
        % Simulate dynamics:
        X(:,ii+1) = rk4(@cweq,dt,X(:,ii),mu,a,u_in);
    end

    % Final states:
    r = X(1:3,end);
    v = X(4:6,end);
    
    % Evaluate the cost:
    u = reshape(u,1,[]);
    pos_error = norm(r - target(1:3));
    y_error = norm(X(2,:));
    z_error = norm(X(3,:));
    vel_error = norm(v - target(4:6));
    val = norm([pos_error,vel_error,norm(u),y_error,z_error]);
end