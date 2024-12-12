function [ desired_state ] = traj_generator(t, state, waypoints)

persistent coef_x coef_y coef_z waypoints0 traj_time d0
if nargin > 2
    d = waypoints(:,2:end) - waypoints(:,1:end-1);
    d0 = 2 * sqrt(d(1,:).^2 + d(2,:).^2 + d(3,:).^2);
    traj_time = [0, cumsum(d0)];
    waypoints0 = waypoints;

    coef_x = getCoef(waypoints0(1,1:end)');
    coef_y = getCoef(waypoints0(2,1:end)');
    coef_z = getCoef(waypoints0(3,1:end)');
    
else
    
    if(t > traj_time(end))
        t = traj_time(end) - 0.0001;
    end
    
    t_index = find(traj_time >= t,1)-1;
    
    if (t_index == 0)
        t_index = 1;
    end
    if(t == 0)
        desired_state.pos = waypoints0(:,1);
        desired_state.vel = 0*waypoints0(:,1);
        desired_state.acc = 0*waypoints0(:,1);
    else
        scale = (t-traj_time(t_index))/d0(t_index);
        
        index = (t_index-1)*8+1:t_index*8;
        
        t0 = polyT(8,0,scale)';
        desired_state.pos = [coef_x(index)'*t0; coef_y(index)'*t0; coef_z(index)'*t0];
        
        t1 = polyT(8,1,scale)';
        desired_state.vel = [coef_x(index)'*t1; coef_y(index)'*t1; coef_z(index)'*t1].*(1/d0(t_index));
        
        t2 = polyT(8,2,scale)';
        desired_state.acc = [coef_x(index)'*t2; coef_y(index)'*t2; coef_z(index)'*t2].*(1/d0(t_index)^2);
    end
    
    desired_state.yaw = 0;
    desired_state.yawdot = 0;
end
end

function [T] = polyT(n, k, t)
T = zeros(n,1);
D = zeros(n,1);

for i=1:n
    D(i) = i-1;
    T(i) = 1;
end
for j=1:k
    for i=1:n
        T(i) = T(i) * D(i);
        if D(i) > 0
            D(i) = D(i) - 1;
        end
    end
end
for i=1:n
    T(i) = T(i) * t^D(i);
end

T = T';

end

function [coef, A, b] = getCoef(waypoints)

n = size(waypoints,1)-1;
b = zeros(1,8*n);
for i=1:n
    b(1,i) = waypoints(i);
    b(1,i+n) = waypoints(i+1);
end
A=zeros(8*n,8*n);
for i=1:n
    A(i,((i-1)*8)+1:i*8) = polyT(8,0,0);
end
for i=1:n
    A(i+n,((i-1)*8)+1:i*8) = polyT(8,0,1);
end
for k=1:3
    A(2*n+k,1:8) = polyT(8,k,0);
end
for k=1:3
    A(2*n+3+k,(end-7):end) = polyT(8,k,1);
end
for i=2:n
    for k=1:6
        A(2*n+6+(i-2)*6+k, (i-2)*8+1:((i-2)*8+n*n)) = [polyT(8,k,1) -polyT(8,k,0)];
    end
end

coef = A\b';

end

