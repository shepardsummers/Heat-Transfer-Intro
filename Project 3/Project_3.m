clear, clc;
%% Defination of Parameters
% Domain Related
N_x = 200;
N_y = N_x;
dx = 1;
dy = 1;
% LBM Related
Ksi = [0 1 0 -1  0 1  -1  -1  1;...
     0 0 1  0 -1 1   1  -1 -1]; % Lattice Velocities for D2Q9
w = [0 1/6 1/6 1/6 1/6 1/12 1/12 1/12 1/12]; % Weights for D2Q9
c_s = 1/sqrt(3); % Speed of Sound for D2Q9
Tau = 0.8; % Relaxation Time
%% Curved Stuff
R = N_x/5;
D = R*2;
N_x_circ_center = round(N_x/2);
N_y_circ_center = round(N_y/2);
center_x = dx * (N_x_circ_center - 1);
center_y = dy * (N_y_circ_center - 1);


for i = 1:N_x
    x(i)=dx*(i-1);
end
for j = 1:N_y
    y(j)=dy*(j-1);
end
% Domain ID
% Domain=0; Solid
% Domain=1; Fluid
for j = 1:N_y
    for i = 1:N_x
        if test_circle(x(i), y(j), R, center_x, center_y)
            Domain_ID(j,i) = 0;
        else 
            Domain_ID(j,i) = 1;
        end
    end
end
%figure
%contourf(Domain_ID, 30)
%axis equal tight

% Zone ID
% Zone ID=0 --> dead zone
% Zone_ID=1 --> b nodes
% Zone_ID=2 --> f nodes
% Zone_ID=3 --> other nodes in fluid domain
for j = 1:N_y
    for i = 1:N_x
        if Domain_ID(j,i) == 0 % solid domain
            if Domain_ID(j,i+1) == 1 || Domain_ID(j,i-1) == 1 || Domain_ID(j+1,i) == 1 || Domain_ID(j-1,i) == 1 ...
                    || Domain_ID(j+1,i+1) == 1 || Domain_ID(j-1,i+1) == 1 || Domain_ID(j+1,i-1) == 1 || Domain_ID(j-1,i-1) == 1
                Zone_ID(j,i) = 1;
            else
                Zone_ID(j,i) = 0;
            end
        else % fluid domain
            if i == 1 || i == N_x || j == 1 || j == N_y
                Zone_ID(j,i) = 3;
            else
                if Domain_ID(j,i+1) == 0 || Domain_ID(j,i-1) == 0 || Domain_ID(j+1,i) == 0 || Domain_ID(j-1,i) == 0 ...
                    || Domain_ID(j+1,i+1) == 0 || Domain_ID(j-1,i+1) == 0 || Domain_ID(j+1,i-1) == 0 || Domain_ID(j-1,i-1) == 0
                    Zone_ID(j,i) = 2;
                else
                    Zone_ID(j,i) = 3;
                end
            end
            
        end
    end
end
figure
contourf(Zone_ID, 30)
axis equal tight
%% Initialization
Rho_ref=2;
R = 8.314;
Rho=ones(1,N_y,N_x)*Rho_ref;
T_H = 1/3;
T_L = 1/2;
T_1 = 1/2;
T_2 = 1/3;


T=zeros(1,N_y,N_x)*T_L;

g=zeros(9,N_y,N_x); % PDF for all 9 directions and at all locations
for j=1:N_y
    for i=1:N_x
        g(:,j,i) = eqm_t_d2q9(squeeze(Rho(1,j,i)),squeeze(T(1,j,i)));
    end
end
g_new=g;
g_eq=g;

%% Timer
Timer = 15000;
tic;
x_circ = center_x;
y_circ = center_y;
%% Solving
for t=1:Timer
% Streaming/Boundary Conditions
    for j=1:N_y
        for i=1:N_x
            if Zone_ID(j,i) == 0 % dead zone
                % Do nothin
            elseif Zone_ID(j,i) == 1 % b node
                % Do nothin
            elseif Zone_ID(j,i) == 2 % f node
                
                g_new(1,j,i) = g(1,j,i);
                 
                % GZS curved boundary scheme goes here
                if Zone_ID(j,i+1) == 1 % Direction 2
                    g_new(4,j,i) = GZS_t_scheme(x(i+1),y(j),x(i),y(j),R,x_circ,y_circ,T(:,j,i),T(:,j,i-1),w(2),Rho(1,j,i),g(2,j,i),g_eq(2,j,i),g(2,j,i-1),g_eq(2,j,i-1),Tau);
                else
                    g_new(4,j,i) = g(4,j,i+1);
                end
                if Zone_ID(j-1,i) == 1 % Direction 3
                    g_new(5,j,i) = GZS_t_scheme(x(i),y(j-1),x(i),y(j),R,x_circ,y_circ,T(:,j,i),T(:,j+1,i),w(3),Rho(1,j,i),g(3,j,i),g_eq(3,j,i),g(3,j+1,i),g_eq(3,j+1,i),Tau);
                else
                    g_new(5,j,i) = g(5,j-1,i);
                end
                if Zone_ID(j,i-1) == 1 % Direction 4
                    g_new(2,j,i) = GZS_t_scheme(x(i-1),y(j),x(i),y(j),R,x_circ,y_circ,T(:,j,i),T(:,j,i+1),w(4),Rho(1,j,i),g(4,j,i),g_eq(4,j,i),g(4,j,i+1),g_eq(4,j,i+1),Tau);
                else
                    g_new(2,j,i) = g(2,j,i-1);
                end
                if Zone_ID(j+1,i) == 1 % Direction 5
                    g_new(3,j,i) = GZS_t_scheme(x(i),y(j+1),x(i),y(j),R,x_circ,y_circ,T(:,j,i),T(:,j-1,i),w(5),Rho(1,j,i),g(5,j,i),g_eq(5,j,i),g(5,j-1,i),g_eq(5,j-1,i),Tau);
                else
                    g_new(3,j,i) = g(3,j+1,i);
                end
                if Zone_ID(j-1,i+1) == 1 % Direction 6
                    g_new(8,j,i) = GZS_t_scheme(x(i+1),y(j-1),x(i),y(j),R,x_circ,y_circ,T(:,j,i),T(:,j+1,i-1),w(6),Rho(1,j,i),g(6,j,i),g_eq(6,j,i),g(6,j+1,i-1),g_eq(6,j+1,i-1),Tau);
                else
                    g_new(8,j,i) = g(8,j-1,i+1);
                end
                if Zone_ID(j-1,i-1) == 1 % Direction 7
                    g_new(9,j,i) = GZS_t_scheme(x(i-1),y(j-1),x(i),y(j),R,x_circ,y_circ,T(:,j,i),T(:,j+1,i+1),w(7),Rho(1,j,i),g(7,j,i),g_eq(7,j,i),g(7,j+1,i+1),g_eq(7,j+1,i+1),Tau);
                else
                    g_new(9,j,i) = g(9,j-1,i-1);
                end         
                if Zone_ID(j+1,i-1) == 1 % Direction 8
                    g_new(6,j,i) = GZS_t_scheme(x(i-1),y(j+1),x(i),y(j),R,x_circ,y_circ,T(:,j,i),T(:,j-1,i+1),w(8),Rho(1,j,i),g(8,j,i),g_eq(8,j,i),g(8,j-1,i+1),g_eq(8,j-1,i+1),Tau);
                else
                    g_new(6,j,i) = g(6,j+1,i-1);
                end
                if Zone_ID(j+1,i+1) == 1 % Direction 9
                    g_new(7,j,i) = GZS_t_scheme(x(i+1),y(j+1),x(i),y(j),R,x_circ,y_circ,T(:,j,i),T(:,j-1,i-1),w(9),Rho(1,j,i),g(9,j,i),g_eq(9,j,i),g(9,j-1,i-1),g_eq(9,j-1,i-1),Tau);
                else
                    g_new(7,j,i) = g(7,j+1,i+1);
                end
            else
                if j == 1 % This is the top boundary nodes
                    if i ==1 % Top-Left corner node
                        g_new(1,j,i) = g(1,j,i);
                        g_new(3,j,i) = g(3,j+1,i);
                        g_new(4,j,i) = g(4,j,i+1);
                        g_new(7,j,i) = g(7,j+1,i+1);
        
                        % Unknown
                        T_w = T_1;
                        g_new(2,j,i) = g_new(4,j,i);
                        g_new(5,j,i) = g_new(3,j,i);
                        g_new(9,j,i) = g_new(7,j,i);
                        g_new(6,j,i) = (Rho(1,j,i)*R*T_w - g_new(1,j,i) - 2*(g_new(3,j,i) + g_new(4,j,i) + g_new(7,j,i)))/2;
                        g_new(8,j,i) = g_new(6,j,i);
                    elseif i == N_x % Top-right corner node
                        g_new(1,j,i) = g(1,j,i);
                        g_new(2,j,i) = g(2,j,i-1);
                        g_new(3,j,i) = g(3,j+1,i);
                        g_new(6,j,i) = g(6,j+1,i-1);
        
                        % Unknown
                        T_w = T_1;
                        g_new(4,j,i) = g_new(2,j,i);
                        g_new(5,j,i) = g_new(3,j,i);
                        g_new(8,j,i) = g_new(6,j,i);
                        g_new(7,j,i) = (Rho(1,j,i)*R*T_w - g_new(1,j,i) - 2*(g_new(2,j,i) + g_new(3,j,i) + g_new(6,j,i)))/2;
                        g_new(9,j,i) = g_new(7,j,i);
                    else % All other nodes on the top boundary
                        g_new(1,j,i) = g(1,j,i);
                        g_new(2,j,i) = g(2,j,i-1);
                        g_new(3,j,i) = g(3,j+1,i);
                        g_new(4,j,i) = g(4,j,i+1);
                        g_new(6,j,i) = g(6,j+1,i-1);
                        g_new(7,j,i) = g(7,j+1,i+1);
        
                        % Unknown
                        g_new(5,j,i) = g_new(3,j,i);
                        T_w = T_1;
                        g_new(8,j,i) = (Rho(1,j,i)*R*T_w - g_new(1,j,i) - g_new(2,j,i) - 2*g_new(3,j,i) - g_new(4,j,i) - g_new(6,j,i) - g_new(7,j,i))/2;
                        g_new(9,j,i) = g_new(8,j,i);
                    end
                elseif j == N_y % This is the bottom boundary nodes
                    if i ==1 % Bottom-Left corner node
                        g_new(1,j,i) = g(1,j,i);
                        g_new(4,j,i) = g(4,j,i+1);
                        g_new(5,j,i) = g(5,j-1,i);
                        g_new(8,j,i) = g(8,j-1,i+1);
                        
                        % Unknown
                        T_w = T_2;
                        g_new(2,j,i) = g_new(4,j,i);
                        g_new(3,j,i) = g_new(5,j,i);
                        g_new(6,j,i) = g_new(8,j,i);
                        g_new(7,j,i) = (Rho(1,j,i)*R*T_w - g_new(1,j,i) - 2*(g_new(4,j,i) + g_new(5,j,i) + g_new(8,j,i)))/2;
                        g_new(9,j,i) = g_new(7,j,i);
                    elseif i == N_x % Bottom-right corner node
                        g_new(1,j,i) = g(1,j,i);
                        g_new(2,j,i) = g(2,j,i-1);
                        g_new(5,j,i) = g(5,j-1,i);
                        g_new(9,j,i) = g(9,j-1,i-1);
        
                        % Unknown
                        slope = (T(1,j-1,i) - T(1,j,i))/((j-1)-j);
                        T_w = -(1/8)*slope + (4/5);
                        g_new(3,j,i) = g_new(5,j,i);
                        g_new(4,j,i) = g_new(2,j,i);
                        g_new(7,j,i) = g_new(9,j,i);
                        g_new(6,j,i) = (Rho(1,j,i)*R*T_w - g_new(1,j,i) - 2*(g_new(2,j,i) + g_new(5,j,i) + g_new(9,j,i)))/2;
                        g_new(8,j,i) = g_new(6,j,i);
                    else % All other nodes on the bottom boundary
                        g_new(1,j,i) = g(1,j,i);
                        g_new(2,j,i) = g(2,j,i-1);
                        g_new(4,j,i) = g(4,j,i+1);
                        g_new(5,j,i) = g(5,j-1,i);
                        g_new(8,j,i) = g(8,j-1,i+1);
                        g_new(9,j,i) = g(9,j-1,i-1);
        
                        % Unknown
                        g_new(3,j,i) = g_new(5,j,i);
                        slope = (T(1,j-1,i) - T(1,j,i))/((j-1)-j);
                        T_w = -(1/8)*slope + (4/5);
                        g_new(6,j,i) = (Rho(1,j,i)*R*T_w - g_new(1,j,i) - g_new(2,j,i) - 2*g_new(5,j,i) - g_new(4,j,i) - g_new(8,j,i) - g_new(9,j,i))/2;
                        g_new(7,j,i) = g_new(6,j,i);
                    end
                elseif i == 1 % This is the left boundary nodes
                    g_new(1,j,i) = g(1,j,i);
                    g_new(3,j,i) = g(3,j+1,i);
                    g_new(4,j,i) = g(4,j,i+1);
                    g_new(5,j,i) = g(5,j-1,i);
                    g_new(7,j,i) = g(7,j+1,i+1);
                    g_new(8,j,i) = g(8,j-1,i+1);
        
                    % Unknown
                    T_w = T_2;
                    g_new(2,j,i) = g_new(4,j,i);
                    g_new(6,j,i) = (Rho(1,j,i)*R*T_w - g_new(1,j,i) - g_new(3,j,i) - 2*g_new(4,j,i) - g_new(5,j,i) - g_new(7,j,i) - g_new(8,j,i))/2;
                    g_new(9,j,i) = g_new(6,j,i);
                elseif i == N_x % This is the right boundary nodes
                    g_new(1,j,i) = g(1,j,i);
                    g_new(2,j,i) = g(2,j,i-1);
                    g_new(3,j,i) = g(3,j+1,i);
                    g_new(5,j,i) = g(5,j-1,i);
                    g_new(6,j,i) = g(6,j+1,i-1);
                    g_new(9,j,i) = g(9,j-1,i-1);
                
                    % Unknown
                    T_w = T(1,j,i-1);
                    g_new(4,j,i) = g_new(2,j,i);
                    g_new(7,j,i) = (Rho(1,j,i)*R*T_w - g_new(1,j,i) - g_new(3,j,i) - 2*g_new(2,j,i) - g_new(5,j,i) - g_new(6,j,i) - g_new(9,j,i))/2;
                    g_new(8,j,i) = g_new(7,j,i);
                else  % All interior nodes
                    g_new(1,j,i) = g(1,j,i);
                    g_new(2,j,i) = g(2,j,i-1);
                    g_new(3,j,i) = g(3,j+1,i);
                    g_new(4,j,i) = g(4,j,i+1);
                    g_new(5,j,i) = g(5,j-1,i);
                    g_new(6,j,i) = g(6,j+1,i-1);
                    g_new(7,j,i) = g(7,j+1,i+1);
                    g_new(8,j,i) = g(8,j-1,i+1);
                    g_new(9,j,i) = g(9,j-1,i-1);
                end
            end
    end
end

% Collision
% T calculation
T = sum(g_new, 1)/R./Rho;

% g_eq calculation
for j=1:N_y
    for i=1:N_x
        g_eq(:,j,i) = eqm_t_d2q9(squeeze(Rho(1,j,i)),squeeze(T(1,j,i)));
    end
end
g_eq_d2q9 = w'.*Rho*R.*T;

% BGK Collision & Update
g = g_new - (g_new - g_eq)/Tau;

fprintf("Itt: %i\n", t);

end
Runtime=toc;
%% Post-Processing/Visualizaation

figure
x_h = (1:N_x)/N_x;
T_h = squeeze((T(1,round(N_y/2),:) - T_2) / (1 - T_2));
plot(x_h, T_h, 'red');
hold on
load 'Project3_Benchmark Data.mat'
plot(x_benchmark, T_benchmark_hori, 'blue');
hold off

figure
y_v = (1:N_y)/N_y;
T_v = flip(squeeze(T(1,:,round(N_x/2)) - T_2)/(1 - T_2));
plot(y_v, T_v, 'red');
hold on
plot(y_benchmark, T_benchmark_vert, 'blue');
hold off

figure
T_nd = (T - T_2)/(1 - T_2);
contourf(flipud(squeeze(T_nd)),30)
axis equal tight