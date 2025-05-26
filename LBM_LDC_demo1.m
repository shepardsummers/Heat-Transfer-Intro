clear, clc;
%% Defination of Parameters
% Domain Related
N_x=129;
N_y=31;
dx=1;
dy=1;
% LBM Related
Ksi=[0 1 0 -1  0 1  -1  -1  1;...
     0 0 1  0 -1 1   1  -1 -1]; % Lattice Velocities for D2Q9
w=[0 1/6 1/6 1/6 1/6 1/12 1/12 1/12 1/12]; % Weights for D2Q9
c_s=1/sqrt(3); % Speed of Sound for D2Q9
Tau=0.8; % Relaxation Time
%% Initialization
Rho_ref=2;
R = 8.314;
Rho=ones(1,N_y,N_x)*Rho_ref;
T_H = 10;
T_L = 2;

T=zeros(1,N_y,N_x)*T_L;

g=zeros(9,N_y,N_x); % PDF for all 9 directions and at all locations
for j=1:N_y
    for i=1:N_x
        g(:,j,i) = eqm_t_d2q9(squeeze(Rho(1,j,i)),squeeze(T(1,j,i)));
    end
end
g_new=g;
g_eq=g;

Timer = 500;
tic;
%% Solving
for t=1:Timer
% Streaming/Boundary Conditions
for j=1:N_y
    for i=1:N_x
        if j == 1 % This is the top boundary nodes
            if i ==1 % Top-Left corner node
                g_new(1,j,i) = g(1,j,i);
                g_new(3,j,i) = g(3,j+1,i);
                g_new(4,j,i) = g(4,j,i+1);
                g_new(7,j,i) = g(7,j+1,i+1);

                % Unknown
                T_w = T_H;
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
                T_w = T_L;
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
                T_w = T(1,j+1,i);
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
                T_w = T_H;
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
                T_w = T_L;
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
                T_w = T(1,j-1,i);
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
            T_w = T_H;
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
            T_w = T_L;
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

% Collision
% T calculation
T = sum(g_new, 1)/R./Rho;

% % g_eq calculation
% for j=1:N_y
%     for i=1:N_x
%         g_eq(:,j,i) = eqm_t_d2q9(squeeze(Rho(1,j,i)),squeeze(T(1,j,i)));
%     end
% end
g_eq_d2q9 = w'.*Rho*R.*T;

% BGK Collision & Update
g = g_new - (g_new - g_eq)/Tau;

fprintf("Itt: %i\n", t);

end
Runtime=toc;
%% Post-Processing/Visualizaation
figure
contourf(flipud(squeeze(T)),30)
axis equal tight

figure
plot(1:N_x, squeeze(T(1,round(N_y/2),:)), 'red');
hold on
plot(1:N_x, (T_L - T_H)/(N_x - 1)*(1:N_x)+10, "blue");
