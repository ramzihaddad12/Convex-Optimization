function [Csts,obj,ops,P,V,A,J] = convprob(myCsts,OneFinalPos,TwoFinalPos,ThreeFinalPos,FourFinalPos)
%% No Obstacle Avoidance constraint

% Number of Drones
N = 4;  % Drones

% Number of discrete time steps
K = 50; % Steps

% Discrete Time Step Resolution
h = 1;  % Seconds

% Minimum Seperating Radius
R = 10;

% Gravity Vector
g = [0,0,-9.81];

% Initial Positions
OneInitPos   = [0,0,0];
TwoInitPos   = [0,0,2000];
ThreeInitPos = [0,0,-1000];
FourInitPos  = [0,2000,1000];
% 
% OneInitPos   = [0,0,1000];
% TwoInitPos   = [0,0,-1000];
% ThreeInitPos = [0,1000,0];
% FourInitPos  = [0,-1000,0];

% OneInitPos   = [100,100,0];
% TwoInitPos   = [-100,100,0];
% ThreeInitPos = [100,-100,0];
% FourInitPos  = [-100,-100,0];

% Initial Velocities
% OneInitVel   = [0,0,0];
% TwoInitVel   = [0,0,0];
% ThreeInitVel = [0,0,0];
% FourInitVel  = [0,0,0];
OneInitVel   = [10,10,10];
TwoInitVel   = [10,20,30];
ThreeInitVel = [30,30,30];
FourInitVel  = [30,40,50];

% Final Positions
% OneFinalPos   = [1000*rand(),1000*rand(),1000*rand()];
% TwoFinalPos   = [1000*rand(),1000*rand(),1000*rand()];
% ThreeFinalPos = [1000*rand(),1000*rand(),1000*rand()];
% FourFinalPos  = [1000*rand(),1000*rand(),1000*rand()];

% OneFinalPos   = [4000,4000,3000];
% TwoFinalPos   = [4000,4000,3000];
% ThreeFinalPos = [4000,4000,3000];
% FourFinalPos  = [4000,4000,3000];

% Final Velocities
OneFinalVel   = [0,0,0];
TwoFinalVel   = [0,0,0];
ThreeFinalVel = [0,0,0];
FourFinalVel  = [0,0,0];

% Min Positions
P_min = [-1000,-1000,-2000];

% Max Position Values
P_max = [6000,5000,4000];

% Min Velocity
V_min = [-200,-200,-200];

% % Max Velocity
V_max = [200,200,200];

% Min Acceleration
A_min = [-20,-20,-20];

% Max Acceleration
A_max = [20,20,20];

% Jerk Min
J_min = [-5,-5,-5];

% Jerk Max
J_max = [5,5,5];

%% Problem Formulation

A = sdpvar(3,K,N);
V = sdpvar(3,K,N);
P = sdpvar(3,K,N);
J = sdpvar(3,K-1,N);

Csts=myCsts;

% Setting the Velocity and Position Vectors' Equations
for i = 1:N
    for k = 1:K-1
        Csts = Csts + [V(:,k+1,i) == V(:,k,i) + h*A(:,k,i)];
        Csts = Csts + [P(:,k+1,i) == P(:,k,i) + h*V(:,k,i) + (0.5*h*h)*A(:,k,i)];
    end
end

% Setting up the Jerk Vector
for i = 1:N
    Csts = Csts + [J(:,:,i) == (A(:,2:K,i) - A(:,1:K-1,i))/h];
end

for l = 1:3 %3D
    for i = 1:N
        for k = 1:K
            % Velocity Bounds
            Csts = Csts + [V(l,k,i) <= V_max(l)];
            Csts = Csts + [V(l,k,i) >= V_min(l)];
            % Acceleration Bounds
            Csts = Csts + [A(l,k,i) <= A_max(l)];
            Csts = Csts + [A(l,k,i) >= A_min(l)];
            % Position Constraints
            Csts = Csts + [P(l,k,i) <= P_max(l)];
            Csts = Csts + [P(l,k,i) >= P_min(l)];
            if(k<K)
                % Jerk Constraints
                Csts = Csts + [J(l,k,i) <= J_max(l)];
                Csts = Csts + [J(l,k,i) >= J_min(l)];
            end
        end
    end
end

% Initialization Position Constraints
Csts = Csts + [P(:,1,1) == OneInitPos'];
Csts = Csts + [P(:,1,2) == TwoInitPos'];
Csts = Csts + [P(:,1,3) == ThreeInitPos'];
Csts = Csts + [P(:,1,4) == FourInitPos'];

% Final Position Constraints
Csts = Csts + [P(:,K,1) == OneFinalPos'];
Csts = Csts + [P(:,K,2) == TwoFinalPos'];
Csts = Csts + [P(:,K,3) == ThreeFinalPos'];
Csts = Csts + [P(:,K,4) == FourFinalPos'];

% Initial Velocity Constraints
Csts = Csts + [V(:,1,1) == OneInitVel'];
Csts = Csts + [V(:,1,2) == TwoInitVel'];
Csts = Csts + [V(:,1,3) == ThreeInitVel'];
Csts = Csts + [V(:,1,4) == FourInitVel'];

% Final Velocity Constraints
Csts = Csts + [V(:,K,1) == OneFinalVel'];
Csts = Csts + [V(:,K,2) == TwoFinalVel'];
Csts = Csts + [V(:,K,3) == ThreeFinalVel'];
Csts = Csts + [V(:,K,4) == FourFinalVel'];
%%%%%%%%%%%%%%%

obj = 0;
for i = 1:N
    for k = 1:K
        obj = obj + (norm(A(:,k,i)+g'))^2;
    end
end

ops = sdpsettings('solver','sdpt3','verbose',1);
optimize(Csts,obj,ops);
pq = value(P);

end

