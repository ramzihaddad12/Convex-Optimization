%% Inserting linearized non-convex constraint and non-convex constraint
clc; clear all;% close all;

% OneFinalPos   = [1000*rand(),1000*rand(),1000*rand()];
% TwoFinalPos   = [1000*rand(),1000*rand(),1000*rand()];
% ThreeFinalPos = [1000*rand(),1000*rand(),1000*rand()];
% FourFinalPos  = [1000*rand(),1000*rand(),1000*rand()];

% OneFinalPos   = [0,0,-1000];
% TwoFinalPos   = [0,0,1000];
% ThreeFinalPos = [0,-1000,0];
% FourFinalPos  = [0,1000,0];
% OneFinalPos   = [1000,1000,1000];
% TwoFinalPos   = [1000,1000,1000];
% ThreeFinalPos = [1000,1000,1000];
% FourFinalPos  = [1000,1000,1000];
OneFinalPos   = [5000,4000,3000];
TwoFinalPos   = [5100,4000,3000];
ThreeFinalPos = [4900,4000,3000];
FourFinalPos  = [5000,4000,2900];


CenterObs=[3500,3000,2000];
Rball=1000;


N = 4;  % Drones

% Number of discrete time steps
K = 50; % Steps

% Discrete Time Step Resolution
h = 1;  % Seconds

% Minimum Seperating Radius
R = 10;

% Gravity Vector
g = [0,0,-9.81];
%Tolerance
eps=1;

[Csts,obj0,ops,P,V,A,J] = convprob([],OneFinalPos,TwoFinalPos,ThreeFinalPos,FourFinalPos);
%% BIGGIE LOOP
over=false;
done=false;
f0_prev=value(obj0);
fff=0;

while(~over || ~done || ~finish)
    fff=fff+1;
pq = value(P);

yalmip('clear')%%%%%%%%%%%%
% A = sdpvar(3,K,N);
% V = sdpvar(3,K,N);
% P = sdpvar(3,K,N);
% J = sdpvar(3,K-1,N);
% Graphing(pq);
Csts_lin=[];

for i=1:N-1
    for j=i+1:N
        for k=1:K
        alpha=pq(:,k,i)-pq(:,k,j);
        beta=norm(alpha);
        eta=alpha/beta;
        Csts_lin=Csts_lin+[(beta+(eta')*((P(:,k,i)-P(:,k,j))-alpha))>=R];
        

        end
        
    end
end

for i=1:N
    for k=1:K
    alpha=pq(:,k,i)- CenterObs';
    beta=norm(alpha);
    eta=alpha/beta;
    Csts_lin = Csts_lin+[(beta+(eta')*((P(:,k,i)-CenterObs')-alpha))>=(Rball+R)];

    end


end
%Solving with linear constraint
% [Csts_new,obj,ops,P,V,A,J] = convprob(Csts_lin,OneFinalPos,TwoFinalPos,ThreeFinalPos,FourFinalPos);
% pq = value(P)


[Csts_new,obj,ops,P,V,A,J] = convprob(Csts_lin,OneFinalPos,TwoFinalPos,ThreeFinalPos,FourFinalPos);
pq = value(P)
% Non convex constraint and convergence check

done=true;
for i=1:N-1
    for j=i+1:N
        for k=1:K
            if(norm(pq(:,k,i)-pq(:,k,j))<=R)
                done=false;
                break
            end
        end
        if ~done
            break;
        end
    end
if ~done
    break;
end
end

finish=true;
for i=1:N
    for k=1:K
        if(norm(pq(:,k,i)-(CenterObs)')<=Rball+R)
            finish=false;
            break
        end
    end
    if ~finish
        break;
    end
end

over=true
if abs(value(obj)-f0_prev)>eps
    over=false;
end
f0_prev=value(obj);

end
pq = value(P)
Graphing1(pq);