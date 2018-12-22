% *************************************************************************
% Econ 899 Homework 7
%  Shiyan Wei
%**************************************************************************
clear
clc
close

%% Question 2 
% ------------- parameter initiated ------------------
global alpha beta

beta  = 0.99;
alpha = 0.36;


nk = 20;
nK = 10; 

klb = 0.001;
kub  = 0.5;
k = linspace(klb,kub,nk)';

Klb = 0.15;
Kub  = 0.25;
K = linspace(Klb,Kub,nK)';

kk = linspace(klb,kub,1000)';



% create a cartesian product of k and K
[k_tem, K_tem] = meshgrid(k,K);
k_space = [k_tem(:),K_tem(:)];

%--------------- estimate value function ----------------

% initial guess A =1 B =1
A0 =1;
B0 =1;
x = log(k_space(:,1));
v0 = A0 + B0 * x;
v1 = zeros(size(k_space,1),1);
dec_k = zeros(size(k_space,1),1);

Iter = 0;
tol = 10;

while tol>0.4
for i = 1:size(k_space,1)    
 % bundling parameter
 inputPar.k = k_space(i,1); 
 inputPar.K = k_space(i,2);
 inputPar.A = A0;
 inputPar.B = B0;
 
 % calculate the flot utility
 w = utility(inputPar,kk);
 [v1(i), dec_k(i)]=max(w,[],2);
 % using the fminunc to find the minimum function 
end

X = [ones(size(x)),x];
b = mvregress(X,v1);
 % updateing A and B
 A0 = b(1);
 B0 = b(2);
 
 tol = max(abs(v1-v0));
 Iter = Iter +1;
 
 v0 = A0 + B0 * x;
 
%  fprintf('The current iteration is %d, the diff is %.3f \n',Iter,tol);
end


%----------------- Bilinear Interpolation -------------------------
k_I = kk; % k interpolation point
K_I = linspace(Klb,Kub,1000)'; % K interpolation point

[k_I_mesh,K_I_mesh] = meshgrid(k_I,K_I);
k_choice = kk(dec_k);

% reshape the target function
v1_resh = reshape(v1,nK, nk);
dec_k_resh = reshape(k_choice,nK, nk);

% Interpolation
v1_I = interp2(k_tem,K_tem,v1_resh,k_I_mesh,K_I_mesh,'linear'); % interpolation over value function
dec_k_I = interp2(k_tem,K_tem,dec_k_resh,k_I_mesh,K_I_mesh,'linear'); % interpolation over value function

% plot the orginial and interpolate value function figure
 figure (1)
plot3(k_tem,K_tem,v1_resh,'o');
hold on
mesh(k_I_mesh,K_I_mesh,v1_I);
hold off
xlabel('k')
ylabel('K')
zlabel('v(k:K)')
title('Bilinear Interpolation result for value function')

% plot the orginial and interpolate value function figure
 figure (2)
plot3(k_tem,K_tem,dec_k_resh,'o');
hold on
mesh(k_I_mesh,K_I_mesh,dec_k_I);
hold off
xlabel('k')
ylabel('K')
zlabel('kk(k:K)')
title('Bilinear Interpolation result for decision function')

% ----------------- Cubic Splines Interpolation --------------------
% Interpolation
v1_Ic = interp2(k_tem,K_tem,v1_resh,k_I_mesh,K_I_mesh,'cubic'); % interpolation over value function
dec_k_Ic = interp2(k_tem,K_tem,dec_k_resh,k_I_mesh,K_I_mesh,'cubic'); % interpolation over value function

% plot the orginial and interpolate value function figure
 figure (3)
plot3(k_tem,K_tem,v1_resh,'o');
hold on
mesh(k_I_mesh,K_I_mesh,v1_Ic);
hold off
xlabel('Individual Asset Holding')
ylabel('Aggregate Asset Holding')
zlabel('v(k:K)')
title('Cubic Spline Interpolation result for value function')

% plot the orginial and interpolate value function figure
 figure (4)
plot3(k_tem,K_tem,dec_k_resh,'o');
hold on
mesh(k_I_mesh,K_I_mesh,dec_k_Ic);
hold off
xlabel('Individual Asset Holding')
ylabel('Aggregate Asset Holding')
zlabel('kk(k:K)')
title('Cubic Spline Interpolation result for decision function')



%% Question 3
% -------------- Calculate the computational solution when k = K --------


% v0 = 1 + 1 * x ; % start with initial guess A = 1 B = 1
% initial guess A =1 B =1
A0 =1;
B0 =1;
x = log(k_space(:,1));
v0_2 = A0 + B0 * x;
v1_2 = zeros(size(k_space,1),1);
dec_k = zeros(size(k_space,1),1);

Iter = 0;
tol = 10;

while tol>0.35
for i = 1:size(k_space,1)    
 % bundling parameter
 inputPar.k = k_space(i,1); 
 inputPar.K = k_space(i,2);
 inputPar.A = A0;
 inputPar.B = B0;
 
 % calculate the flot utility
 w = utility2(inputPar,kk);
 [v1_2(i), dec_k(i)]=max(w,[],2);
 % using the fminunc to find the minimum function 
end

X = [ones(size(x)),x];
b = mvregress(X,v1_2);
 % updateing A and B
 A0 = b(1);
 B0 = b(2);
 
 tol = max(abs(v1_2-v0_2));
 Iter = Iter +1;
 
 v0_2 = A0 + B0 * x;
 
%  fprintf('The current iteration is %d, the diff is %.3f \n',Iter,tol);
end

% reshape the target function
v0_resh = reshape(v0_2,nK, nk);
dec_k_resh = reshape(dec_k,nK, nk);

% Interpolation
v0_I = interp2(k_tem,K_tem,v0_resh,k_I_mesh,K_I_mesh,'linear'); % interpolation over value function
dec_k_I = interp2(k_tem,K_tem,dec_k_resh,k_I_mesh,K_I_mesh,'linear'); % interpolation over value function

% Interpolation
v0_Ic = interp2(k_tem,K_tem,v0_resh,k_I_mesh,K_I_mesh,'cubic'); % interpolation over value function
dec_k_Ic = interp2(k_tem,K_tem,dec_k_resh,k_I_mesh,K_I_mesh,'cubic'); % interpolation over value function
%--------- Steady state calculations ------------------
% calculate the steady state
Kss=(alpha*beta)^(1/(1-alpha));
%Closed form solution
A=(beta*alpha*log(beta*alpha)/(1-alpha*beta)+log((1-alpha*beta)))/(1-beta);
B=alpha/(1-alpha*beta);

% calculate value function given Kss
valueclosedform=A+B*log(k_I);

% find the interpolation point with K_I = K_ss
inc = (Kub-Klb) / 1000 ;
Kss_loc = ceil((Kss - Klb)/inc)+1;

% valuefunction with K_I = K_ss
v_ss_I = v0_I(Kss_loc,:)';
diff =  max(abs(v_ss_I - valueclosedform));
v_ss_I = v_ss_I - diff;


% valuefunction with K_I = K_ss
v_ss_Ic = v0_Ic(Kss_loc,:)';
diffc =  max(abs(v_ss_Ic - valueclosedform));
v_ss_Ic = v_ss_Ic - diffc;

% plot the value function in steady state
figure (5)
plot(k_I,valueclosedform,k_I,v_ss_I,'--',k_I,v_ss_Ic,'g');
xlabel('Asset Holding')
ylabel('Value function')
title('Cubic Spline Interpolation result for decision function')
legend('Closed form solution','Bilinear interpolation','Cubic Spline interpolation'...
    ,'Location','southeast')