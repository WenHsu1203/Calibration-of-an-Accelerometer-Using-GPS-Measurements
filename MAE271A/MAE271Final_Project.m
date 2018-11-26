
clear
close all

%% Parameters
N = 3000;                  %Nave (should be large enough) 
tsA = 1/200;              %Accelerometer sampling time
ta = 0:tsA:30;
Na = length(ta);
tsG = 1/5;                %GPS sampling time
tg = 0:tsG:30;
Ng = length(tg);
w_freq = 0.1;             %Frequency of true acceleration
amp = 10;                 %Amplititude of true acceleration
el = zeros(3*N , Ng);     %Actual error for realization l
rl = zeros(2*N , Ng);     %Residual for realization l
est_x = zeros(3*N , Ng);  %estimate value of the state (deltaX+Xc)


for i = 1:N
    %% Derive Model
    w = 0.02*randn(1,Na); %Gaussian noise with 0 mean amd variance 0.0004
    b = 0.1*randn();      %Bias ~ N(0,0.01)
    p0 = 10*randn();      %x0 ~ N(0,10^2)
    v0 = randn() + 100;   %v0 ~ N(100,1^2)
    eta = [randn(1,Ng) ; 0.04*randn(1,Ng)]; %[etaX ; etaV] ~ N([0;0],[1 0;0 0.04^2]) Note: Ng not Na
    %True Model
    a = amp*sin(w_freq*ta);                                         %True acceleration
    v = v0 + amp/w_freq - amp/w_freq*cos(w_freq*ta);                %True velocity
    p = p0 + (v0 + amp/w_freq)*ta - amp/w_freq^2*sin(w_freq*ta);    %True position
    %Accelerometer Model
    ac = a + b + w;
    vc = zeros(1,Na);
    vc(1) = 100;           %Vc(0) = mean of v0
    pc = zeros(1,Na);
    pc(1) = 0;             %Pc(0) = mean of p0
    %Euler integration formula
    ve = zeros(1,Na);
    ve(1) = v0;
    pe = zeros(1,Na);
    pe(1) = p0;
    
    for j = 2:Na
        %Derive accelerometer model
        vc(j) = vc(j-1) + ac(j-1) * tsA;
        pc(j) = pc(j-1) + vc(j-1) * tsA + ac(j-1) * tsA^2/2;
        %Derive Euler integration model
        ve(j) = ve(j-1) + a(j-1) * tsA;
        pe(j) = pe(j-1) + ve(j-1) * tsA + a(j-1) * tsA^2/2;
    end
    %To get dynamical model independent of acceleration,
    %subtract Euler Model by accelerometer model
    delta_Xe = [(pe - pc) ; (ve - vc) ; b * ones(1,Na)];
    
    %to get the measure equations,
    %subtract measure model by accelerometer model
    delta_Z = [(p(1:40:end) - pc(1:40:end)) ; (v(1:40:end) - vc(1:40:end))] + eta; %Note difference of sampling rate 
   
    %% Kalman Filter
    %Parameters
    V = [1 0 ; 0 0.04^2];
    phi = [1 tsA -tsA^2/2 ; 0 1 -tsA ; 0 0 1];
    gamma = [-tsA^2/2 ; -tsA ; 0];
    H = [1 0 0 ; 0 1 0];
    W = 0.04;
    delta_Xm = zeros(3,Na);   % a priori and delta_Xm0 = [0 0 0]'
    delta_Xhat = zeros(3,Na); % a posteriori
    M0 = [10^2 0 0 ; 0 1^2 0 ; 0 0 0.1^2];
    K0 = M0*H'/(H*M0*H' + V);
    P0 = inv(inv(M0) + H'*inv(V)*H);
    delta_Xhat(:,1) = delta_Xm(:,1) + K0*(delta_Z(:,1) - H*delta_Xm(:,1)); %Kalman algorithm No.3
    %Now we have P0, Xhat0 => we can calculate X1m,M1...
    %start iterative process
    M = M0;
    P = P0;
    P_track = P0;
    for k = 2:Na
        delta_Xm(:,k) = phi * delta_Xhat(:,k-1);    %Kalman algorithm No.1
        if (mod(k,40) == 1)                         %When we have GPS measurement 
            u = ceil(k/40);                         %To determine which Z to use
            M = phi * P * phi' + gamma * W * gamma';%Kalman algorithm No.2
            K = M*H'/(H*M*H' + V);
            P = inv(inv(M) + H'*inv(V)*H);
            delta_Xhat(:,k) = delta_Xm(:,k) + K*(delta_Z(:,u)-H*delta_Xm(:,k));
        else                                        %When no measurement
            M = phi * P * phi' + gamma * W * gamma';
            P = M;
            delta_Xhat(:,k) = delta_Xm(:,k);
        end
        P_track = cat(2,P_track,P);
    end
    %Monte Carlo simulation
    delta_x = [(p-pc) ; (v-vc) ; b*ones(1,Na)];                  %True state - accelerometer state
    priori_E = delta_x - delta_Xm;
    posteriori_E = delta_x(:,1:40:end) - delta_Xhat(:,1:40:end);
    el((3*i-2):3*i , :) = posteriori_E;
    rl(2*i-1:2*i , :) = H*priori_E(:,1:40:end)-eta;
    estimate_state = delta_Xhat + [pc ; vc ; zeros(1,Na)];       %estimate p, v
    st_x((3*i-2):3*i , :) = estimate_state(:,1:40:end);
    
end

e_ave = [sum(el(1:3:end,:));sum(el(2:3:end,:));sum(el(3:3:end,:))] / N;
pl = zeros(3,3*Ng);
p_total = zeros(3,3*Ng);
orth = zeros(3,3*Ng);
orth_total = zeros(3,3*Ng);
riXrm = zeros(2,2*Ng);
riXrm_total = zeros(2,2*Ng);

for m = 1:N
    for n = 1:Ng
        pl(:,3*n-2:3*n) = (el(3*m-2:3*m,n)-e_ave(:,n)) * (el(3*m-2:3*m,n)-e_ave(:,n))';     %check P_ave-P ~0
        orth(:,3*n-2:3*n) = (el(3*m-2:3*m,n)-e_ave(:,n)) * st_x(3*m-2:3*m,n)';              %check orthogonality of the error in estimates with the estimate
        riXrm(:, 2*n-1:2*n) = rl(2*m-1:2*m,ceil((n+100)/2))*rl(2*m-1:2*m,ceil((n+50)/2))';  %check residuals~0
    end
    p_total = p_total + pl;
    orth_total = orth_total + orth;
    riXrm_total = riXrm_total + riXrm;
end
p_ave = p_total / (N-1);
orth_ave = orth_total / N;
riXrm_ave = riXrm_total / N;

%% drawing
run('MAE271Final_Project_Plot');



