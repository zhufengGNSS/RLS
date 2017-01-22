%% Example Recursive Least Squares %%
clc
close all
clear all
%%
Gs = tf([0.5],[1 1 1]);                     % Actual system
Ts = 0.2;                                   %Sampling Time
Gz = c2d(Gs,Ts);
%%
u = @(t) 2*exp(-0.1*t)*sin(2*pi/2.5*t);     % Input
t = 0:0.2:50;                               % Time of Operation
% for i=1:length(t)
% U(i) = u(t(i)); % Input
% end
U = wgn(1,length(t),2);
%% Calculations of Orders
num = get(Gz,'num');
den = get(Gz,'den');
delay = get(Gz,'iodelay');
Az = cell2mat(den);                  
Az = Az./Az(1);                             % Contains [1 a1 a2 ... ana]
if delay >0
Bz = [zeros(1,delay),cell2mat(num)]./Az(1); % Contains [zeros(d) b0 b1 b2 ... bnb]
else
Bz = cell2mat(num)./Az(1); % Contains [zeros(d) b0 b1 b2 ... bnb]  
end

% To get delay "d ==  num of first zeros in Bz"
ind = find(Bz == 0);
test = isempty(ind);
if test == 1
        d = 0;
elseif ind(1) == 1
        d = 1;
elseif ind(1)~= 1 
        d = 0;
end
if length(ind)>1
for i=1:length(ind)-1
    if ind(i+1)-ind(i) == 1
       d = i+1;
    else
        break
    end
end
end

B = Bz(d+1:end);                          % Contains [b0 b1 b2 ... bnb] 
A = Az(2:end);                            % Contains [a1 a2 ... ana]

nb=length(B)-1;
na=length(A);
nu = na+nb+1;


%% Actual Output of the System
Y = lsim(Gz,U,t);

%% Estimation (simulated as online, where i represents time)

Theta = zeros(nu,1);                      % Initial Parameters
P = 10^6 * eye(nu,nu);                    % Initial Covariance Matrix
Phi = zeros(1,nu);                        % Initial phi
for i = 1 : length(U)
[a(:,i),b(:,i),P,Theta,Phi,K(:,i)] = RecursiveLeastSquares(U(1:i),Y(1:i),d,nb,na,P,Theta,Phi,i);
end

%% Plots

% Convergence of parameters
Sa = size(a);
Sb = size(b);
colors = ['b','r','g','k'];

figure
hold on
for m = 1:Sa(1)
    plot(t,a(m,:),'color',colors(m),'LineWidth',1.5)
    plot(t,A(m)*ones(length(t),1),'--','color','k')
    grid on
    title('Convergance of paramters a_i')
    xlabel('time (sec.)')
    ylabel('Paramter')
end
legend('a_1','','a_2','','a_3','a_4','a_5')
figure
hold on
for m = 1:Sb(1)
    plot(t,b(m,:),'color',colors(m),'LineWidth',1.5)
    plot(t,B(m)*ones(length(t),1),'--','color','k')
    grid on
    title('Convergance of paramters b_i')
    xlabel('time (sec.)')
    ylabel('Paramter')
end
legend('b_0','','b_1','','b_2','b_3','b_4','b_5')

% Check the estimated T.F.
y = lsim(tf(b(:,end)',[1 a(:,end)'],Ts,'iodelay',d,'variable','z^-1'),U,t);

% Compare the actual and estimated systems outputs
figure
grid on 
hold on
plot(t,Y(1:end),'--','LineWidth',1.5)
plot(t,y(1:end),'o','LineWidth',1.5,'color','r')
xlabel('time (sec.)')
ylabel('Amplitude')
title('2^n^d Order Estimated System')
legend('System Output','2^n^d Order Estimated Output') 
