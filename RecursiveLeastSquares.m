function [a,b,P,Theta,phi,K] = RecursiveLeastSquares(U,Y,d,nb,na,P,Theta,phi,n)
% This function Identify the system parameters for a known system input and
% output.
% The system Transfer Function is as in the following form:
%
%         z^-d * (bo + b1*z^-1 + b2*z^-2 + ... + b_nb*z^-nb)
% G(z) =  ---------------------------------------------------
%                1 + a1*z^-1 + ... + a_na*z^-na
%
% OUTPUTS:
% a and b are vectors of the system estimated parameters.
% INPUTS:
% u : is the system input raw vector, y is the system output raw vector.
% d : is the delay.
% nb: is the number of zeros of the equired system model.
% na: is the number of poles of the equired system model.
% k : is the instant of time at which parameters are to be calculated.
% P,Theta,phi: are the last calculated ones.

nu = na+nb+1;                           % Number of unknowns


for j = 1:nu
        if j <= na % terms of y
            if (n-j)<=0
                phi(n,j) = 0;
            else
                phi(n,j) = -Y(n-j);
            end
        else       % terms of u
            if (n-d-(j-(na+1)))<=0
                phi(n,j) = 0;
            else
                phi(n,j) = U(n-d-(j-(na+1)));
            end
        end
end

                    % Estimation   
           
           K = P*phi(n,:)'*inv(1+phi(n,:)*P*phi(n,:)');
           Theta = Theta+K*(Y(n)-phi(n,:)*Theta);
           P = P-P*phi(n,:)'*inv(1+phi(n,:)*P*phi(n,:)')*phi(n,:)*P;
                    % Estimated System Parameters
           a = Theta(1:na);
           b = Theta(na+1:end);
end