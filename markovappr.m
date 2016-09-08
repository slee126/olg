function [Tran,s] = markovappr(lambda, sigma, m, N)
% the simple case of approximating first-order 
% autoregressive process with Markov chain


% same as previous except accounting for variance of 20 years old are
% different

% y_t = lambda * y_(t-1) + u_t

% u_t is a Gaussian white noise process with standard deviation sigma.

% m determines the width of discretized state space, Tauchen uses m=3
% ymax=m*vary,ymin=-m*vary, ymax and ymin are two boundary points

% N is the number of possible states chosen to approximate
% the y_t process, usually N=9 should be fine

% Tran is the transition matrix of the Markov chain

% s is the discretized state space of y_t

% alambda is the theoretical first order autoregression coefficient 
% for Markov chain

% asigma is the theoretical standard deviation for Markov chain Y_t

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% discretize the state space

stvy = sqrt(sigma^2/(1-lambda^2)); % standard deviation of y_t
ymax = m*stvy;                     % upper boundary of state space
ymin = -ymax;                      % lower boundary of state space

% w = (ymax-ymin)/(N-1);             % length of interval 
% 
% s = ymin:w:ymax;                   % the discretized state space
s = linspace(ymin, ymax, N);
w = s(2) - s(1);
% calculate the transition matrix

for j=1:N;
   
   for k=2:N-1;
      Tran(j,k)= normcdf(s(k)-lambda*s(j)+w/2,0,sigma)...
         - normcdf(s(k)-lambda*s(j)-w/2,0,sigma);
      
   end
   
   Tran(j,1) = normcdf(s(1)-lambda*s(j)+w/2,0,sigma);
   Tran(j,N) = 1 - normcdf(s(N)-lambda*s(j)-w/2,0,sigma);
   
end

if sum(Tran') ~= ones(1,N)
   
   str = find(Tran'-ones(1,N));  % find rows not adding up to one
   disp('error in transition matrix');
   disp(['rows ',num2str(str),' does not sum to one']);
   
end


% calculate the invariant distribution of Markov chain

Trans= Tran';
