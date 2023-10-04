function [xMem,stats] = AndersonAcceleration(xMem,typ,reg)
% see GLOBALLY CONVERGENT TYPE-I ANDERSON ACCELERATION FOR NONSMOOTH FIXED-POINT ITERATION
% JUNZI ZHANG, BRENDAN O'DONOGHUE,AND STEPHEN BOYD
% https://web.stanford.edu/~boyd/papers/pdf/scs_2.0_v_global.pdf

if nargin<2
    typ = 1;
end
if nargin<3
    reg = 0;
end

G = (-1)*diff(xMem,1,2); % g_i=x_i -  f(x_i)
Y = diff(G,1,2);
dim=size(Y,2);

% under the condition that xMem = [x1 x2=f(x1) x3=f(x2) x4=f(x3)], 
% S_k from the SCS acceleration documentation satisfies s_i=-g_i,
% in the notation we use here: S=-G(:,1:end-1)

switch typ 
    case 2
        %% type II
        % x_k+1 = x_k - B_k*g_k; 
        % B_k = I + (S_k - Y_k)*(Y_k^T * Y_k)^inv *Y_k^T
        % x_k+1 = x_k - g_k - (S_k - Y_k)*(Y_k^T * Y_k)^inv *Y_k^T * g
        % x_k+1 = x_k - (x_k-f(x_k)) - (S_k - Y_k)*(Y_k^T * Y_k)^inv *Y_k^T * g_k 
        %       = f(x_k) - ... 
        xMem(:,1) = xMem(:,end) - (-G(:,1:end-1) - Y) * ((Y'*Y + reg*eye(dim))\(Y' * G(:,end)));
        if nargout>1
            stats.condNum= cond((Y'*Y + reg*eye(dim)));
            stats.gammaNorm= norm(((Y'*Y + reg*eye(dim))\(Y' * G(:,end))));
        end
    case 1
        %% type I
        %  x_k+1 = x_k - H_k^inv * g_k; 
        % H_k^inv  = I + (S_k - Y_k)*(S_k^T * Y_k)^inv *S_k^T
        %  x_k+1 = x_k - (I + (S_k - Y_k)*(S_k^T * Y_k)^inv *S_k^T) * g_k = 
        %        = x_k - g_k - (S_k - Y_k)*(S_k^T * Y_k)^inv *S_k^T *g_k = 
        %        =  f(x_k) - (...)*g_k
        xMem(:,1) = xMem(:,end) - (-G(:,1:end-1) - Y) * ((-G(:,1:end-1)' * Y + reg*eye(dim))\(-G(:,1:end-1)' *G(:,end)));
         if nargout>1
            stats.condNum= cond(-G(:,1:end-1)' * Y + reg*eye(dim));
            stats.gammaNorm= norm(((-G(:,1:end-1)' * Y + reg*eye(dim))\(-G(:,1:end-1)' *G(:,end))));
        end
    otherwise 
        error('Anderson acceleration type should be 1 or 2')
end