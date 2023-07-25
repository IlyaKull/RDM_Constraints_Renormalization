function [K] = dontUseKron(X,D,d,order)
% produces  kron(eye(d),X) [[[order='IX']]
% and kron(X,eye(d)) [[[order='XI']]
% wihout using kron (kron for sdpvars is very slow for large variables).

% 
assert(D==size(X,1));
assert(D==size(X,2));
switch order
    case 'IX'
        
        % block copies of X instead of kron(I,X) 
        X=repmat({X},d,1);
        K=blkdiag(X{:});
        
    case 'XI'

        % intelace columns instead of kron(X,I)
%         W=whos('X');
%         switch W.class
%             case 'double'
%                 K = zeros(D*d);
%             case 'sdpvar'
%                 K = sdpvar(D*d);
%         end
        z=zeros(D,1);
        K=[];
        for j = 1:D
            col = X(:,j);
            for dd=1:d-1 , col = [col,z]; end
            col=col.';   
            col=col(:);
            cols=col;
            for dd=1:d-1
                cols= [cols,circshift(col,dd)]; 
            end
            K=[K,cols];
        end

    otherwise 
        error('specify order as "IX" or "XI" ' )

end
