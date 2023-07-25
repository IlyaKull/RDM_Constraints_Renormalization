function A=tensor(varargin)
%computes the tensor product of the input
A=1;
for k=1:nargin
    A=kron(A,varargin{k});
end

end