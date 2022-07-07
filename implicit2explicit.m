function A = implicit2explicit(Afun, m, n)
%IMPLICIT2EXPLICIT takes linear function A(x) and builts corresponding matrix
%   makes an explicit matrix using the linear function
%   in the function handle "Afun", where the domain is R^n
%   and the range is in R^m
% Usage: implicit2explicit(Afun,m,n)
%   If n = [n1,n2], the domain is the space of n1 x n2 matrices
%
% Stephen Becker, stephen.beckr@gmail.com, 2009

A = zeros(m,prod(n));
if numel(n) == 1
    e = zeros(n,1);
else
    if numel(n) ~= 2, error('bad value for size of domain'); end
    e = zeros(n(1),n(2));
end
for j = 1:prod(n)
    e(j) = 1;
    A(:,j) = Afun(e);
    e(j) = 0;
end
%{
Feb 2013
another way:
rr = m*n;
inpt=mat2cell( eye(rr^2), ones(rr^2,1) );
LINMAP  = cellfun( @(in) linmap_vec(in'), inpt,'UniformOutput',false );
LINMAP  = cell2mat(LINMAP');

March 2013
Also try arrayfun
 or even better
 bsxfun
 http://yagtom.googlecode.com/svn/trunk/html/speedup.html#29
 applies to every column... (or any singleton dimension)
 It's supposed to be super fast
%}
