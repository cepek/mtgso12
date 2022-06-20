# random ill-conditioned matrix
#
#  Based on https://www.cis.upenn.edu/~cis610/Gram-Schmidt-Bjorck.pdf
#
N = randi([9 20]);
M = N + randi([1 N]);

T = rand(M, N);

[U S V] = svd(T);

q = 1;
for i=1:N
  S(i,i) = q;
  q = q/10;
end

A = U*S*V;
b = A*ones(N,1);

B = [A -b; eye(N,N) zeros(N,1)];


size_B = size(B);
size_A = size(A);
dim_b = length(b);

save randmat.dat size_B B size_A A dim_b b
