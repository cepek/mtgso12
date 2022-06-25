# random ill-conditioned matrix
#
#  Based on https://www.cis.upenn.edu/~cis610/Gram-Schmidt-Bjorck.pdf
#
for K = 70:50:1100

N = randi([K cast(K*1.32,"int32")]);
M = N + randi([1 cast(N/2,"int32")]);

T = rand(M, N);

[U S V] = svd(T);

for i=1:N
  S(i,i) = 1;
end
S(1,1) = 1e-8;

A = U*S*V;
b = A*ones(N,1);

B = [A -b; eye(N,N) zeros(N,1)];


size_B = size(B);
size_A = size(A);
dim_b = length(b);

output_file = strcat(num2str(N,"%04d"), num2str(M,"%04d"), "-randm-a.dat");
save(output_file, "size_B", "B", "size_A", "A", "dim_b", "b")
display(output_file);

endfor
