N = 512
M = 16

A = rand(N,N);

V = rand(M,N);
U = [eye(M,M);zeros(N-M,M)];

Ainv = A^(-1);
Anew = A + U*V;
Ainv_new = Anew^(-1);

Ainv_WF = Ainv - Ainv*U*(eye(M,M) + V*Ainv*U)^(-1)*V*Ainv;

max(max(abs(Ainv_WF - Ainv_new)))