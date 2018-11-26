M = [2 0;2 2;2 4; 2 6];
N = [0 1 2 3; 1 1 1 1];
Ma = [0 2; 2 4; 4 6;6 8];
Nb = [1 0 -1 -2;0 1 2 3];
M1 = inv(ctranspose(M)*M)*ctranspose(M);
N1 = ctranspose(N)*inv(N*ctranspose(N));
M2 = inv(ctranspose(Ma)*Ma)*ctranspose(Ma);
N2 = ctranspose(Nb)*inv(Nb*ctranspose(Nb));
H = [2 4 6 8;4 6 8 10; 6 8 10 12; 8 10 12 14];
M2*H*N2