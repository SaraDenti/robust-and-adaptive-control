clc
clear all
close all

g=9.81;

%dinamica disturbi d2 e d3
omega = 1;
alpha = 1;

% matrici

A = [0 1 0 0; 0 -1 0 0; 0 0 0 1;0 1 g 0];
B = [0; 1; 0; -1];
P = [0 0 0; 1 0 0; 0 0 0; -1 0 0];

C=[1 0 0 0; 0 0 1 0]; % matrice dell'uscita y
Cs = [1 0 0 0]; %matrice errore e

Q=[0 -1 0];

S = [0 0 0; 0 0 1; 0 -omega^2 0]

%% Calcolo matrice K stabilizzabile

% Vettore poli desiderati
v = [-1,-sqrt(9.81),-1.2,-0.8];

% Poli desiderati per A+BK (distinti per evitare errore place)
K = place(A,-B,v);

%% Soluzione delle FBI equations

% Questa funzione ritorna la soluzione X dell'equazione AXB + CXD = E
function [X] = generalized_sylvester(A,B,C,D,E)
n = length(A);
m = length(B);

% All inputs must be numeric
% A and C must be square matrices of the same size 
% B and D must be square matrices of the same size
% E must be a matrix with number of rows being the same as A and C,
% and number of columns being the same as B and D
required_classes = {'numeric'};
required_attributes_A_C = {'nonempty','finite','nonnan','square','size',[n,n]};
required_attributes_B_D = {'nonempty','finite','nonnan','square','size',[m,m]};
required_attributes_E = {'nonempty','finite','nonnan','size',[n,m]};

% Validate the properties described above
validateattributes(A,required_classes,required_attributes_A_C,'','A',1);
validateattributes(B,required_classes,required_attributes_B_D,'','B',2);
validateattributes(C,required_classes,required_attributes_A_C,'','C',3);
validateattributes(D,required_classes,required_attributes_B_D,'','D',4);
validateattributes(E,required_classes,required_attributes_E,'','E',5);

% NOTE ON OBTAINING THE SOLUTION: Here I use a vectorization approach to easily determine the solution. However, this approach requires inverting a nm times nm matrix, which could be problematic for large problems. To avoid this, a generalized eigenvalue/eigenvector approach should be implemented in the future.
% To do this, add a check for regularity and disjoint generalized spectra of the matrix penxils (A + x C), (D - x B), then apply generalized eigenvectors to calculate X element by element.

% Calculate the vectorization of the solution X by using the kronecker product (for any matrices A, B, C, vec(A*B*C) = kron(transpose(C),A)*vec(B))
kronecker_vectorization_factor = kron(transpose(B),A) + kron(transpose(D),C);

% Before trying to calculate the inverse, check that it is invertible. If it is not invertible, then there does not exist a unique solution to the generalized Sylvester equation (there is either no solution, or there is an infinite set solutions)
validateattributes(det(kronecker_vectorization_factor),required_classes,{'nonempty','finite','nonnan','nonzero'},'','the determinant of the Kronecker vectorization factor (required to be invertible for the existence of a unique solution to the generalized Sylvester equation)');
vecX = kronecker_vectorization_factor\E(:);

% Check that vec(X) is acceptable before reshaping (numeric, nonempty, finite, nonnan, and has size [n*m,1])
required_attributes_vecX = {'nonempty','finite','nonnan','size',[n*m,1]};
validateattributes(vecX,required_classes,required_attributes_vecX,'','the solution to the generalized Sylvester equation, X,');

%Reshape vec(X) into the n times m matrix X
X = reshape(vecX,n,m);
end

%AXB + CXD = E
A_ext = [-A, -B; -Cs, zeros(size(Cs,1),size(B,2))];
B_ext = eye(size(S));
C_ext = [eye(size(A,1)) zeros(size(A,1), size(B,2)); zeros(size(Cs,1),size(A,1)+size(B,2))];
D_ext = S;
E_ext = [P; Q];

X = generalized_sylvester(A_ext, B_ext, C_ext, D_ext, E_ext);

% Estraggo Pi e Gamma
Pi = X(1:size(A,1), :);
Gamma = X(size(A,1)+1:end, :);

L = Gamma - K * Pi;




