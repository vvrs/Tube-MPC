


A = [1 1; 0 1];
B = [0.5; 1];

Q = eye(2);
R = 0.01;

% state bounds {x | [0 1]<=2}

X_bounds = [-8 2; -8 -3; 4 -3; 4 2 ];

U_bounds = [1; -1];

W_bounds = [0.1 0.1; 0.1 -0.1; -0.1 0.1; -0.1 0.1];

X = Polyhedron(X_bounds);
U = Polyhedron(U_bounds);
W = Polyhedron(W_bounds);

x_0 = [-5; -2];
N = 9;

[x_tmpc,u_tmpc,x_seq,u_seq,Z,Xf] = TubeMPC(A,B,Q,R,X,U,W,N,x_0);
xz = Xf+Z;
figure(1)
plot3(x_seq(1,:),x_seq(2,:),0:N,'-*');
hold on 
plot3(x_tmpc(1,:),x_tmpc(2,:),0:N,'-o');
hold on 

