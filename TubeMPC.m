function [x_tmpc,u_tmpc,x_seq,u_seq,Z,Xf] = TubeMPC(A,B,Q,R,X,U,W,N,x_init)

    global Wc;
    Wc = W;
    [nx,nu] = size(B);
    [K_tmp, P] = dlqr(A, B, Q, R);
    K = -K_tmp;
    Ak = (A+B*K);
    n_opt = nx*(N+1)+nu*N; % state s = [x(0),x(1),......,x(N+1),u(1),u(2),....,u(N)]

    % compute Z
    Z = (W+Ak*W+Ak^2*W+Ak^3*W);

    X_robust = X - Z;
    U_robust = U - K*Z;

    % compute Final set Xf
    Xf = compute_terminalConstraintSet(X,U,Z,K,Ak);
    X_init = x_init + Z;

    solved = 0;
    H = compute_cost(Q,R,P,N);

    [Aineq,bineq,nc] = compute_inequalities(X_robust,U_robust,N);

    n.nx = nx;
    n.nu = nu;
    n.n_opt = n_opt;
    n.nc = nc;

    [Aineq,bineq] = addIneqConstr(Aineq,bineq,Xf,N+1,n);
    [Aineq,bineq] = addIneqConstr(Aineq,bineq,X_init,1,n);

    [Aeq,beq] = compute_equalities(A,B,x_init,n_opt,N);

    Aeq = Aeq(nx+1:size(Aeq)*[1; 0], 1:size(Aeq)*[0; 1]);
    beq = zeros(size(Aeq)*[1;0], 1);
    options = optimoptions('quadprog', 'Display', 'none');
    itr = 0;
    while(solved ~=1 )
        [optim, ~, exitflag] = quadprog(H, [], Aineq, bineq, Aeq, beq, [], [], [], options);
        
        solved = (exitflag==1);
        itr = itr + 1;
        if (itr>10)
            error('Not feasible');
        end
    end
    x_seq = reshape(optim(1:nx*(N+1)), nx, N+1);
    u_seq = reshape(optim(nx*(N+1)+1:n_opt), nu, N);

    [x_tmpc,u_tmpc] = controller(A,B,x_init,x_seq,u_seq,K,N);

    % add initial constraint to inequality constraints
    % solve optimal control (nominal)
    % use nominal trajectory and ancillary controller to get robust trajectory
    %   - state update
end

function [Xf] = compute_terminalConstraintSet(X,U,Z,K,Ak)
    % MPIset is computed only once in the constructor;
    
    Xf = Xpi(0,X,U,K,Ak);
    i= 0;
    while(1) % 
        i = i + 1;
        Xf_tmp = and(Xf, Xpi(i,X,U,K,Ak));
        if Xf_tmp == Xf
            break;
        else
            Xf = Xf_tmp;
        end
    end
    Xf = Xf - Z;
end


function [Xret] = Xpi(i,X,U,K,Ak)
    [F, G, ~] = Polyhedron2Matrix(X, U);
    Fpi = (F+G*K)*Ak^i;
    Xret = Polyhedron(Fpi, ones(size(Fpi,1), 1));
end
function H = compute_cost(Q,R,P,N)
    % compute H
    Q_ = [];
    R_ = [];
    for itr=1:N
        Q_ = blkdiag(Q_, Q);
        R_ = blkdiag(R_, R);
    end
    H = blkdiag(Q_, P, R_);
end

function [Aeq,beq] = compute_equalities(A,B,x_init,n_opt,N)

    % compute Aeq and beq

    [nx,nu] = size(B);
    A_ = [];
    B_ = [];
    for itr=1:N
        A_ = blkdiag(A_, A);
        B_ = blkdiag(B_, B);
    end         


    %Note: [x(0)...x(N)]^T = C_bar*[x(0)...x(N), u(0)...u(N-1)] + beq

    C_bar = [zeros(nx, n_opt);
             A_, zeros(nx*N, nx), B_]; 
    
    Aeq = eye(nx*(N+1), n_opt) - C_bar;
    beq = [x_init; zeros(size(Aeq)*[1; 0]-nx, 1)]; 
    
end

function [Aineq,bineq,nc] = compute_inequalities( Xc, Uc,N)

    % Compute inequality constraints from the given bounds defined as polyhedrons (MPT3)
    nx = size(Xc.A,2);
    nu = size(Uc.A,2);
    [F, G, nc] = Polyhedron2Matrix(Xc, Uc);
    
    F_block = [];
    G_block = [];
    for itr = 1:N
        G_block = blkdiag(G_block, G);
    end
    for itr = 1:N+1
        F_block = blkdiag(F_block, F);
    end
    Aineq = [F_block, [G_block; zeros(nc, nu*N)]];
    nc_total = size(Aineq,1);
    bineq = ones(nc_total, 1);
end


function [x_tmpc,u_tmpc] = controller(A,B,x_init,x_seq,u_seq,K,N)
    global Wc;
    x = x_init;
    x_tmpc = x;
    u_tmpc = [];
    time = 1;
    for i=1:N
        u = u_seq(:, time) + K*(x-x_seq(:, time));
        x = updateState(A,B,x, u, Wc);
        x_tmpc = [x_tmpc, x];
        u_tmpc = [u_tmpc, u];
        time = time + 1;
    end
    
end

function x = updateState(A,B,x, u, W)

    [nx,~] = size(B);
    w_min = min(W.V,[],1)';
    w_max = max(W.V,[],1)';
    w = w_min + (w_max - w_min).*rand(nx, 1);
    x = A*x+B*u+w;
end

function [F, G, nc] = Polyhedron2Matrix(X, U)
    
    nx = size(X.A,2);
    nu = size(U.A,2);

    F_tmp = X.A./repmat(X.b, 1, size(X.A)*[0; 1]);
    G_tmp = U.A./repmat(U.b, 1, size(U.A)*[0; 1]);

    if numel(F_tmp)==0
        F_tmp = zeros(0, nx);
    end
    if numel(G_tmp)==0
        G_tmp = zeros(0, nu);
    end
    F = [F_tmp; zeros(size(G_tmp,1), nx)];
    G = [zeros(size(F_tmp,1), nu); G_tmp];
    nc = size(F,1);
end

function [Aineq,bineq] = addIneqConstr(Aineq,bineq,Xadd, k_add, n)
    nx = n.nx;
    nu = n.nu;
    n_opt = n.n_opt;
    nc = n.nc;

    % add a new constraint at time step k 
    if Xadd.contains(zeros(2, 1)) % If Xadd contains the origin, the contraint can be expressed as C1*x<=1
        [F_add, ~, nc_add] = Polyhedron2Matrix(Xadd, Polyhedron());
        bineq = [bineq; ones(nc_add, 1)];
        
    else % in other cases, expressed in a general affine form C1*x<=C2
        F_add = Xadd.A;
        nc_add = size(F_add)*[1; 0];
        bineq = [bineq; Xadd.b];
    end
    
    block_add = zeros(nc_add, n_opt);
    block_add(:, (k_add-1)*nx+1:k_add*nx) = F_add;
    
    Aineq = [Aineq;
        block_add];
    nc = nc + nc_add;
end