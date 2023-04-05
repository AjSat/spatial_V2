%clearvars

clear K_con k_con k_con_inp Soft
import casadi.*

cs = SX;
csX = @SX;

model = load('iiwa_model.mat');
model = model.model;
%n = 7;  

%model = autoTree(n, 1, pi/2, 1);
constraint_per_ee = 6;

for i = 1:model.NB
    model.Xtree{1,i} = sparsify(model.Xtree{1,i}, 1e-10);
end

n = model.NB;
q = cs.sym('q', n, 1);
qd = cs.sym('qd', n, 1);
tau = cs.sym('tau', n, 1);
k_con_sym = cs.sym('tau', constraint_per_ee, 1);
k_con_inp{model.NB} = k_con_sym;

taurand = randn(n,1);
qrand = randn(n, 1);
qdrand = randn(n,1);
%taurand = ID(model, qrand, qdrand, q*0);

K_uncon{n} = [];
k_uncon{n} = [];
% FD computed using PV
K_con{model.NB} = randn(constraint_per_ee,6);
k_con{model.NB} = SX(constraint_per_ee,1); %zeros(6,1); %randn(6,1);
[qdd, nu, a_ee, Xee] = PV_tree(model, q, qd, tau, {}, K_con, k_con);
PV_tree_fun = Function('f_rob_dyn', {q, qd, tau, k_con_sym}, {qdd, nu});
[qdd, nu, a_ee, Xee] = PV_early_fast(model, q, qd, tau, {}, K_con, k_con_inp);
PV_tree_early_fun = Function('f_rob_dyn', {q, qd, tau, k_con_sym}, {qdd, nu});

[H,C,Vs,avp, a_grav_links, J, Xa, Xup] = HandC(model, q, qd);
% H_fun = Function('H_fun', {q}, {H, C});
% H_fun.n_instructions
J = jacobian(K_con{n}*Vs{n}, qd);
k_con_cor = K_con{n}*(avp{n} - a_grav_links{n});
H_chol = cholesky(H);
Y = forward_sub(H_chol, J');
inv_OSIM = casadi_symmetric(Y'*Y);
OSIM_chol = cholesky(inv_OSIM);
qdd_unconstrained = back_sub(H_chol', forward_sub(H_chol,tau - C));
if constraint_per_ee > 0
    lambda = back_sub(OSIM_chol', forward_sub(OSIM_chol,J*qdd_unconstrained + k_con_cor - k_con_sym));
    % lambda = back_sub(OSIM_chol', forward_sub(OSIM_chol,K_con{n}*Xa{n}*(J*qdd_unconstrained) + k_con_cor - k_con_sym));
    qdd = qdd_unconstrained - back_sub(H_chol', Y*lambda);
else
    lambda = 0;
    qdd = qdd_unconstrained;
end
LTL_fun = Function('f_rob_dyn', {q, qd, tau, k_con_sym}, {qdd, lambda});

Soft{n}.Ki = K_con{model.NB}; %SX.sym('K_con', 6, 6);
Soft{n}.ki = k_con{model.NB}; %SX.sym('K_con', 6, 1);
K_con_soft{model.NB} = {};
%k_con{model.NB} = {};
Soft{n}.Ri = 1e8;
[qdd, nu] = PV_tree(model, q, qd, tau, {}, K_con_soft, k_con, Soft);
PV_tree_soft_fun = Function('f_rob_dyn', {q, qd, tau, k_con_sym}, {qdd, nu});

% Mujoco-type soft function
[H,C,Vs,avp, a_grav_links, J, Xa, Xup] = HandC(model, q, qd);
J = jacobian(K_con{n}*Vs{n}, qd);
k_con_cor = K_con{n}*(avp{n} - a_grav_links{n});
H_nu = H + J'*Soft{n}.Ri*J;
H_chol = cholesky(H_nu);
qdd = back_sub(H_chol', forward_sub(H_chol, J'*Soft{n}.Ri*(-k_con_cor + k_con_sym) + tau - C));
LTL_soft_fun = Function('f_rob_dyn', {q, qd, tau, k_con_sym}, {qdd});


[qdd2, nu2] = PV_tree_fun(qrand, qdrand, taurand, k_con{model.NB});
[qdd3, ~] = PV_tree_early_fun( qrand, qdrand, taurand, k_con{model.NB});
[qdd4, nu4] = LTL_fun( qrand, qdrand, taurand, k_con{model.NB});
[qdd5, nu5] = PV_tree_soft_fun( qrand, qdrand, taurand, k_con{model.NB});
[qdd6] = LTL_soft_fun( qrand, qdrand, taurand, k_con{model.NB});

PV_tree_fun.generate(strcat(strcat('PV_chain_', num2str(model.NB)), '_6con.c'), struct('with_header', true))
PV_tree_early_fun.generate(strcat(strcat('PV_early_chain_', num2str(model.NB)), '6con.c'), struct('with_header', true))
LTL_fun.generate(strcat(strcat('LTL_chain_', num2str(model.NB)), '_6con.c'), struct('with_header', true))
PV_tree_soft_fun.generate(strcat(strcat('PV_soft_', num2str(model.NB)), '_6con.c'), struct('with_header', true))
LTL_soft_fun.generate(strcat(strcat('LTL_soft_', num2str(model.NB)), '_6con.c'), struct('with_header', true))

sqrt(sumsqr(qdd2 - qdd3))
sqrt(sumsqr(qdd2 - qdd4))
sqrt(sumsqr(qdd2 - qdd5))
sqrt(sumsqr(qdd6 - qdd2))
PV_tree_fun.n_instructions - PV_tree_fun.nnz_in - PV_tree_fun.nnz_out
PV_tree_early_fun.n_instructions - PV_tree_early_fun.nnz_in - PV_tree_early_fun.nnz_out
LTL_fun.n_instructions - LTL_fun.nnz_in - LTL_fun.nnz_out
PV_tree_soft_fun.n_instructions - PV_tree_soft_fun.nnz_in - PV_tree_soft_fun.nnz_out
LTL_soft_fun.n_instructions - LTL_soft_fun.nnz_in - LTL_soft_fun.nnz_out
