clearvars

import casadi.*

cs = SX;
csX = @SX;

model = load('talos_model.mat');
model = model.model;
% 

n = model.NB-1;
q = cs.sym('q', n, 1);
qd = cs.sym('qd', n, 1);
tau = cs.sym('tau', n, 1);
x_fb = cs.sym('xfb', 13, 1);

taurand = randn(n,1);
qrand = randn(n, 1);
qdrand = randn(n,1);
temp_rand = randn(4,1);
x_fb_rand = [temp_rand/sqrt(sumsqr(temp_rand)); randn(9,1)];


% FD computed using PV
constraint_per_ee = 6;
k_con_sym = [];
K_con{model.NB} = [];
k_con{model.NB} = [];
K_con{43} = cs.eye(6);%;
k_con{43} = cs.sym('k', constraint_per_ee, 1); k_con_sym = [k_con_sym; k_con{43}];
K_con{7} = cs.eye(6);%;
k_con{7} = cs.sym('k', constraint_per_ee, 1); k_con_sym = [k_con_sym; k_con{7}];
K_con{26} = cs.eye(6);%;
k_con{26} = cs.sym('k', constraint_per_ee, 1); k_con_sym = [k_con_sym; k_con{26}];
K_con{13} = cs.eye(6);%;
k_con{13} = cs.sym('k', constraint_per_ee, 1); k_con_sym = [k_con_sym; k_con{13}];

[qdd, nu, xd_fb, ~, ~, ~, ~, ~, ~, a_ee] = PV_tr_fb(model, x_fb, q, qd, tau, {}, K_con, k_con);
PV_tree_fun = Function('f_rob_dyn', {x_fb, q, qd, tau, k_con_sym}, {qdd, xd_fb});
[qdd, nu, xd_fb, a_ee] = PV_early_fast_fb(model, x_fb, q, qd, tau, {}, K_con, k_con);
PV_early_fun = Function('f_rob_dyn', {x_fb, q, qd, tau, k_con_sym}, {qdd, xd_fb});

soft_gauss_weight = SX.sym('weight', 1, 1);

% comparing with LTL
[H,C,~, Vs, avp, a_grav_links] = HandC_fb(model, x_fb, q, qd);
J = [];
feet_indices = [];
for i = 1:model.NB
    if size(K_con{i}, 1) > 0
        feet_indices = [feet_indices, i];
    end
end
k_con_all = [];
for ind = feet_indices
  J = [J; jacobian(K_con{ind}*Vs{ind}, vertcat(x_fb(8:13), qd))];
  k_con_all = [k_con_all; k_con{ind}];
end
%J = SX.sym('J', J.sparsity());
J = fliplr(J);

k_con_cor = [];
for ind = feet_indices
  k_con_cor = [k_con_cor; K_con{ind}*(avp{ind} - a_grav_links{ind})];
end

H = casadi_symmetric(H);
H = flipud(fliplr(H));
H_chol = cholesky(H);

Y = forward_sub(H_chol, J');
inv_OSIM = casadi_symmetric(Y'*Y);
% HinvJ = back_sub(H_chol', forward_sub(H_chol,J'));
%inv_OSIM = casadi_symmetric(J*HinvJ);
OSIM_chol = cholesky(inv_OSIM);
C = flipud(C);
tauLTL = flipud(tau);
qdd_unconstrained = back_sub(H_chol', forward_sub(H_chol,[tauLTL; zeros(6,1)] - C));
lambda = back_sub(OSIM_chol', forward_sub(OSIM_chol,J*qdd_unconstrained + k_con_cor - k_con_all));
qdd = flipud(qdd_unconstrained - back_sub(H_chol', Y*lambda));

LTL_fun = Function('f_rob_dyn', {x_fb, q, qd, tau, k_con_sym}, {qdd, J, k_con_cor});

% Mujoco type with LTL
H_nu = H + J'*soft_gauss_weight*J;
H_chol = cholesky(H_nu);
qdd = back_sub(H_chol', forward_sub(H_chol, J'*soft_gauss_weight*(-k_con_cor + k_con_sym) + [tauLTL; zeros(6,1)] - C));
qdd = flipud(qdd);
LTL_soft_fun = Function('f_rob_dyn', {x_fb, q, qd, tau, k_con_sym, soft_gauss_weight}, {qdd});

% The same for Mujoco type constraint now
K_con_soft{model.NB} = [];
Soft{model.NB} = [];
for i = 1:model.NB
    if size(K_con{i},1) > 0
    Soft{i}.Ki = K_con{i}; %SX.sym('K_con', 6, 6);
    Soft{i}.ki = k_con{i};
    Soft{i}.Ri = soft_gauss_weight;
    end
end
[qdd, nu, xd_fb, ~, ~, ~, ~, ~, ~, a_ee]= PV_tr_fb(model, x_fb, q, qd, tau, {}, K_con_soft, k_con, Soft);
PV_tree_soft_fun = Function('f_rob_dyn', {x_fb, q, qd, tau, k_con_sym, soft_gauss_weight}, {qdd, xd_fb});% qdd4_uncon = full(DM(flipud(qdd_unconstrained)));

early_q_res = zeros(1000,10);
LTL_q_res = zeros(1000,10);
soft_q_res = zeros(1000,10);
LTL_soft_q_res = zeros(1000,10);

PV_con_res = zeros(1000, 10);
early_con_res = zeros(1000,10);
LTL_con_res = zeros(1000,10);
soft_con_res = zeros(1000,10);
LTL_soft_con_res = zeros(1000,10);

%%
import casadi.*

for j = 1:1
    rng(100);
    weight = 10^j;
    j
for i = 1:1
taurand = randn(n,1);
qrand = randn(n, 1);
qdrand = randn(n,1);
temp_rand = randn(4,1);
x_fb_rand = [temp_rand/sqrt(sumsqr(temp_rand)); randn(9,1)];

n_con = length(feet_indices)*6;
[qdd2, xd_fb2] = PV_tree_fun(x_fb_rand, qrand, qdrand, taurand, zeros(n_con,1)); 
[qdd3, xd_fb3] = PV_early_fun(x_fb_rand, qrand, qdrand, taurand, zeros(n_con,1)); 
qn = x_fb_rand(1:4);
r = x_fb_rand(5:7);
Xup_fb = plux(rq(qn), r);
x_fb_fb = [x_fb_rand(1:7); Xup_fb*x_fb_rand(8:13)];
[qdd4, J, k_con_cor] = LTL_fun( x_fb_fb, qrand, qdrand, taurand, zeros(n_con,1));
[qdd5, xd_fb5] = PV_tree_soft_fun( x_fb_rand, qrand, qdrand, taurand, zeros(n_con,1), weight);
[qdd6] = LTL_soft_fun( x_fb_fb, qrand, qdrand, taurand, zeros(n_con,1), weight);

early_q_res(i,j) = full(DM(sqrt(sumsqr(qdd2 - qdd3))));
LTL_q_res(i,j) = full(DM(sqrt(sumsqr(qdd4(7:end) - qdd2(2:end)))));
soft_q_res(i,j) = full(DM(sqrt(sumsqr(qdd2 - qdd5))));
LTL_soft_q_res(i,j) = full(DM(sqrt(sumsqr(qdd6(7:end) - qdd2(2:end)))));

PV_con_res(i,j) = full(DM(sqrt(sumsqr(fliplr(J)*[Xup_fb*xd_fb2(7:12); qdd2(2:end)] + k_con_cor))));
early_con_res(i, j) = full(DM(sqrt(sumsqr(fliplr(J)*[Xup_fb*xd_fb3(7:12); qdd3(2:end)] + k_con_cor))));
LTL_con_res(i, j) = full(DM(sqrt(sumsqr(fliplr(J)*qdd4 + k_con_cor))));
soft_con_res(i, j) = full(DM(sqrt(sumsqr(fliplr(J)*[Xup_fb*xd_fb5(7:12); qdd5(2:end)] + k_con_cor))));
LTL_soft_con_res(i, j) = full(DM(sqrt(sumsqr(fliplr(J)*qdd6 + k_con_cor))));


PV_tree_fun.n_instructions - PV_tree_fun.nnz_in - PV_tree_fun.nnz_out
PV_early_fun.n_instructions - PV_early_fun.nnz_in - PV_early_fun.nnz_out
LTL_fun.n_instructions - LTL_fun.nnz_in - LTL_fun.nnz_out
PV_tree_soft_fun.n_instructions - PV_tree_soft_fun.nnz_in - PV_tree_soft_fun.nnz_out
LTL_soft_fun.n_instructions - LTL_soft_fun.nnz_in - LTL_soft_fun.nnz_out
%assert(full(DM(sumsqr(qdd2(2:end) - qdd))) <= 1e-12);

end
end

% Uncomment below to code-generate C files
% PV_tree_fun.generate('PV_talos_24con.c', struct('with_header', true))
% PV_early_fun.generate('PV_early_talos_24con.c', struct('with_header', true))
% LTL_fun.generate('LTL_talos_24con.c', struct('with_header', true))
% PV_tree_soft_fun.generate('PV_tree_soft_24con.c', struct('with_header', true))
% LTL_soft_fun.generate('LTL_soft_24con.c', struct('with_header', true))


