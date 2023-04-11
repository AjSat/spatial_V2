clearvars

import casadi.*

cs = SX;
csX = @SX;

model = load('../robot_models/talos_model.mat');
model = model.model;


%% Create functions by passing symbolic variables

n = model.NB - 1;
import casadi.*
clear K_con k_con
q = cs.sym('q', n, 1);
qd = cs.sym('qd', n, 1);
tau = cs.sym('tau', n, 1);
x_fb = cs.sym('x_fb', 13, 1);
feet_indices = [7, 13, 26, 40];
feet_indices = [40];
K_con{model.NB} = [];
k_con{model.NB} = [];
K_con_all = [];
k_con_all = [];


import casadi.*

%storing the gca for efp algorithm for the atlas robot
gca = zeros(4,4);
gca(1,2) = 1;
gca(1,3) = 1;
gca(1,4) = 1;
gca(2,3) = 1;
gca(2,4) = 1;
gca(3,4) = 15;



constraints_per_ee = 6;
feet_indices = [7, 13]; %7, 13
K_con_all = [];
k_con_all = [];
constraints_per_ee = 6;
for ind = feet_indices
   K_con{ind} =  cs.sym(strcat('k', num2str(ind)), constraints_per_ee ,6);% cs.eye(6); %[csX(3, 3), cs.eye(3)]; %cs.sym(strcat('K', num2str(ind)), 3, 6); %
   k_con{ind} = cs.sym(strcat('k', num2str(ind)), constraints_per_ee ,1); %csX(6,1);%
   K_con_all = [K_con_all; K_con{ind}];
   k_con_all = [k_con_all; k_con{ind}];
end

feet_indices2 = [28, 43];
feet_indices2 = [28, 43];
constraints_per_ee = 6;
for ind = feet_indices2
   K_con{ind} =  cs.sym(strcat('k', num2str(ind)), constraints_per_ee ,6);% ; %[csX(3, 3), cs.eye(3)]; %cs.sym(strcat('K', num2str(ind)), 3, 6); %
   k_con{ind} = cs.sym(strcat('k', num2str(ind)), constraints_per_ee ,1); %csX(6,1);%
   K_con_all = [K_con_all; K_con{ind}];
   k_con_all = [k_con_all; k_con{ind}];
end

% constraints = SX(24,24);
% for i = 1:4
%     constraints((i-1)*constraints_per_ee+1:i*constraints_per_ee, (i-1)*6+1:i*6) = K_con{feet_indices(i)};%rand(constraints_per_ee, 6);
% end

% kjr_osim = constraints*kjr_osim*constraints';
% 
% kjr_osim_chol = cholesky(kjr_osim);
% KJR_osim = Function('f_rob_dyn', {x_fb, q}, {kjr_osim});
% KJR_osim.n_instructions
%KJR_osim.generate('kjr_osim_24con.c', struct('with_header', true))


[Omega2, IA, KA] = EFP(model, q, x_fb, K_con, gca(1:end, 1:end));
efp_osim = [];

for i = fliplr([feet_indices, feet_indices2])
   efp_osim_row = []; 
   for j = fliplr([feet_indices, feet_indices2])
       efp_osim_row = [efp_osim_row, Omega2{i, j}];
   end
   efp_osim = [efp_osim; efp_osim_row];
end
% for i = 1:m/6
%     for j = 1:m/6
%         efp_osim((i-1)*constraints_per_ee + 1 : i*constraints_per_ee, (j-1)*constraints_per_ee + 1 : j*(constraints_per_ee)) = Omega2{feet_indices_r(i), feet_indices_r(j)};
%     end
% end

K_all = [];
for i = [feet_indices, feet_indices2]
   K_all = [K_all; KA{1, i}]; 
end

qrand = randn(44,1);
xfb_rand = randn(13,1);
Krand = randn(size(K_con_all));
krand = randn(size(k_con_all));

disp('EFP stuff')
efp_osim_fun = Function('f_rob_dyn', {x_fb, q}, {IA{1}});
efp_osim_fun.n_instructions - efp_osim_fun.nnz_in - efp_osim_fun.nnz_out
efp_osim_fun = Function('f_rob_dyn', {x_fb, q}, {K_all, IA{1}});
efp_osim_fun.n_instructions - efp_osim_fun.nnz_in - efp_osim_fun.nnz_out
efp_osim_fun = Function('f_rob_dyn', {x_fb, q}, {efp_osim, K_all, IA{1}});
efp_osim_fun.n_instructions - efp_osim_fun.nnz_in - efp_osim_fun.nnz_out
efp_osim_chol = cholesky(efp_osim);
efp_osim_fun = Function('f_rob_dyn', {x_fb, q, K_con_all, k_con_all}, {efp_osim_chol, efp_osim, K_all, IA{1}});
efp_osim_fun.n_instructions - efp_osim_fun.nnz_in - efp_osim_fun.nnz_out

disp('PV-OSIM stuff')

[pv_invosim, IA, KA, LA] = OSIM_fb(model, x_fb, q, K_con, []);
pv_invosim = casadi_symmetric(pv_invosim);
pv_invosim_fun = Function('f_pv_invosim', {x_fb, q}, {IA{1}});
pv_invosim_fun.n_instructions - pv_invosim_fun.nnz_in - pv_invosim_fun.nnz_out
pv_invosim_fun = Function('f_pv_invosim', {x_fb, q}, {KA{1}, IA{1}});
pv_invosim_fun.n_instructions - pv_invosim_fun.nnz_in - pv_invosim_fun.nnz_out
pv_invosim_fun = Function('f_pv_invosim', {x_fb, q}, {pv_invosim, KA{1}, IA{1}});
pv_invosim_fun.n_instructions - pv_invosim_fun.nnz_in - pv_invosim_fun.nnz_out
pv_osim_chol = cholesky(pv_invosim);
pv_invosim_fun = Function('f_pv_invosim', {x_fb, q, K_con_all, k_con_all}, {pv_osim_chol, pv_invosim, KA{1}, IA{1}});
pv_invosim_fun.n_instructions - pv_invosim_fun.nnz_in - pv_invosim_fun.nnz_out

[~, efp_invosim_val, ~, ~] = efp_osim_fun(xfb_rand, qrand, Krand, krand);
[~, pv_invosim_val, ~, ~] = pv_invosim_fun(xfb_rand, qrand, Krand, krand);

% Verifying the equality of OSIM using EFPA and PV-OSIM
assert(full(DM(sqrt(sumsqr(efp_invosim_val - pv_invosim_val)))) < 1e-10)

disp('PV-OSIM fast stuff')
[~, ~, ~, ~, ~, ~, LA, IA, KA] = PV_tr_fb(model, x_fb, q, qd, tau, [], K_con, k_con); 
Lambda_b = cholesky(LA{1});
Lcrossk = back_sub(Lambda_b', forward_sub(Lambda_b, KA{1}));
H_inv = inv(IA{1}); H_inv = casadi_symmetric(H_inv);
H_plus_cross_K = cholesky(IA{1} + casadi_symmetric(KA{1}'*Lcrossk));
%full_osim = Lambda_b + (Lcrossk * casadi_symmetric( H_inv + KA{1}'*Lcrossk))*Lcrossk';
%full_osim = casadi_symmetric(full_osim);
pv_invosim_fun = Function('f_pv_invosim', {x_fb, q}, {IA{1}});
pv_invosim_fun.n_instructions - pv_invosim_fun.nnz_in - pv_invosim_fun.nnz_out
pv_invosim_fun = Function('f_pv_invosim', {x_fb, q}, {KA{1}, IA{1}});
pv_invosim_fun.n_instructions - pv_invosim_fun.nnz_in - pv_invosim_fun.nnz_out
pv_invosim_fun = Function('f_pv_invosim', {x_fb, q}, {LA{1}, KA{1}, IA{1}});
pv_invosim_fun.n_instructions - pv_invosim_fun.nnz_in - pv_invosim_fun.nnz_out
pv_invosim_fun = Function('f_pv_invosim', {x_fb, q, K_con_all, k_con_all}, {Lambda_b, Lcrossk, H_plus_cross_K, LA{1}, KA{1}, IA{1}});
pv_invosim_fun.n_instructions - pv_invosim_fun.nnz_in - pv_invosim_fun.nnz_out
[Lambda_b_val, Lcrossk_val, H_plus_cross_K_val, LA1_val, KA1_val, IA1_val] = pv_invosim_fun(xfb_rand, qrand, Krand, krand);

% Verifying the correctness of PV-OSIM-fast
lambda_sol = back_sub(Lambda_b_val', forward_sub(Lambda_b_val, krand)) - Lcrossk_val*(back_sub(H_plus_cross_K_val', forward_sub(H_plus_cross_K_val, (Lcrossk_val'*krand))));
lambda_sol_from_pv_osim = inv(pv_invosim_val)*krand;
disp("Error between PV-OSIM and PV-OSIM-fast (not zero because of numerical reasons)")
full(DM(sqrt(sumsqr(lambda_sol - lambda_sol_from_pv_osim)/sumsqr(lambda_sol))))

% LTL-OSIM algorithm
[H,C,~, Vs] = HandC_fb(model, x_fb, q, qd);

J = [];
for ind = fliplr([feet_indices, feet_indices2])
  J = [J; jacobian(K_con{ind}*Vs{ind}, vertcat(x_fb(8:13), qd))];
end

%J = SX.sym('J', J.sparsity());
J = fliplr(J);

H = casadi_symmetric(H);
H = flipud(fliplr(H));
H_chol = cholesky(H);

Y = forward_sub(H_chol, J');
inv_OSIM = casadi_symmetric(Y'*Y);
% HinvJ = back_sub(H_chol', forward_sub(H_chol,J'));
%inv_OSIM = casadi_symmetric(J*HinvJ);
OSIM_chol = cholesky(inv_OSIM);

invOSIM_fun =  Function('f_rob_dyn', {q, x_fb}, {H});
invOSIM_fun.n_instructions - invOSIM_fun.nnz_in - invOSIM_fun.nnz_out
invOSIM_fun =  Function('f_rob_dyn', {q, x_fb}, {H_chol, H});
invOSIM_fun.n_instructions - invOSIM_fun.nnz_in - invOSIM_fun.nnz_out
invOSIM_fun =  Function('f_rob_dyn', {q, x_fb}, {J, H, H_chol});
invOSIM_fun.n_instructions - invOSIM_fun.nnz_in - invOSIM_fun.nnz_out
invOSIM_fun =  Function('f_rob_dyn', {q, x_fb}, {Y, J, H, H_chol});
invOSIM_fun.n_instructions - invOSIM_fun.nnz_in - invOSIM_fun.nnz_out
invOSIM_fun =  Function('f_rob_dyn', {q, x_fb}, {inv_OSIM, Y, J, H, H_chol});
invOSIM_fun.n_instructions - invOSIM_fun.nnz_in - invOSIM_fun.nnz_out
OSIM_chol = cholesky(inv_OSIM);
invOSIM_fun =  Function('f_rob_dyn', {x_fb, q, K_con_all, k_con_all}, {OSIM_chol, inv_OSIM, Y, J, H, H_chol});
invOSIM_fun.n_instructions - invOSIM_fun.nnz_in - invOSIM_fun.nnz_out

[~, LTL_invosim_val, Y_val, J_val, M_val, M_chol_val] = invOSIM_fun(xfb_rand, qrand, Krand, krand);

% uncomment below to see sparsity patterns
%spy(M_chol_val)
%figure, spy(Y_val)

assert(full(DM(sqrt(sumsqr(LTL_invosim_val - pv_invosim_val)))) < 1e-10)



