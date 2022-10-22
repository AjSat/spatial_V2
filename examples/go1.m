clearvars

import casadi.*

cs = SX;
csX = @SX;

% Load the floating-base robot model for the PV solver
model = load('../robot_models/go1_altered_prefb.mat');
model = model.model;

% Load robot model suitable for Featherstone's floating-base ABA algorithm
floating_model = load('../robot_models/go1_altered_fb.mat');
floating_model = floating_model.floating_model;


% Verify that PV floating base and Featherstone's algorithm agree for the
% unconstrained case
qrand = rand(12,1);
taurand = rand(12,1);
qdrand = rand(12,1);
x_base = [1,0,0,0,0,0,0,0,0,0,0,0,0]';
[xdfb, qdd] = FDfb(floating_model, x_base, qrand, qdrand, taurand);

K_con{13} = [];
k_con{13} = [];
[qdd2, ~, xd_fb] = PV_tr_fb(model, x_base, qrand, qdrand, taurand, {}, K_con, k_con);
assert(full(DM(sumsqr(qdd2(2:end) - qdd))) <= 1e-12);

% Add a 3D constraint on the all the four feet of the quadruped
K_con{13} = [cs.eye(3), csX(3,3)];
k_con{13} = csX(3,1);
K_con{10} = [cs.eye(3), csX(3,3)];
k_con{10} = csX(3,1);
K_con{7} = [cs.eye(3), csX(3,3)];
k_con{7} = csX(3,1);
K_con{4} = [cs.eye(3), csX(3,3)];
k_con{4} = csX(3,1);
q = csX([0,0,0,0,0,0,0,0,0,0,0,0]');
[qdd3, ~, xd_fb] = PV_tr_fb(model, csX([0,0,0,1,0,0,0,0,0,0,0,0,0]'), q, zeros(12,1), taurand, {}, K_con, k_con);


% Add the base constraints as soft constraints
clear K_con k_con
Soft{13} = struct;
Soft{13}.Ri = 1e12;
Soft{13}.Ki = [cs.eye(3), csX(3,3)];
Soft{13}.ki = csX(3,1);
Soft{10}.Ri = 1e12;
Soft{10}.Ki = [cs.eye(3), csX(3,3)];
Soft{10}.ki = csX(3,1);
Soft{7}.Ri = 1e12;
Soft{7}.Ki = [cs.eye(3), csX(3,3)];
Soft{7}.ki = csX(3,1);
Soft{4}.Ri = 1e12;
Soft{4}.Ki = [cs.eye(3), csX(3,3)];
Soft{4}.ki = csX(3,1);
K_con{13} = [];
k_con{13} = [];
[qdd4, ~, xd_fb4] = PV_tr_fb(model, csX([0,0,0,1,0,0,0,0,0,0,0,0,0]'), csX(zeros(12,1)), zeros(12,1), taurand, {}, K_con, k_con, Soft);

%% Create functions by passing symbolic variables

import casadi.*
clear K_con k_con
q = cs.sym('q', 12, 1);
qd = cs.sym('qd',12, 1);
tau = cs.sym('tau', 12, 1);
x_fb = cs.sym('x_fb', 13, 1);
feet_indices = [4, 7, 10, 13];
K_con_all = [];
k_con_all = [];
for ind = feet_indices
   K_con{ind} = cs.sym(strcat('K', num2str(ind)), 3, 6); %[csX(3, 3), cs.eye(3)]; %
   k_con{ind} = cs.sym(strcat('k', num2str(ind)), 3 ,1);
   K_con_all = horzcat(K_con_all, K_con{ind});
   k_con_all = [k_con_all; k_con{ind}];
end

[qdd, nu, xd_fb, Xs, Vs, As] = PV_tr_fb(model, x_fb, q, qd, tau, {}, K_con, k_con);

PV_jac = jacobian(qdd, x_fb);

PV_fun = Function('f_rob_dyn', {x_fb, q, qd, tau, K_con_all, k_con_all }, {qdd, nu, xd_fb, Xs, Vs, As});
% PV_fun.save('go1_all_con.casadi')

% Uncomment below to C-code generat the constrained dynamics solver
% PV_fun = Function('f_rob_dyn', {x_fb, q, qd, tau, K_con_all, k_con_all }, {qdd, nu, xd_fb});%, tau_ctrl, [seg_accs{30}; seg_accs{20}],[Xa{30}; Xa{20}], [v{30}; v{20}]});%,  q_acc_X_fun, tau_X_fun});
%PV_fun.generate('PV2_go1.c', struct('with_header', true))
PV_fun.n_instructions
%%
% create symbolic variables for soft PV

clear K_con k_con Soft
Soft{13} = struct;
K_con_all = [];
k_con_all = [];
for ind = feet_indices
    Soft{ind}.Ki = [csX(3, 3), cs.eye(3)]; %cs.sym(strcat('Ks', num2str(ind)), 3, 6);
    Soft{ind}.ki = cs.sym(strcat('ks', num2str(ind)), 3 ,1);
    Soft{ind}.Ri = 1e12;
    K_con_all = [K_con_all, Soft{ind}.Ki];
    k_con_all = [k_con_all; Soft{ind}.ki];
end

K_con{13} = [];
k_con{13} = [];
[qdd, nu, xd_fb] = PV_tr_fb(model, x_fb, q, qd, tau, {}, K_con, k_con, Soft);
PV_fun_soft = Function('f_rob_dyn', {x_fb, q, qd, tau, k_con_all}, {qdd, nu, xd_fb});

% uncomment below to C-code generate the soft constrained dynamics solver
% PV_fun_soft = Function('f_rob_dyn', {x_fb, q, qd, tau, K_con_all, k_con_all}, {qdd, nu, xd_fb});%, tau_ctrl, [seg_accs{30}; seg_accs{20}],[Xa{30}; Xa{20}], [v{30}; v{20}]});%,  q_acc_X_fun, tau_X_fun});
% PV_fun_soft.generate('PV2_soft_go1.c', struct('with_header', true))

%% Verifying the output of hard-constrained dynamics solver

import casadi.*

clear K_con k_con
q = rand(12, 1);
qd = rand(12, 1);
tau = rand(12, 1);
x_fb = rand(13, 1);
x_fb(1:4) = x_fb(1:4) / norm(x_fb(1:4));

feet_indices = [4, 7, 10, 13];
K_con_all = [];
k_con_all = [];

% random constraint on the translation component of the feet
for ind = feet_indices
   K_con{ind} = [zeros(3,3), eye(3,3)]; %
   k_con{ind} = rand(3,1);
end

[qdd, nu, xd_fb, Xs, Vs, As] = PV_tr_fb(model, x_fb, q, qd, tau, {}, K_con, k_con);

% Verify that the resulting accelerations agree with the constraint
for ind = feet_indices
   aee = As(6*(ind - 1)+1 : 6*(ind));
   assert(full(DM(norm(aee(4:6) - k_con{ind}))) <= 1e-4);
end