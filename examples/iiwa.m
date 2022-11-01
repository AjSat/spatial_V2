% PV constrained dynamics example for an Iiwa robot / a chain of arbitrary size
% with an end-effector constraint

clearvars

import casadi.*

% choosing SX versions of Casadi symbolic variables
cs = SX;
csX = @SX;

clear K_con k_con Soft
% Run below for using iiwa
model = load('../robot_models/iiwa_model.mat');
model = model.model;
n = model.NB;

%% Uncomment the two lines below to instead execute the solver on a chain of arbitrary size

% n = 31;
% model = autoTree(n, 1, 1, 1);

%% Build CasADi object of the solver

q = cs.sym('q', n, 1);
qd = cs.sym('qd', n, 1);
tau = cs.sym('tau', n, 1);

m = 6; % constraint dimension

% constraint matrix
K_con{n} = rand(m,6); %SX.eye(6); %[SX(3,3), SX.eye(3)]; %SX.sym('K_con', 6, 6);

% desired constraint accelerations
k_con{n} = SX.sym('k_con', m, 1); %
% Soft{n}.Ki = SX.eye(6); %SX.sym('K_con', 6, 6);
% Soft{n}.ki = SX(6,1); s%SX.sym('K_con', 6, 1);
% Soft{n}.Ri = 1e12;

% q = [0 1 2 3 4 5 6]';
% qd = zeros(7,1);
% tau = q;
% Soft{7}.Ki = cs.eye(6);
% Soft{7}.ki = csX(6,1);

[qdd, nu] = PV_tree(model, q, qd, tau, {}, K_con, k_con);
%[qdd, nu] = PV_tree(model, q, qd, tau, {}, K_con, k_con, Soft);

PV_fun = Function('f_rob_dyn', {q, qd, tau}, {qdd, nu});

PV_fun.n_instructions % print the number of instructions in the CasADi function

% PV_fun = Function('f_rob_dyn', {q, qd, tau, Soft{n}.Ki, Soft{n}.ki}, {qdd, nu});%, tau_ctrl, seg_accs{7},Xa{7}, v{7},  q_acc_X_fun, tau_X_fun});
% PV_fun = Function('f_rob_dyn', {q, qd, tau, Soft{n}.Ki, Soft{n}.ki}, {qdd, nu});

% uncomment the line below for C-code generation of the dynamics solver
%PV_fun.generate(strcat(strcat('iiwa_6con_', num2str(n)), '.c'), struct('with_header', true))

%% Testing
import casadi.*

% random initialization of robot state and torque input
qrand = rand(n,1);
qdrand = rand(n,1);
taurand = rand(n,1);
%taurand = inverseDynamics(robot, qrand, qdrand, zeros(7,1))*0.9;

% Implement forward dynamics using Featherstone's ABA solver
qdd_fd = FDab(model, qrand, qdrand, taurand);

% Evaluate constrained dynamics without any constraints (reduces to ABA)
K_con = {};
k_con = {};
K_con{n} = [];
k_con{n} = [];
[qdd_res, nu_res] = PV_tree(model, qrand, qdrand, taurand, {}, K_con, k_con);

% Verify that ABA and PV agree in the absence of constraints
assert(norm(full(DM(qdd_fd)) - full(DM(qdd_res))) < 1e-10)

% constrained case, requiring random end-effector acceleration
K_con{n} = SX.eye(6);
k_con{n} = rand(6,1);
[qdd_res, nu_res, a_ee, Xee] = PV_tree(model, qrand, qdrand, taurand, {}, K_con, k_con);
% verify that CD acceleration is consistent with the constraints
% Error not close to machine zero because of regularization
assert(norm(full(DM(Xee*a_ee - k_con{n}))) <= 1e-4);

% verifying joint accelerations of PV solver matches with Featherstone's
% ABA when the computed constraint force is added as an external force
f_ee = Xee' * (K_con{n}'*nu_res);
f_ext{n} = -f_ee;
qdd_fd = FDab(model, qrand, qdrand, taurand, f_ext);

assert(norm(full(DM(qdd_fd)) - full(DM(qdd_res))) < 1e-6)

% Verifying the PV_early output
[qdd_res2, nu_res2, a_ee2] = PV_tree_early(model, qrand, qdrand, taurand, {}, K_con, k_con);
full(DM(sqrt(sumsqr(qdd_res - qdd_res2))))
full(DM(sqrt(sumsqr(a_ee - a_ee2))))

%% Testing soft Gauss'

import casadi.*
% random initialization of robot state and torque input
qrand = rand(n,1);
qdrand = rand(n,1);
taurand = rand(n,1);
%taurand = inverseDynamics(robot, qrand, qdrand, zeros(7,1))*0.9;

% Implement forward dynamics using Featherstone's ABA solver
qdd_fd = FDab(model, qrand, qdrand, taurand);

% Evaluate constrained dynamics without any constraints (reduces to ABA)
K_con{n} = [];
k_con{n} = [];
[qdd_res, nu_res] = PV_tree(model, qrand, qdrand, taurand, {}, K_con, k_con);

% Verify that ABA and PV agree in the absence of constraints
assert(norm(full(DM(qdd_fd)) - full(DM(qdd_res))) < 1e-10)

% soft-constrained case, requiring random end-effector acceleration
Soft{n}.Ki = SX.eye(6); %SX.sym('K_con', 6, 6);
Soft{n}.ki = rand(6,1); %SX(6,1); %SX.sym('K_con', 6, 1);
Soft{n}.Ri = 1e12; % high quadratic penalty weight on soft-constraints
[qdd_res, nu_res, a_ee, Xee] = PV_tree(model, qrand, qdrand, taurand, {}, K_con, k_con, Soft);
a_ee

% verify that CD acceleration is consistent with the constraints
% not close to machine zero because of regularization
assert(norm(full(DM(Xee*a_ee - Soft{n}.ki))) <= 1e-4);
