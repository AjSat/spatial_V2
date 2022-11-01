clearvars

import casadi.*

cs = SX;
csX = @SX;

model = load('atlas_model.mat');
model = model.model;
% 
floating_model = load('atlas_fb_model.mat');
floating_model = floating_model.fb_model;

n = model.NB - 1;
taurand = rand(n,1);
qrand = rand(n, 1);
qdrand = rand(n,1);
x_fb = [1,0,0,0,0,0,0,0,0,0,0,0,0]';

% Verifying that Featherstone's FD and PV agree in the absence of
% constraints
% floating base FB computed using featherstone's function
[xdfb, qdd] = FDfb(floating_model, [1,0,0,0,0,0,0,0,0,0,0,0,0]', qrand, zeros(n,1), taurand);
% FD computed using PV
K_con{model.NB} = [];
k_con{model.NB} = [];
[qdd2, ~, xd_fb] = PV_tr_fb(model, [1,0,0,0,0,0,0,0,0,0,0,0,0]', qrand, zeros(n,1), taurand, {}, K_con, k_con);
assert(full(DM(sumsqr(qdd2(2:end) - qdd))) <= 1e-12);

%% Compute the inverse OSIM using EFP and KRJ

import casadi.*

%storing the gca for efp algorithm for the atlas robot
gca = zeros(4,4);
gca(1,2) = 4;
gca(1,3) = 1;
gca(1,4) = 1;
gca(2,3) = 1;
gca(2,4) = 1;
gca(3,4) = 1;

K_con = {};
feet_indices = [11, 19, 25, 31];
feet_indices = [11, 19, 25, 31];
K_con = {};

% initialize constraints to identity to allow propagation in KRJ and EFP
for ind = feet_indices
   K_con{ind} = cs.eye(6); %[csX(3, 3), cs.eye(3)]; %cs.sym(strcat('K', num2str(ind)), 3, 6);
end

Omega = KRJ(model, qrand, x_fb, K_con);
feet_indices_r = fliplr(feet_indices);
m = length(feet_indices)*6;
kjr_invosim = SX(m,m);
for i = 1:m/6
    for j = 1:m/6
        kjr_invosim((i-1)*6 + 1 : i*6, (j-1)*6 + 1 : j*6) = Omega{feet_indices_r(i), feet_indices_r(j)};
    end
end

% Randomly generate constraint matrices
constraints_per_ee = 6;
constraints = zeros(constraints_per_ee*length(feet_indices), constraints_per_ee*length(feet_indices));
for i = 1:length(feet_indices)
    constraints((i-1)*constraints_per_ee+1:i*constraints_per_ee, (i-1)*6+1:i*6) = rand(constraints_per_ee, 6);
end


% Compute the OSIM after projecting onto the constraint space
kjr_invosim = constraints*kjr_invosim*constraints';


% Compute the OSIM using EFP
Omega2 = EFP(model, qrand, x_fb, K_con, gca(1:end, 1:end));
efp_invosim = SX(m,m);
for i = 1:m/6
    for j = 1:m/6
        efp_invosim((i-1)*6 + 1 : i*6, (j-1)*6 + 1 : j*6) = Omega2{feet_indices_r(i), feet_indices_r(j)};
    end
end

efp_invosim = constraints*efp_invosim*constraints';

% Initialize the constraint matrices with the randomly generated constraint
% matrices for PV
for i = 1:length(feet_indices)
    K_con{feet_indices_r(i)} =  (constraints((i-1)*constraints_per_ee+1:i*constraints_per_ee, (i-1)*6+1:i*6));
end

% Compute inv_OSIM using PV-OSIM
pv_invosim = OSIM_fb(model, x_fb, qrand, K_con, k_con); 
    

% Verify that EFP and KRJ algorithms give the same solution
assert(full(DM(sqrt(sumsqr(efp_invosim - kjr_invosim)))) < 1e-10)


% Verify that PV-OSIM and KRJ algorithms give the same solution
assert(full(DM(sqrt(sumsqr(pv_invosim - kjr_invosim)))) < 1e-10)
