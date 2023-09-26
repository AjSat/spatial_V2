% long tree that favours low-complexity algorithms like EFPA

%clearvars;

clear K_con gca k_con
import casadi.*

n = 10;
model = autoTree(n, 1, 1, 1);

% the tips of the tree

tip_dof = 7; 
model_tip = autoTree(tip_dof, 1, pi/3, 1);

% number of radiating tips

num_tips = 4;

ee_links = [];

gca = ones(num_tips*2,num_tips*2);
gca(num_tips+1:num_tips*2, num_tips + 1:num_tips*2) = n*ones(num_tips,num_tips);

    
point_of_attachments = [n, 1];

for k = point_of_attachments
    for i = 1:num_tips
        model.parent = [model.parent, model.NB + model_tip.parent];
        model.parent(model.NB+1) = k;
        model.NB = model.NB + tip_dof;
        if k == 1
            K_con{model.NB} = SX.eye(6);
        else
            K_con{model.NB} = SX.eye(6);
        end
        ee_links = [ee_links, model.NB];
        for j = 1:tip_dof
            model.Xtree{1, length(model.Xtree) + 1} = model_tip.Xtree{1, j};
            model.jtype{1, length(model.jtype) + 1} = model_tip.jtype{1, j};
            model.I{1, length(model.I) + 1} = model_tip.I{1, j};
            model.appearance.body{1, length(model.appearance.body) + 1} = model_tip.appearance.body{1, j};
        end
    end
end


for i = 1:model.NB
    model.I{1, i} = sparsify(SX(model.I{1,i}), 1e-12); 
    model.I{1,i} = casadi_symmetric(model.I{1,i});
    model.Xtree{1,i}(1:3,1:3) = SX.eye(3);
    model.Xtree{1,i}(4:6,4:6) = SX.eye(3);
end

model.appearance.base{1} = 'cyl';
model.appearance.base{2} = [0,0,0;1,0,0];
model.appearance.base{3} = 0.0005;


%showmotion(model, [0:1:4], randn(model.NB, 5))

% Compute inv_OSIM using PV-OSIM

q = randn(model.NB, 1);
x_fb0 = [1,0,0,0,0,0,0,0,0,0,0,0,0]';
qrand = randn(n + tip_dof*num_tips*2-1, 1);
q = SX.sym('q', model.NB-1, 1);
qd = SX.sym('q', model.NB-1, 1);
x_fb = SX.sym('x_fb', 13, 1);
pv_invosim = OSIM_fb(model, x_fb, q, K_con, []); 
pv_invosim2 = LOSIM(model, x_fb, q, K_con); 

disp('PV-OSIM stuff')
[pv_invosim, IA, KA, LA] = OSIM_fb(model, x_fb, q, K_con, []);
pv_invosim = casadi_symmetric(pv_invosim);
% pv_invosim_fun = Function('f_pv_invosim', {x_fb, q}, {IA{1}});
% pv_invosim_fun.n_instructions - pv_invosim_fun.nnz_in - pv_invosim_fun.nnz_out;
% pv_invosim_fun = Function('f_pv_invosim', {x_fb, q}, {KA{1}, IA{1}});
% pv_invosim_fun.n_instructions - pv_invosim_fun.nnz_in - pv_invosim_fun.nnz_out;
pv_invosim_fun = Function('f_pv_invosim', {x_fb, q}, {pv_invosim});
PVcount = pv_invosim_fun.n_instructions - pv_invosim_fun.nnz_in - pv_invosim_fun.nnz_out;
% pv_osim_chol = cholesky(pv_invosim);
% pv_invosim_fun_chol = Function('f_pv_invosim', {x_fb, q}, {pv_osim_chol, pv_invosim, KA{1}, IA{1}});
% pv_invosim_fun_chol.n_instructions - pv_invosim_fun_chol.nnz_in - pv_invosim_fun_chol.nnz_out

disp('LOSIM stuff')
[pv_invosim2, IA, KA] = LOSIM(model, x_fb, q, K_con); 
% pv_invosim2_fun = Function('f_pv_invosim2', {x_fb, q}, {IA{1}});
% LOSIMcount = pv_invosim2_fun.n_instructions - pv_invosim2_fun.nnz_in - pv_invosim2_fun.nnz_out;
% pv_invosim2_fun = Function('f_pv_invosim2', {x_fb, q}, {IA{1}, KA});
% LOSIMcount = pv_invosim2_fun.n_instructions - pv_invosim2_fun.nnz_in - pv_invosim2_fun.nnz_out;
pv_invosim2_fun = Function('f_pv_invosim2', {x_fb, q}, {pv_invosim2});
LOSIMcount = pv_invosim2_fun.n_instructions - pv_invosim2_fun.nnz_in - pv_invosim2_fun.nnz_out;



% pv_invosim_fun = Function('f_pv_invosim', {x_fb, q}, {pv_invosim});
% PVcount = pv_invosim_fun.n_instructions - pv_invosim_fun.nnz_in - pv_invosim_fun.nnz_out;
osim1 = pv_invosim_fun([1, zeros(1, 12)], qrand);
clear pv_invosim_fun pv_invosim

% pv_invosim2_fun = Function('f_pv_invosim2', {x_fb, q}, {pv_invosim2});
% PVrcount = pv_invosim2_fun.n_instructions - pv_invosim2_fun.nnz_in - pv_invosim2_fun.nnz_out;
osim2 = pv_invosim2_fun([1, zeros(1, 12)], qrand);
clear pv_invosim2_fun pv_invosim2

% pv_invosim2 = full(DM(pv_invosim2));
% pv_invosim = full(DM(pv_invosim));

gca = fliplr(flipud(gca));

Omega2 = EFP(model, q, x_fb, K_con, gca(1:end, 1:end));
constraints_per_ee = 6;
efp_osim = SX(num_tips*constraints_per_ee*2,num_tips*constraints_per_ee*2);
for i = 1:num_tips*2
    for j = 1:num_tips*2
        efp_osim((i-1)*constraints_per_ee + 1 : i*constraints_per_ee, (j-1)*constraints_per_ee + 1 : j*(constraints_per_ee)) = Omega2{ee_links(i), ee_links(j)};
    end
end
% 
% %efp_osim = constraints*efp_osim*constraints';
% %efp_osim_chol = cholesky(efp_osim);
efp_osim_fun = Function('f_rob_dyn', {x_fb, q}, {efp_osim});
EFPcount = efp_osim_fun.n_instructions - efp_osim_fun.nnz_in - efp_osim_fun.nnz_out;
osim3 = efp_osim_fun(x_fb0, qrand);
clear efp_invosim_fun efp_osim Omega2

LOSIMcount
PVcount
PVcount/LOSIMcount

% n
% num_tips
% 
% 
% % [H,C,X_up_all, Vs] = HandC_fb(model, x_fb, q, qd);
% % 
% % J = [];
% % for ind = ee_links
% %   J = [J; jacobian(K_con{ind}*Vs{ind}, vertcat(x_fb(8:13), qd))];
% % end
% % 
% % J = fliplr(J);
% % 
% % H = casadi_symmetric(H);
% % H = flipud(fliplr(H));
% % H_chol = cholesky(H);
% % 
% % Y = forward_sub(H_chol, J');
% % inv_OSIM = casadi_symmetric(Y'*Y);
% % LTL_invosim_fun = Function('f_LTL_invosim', {x_fb, q}, {inv_OSIM});
% % LTL_invosim_fun.n_instructions - LTL_invosim_fun.nnz_in - LTL_invosim_fun.nnz_out
% 
% PVcount
% PVrcount
PV_OSIM(num_tips, n) = PVcount;
PV_OSIMr(num_tips, n) = LOSIMcount;
EFPAr(num_tips, n) = EFPcount;
% 
losim_error = sqrt(sumsqr(vec(osim1 - osim2)))
% efpa_error = sqrt(sumsqr(vec(osim1 - osim3)))
% 
% 
% 
% 
% 
