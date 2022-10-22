% PV dynamics for an Iiwa robot
import casadi.*
clear
robot = importrobot('iiwa2.urdf');
Tf_base = eye(4);

setFixedTransform(robot.Bodies{1,1}.Joint, Tf_base);
    
%robot = importrobot('ur10_robot.urdf');

%to check if the ill-conditioning due to low inertia of later links is
 
model = urdf_to_featherstone(robot); %only for KUKA at the moment
model.joint_frictions = [1, 1, 1, 1, 0, 0, 0];

first_joint_no = 2; %3 for ur_10, 2 for iiwa

opti = casadi.Opti();

robot.DataFormat = 'column';
robot.Gravity = [0; 0; -9.81];
g = 9.81;
rob_casadi_model.robot = robot;
rob_casadi_model.model = model;

q = SX.sym('q', 7, 1);
qd = SX.sym('qd', 7, 1);
tau = SX.sym('tau', 7, 1);
x_acc_command = SX.sym('x_acc', 6, 1);
f_ext = {};
f_ext{7} = SX.sym('f_ext', 6, 1);
[qdd, nu, tau_ctrl, seg_accs, Xa, v, q_acc_X_fun, tau_X_fun] = PV_deprec(model, q, qd, tau, [], [ eye(6,6)], x_acc_command);
rob_dyn = [qd; qdd];

PV_fun = Function('f_rob_dyn', {q, qd, tau, x_acc_command, f_ext{7}}, {qdd, nu, tau_ctrl, seg_accs{7},Xa{7}, v{7},  q_acc_X_fun, tau_X_fun});
%PV_fun = Function('f_rob_dyn', {q, qd, tau, x_acc_command, f_ext{7}}, {qdd, nu});
%PV_fun.generate('PV_test.c', struct('main', true));
n_instructions = PV_fun.n_instructions
% Benchmark versus the forward dynamics algorithm
qrand = rand(7,1);
qdrand = rand(7,1)*0;

taurand = rand(7,1)*0;
%taurand = inverseDynamics(robot, qrand, qdrand, zeros(7,1))*0.9;
qdd_fd = forwardDynamics(robot, qrand, qdrand, taurand)

[qdd_res, nu_res, tau_ctrl_res, ee_acc_bf, Xa, vb] = PV_fun(qrand, qdrand, taurand, zeros(6,1), zeros(6,1))
% [qdd_res, nu_res, tau_ctrl_res, ee_acc_bf, Xa, ee_f_fun, q_acc_X_fun, tau_X_fun] = PV_fun(qrand, qdrand, taurand, [0; zeros(5,1)], zeros(6,1))

%end effector acceleration world frame
nu_res
ee_acc_wf = Xa\ee_acc_bf
Xa(1:3,1:3)'*ee_acc_bf(4:6)
tau_id = inverseDynamics(robot, qrand, qdrand, full(qdd_res))'
qdd_verify = forwardDynamics(robot, qrand, qdrand, full(tau_ctrl_res))'
full(qdd_res)'
tau_ctrl_res = full(tau_ctrl_res)'

% qdd = SX.sym('qdd', 7, 1);
% tau = IDwf(model, q, qd, qdd);
% ID_fun = Function('ID_wf', {q, qd, qdd}, {tau});

%% Plot the simulation of the KUKA robot executing the policy

q_traj = qrand;
qd_traj = qdrand;
[qdd_res, nu_res, tau_ctrl_res, ee_acc_bf, Xa, vb] = PV_fun(qrand, qdrand, taurand, zeros(6,1), zeros(6,1));
dt = 0.01;
q_traj_next = qrand;
qd_traj_next = qdrand;
vb = zeros(6,1);

xdesired = getTransform(robot, qrand, robot.BodyNames{8});
xdesired = xdesired(1:3,4);
vdesired = zeros(3,1);
vhist = zeros(3,1);

motion_con = zeros(6,1);
motion_con(4) = 0.0;
m_corr = motion_con;
for i = 1:1000
    
    % RK4
%     xnow = getTransform(robot, full(q_traj_next), robot.BodyNames{9});
%     verr = vdesired - xnow(1:3,1:3)*vb(4:6);
%     xnow = xnow(1:3,4);
%     m_corr = motion_con + 100*(xdesired - xnow) + 20*verr;
%     qd1 = qd_traj(:,end);
%     [qdd1, ~, ~, ~,~,vb] = PV_fun(q_traj_next, qd_traj_next,  0.1*qd1, m_corr, zeros(6,1));
%     
%     qp2 = q_traj_next + qd1*dt/2;
%     qd2 = qd1 + qdd1*dt/2;
%     xnow = getTransform(robot, full(q_traj_next), robot.BodyNames{9});
%     verr = vdesired - xnow(1:3,1:3)*vb(4:6);
%     xnow = xnow(1:3,4);
%     m_corr = motion_con + 100*(xdesired - xnow) + 20*verr;
%     [qdd2, ~, ~, ~,~,vb] = PV_fun(qp2, qd2, - 0.1*qd2, m_corr, zeros(6,1));
%     
%     qp3 = q_traj_next + qd2*dt/2;
%     qd3 = qd1 + qdd2*dt/2;
%     xnow = getTransform(robot, full(qp3), robot.BodyNames{9});
%     verr = vdesired - xnow(1:3,1:3)*vb(4:6);
%     xnow = xnow(1:3,4);
%     m_corr = motion_con + 100*(xdesired - xnow) + 20*verr;
%     [qdd3, ~, ~, ~,~,vb] = PV_fun(qp3, qd3, - 0.1*qd3, m_corr, zeros(6,1));
%     
%     qp4 = q_traj_next + qd3*dt;
%     qd4 = qd1 + qdd3*dt;
%     xnow = getTransform(robot, full(qp3), robot.BodyNames{9});
%     verr = vdesired - xnow(1:3,1:3)*vb(4:6);
%     xnow = xnow(1:3,4);
%     m_corr = motion_con + 100*(xdesired - xnow) + 20*verr;
%     [qdd4, ~, ~, ~,~,vb] = PV_fun(qp4, qd4, - 0.1*qd4, m_corr, zeros(6,1));
%     
%     q_traj_next = q_traj_next + 1/6*(qd1 + 2*qd2 + 2*qd3 + qd4)*dt;
%     qd_traj_next = 0.99*qd_traj_next + 1/6*(qdd1 + 2*qdd2 + 2*qdd3 + qdd4)*dt;
    
    
    
     xdes = xdesired;
     xdes(1) = xdes(1) + (i*dt)^2*0.5*motion_con(4);
     vdes = vdesired;
     vdes(1) = (i*dt)*motion_con(4);
    
    qd_traj(:,end) = qd_traj(:,end);
    q_traj_next = q_traj(:,end) + qd_traj(:,end)*dt + 0.5*dt^2*full(qdd_res);
    qd_traj_next = qd_traj(:,end) + dt*full(qdd_res);
    %[~, ~, ~, ~, ~, vb] = PV_fun(q_traj_next, qd_traj_next, - 0.0*qd_traj_next, m_corr, zeros(6,1));
    xnow = getTransform(robot, full(q_traj_next), robot.BodyNames{8});
    verr = vdes- xnow(1:3,1:3)*vb(4:6);
    vel = xnow(1:3,1:3)*vb(4:6)
    xnow = xnow(1:3,4);
    m_corr(4:6) = full(motion_con(4:6) + 100*(xdes - xnow) + 20*verr);


    error = norm(xdes - xnow)
    vel_err_norm = norm(verr)
    %tau_id = inverseDynamics(robot, full(q_traj_next), qdrand*0, zeros(7,1));
    [qdd_res, nu_res, tau_ctrl_res, ee_acc_bf, Xa, vb] = PV_fun(q_traj_next, qd_traj_next,  - 0.01*qd_traj_next, m_corr, zeros(6,1));
    
    q_traj = [q_traj, q_traj_next];
    qd_traj = [qd_traj, qd_traj_next];
    vhist = [vhist, -verr + vdesired];
end

%visualize_motion(robot, q_traj', 1000);

fk_end = getTransform(robot, qrand, robot.BodyNames{8})
fk_end2 = getTransform(robot, full(q_traj(:,end)), robot.BodyNames{8})

error = norm(fk_end(1:3, 4) - fk_end2(1:3, 4))
plot(full(q_traj'))
figure, plot(full(vhist'))
% Show the linearity of the dynamics w.r.t to all the motion acceleration
% constraints and the external wrenches.


% Compare the computation times of PV with the textbook constrained
% dynamics.