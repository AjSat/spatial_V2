% Compute the inv_OSIM for a floating-base robot

% Compute the Omega matrices at all branching points. Compute Projection matrix
% at all branching points. Then compute Kn for
% all the constraints and its ancestor branching points. Then populate the OSIM
% matrix one-by-one.

function  [inv_OSIM, IA, KA1] = LOSIM( model, x_fb, q, K_con)

import casadi.*;

model = branching_structure(model, K_con);

if strcmp(class(q), 'casadi.MX')
    cs = MX;
    csX = @MX;
else
    cs = SX;
    csX = @SX;
end


q = [0; q];

n = size(q, 1);

qn = x_fb(1:4);				% unit quaternion fixed-->f.b.
r = x_fb(5:7);				% position of f.b. origin
Xup{1} = plux( rq(qn), r );		% xform fixed --> f.b. coords


IA{1} = model.I{1};
KA{model.NB, model.NB} = [];
LA{model.NB, model.NB} = [];

for i = 2:model.NB
    if strcmp(model.jtype{i}, 'R3')
        [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) , model.axis{i});
    else
        [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i));
    end
    Xup{i} = XJ * model.Xtree{i};
    IA{i} = model.I{i};
end

ee_links = [];
for i = 1:model.NB
    m_i = size(K_con{i}, 1);
    if m_i > 6
        KA{i,i} = cs.eye(6);
        LA{i, i} = csX(6,6);
    else
        KA{i,i} = K_con{i};
        LA{i, i} = csX(m_i, m_i);
    end
    if m_i > 0
        ee_links  = [ee_links, i];
    end
end


parent = model.parent;
for i = 1:length(parent)
    if parent(i) == 0
        Xa{i} = Xup{i};
    else
        Xa{i} = Xup{i} * Xa{parent(i)};
    end
end

% First backward sweep
for i = model.NB:-1:2
    
    no_constraints = 0;
    U{i} = IA{i} * S{i};
    d{i} = S{i}' * U{i};
    
    Ia = IA{i} - U{i}/d{i}*U{i}';
    if strcmp(class(cs), 'casadi.SX')
        Ia = casadi_symmetric(Ia);
    end
    
    IA{model.parent(i)} = IA{model.parent(i)} + Xup{i}' * Ia * Xup{i};
    if strcmp(class(cs), 'casadi.SX')
        IA{model.parent(i)} = casadi_symmetric(IA{model.parent(i)});
    end
    
    if  ismember(i, model.repeating_links)
        child_link = i;
        if length(model.descendants{i}) > 1
            KA{i, i} = cs.eye(6);
            LA{i, i} = csX(6,6);
            child_link = i;
        end
    elseif    ~isempty(model.descendants{i})
        child_link = model.descendants{i}{1};
    else
        no_constraints = 1;
    end
    
    if no_constraints == 0
        child_link = model.Ndes{i};
        FS = KA{i, child_link}*S{i};
        KA{model.parent(i), child_link} =  (KA{i, child_link} - FS/d{i}*U{i}')*Xup{i};
        LA{model.parent(i), child_link} = casadi_symmetric(LA{i, child_link}) + FS/d{i}*(FS');
    end
    
end

% Compute Omega for the floating base link
LA{1,1} = casadi_symmetric(inv(IA{1}));

KLji{model.NB, model.NB} = [];
for i = model.repeating_links
    if i > 1
        i_anc = model.Nanc{i};
        KLji{i_anc, i} = KA{i_anc,i}*LA{1, i_anc};
        LA{1, i} = casadi_symmetric(casadi_symmetric(KLji{i_anc, i}*KA{i_anc,i}') + LA{i_anc, i});
    end
end



for ee = ee_links
    ee_now = ee;
    while model.Nanc{ee_now} ~= 0
        anc = model.Nanc{ee_now};
        if isempty(KA{anc, ee})
            KA{anc, ee} = KA{ee_now, ee}*KA{anc, ee_now};
        end
        ee_now = model.Nanc{ee_now};
    end
end

%     for i = model.repeating_links
%         for j = flatten_list(model.descendants{i})
%             if isempty(KLji{i,j})
%                 KLji{i, j} = KA{i,j}*LA{1, i};
%             end
%         end
%     end

%     for i = ee_links
%         if ~isempty(KLji{model.Nanc{i}, i})
%             KLji{model.Nanc{i}, i} = K_con{i}*KLji{model.Nanc{i}, i};
%         end
%     end

Omega{model.NB, model.NB} = [];

% assembling the OSIM matrix
for i = ee_links
    for j = ee_links
        if j > i
            gca = model.gca(i,j);
            if isempty(KLji{gca,j})
                KLji{gca, j} = KA{gca,j}*LA{1, gca};
            end
            Omega{j, i} = KLji{gca, j}*KA{gca, i}';
            if size(K_con{j}, 1) > 6
                Omega{j,i} = K_con{j}*Omega{j,i};
            end
            if size(K_con{i}, 1) > 6
                Omega{j,i} = Omega{j,i}*K_con{i}';
            end
            Omega{i, j} = Omega{j,i}';
        end
    end
    
    if size(K_con{i}, 1) <= 6
        Omega{i,i} = casadi_symmetric(LA{1,i}); % + KAjLAi*KA{i, j_descendants}');
    else
        Omega{i,i} = casadi_symmetric(K_con{i}*LA{1,i}*K_con{i}');
    end
end


%     % assembling the OSIM matrix
%     for i = model.repeating_links
%         for j = 1:length(model.descendants{i})
%             for j_descendants = model.descendants{i}{j}
% %                 fprintf("j_descendant = %d, i = %d\n", j_descendants, i)
%                 if isempty(KLji{i,j_descendants})
%                     KLji{i, j_descendants} = KA{i,j_descendants}*LA{1, i};
%                 end
% %                 KAjLAi = KA{i, j_descendants}*LA{1,i};
%                 for k = j+1:length(model.descendants{i})
%                     for k_descendants = model.descendants{i}{k}
%                         Omega{j_descendants,k_descendants} = KLji{i, j_descendants}*KA{i, k_descendants}';
%                         if size(K_con{j_descendants}, 1) > 6
%                              Omega{j_descendants,k_descendants} = K_con{j_descendants}*Omega{j_descendants,k_descendants};
%                         end
%                         if size(K_con{k_descendants}, 1) > 6
%                              Omega{j_descendants,k_descendants} = Omega{j_descendants,k_descendants}*K_con{k_descendants}';
%                         end
%                         Omega{k_descendants, j_descendants} = Omega{j_descendants,k_descendants}';
%                     end
%                 end
%                 if length(model.descendants{i}{j}) == 1
%                     if size(K_con{j_descendants}, 1) <= 6
%                         Omega{j_descendants,j_descendants} = casadi_symmetric(LA{1,j_descendants}); % + KAjLAi*KA{i, j_descendants}');
%                     else
%                         Omega{j_descendants,j_descendants} = casadi_symmetric(K_con{j_descendants}*LA{1,j_descendants}*K_con{j_descendants}');
%                     end
%                 end
%             end
%         end
%     end

% get the order in which to construct the invosim
constraints_supported{model.NB} = [];
for i = model.NB:-1:2
    if size(K_con{i}, 1) > 0
        constraints_supported{i} = [i, constraints_supported{i}];
    end
    constraints_supported{model.parent(i)} = [constraints_supported{model.parent(i)}, constraints_supported{i}];
end


%     gca = ones(4,4)*4;
%     gca(3,:) = 1;
%     gca(:,3) = 1;

%     gca = ones(12,12);
%     gca(1:2,1:2) = 4;
%     gca(3:7, 3:7) = 33;
%     gca(8:12, 8:12) = 57;
%     for i = 1:length(ee_links)
%         i_link = ee_links(i);
%         for j = i:length(ee_links)
%             j_link = ee_links(j);
%             if j == i
%                 Omega{i_link, i_link} = casadi_symmetric(LA{gca(i,j),j_link} + KLji{gca(i,j), i_link}*KA{gca(i,j), j_link}');
%             else
%                 Omega{i_link, j_link} = Omega{ee_links(i), ee_links(j)};
%                 Omega{j_link, i_link} = Omega{i_link, j_link}';
%             end
%         end
%     end


ee_links = constraints_supported{1}; %fliplr(ee_links);
inv_OSIM = [];
for i = 1:length(ee_links)
    inv_OSIMrow = [];
    for j = 1:length(ee_links)
        inv_OSIMrow = [inv_OSIMrow, Omega{ee_links(i), ee_links(j)}];
    end
    inv_OSIM = [inv_OSIM; inv_OSIMrow];
end
inv_OSIM = casadi_symmetric(inv_OSIM);

KA1 = collect_terms1(KA, model.NB);

end

function x = collect_terms (a, n)

x = [];
for i = 1:n
    if ~isempty(a{i})
        x = [x; a{i}];
    end
end
end

function x = collect_terms1(a, n)
x = [];
for i = 1:n
    if ~isempty(a{1,i})
        x = [x; vec(a{1,i})];
    end
end
end

% creates the branching structure. Generates a list of nodes where the
% kinematic tree branches out. Also, for every node has a list of lists of
% branching nodes. Each nested list corresponding to all the branching
% nodes from a specific branch rising from that node.

function model = branching_structure(model, K_con)

% go over the model.parent object. Run a counter for the number of
% times a link is another link's parent.

model.children{model.NB} = {};
body_child_counter = zeros(1, model.NB);
for i = 2:model.NB
    body_child_counter(model.parent(i)) = body_child_counter(model.parent(i)) + 1;
    model.children{model.parent(i)} = [model.children{model.parent(i)}, i];
end

% Then compile all the repeating links as the branching links.
repeating_links = [];
leaf_nodes = [];
%     for i = 1:model.NB
%         if body_child_counter(i) > 1
%             repeating_links = [repeating_links, i];
%         end
%         if body_child_counter(i) == 0
%             leaf_nodes = [leaf_nodes, i];
%         end
%     end

%     model.repeating_links = repeating_links;
%     model.leaf_nodes = leaf_nodes;

%reversed_parent_list = flip(model.parent);
%duplicate_parent_list = reversed_parent_list;

model.descendants{model.NB} = {};
model.child_branching_points{model.NB} = [];

for i = model.NB:-1:2
    
    if ~isempty(K_con{i})
        model.descendants{i}{length(model.descendants{i}) + 1} = i;
    end
    
    if ~isempty(model.descendants{i})
        
        if length(model.descendants{i}) > 1 || size(K_con{i}, 1) > 0
            
            repeating_links = [i, repeating_links];
            
            if i > 1
                model.child_branching_points{model.parent(i)} = [i, model.child_branching_points{model.parent(i)}];
            end
            
        end
        
        if i > 1
            model.descendants{model.parent(i)}{length(model.descendants{model.parent(i)}) + 1} = flatten_list(model.descendants{i});
            if ~ismember(i, repeating_links)
                model.child_branching_points{model.parent(i)} = [model.child_branching_points{i}, model.child_branching_points{model.parent(i)}];
            end
        end
    end
end

model.repeating_links = [1, repeating_links];

model.Nanc{model.NB} = {};

%     for i = model.NB:-1:1
%        for j = 1:length(model.descendants{i})
%            if length(model.descendants{i}{j}) == 1
%                model.Nanc{model.descendants{i}{j}} = i;
%            end
%        end
%     end
for i = model.repeating_links
    j = i;
    while model.parent(j) > 0
        j = model.parent(j);
        if ismember(j, model.repeating_links)
            model.Nanc{i} = j;
            break;
        end
    end
end

model.Ndes{model.NB} = [];
for i = model.NB:-1:2
    if ~isempty(model.Nanc{i})
        model.Ndes{i} = i;
    end
    model.Ndes{model.parent(i)} = model.Ndes{i};
end

% compute gca
gca = zeros(model.NB, model.NB);
for i = 1:model.NB
    for j = i:model.NB
        i_parent = i;
        j_parent = j;
        while(i_parent ~= j_parent)
            if i_parent > j_parent
                i_parent = model.parent(i_parent);
            else
                j_parent = model.parent(j_parent);
            end
        end
        gca(i,j) = i_parent;
        gca(j,i) = j_parent;
    end
end
model.gca = gca;

end

function flat_list =  flatten_list(cell_array)
flat_list = [];
for i = 1:length(cell_array)
    flat_list = [flat_list, cell_array{i}];
end
end

