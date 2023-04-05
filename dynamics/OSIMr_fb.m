% Compute the inv_OSIM for a floating-base robot

% Compute the Omega matrices at all branching points. Compute Projection matrix
% at all branching points. Then compute Kn for
% all the constraints and its ancestor branching points. Then populate the OSIM
% matrix one-by-one.

function  [inv_OSIM] = OSIMr_fb( model, x_fb, q, K_con, k_con)

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
    [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );
    Xup{i} = XJ * model.Xtree{i};
    IA{i} = model.I{i};   
end

ee_links = [];
for i = 1:model.NB
    KA{i,i} = K_con{i};
    m_i = size(K_con{i}, 1);
    if m_i > 0
        ee_links  = [ee_links, i];
    end
    LA{i, i} = csX(m_i, m_i);    
end

parent = model.parent;
    for i = 1:length(parent)
        if parent(i) == 0
            Xa{i} = Xup{i};
        else
            Xa{i} = Xup{i} * Xa{parent(i)};
        end
    end

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
        
        if isempty(model.child_branching_points{i}) && ~isempty(model.descendants{i}) && ~ismember(i, model.repeating_links)
            child_link = model.descendants{i}{1};
        elseif ismember(i, model.repeating_links)
            KA{i, i} = cs.eye(6);
            LA{i, i} = csX(6,6);
            child_link = i;
        elseif ~isempty(model.child_branching_points{i})
            child_link = model.child_branching_points{i};
        else
            no_constraints = 1;
        end
        
        if no_constraints == 0
            FS = KA{i, child_link}*S{i};
            KA{model.parent(i), child_link} =  (KA{i, child_link} - FS/d{i}*U{i}')*Xup{i};
            %zero_off_diag_LM = csX(size(LA{model.parent(i)}, 1), size(KA{i}, 1));
            LA{model.parent(i), child_link} = LA{i, child_link} + FS/d{i}*(FS');
            if strcmp(class(cs), 'casadi.SX')
                LA{model.parent(i), child_link} = casadi_symmetric(LA{model.parent(i), child_link});
            end
            
        end
            
%         if size( KA{i}, 1) > 0
%             %PV
%             FS = KA{i}*S{i};
%             KA{model.parent(i)} = [KA{model.parent(i)};  (KA{i} - FS/d{i}*U{i}')*Xup{i}];
%             
%             LA{model.parent(i)} = [LA{model.parent(i)}, zero_off_diag_LM; zero_off_diag_LM', LA{i} + FS/d{i}*(FS')];
%             
%         end
        
    end
    
    LA{1,1} = casadi_symmetric(inv(IA{1}));
    
    for i = model.repeating_links
        if ~isempty(model.child_branching_points)
            for j = model.child_branching_points{i}
                LA{1, j} = casadi_symmetric(KA{i,j}*LA{1,i}*KA{i,j}' + LA{i, j});
            end
        end
    end
    
    Omega{model.NB, model.NB} = [];
    % assembling the OSIM matrix
    
    % assembling the diagonal terms
    for i = model.repeating_links
        for j = 1:length(model.descendants{i})
            if length(model.descendants{i}{j}) == 1
                ee_link = model.descendants{i}{j};
                Omega{ee_link, ee_link} = casadi_symmetric(LA{i,ee_link} + KA{i, ee_link}*LA{1,i}*KA{i, ee_link}');
            end
        end
    end
    
    % compute K{i, ES(i)} for all repeating links
    for i = fliplr(model.repeating_links)
       if ~isempty(model.child_branching_points{i})
           for j = model.child_branching_points{i}
               for k = flatten_list(model.descendants{j})
                   KA{i,k} = KA{j,k}*KA{i,j};
               end
           end
       end
    end
    
    % assembling all the sub-diagonal terms
    for i = fliplr(model.repeating_links)
        for j = 1:length(model.descendants{i})-1
            for j_descendants = model.descendants{i}{j}
                for k = j+1:length(model.descendants{i})
                    for k_descendants = model.descendants{i}{k}
                        Omega{j_descendants,k_descendants} = KA{i, j_descendants}*LA{1,i}*KA{i, k_descendants}';
                        Omega{k_descendants, j_descendants} = Omega{j_descendants,k_descendants}';
                    end
                end
            end
        end
    end
    
    ee_links = fliplr(ee_links);
    inv_OSIM = [];
    for i = 1:length(ee_links)
        inv_OSIMrow = [];
        for j = 1:length(ee_links)
            inv_OSIMrow = [inv_OSIMrow, Omega{ee_links(i), ee_links(j)}];
        end
        inv_OSIM = [inv_OSIM; inv_OSIMrow];
    end
    if strcmp(class(cs), 'casadi.SX')
        inv_OSIM = casadi_symmetric(inv_OSIM);
    end
    %inv_OSIM = Omega;
% if isempty(LA{1})
%     inv_OSIM = [];
% else
%     inv_OSIM = LA{1} + KA{1}*inv(IA{1})*KA{1}';
%     inv_OSIM = casadi_symmetric(inv_OSIM);
    %invosim_fun = Function('f_osim', {q(2:end)}, {casadi_symmetric(inv(full_osim))});
    
    % Uncomment to create a function
%     invosim_fun = Function('f_osim', {q(2:end)}, {cholesky(inv_OSIM)});
%     invosim_fun.n_instructions
    %invosim_fun.generate(strcat(strcat('osiminv_atlas_18D_', num2str(size(LA{1},1))), '.c'), struct('with_header', true))
%end

end

function x = collect_terms (a, n)
    
    x = [];
    for i = 1:n
        x = [x; a{i}];
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
    
    for i = model.NB:-1:1
        
        if ~isempty(K_con{i})
            model.descendants{i}{length(model.descendants{i}) + 1} = i;
        end
        
        if ~isempty(model.descendants{i})
            
            if length(model.descendants{i}) > 1
            
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
    
    model.repeating_links = repeating_links;
    
end
    
function flat_list =  flatten_list(cell_array)
    flat_list = [];
    for i = 1:length(cell_array)
        flat_list = [flat_list, cell_array{i}];
    end
end

