function horz_catenate = cell_to_mat(a)
    
    import casadi.*
    
    horz_catenate = a{1};
    for i = 2:length(a)
        horz_catenate = [horz_catenate, a{i}];
    end
end
        
        