function m_obj = py2m_obj(py_obj)    
    props = properties(py_obj);
    for iprop = 1:length(props)
        thisprop = props{iprop};
        thisprop_value = py_obj.(thisprop);
        if thisprop_value ~= py.NoneType
            error('Can only convert python objects with all member values set to None')
        end
    end
    n_props = length(props);
    empty_cell = cell(n_props, 1);
    m_obj = cell2struct(empty_cell, props, 1);   
end