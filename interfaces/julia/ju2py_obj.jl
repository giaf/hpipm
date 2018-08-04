function ju2py_obj(ju_obj, py_obj)
    py_attr = keys(py_obj)
    n_attr = length(py_attr)
    for i = 1:n_attr
        py_thisattr = py_obj[py_attr[i]]
        if py_obj[py_attr[i]] != nothing
            error("Cannot assign to non-None attribute $(py_obj[py_attr[i]]) in Python object.")
        end
        eval(parse(string("py_attr_name = py_qp_data[:", String(py_attr[i]), "]")))
        ju_key_name = String(py_attr_name)
        if ju_obj[ju_key_name] != nothing
            ju_eltype = eltype(qp_data[ju_key_name])
            if ju_eltype == Array{Int64, 1} || ju_eltype == Array{Float64,1} || ju_eltype == Array{Int64,2} || ju_eltype == Array{Float64, 2}
            elseif
            elseif
            end
                        
        end
    end
    
end
