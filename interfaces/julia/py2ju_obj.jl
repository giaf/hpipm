function py2ju_obj(py_obj)
    prop = keys(py_obj)
    n_prop = length(prop)
    ju_obj = Dict()
    for i = 1:n_prop
        # skip base class attributes
        if !contains(String(prop[i]), "__") 
            if py_obj[prop[i]] != nothing
                error("Cannot convert Python object with non-None attribute value $(prop[i]).")
            end
            entry_name = String(prop[i])
            ju_obj[entry_name] = nothing
        end
    end
    return ju_obj
end
