function hpipm_data()
    py_qp_data = hp.hpipm_data()
    ju_qp_data = py2ju_obj(py_qp_data)
    return ju_qp_data
end

function hpipm_dims()
    py_qp_dims = hp.hpipm_dims()
    ju_qp_dims = py2ju_obj(py_qp_dims)
    return ju_qp_dims
end
