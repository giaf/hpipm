function solver = hpipm_solver(m_qp_dims, m_qp_data)
    import hpipm_matlab.*
    hp = py.importlib.import_module('hpipm_python');
    py_qp_dims = hp.hpipm_dims();
    py_qp_dims = m2py_obj(m_qp_dims, py_qp_dims);
    
    py_qp_data = hp.hpipm_data();
    py_qp_data = m2py_obj(m_qp_data, py_qp_data);
    solver = hp.hpipm_solver(py_qp_dims, py_qp_data);
    
end