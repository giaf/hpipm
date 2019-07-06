SRCS =
SRCS += ocp_qp_dim_create.c
SRCS += ocp_qp_dim_destroy.c
SRCS += ocp_qp_dim_set.c
SRCS += ocp_qp_dim_print.c
SRCS += ocp_qp_dim_codegen.c
SRCS += ocp_qp_create.c
SRCS += ocp_qp_destroy.c
SRCS += ocp_qp_set.c
SRCS += ocp_qp_print.c
SRCS += ocp_qp_codegen.c
SRCS += ocp_qp_sol_create.c
SRCS += ocp_qp_sol_destroy.c
SRCS += ocp_qp_sol_get.c
SRCS += ocp_qp_sol_print.c
SRCS += ocp_qp_solver_arg_create.c
SRCS += ocp_qp_solver_arg_destroy.c
SRCS += ocp_qp_solver_arg_set.c
SRCS += ocp_qp_solver_arg_codegen.c
SRCS += ocp_qp_solver_create.c
SRCS += ocp_qp_solver_destroy.c
SRCS += ocp_qp_solver_solve.c

compile_mex_with_octave: $(SRCS)
	octave-cli compile_mex_ocp_qp.m

clean:
	rm -f *.so
	rm -f *.mex
	rm -f *.mexa64
	rm -f *.o
	rm -f cflags_octave.txt
