_Pragma("once")

#include <gsl/gsl_matrix.h>

struct gsl_quad_tensor
{
    int i,j,k,l;
    gsl_matrix *** element;
};



gsl_quad_tensor *gsl_quad_tensor_alloc(int i, int j, int k, int l);
gsl_quad_tensor *gsl_quad_tensor_calloc(int i, int j, int k, int l);
void gsl_quad_tensor_free(gsl_quad_tensor *q);
double gsl_quad_tensor_get(const gsl_quad_tensor *q, int i, int j, int k, int l);
void gsl_quad_tensor_get_matrix_copy(gsl_matrix *m, gsl_quad_tensor *q, int i, int j);
void gsl_quad_tensor_set(gsl_quad_tensor *q, int i, int j, int k, int l, double value);
void gsl_quad_tensor_set_matrix(gsl_quad_tensor *q, int i, int j, int k, int l, const gsl_matrix *m);
void gsl_quad_tensor_add(gsl_quad_tensor *dst, gsl_quad_tensor *src);
void gsl_quad_tensor_sub(gsl_quad_tensor *dst, gsl_quad_tensor *src);
void gsl_quad_tensor_scale(gsl_quad_tensor *dst, double x);
void gsl_eigen_Lowdin_diag(gsl_matrix *Fock, gsl_matrix *S, gsl_vector *eigen, gsl_matrix *eigenvector);

void gsl_matrix_print(gsl_matrix * m,int precision=6,int width=14);
void gsl_matrix_print(gsl_matrix * m,int rows,int cols ,int precision=6,int width=14);
