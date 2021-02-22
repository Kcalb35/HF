#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <cmath>

#include "../include/gslextra.h"

gsl_quad_tensor * gsl_quad_tensor_alloc(int i, int j, int k, int l)
{
    gsl_quad_tensor *result = new gsl_quad_tensor;
    result->i = i;
    result->j = j;
    result->k = k;
    result->l = l;

    result->element = (gsl_matrix ***)malloc(i * sizeof(gsl_matrix **));
    for (int a = 0; a < i; a++)
    {
        *(result->element + a) = (gsl_matrix **)malloc(j * sizeof(gsl_matrix *));
        for (int b = 0; b < j; b++)
        {
            *(*(result->element + a) + b) = gsl_matrix_alloc(k, l);
        }
    }
    return result;
}

gsl_quad_tensor *gsl_quad_tensor_calloc(int i, int j, int k, int l)
{
    gsl_quad_tensor *result = new gsl_quad_tensor;
    result->i = i;
    result->j = j;
    result->k = k;
    result->l = l;

    result->element = (gsl_matrix ***)malloc(i * sizeof(gsl_matrix **));
    for (int a = 0; a < i; a++)
    {
        *(result->element + a) = (gsl_matrix **)malloc(j * sizeof(gsl_matrix *));
        for (int b = 0; b < j; b++)
        {
            *(*(result->element + a) + b) = gsl_matrix_calloc(k, l);
        }
    }
    return result;
}

void gsl_quad_tensor_free(gsl_quad_tensor *q)
{
    for (int i = 0; i < q->i; i++)
    {
        for (int j = 0; j < q->j; j++)
        {
            gsl_matrix_free(*(*(q->element + i) + j));
        }
        free(*(q->element + i));
    }
    free(q->element);
    free(q);
}

double gsl_quad_tensor_get(const gsl_quad_tensor *q, int i, int j, int k, int l)
{
    return gsl_matrix_get(*(*(q->element + i) + j), k, l);
}

void gsl_quad_tensor_get_matrix_copy(gsl_matrix *m, gsl_quad_tensor *q, int i, int j)
{
    gsl_matrix_memcpy(m, *(*(q->element + i) + j));
}

void gsl_quad_tensor_set(gsl_quad_tensor *q, int i, int j, int k, int l, double value)
{
    gsl_matrix_set(*(*(q->element + i) + j), k, l, value);
}

void gsl_quad_tensor_set_matrix(gsl_quad_tensor *q, int i, int j, int k, int l, const gsl_matrix *m)
{
    gsl_matrix_memcpy(*(*(q->element + i) + j), m);
}

void gsl_quad_tensor_add(gsl_quad_tensor *dst, gsl_quad_tensor *src)
{
    for (int i = 0; i < dst->i; i++)
    {
        for (int j = 0; j < dst->j; j++)
        {
            gsl_matrix_add(*(*(dst->element + i) + j), *(*(src->element + i) + j));
        }
    }
}

void gsl_quad_tensor_sub(gsl_quad_tensor *dst, gsl_quad_tensor *src)
{
    for (int i = 0; i < dst->i; i++)
    {
        for (int j = 0; j < dst->j; j++)
        {
            gsl_matrix_sub(*(*(dst->element + i) + j), *(*(src->element + i) + j));
        }
    }
}

void gsl_quad_tensor_scale(gsl_quad_tensor *dst, double x)
{
    for (int i = 0; i < dst->i; i++)
    {
        for (int j = 0; j < dst->j; j++)
        {
            gsl_matrix_scale(*(*(dst->element + i) + j), x);
        }
    }
}


/// Lowdin method
/// \param Fock Fork matrix
/// \param S overlap matrix
/// \param eigen eigen Energy order by ascending
/// \param eigenvector
void gsl_eigen_Lowdin_diag(gsl_matrix *Fock, gsl_matrix *S, gsl_vector *eigen, gsl_matrix *eigenvector)
{

    int length = S->size1;

    gsl_matrix * Fock_,* S_minus_half, * U,* S_cpy,*tmp;

    Fock_ = gsl_matrix_calloc(length,length);
    S_minus_half = gsl_matrix_calloc(length,length);
    S_cpy = gsl_matrix_calloc(length,length);
    U = gsl_matrix_calloc(length,length);
    tmp = gsl_matrix_calloc(length,length);

    // eigenvalues of S
    gsl_vector * S_eigen = gsl_vector_calloc(length);

    // prepare workspace
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(length);

    // copy S to S_cpy , because calculating eigenvalues will result in data destruction
    gsl_matrix_memcpy(S_cpy,S);

    // calculate U and eigen vector
    gsl_eigen_symmv(S,S_eigen,U,w);

    // set up S_minus_half
    for (int i = 0; i < length; i++)
    {
        gsl_matrix_set(S_minus_half,i,i,pow(gsl_vector_get(S_eigen,i),-0.5));
    }

    // calculate S minus half
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,U,S_minus_half,0,tmp);
    gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,tmp,U,0,S_minus_half);

    // calculate F'
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,S_minus_half,Fock,0,tmp);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,tmp,S_minus_half,0,Fock_);
    
    // calculate C'
    // C' store in tmp
    // output eigenvalue E by eigen
    gsl_eigen_symmv(Fock_,eigen,tmp,w);
    gsl_eigen_symmv_sort(eigen,tmp,GSL_EIGEN_SORT_VAL_ASC);

    // left multiply C' with S_minus_half to get C
    // output with eigenvector
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,S_minus_half,tmp,0,eigenvector);

    // free memories
    gsl_matrix_free(Fock_);
    gsl_matrix_free(S_minus_half);
    gsl_matrix_free(U);
    gsl_matrix_free(S_cpy);
    gsl_matrix_free(tmp);

    gsl_vector_free(S_eigen);

    gsl_eigen_symmv_free(w);
    
}