#include <R.h>
#include <Rinternals.h>
#include <R_ext/Random.h>
#include <math.h>
#include <string.h>

/* Build a named list SEXP from parallel arrays of keys and values.
   Caller must PROTECT the value SEXPs before passing them. */
static SEXP make_named_list(const char **names, SEXP *values, int n) {
    SEXP list  = PROTECT(allocVector(VECSXP, n));
    SEXP nmsxp = PROTECT(allocVector(STRSXP, n));
    for (int i = 0; i < n; i++) {
        SET_VECTOR_ELT(list, i, values[i]);
        SET_STRING_ELT(nmsxp, i, mkChar(names[i]));
    }
    setAttrib(list, R_NamesSymbol, nmsxp);
    UNPROTECT(2);
    return list;
}

/* X is n×d (R column-major): X[i,j] = xp[j*n + i]
   centers stored d×nc: centers[k*d + j] = feature j of cluster k
   U stored n×nc column-major: U[k*n + i] = u_{i,k}  (matches allocMatrix(n,nc)) */

SEXP fcm_standard_c(SEXP X_sexp, SEXP k_sexp, SEXP m_sexp,
                    SEXP tol_sexp, SEXP max_sexp)
{
    const double *xp  = REAL(X_sexp);
    int n             = nrows(X_sexp);
    int d             = ncols(X_sexp);
    int nc            = asInteger(k_sexp);
    double m          = asReal(m_sexp);
    double tol        = asReal(tol_sexp);
    int maxiter       = asInteger(max_sexp);
    double exp_inv    = 1.0 / (m - 1.0);

    double *U        = R_Calloc((size_t)n * nc, double);
    double *U_old    = R_Calloc((size_t)n * nc, double);
    double *Um       = R_Calloc((size_t)n * nc, double);
    double *centers  = R_Calloc((size_t)d * nc, double);
    double *dist2    = R_Calloc((size_t)n * nc, double);
    double *obj_buf  = R_Calloc((size_t)maxiter, double);
    int    *zero_idx = R_Calloc((size_t)nc, int);

    /* Random initialisation of U, row-normalised */
    GetRNGstate();
    for (int i = 0; i < n; i++) {
        double rs = 0.0;
        for (int k = 0; k < nc; k++) {
            double v = unif_rand();
            U[k*n + i] = v;
            rs += v;
        }
        for (int k = 0; k < nc; k++) U[k*n + i] /= rs;
    }
    PutRNGstate();

    int iter_count = 0;

    for (int iter = 0; iter < maxiter; iter++) {
        iter_count = iter + 1;

        /* Copy U → U_old */
        memcpy(U_old, U, (size_t)n * nc * sizeof(double));

        /* Um = U^m */
        for (int k = 0; k < nc; k++)
            for (int i = 0; i < n; i++)
                Um[k*n + i] = pow(U[k*n + i], m);

        /* Center update: centers[k*d+j] = sum_i(X[i,j]*Um[k*n+i]) / sum_i(Um[k*n+i]) */
        for (int k = 0; k < nc; k++) {
            double sum_um = 0.0;
            for (int i = 0; i < n; i++) sum_um += Um[k*n + i];
            for (int j = 0; j < d; j++) {
                double acc = 0.0;
                for (int i = 0; i < n; i++)
                    acc += xp[j*n + i] * Um[k*n + i];
                centers[k*d + j] = acc / sum_um;
            }
        }

        /* Squared distances: dist2[k*n+i] = ||x_i - center_k||^2 */
        for (int k = 0; k < nc; k++) {
            for (int i = 0; i < n; i++) {
                double s = 0.0;
                for (int j = 0; j < d; j++) {
                    double diff = xp[j*n + i] - centers[k*d + j];
                    s += diff * diff;
                }
                dist2[k*n + i] = s;
            }
        }

        /* Objective: sum(Um * dist2) */
        double obj = 0.0;
        for (int k = 0; k < nc; k++)
            for (int i = 0; i < n; i++)
                obj += Um[k*n + i] * dist2[k*n + i];
        obj_buf[iter] = obj;

        /* Membership update */
        for (int i = 0; i < n; i++) {
            int nz = 0;
            for (int k = 0; k < nc; k++)
                if (dist2[k*n + i] == 0.0) zero_idx[nz++] = k;

            if (nz > 0) {
                for (int k = 0; k < nc; k++) U[k*n + i] = 0.0;
                double eq = 1.0 / nz;
                for (int z = 0; z < nz; z++) U[zero_idx[z]*n + i] = eq;
            } else {
                double sum_w = 0.0;
                for (int k = 0; k < nc; k++) {
                    double w = pow(dist2[k*n + i], -exp_inv);
                    Um[k*n + i] = w;  /* reuse Um as weight buffer */
                    sum_w += w;
                }
                for (int k = 0; k < nc; k++)
                    U[k*n + i] = Um[k*n + i] / sum_w;
            }
        }

        /* Convergence */
        double max_diff = 0.0;
        for (int k = 0; k < nc; k++)
            for (int i = 0; i < n; i++) {
                double diff = fabs(U[k*n + i] - U_old[k*n + i]);
                if (diff > max_diff) max_diff = diff;
            }
        if (max_diff < tol) break;
    }

    /* Build output */
    SEXP centers_sexp = PROTECT(allocMatrix(REALSXP, nc, d));
    double *cp = REAL(centers_sexp);
    for (int j = 0; j < d; j++)
        for (int k = 0; k < nc; k++)
            cp[j*nc + k] = centers[k*d + j];

    SEXP U_sexp = PROTECT(allocMatrix(REALSXP, n, nc));
    memcpy(REAL(U_sexp), U, (size_t)n * nc * sizeof(double));

    SEXP obj_sexp = PROTECT(allocVector(REALSXP, iter_count));
    memcpy(REAL(obj_sexp), obj_buf, (size_t)iter_count * sizeof(double));

    R_Free(U); R_Free(U_old); R_Free(Um); R_Free(centers);
    R_Free(dist2); R_Free(obj_buf); R_Free(zero_idx);

    const char *names[] = {"centers", "U", "obj"};
    SEXP vals[] = {centers_sexp, U_sexp, obj_sexp};
    SEXP result = PROTECT(make_named_list(names, vals, 3));
    UNPROTECT(4);
    return result;
}

SEXP fcm_infinity_c(SEXP X_sexp, SEXP k_sexp, SEXP m_sexp,
                    SEXP tol_sexp, SEXP max_sexp, SEXP rho_sexp)
{
    const double *xp  = REAL(X_sexp);
    int n             = nrows(X_sexp);
    int d             = ncols(X_sexp);
    int nc            = asInteger(k_sexp);
    double m          = asReal(m_sexp);
    double tol        = asReal(tol_sexp);
    int maxiter       = asInteger(max_sexp);
    double rho        = asReal(rho_sexp);
    double exp_inv    = 1.0 / (m - 1.0);

    double *U        = R_Calloc((size_t)n * nc, double);
    double *U_old    = R_Calloc((size_t)n * nc, double);
    double *Um       = R_Calloc((size_t)n * nc, double);
    double *U_inf    = R_Calloc((size_t)n, double);
    double *U_inf_old = R_Calloc((size_t)n, double);
    double *centers  = R_Calloc((size_t)d * nc, double);
    double *dist2    = R_Calloc((size_t)n * nc, double);
    double *row_sums = R_Calloc((size_t)n, double);
    double *min_dists = R_Calloc((size_t)n, double);
    double *d_inf    = R_Calloc((size_t)n, double);
    double *obj_buf  = R_Calloc((size_t)maxiter, double);
    int    *zero_idx = R_Calloc((size_t)nc, int);

    /* Random initialisation of U (regular clusters only), row-normalised */
    GetRNGstate();
    for (int i = 0; i < n; i++) {
        double rs = 0.0;
        for (int k = 0; k < nc; k++) {
            double v = unif_rand();
            U[k*n + i] = v;
            rs += v;
        }
        for (int k = 0; k < nc; k++) U[k*n + i] /= rs;
        U_inf[i] = 0.0;
    }
    PutRNGstate();

    int iter_count = 0;

    for (int iter = 0; iter < maxiter; iter++) {
        iter_count = iter + 1;

        memcpy(U_old,     U,     (size_t)n * nc * sizeof(double));
        memcpy(U_inf_old, U_inf, (size_t)n      * sizeof(double));

        /* Row sums of regular-cluster memberships */
        for (int i = 0; i < n; i++) {
            double rs = 0.0;
            for (int k = 0; k < nc; k++) rs += U[k*n + i];
            row_sums[i] = rs;
        }

        /* Um = (U / row_sum)^m */
        for (int k = 0; k < nc; k++)
            for (int i = 0; i < n; i++)
                Um[k*n + i] = pow(U[k*n + i] / row_sums[i], m);

        /* Center update */
        for (int k = 0; k < nc; k++) {
            double sum_um = 0.0;
            for (int i = 0; i < n; i++) sum_um += Um[k*n + i];
            for (int j = 0; j < d; j++) {
                double acc = 0.0;
                for (int i = 0; i < n; i++)
                    acc += xp[j*n + i] * Um[k*n + i];
                centers[k*d + j] = acc / sum_um;
            }
        }

        /* Squared distances */
        for (int k = 0; k < nc; k++) {
            for (int i = 0; i < n; i++) {
                double s = 0.0;
                for (int j = 0; j < d; j++) {
                    double diff = xp[j*n + i] - centers[k*d + j];
                    s += diff * diff;
                }
                dist2[k*n + i] = s;
            }
        }

        /* min_dists and rho check */
        for (int i = 0; i < n; i++) {
            double md = dist2[0*n + i];
            for (int k = 1; k < nc; k++)
                if (dist2[k*n + i] < md) md = dist2[k*n + i];
            min_dists[i] = sqrt(md);
        }
        for (int i = 0; i < n; i++) {
            if (min_dists[i] > rho) {
                /* Build rho_violation sentinel before freeing scratch */
                SEXP s_type  = PROTECT(mkString("rho_violation"));
                SEXP s_rho   = PROTECT(ScalarReal(rho));
                SEXP s_point = PROTECT(ScalarInteger(i + 1));
                SEXP s_ndist = PROTECT(ScalarReal(min_dists[i]));

                R_Free(U); R_Free(U_old); R_Free(Um); R_Free(U_inf);
                R_Free(U_inf_old); R_Free(centers); R_Free(dist2);
                R_Free(row_sums); R_Free(min_dists); R_Free(d_inf);
                R_Free(obj_buf); R_Free(zero_idx);

                const char *names[] = {"type", "rho", "point", "nearest_dist"};
                SEXP vals[] = {s_type, s_rho, s_point, s_ndist};
                SEXP sentinel = PROTECT(make_named_list(names, vals, 4));
                UNPROTECT(5);
                return sentinel;
            }
        }

        /* d_inf and objective */
        for (int i = 0; i < n; i++) d_inf[i] = rho - min_dists[i];

        double obj = 0.0;
        for (int k = 0; k < nc; k++)
            for (int i = 0; i < n; i++)
                obj += pow(U[k*n + i], m) * dist2[k*n + i];
        for (int i = 0; i < n; i++)
            obj += pow(U_inf[i], m) * (d_inf[i] * d_inf[i]);
        obj_buf[iter] = obj;

        /* Membership update */
        for (int i = 0; i < n; i++) {
            int nz = 0;
            for (int k = 0; k < nc; k++)
                if (dist2[k*n + i] == 0.0) zero_idx[nz++] = k;

            if (nz > 0) {
                for (int k = 0; k < nc; k++) U[k*n + i] = 0.0;
                double eq = 1.0 / nz;
                for (int z = 0; z < nz; z++) U[zero_idx[z]*n + i] = eq;
                U_inf[i] = 0.0;
            } else {
                double sum_w = 0.0;
                for (int k = 0; k < nc; k++) {
                    double w = pow(dist2[k*n + i], -exp_inv);
                    Um[k*n + i] = w;
                    sum_w += w;
                }
                double w_inf = pow(d_inf[i], -2.0 * exp_inv);
                double total = sum_w + w_inf;
                for (int k = 0; k < nc; k++)
                    U[k*n + i] = Um[k*n + i] / total;
                U_inf[i] = w_inf / total;
            }
        }

        /* Convergence */
        double max_diff = 0.0;
        for (int k = 0; k < nc; k++)
            for (int i = 0; i < n; i++) {
                double diff = fabs(U[k*n + i] - U_old[k*n + i]);
                if (diff > max_diff) max_diff = diff;
            }
        for (int i = 0; i < n; i++) {
            double diff = fabs(U_inf[i] - U_inf_old[i]);
            if (diff > max_diff) max_diff = diff;
        }
        if (max_diff < tol) break;
    }

    /* Build output: U is n×(nc+1), last column = U_inf */
    SEXP centers_sexp = PROTECT(allocMatrix(REALSXP, nc, d));
    double *cp = REAL(centers_sexp);
    for (int j = 0; j < d; j++)
        for (int k = 0; k < nc; k++)
            cp[j*nc + k] = centers[k*d + j];

    SEXP U_sexp = PROTECT(allocMatrix(REALSXP, n, nc + 1));
    double *up = REAL(U_sexp);
    memcpy(up, U, (size_t)n * nc * sizeof(double));
    memcpy(up + n * nc, U_inf, (size_t)n * sizeof(double));

    SEXP obj_sexp = PROTECT(allocVector(REALSXP, iter_count));
    memcpy(REAL(obj_sexp), obj_buf, (size_t)iter_count * sizeof(double));

    R_Free(U); R_Free(U_old); R_Free(Um); R_Free(U_inf);
    R_Free(U_inf_old); R_Free(centers); R_Free(dist2);
    R_Free(row_sums); R_Free(min_dists); R_Free(d_inf);
    R_Free(obj_buf); R_Free(zero_idx);

    const char *names[] = {"centers", "U", "obj"};
    SEXP vals[] = {centers_sexp, U_sexp, obj_sexp};
    SEXP result = PROTECT(make_named_list(names, vals, 3));
    UNPROTECT(4);
    return result;
}
