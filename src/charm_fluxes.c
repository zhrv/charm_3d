//
// Created by zhrv on 26.10.17.
//

#include "charm_fluxes.h"


void rim_orig(  double* RI, double* EI, double* PI, double* UI, double* VI, double* WI,
                double RB, double PB, double UB, double VB, double WB,
                double RE, double PE, double UE, double VE, double WE, double gam) {
    int    step;
    double AGAM = (gam - 1.0);
    double BGAM = (2.0 * sqrt(gam / AGAM));
    double CGAM = (1.0 / gam);
    double DGAM = (2.0 / AGAM);
    double EGAM = (AGAM / (gam + 1.0));
    double GGAM = (sqrt(gam * AGAM));
    double HGAM = (AGAM / 2.0);
    double FGAM = (3.0 * gam - 1.0);
    double OGAM = (AGAM / (2.0 * gam));
    double QGAM = (gam + 1.0);
    double PGAM = (QGAM / (2.0 * gam));
    double RGAM = (4.0 * gam);
    double SGAM = (gam * AGAM);
    double TGAM = (QGAM / 2.0);
    double UGAM = (sqrt(AGAM / gam));

    double RF, RS, EF, ES, SBL, SFL, SSL, SEL, D, FS1, F1, ZNB, PKB, ZFB, PKE, ZNE, F2, FS2, ZFE, DP, UBD, RUBD, UF, UED, RUED, US, PPE, PPB, P;
    double eps = CHARM_RIM_EPS;
    double CB = sqrt(gam * PB / RB);
    double CE = sqrt(gam * PE / RE);
    double EB = CB * CB / SGAM;
    double EE = CE * CE / SGAM;
    double RCB = RB * CB;
    double RCE = RE * CE;
    double DU = UB - UE;
    if (DU < -2.0 * (CB + CE) / AGAM) {
        CHARM_GLOBAL_ESSENTIALF ("%s\n", " RIEMANN PROBLEM SOLVER: ATTENTION!!!  VACUUM!!!");
        RF = 0.0;
        RS = 0.0;
        EF = 0.0;
        ES = 0.0;
        SBL = UB - CB;
        SFL = UB + 2.0 * CB / AGAM;
        SSL = UE - 2.0 * CE / AGAM;
        SEL = UE + CE;
    } else {
        P = (PB * RCE + PE * RCB + DU * RCB * RCE) / (RCB + RCE);
        step = 0;
        do {

            if (P < eps) P = eps;

            PPB = P / PB;
            if (PB <= P) {
                PKB = PGAM * PPB + OGAM;
                ZNB = RCB * sqrt(PKB);
                F1 = (P - PB) / ZNB;
                FS1 = (QGAM * PPB + FGAM) / (RGAM * ZNB * PKB);
            }
            else {
                ZFB = CB * exp(log(PPB) * OGAM);
                F1 = DGAM * (ZFB - CB);
                FS1 = ZFB / (gam * P);
            }
            PPE = P / PE;
            if (PE <= P) {
                PKE = PGAM * PPE + OGAM;
                ZNE = RCE * sqrt(PKE);
                F2 = (P - PE) / ZNE;
                FS2 = (QGAM * PPE + FGAM) / (RGAM * ZNE * PKE);
            }
            else {
                ZFE = CE * exp(log(PPE) * OGAM);
                F2 = DGAM * (ZFE - CE);
                FS2 = ZFE / (gam * P);
            }
            DP = (DU - F1 - F2) / (FS1 + FS2);
            P = P + DP;
        } while ((fabs(DU - F1 - F2) > eps) && (++step < CHARM_RIM_NEWTON_STEPS));


        PPB = P / PB;
        PPE = P / PE;

//       ZFB=CB*PPB**OGAM;
//       ZFE=CE*PPE**OGAM;
        ZFB = CB * exp(log(PPB) * OGAM);
        ZFE = CE * exp(log(PPE) * OGAM);
        if (PB <= P) {
            D = UB - sqrt((TGAM * P + HGAM * PB) / RB);
            UBD = UB - D;
            RUBD = RB * UBD;
            RF = RUBD * RUBD / (PB - P + RUBD * UBD);
            UF = D + RUBD / RF;
            EF = P / (AGAM * RF);
            SBL = D;
            SFL = D;
        }
        else {
            EF = ZFB * ZFB / SGAM;
            UF = UB + DGAM * (CB - ZFB);
            RF = P / (AGAM * EF);
            SBL = UB - CB;
            SFL = UF - ZFB;
        }
        if (PE <= P) {
            D = UE + sqrt((TGAM * P + HGAM * PE) / RE);
            UED = UE - D;
            RUED = RE * UED;
            RS = RUED * RUED / (PE - P + RUED * UED);
            US = D + RUED / RS;
            ES = P / (AGAM * RS);
            SEL = D;
            SSL = D;
        }
        else {
            ES = ZFE * ZFE / SGAM;
            US = UE - DGAM * (CE - ZFE);
            RS = P / (AGAM * ES);
            SSL = US + ZFE;
            SEL = UE + CE;
        }
    }
//
// C     compute the interpolation value
    if (SEL<=0.0) {
        *RI= RE;
        *EI= EE;
        *UI= UE;
        *VI= VE;
        *WI= WE;
        *PI= AGAM*(*EI)*(*RI);
        return;
    }

    if (SBL>=0.0) {
        *RI= RB;
        *EI= EB;
        *UI= UB;
        *VI= VB;
        *WI= WB;
        *PI= AGAM*(*EI)*(*RI);
        return;
    }

    if ((SSL>=0.0)&&(SFL<=0.0)) {
        if (US>=0.0) {
            *RI= RF;
            *EI= EF;
            *UI= UF;
            *VI= VB;
            *WI= WB;
        } else {
            *RI= RS;
            *EI= ES;
            *UI= US;
            *VI= VE;
            *WI= WE;
        }
        *PI= AGAM*(*EI)*(*RI);
        return;
    }

    if (SFL>0.0) {
        *UI= (UB+DGAM*GGAM*sqrt(EB))/(1+DGAM);
        *VI= VB;
        *WI= WB;
        *EI= ((*UI)*(*UI))/SGAM;
        *RI= RB*exp(log((*EI)/EB)*(1/AGAM));
    } else {
        *UI= (UE-DGAM*GGAM*sqrt(EE))/(1+DGAM);
        *VI= VE;
        *WI= WE;
        *EI= ((*UI)*(*UI))/SGAM;
        *RI= RE*exp(log((*EI)/EE)*(1/AGAM)) ;
    }

    *PI= AGAM*(*EI)*(*RI);
}

void charm_calc_flux_godunov(p4est_t *p4est, charm_prim_t prim[2], double* qu, double* qv, double* qw, double* qe, double* qc, double n[3])
{
//    int i,j;
//    double ri, ei, pi, uu[3], uv[3];
//    double nt[3][3], vv[2][3], vn[2][3];
//    double r_[2], u_[2], v_[2], w_[2], p_[2];
//
//    for (i = 0; i < 2; i++) {
//        r_[i] = prim[i].r;
//        u_[i] = prim[i].u;
//        v_[i] = prim[i].v;
//        w_[i] = prim[i].w;
//        p_[i] = prim[i].p;
//    }
//
//    nt[0][0] = n[0];
//    nt[0][1] = n[1];
//    nt[0][2] = n[2];
//
//
//    ri = sqrt(n[0]*n[0]+n[1]*n[1]);
//    if (ri > CHARM_EPS) {
//        nt[1][0] = -n[1]/ri;
//        nt[1][1] =  n[0]/ri;
//        nt[1][2] =  0.;
//    }
//    else {
//        ri = sqrt(n[1]*n[1]+n[2]*n[2]);
//        if (ri > CHARM_EPS) {
//            nt[1][0] =  0.;
//            nt[1][1] = -n[2]/ri;
//            nt[1][2] =  n[1]/ri;
//        }
//        else {
//            ri = sqrt(n[0]*n[0]+n[2]*n[2]);
//            if (ri > CHARM_EPS) {
//                nt[1][0] = -n[2]/ri;
//                nt[1][1] =  0.;
//                nt[1][2] =  n[0]/ri;
//            }
//            else {
//                CHARM_ASSERT(1./sqrt(n[1]*n[1]+n[2]*n[2]) > CHARM_EPS);
//            }
//        }
//    }
//    vect_prod(nt[0], nt[1], nt[2]);
//    for (i = 0; i < 2; i++) {
//        vv[i][0] = u_[i];
//        vv[i][1] = v_[i];
//        vv[i][2] = w_[i];
//    }
//    for (i = 0; i < 2; i++) {
//        for (j = 0; j < 3; j++) {
//            vn[i][j] = scalar_prod(vv[i], nt[j]);
//        }
//    }
//    rim_orig(  &ri, &ei, &pi, &(uu[0]), &(uu[1]), &(uu[2]),
//               r_[0], p_[0], vn[0][0], vn[0][1], vn[0][2],
//               r_[1], p_[1], vn[1][0], vn[1][1], vn[1][2], GAM);
//
//    for (i = 0; i < 3; i++) {
//        uv[i] = 0.;
//        for (j = 0; j < 3; j++) {
//            uv[i] += nt[j][i]*uu[j];
//        }
//    }
//    *qr = ri*uu[0];
//    *qu = (*qr)*uv[0]+pi*n[0];
//    *qv = (*qr)*uv[1]+pi*n[1];
//    *qw = (*qr)*uv[2]+pi*n[2];
//    *qe = (ri*(ei+0.5*(uv[0]*uv[0]+uv[1]*uv[1]+uv[2]*uv[2]))+pi)*uu[0];
}

void charm_calc_flux_lf(p4est_t *p4est, charm_prim_t prim[2], double* qu, double* qv, double* qw, double* qe, double qc[], double n[3])
{
    int      i, j;
    double   alpha;
    double   vn[2], **ff/*[5][2]*/, **uu/*[5][2]*/, fc;
    double **q;
    size_t   f_count;
    double   rs[2], rs_, rr, ur, vr, wr, vrn;

    charm_ctx_t *ctx = charm_get_ctx(p4est);
    f_count = 4+ctx->comp->elem_count;

    q  = CHARM_ALLOC(double*, f_count);
    ff = CHARM_ALLOC(double*, f_count);
    uu = CHARM_ALLOC(double*, f_count);
    for (i = 0; i < f_count; i++) {
        ff[i] = CHARM_ALLOC(double, 2);
        uu[i] = CHARM_ALLOC(double, 2);
    }
    q[0] = qu;
    q[1] = qv;
    q[2] = qw;
    q[3] = qe;
    for (i = 4; i < f_count; i++) {
        q[i] = &(qc[i-4]);
    }

    for (i = 0; i < 2; i++) {
        vn[i] = prim[i].u*n[0]+prim[i].v*n[1]+prim[i].w*n[2];

        uu[0][i] = prim[i].r*prim[i].u;
        uu[1][i] = prim[i].r*prim[i].v;
        uu[2][i] = prim[i].r*prim[i].w;
        uu[3][i] = prim[i].r*prim[i].e_tot;

        ff[0][i] = prim[i].r*vn[i]*prim[i].u+prim[i].p*n[0];
        ff[1][i] = prim[i].r*vn[i]*prim[i].v+prim[i].p*n[1];
        ff[2][i] = prim[i].r*vn[i]*prim[i].w+prim[i].p*n[2];
        ff[3][i] = (prim[i].r*prim[i].e_tot+prim[i].p)*vn[i];

        for (j = 4; j < f_count; j++) {
            uu[j][i] = prim[i].r*prim[i].c[j-4];
            ff[j][i] = prim[i].r*vn[i]*prim[i].c[j-4];
        }

    }

    alpha = _MAX_(fabs(vn[0])+prim[0].cz, fabs(vn[1])+prim[1].cz);

    for (i = 0; i < f_count; i++) {
        *(q[i]) = 0.5*( ff[i][1]+ff[i][0]-alpha*(uu[i][1]-uu[i][0]) );
    }

    for (i = 0; i < f_count; i++) {
        CHARM_FREE(ff[i]);
        CHARM_FREE(uu[i]);
    }
    CHARM_FREE(q);
    CHARM_FREE(uu);
    CHARM_FREE(ff);
}



#define F_HLLC_U(UK, FK, SK, SS, PK, RK, VK) (((SS)*((SK)*(UK)-(FK)) + (SK)*( (PK)+(RK)*((SK)-(VK))*((SS)-(VK)) )) / ((SK)-(SS)))
#define F_HLLC_V(UK, FK, SK, SS, PK, RK, VK) (((SS)*((SK)*(UK)-(FK))) / ((SK)-(SS)))
#define F_HLLC_E(UK, FK, SK, SS, PK, RK, VK) (((SS)*((SK)*(UK)-(FK)) + (SK)*( (PK)+(RK)*((SK)-(VK))*((SS)-(VK)) )*(SS)) / ((SK)-(SS)))


static void _charm_calc_flux_hllc_x_1(p4est_t *p4est, charm_prim_t prim[2], double* qu, double* qv, double* qw, double* qe, double qc[])
{
    int             i;
    size_t          c_count = charm_get_comp_count(p4est);
    double          sl, sr, p_star, s_star, p_pvrs, ql, qr, tmp;

    p_pvrs = 0.5*(prim[0].p0+prim[1].p0)-0.5*(prim[1].u-prim[0].u)*0.25*(prim[0].r+prim[1].r)*(prim[0].cz+prim[1].cz);
    p_star = (p_pvrs > 0.) ? p_pvrs : 0.;

    ql = (p_star <= prim[0].p0) ? 1 : sqrt(1.+(prim[0].gam+1.)*(p_star/prim[0].p0-1.)/(2.*prim[0].gam));
    qr = (p_star <= prim[1].p0) ? 1 : sqrt(1.+(prim[1].gam+1.)*(p_star/prim[1].p0-1.)/(2.*prim[1].gam));

    sl = prim[0].u-prim[0].cz*ql;
    sr = prim[1].u+prim[1].cz*qr;

    if (sl>sr) {
        tmp = sl;
        sl = sr;
        sr = tmp;
    }

    s_star  = prim[1].p0-prim[0].p0;
    s_star += prim[0].r*prim[0].u*(sl-prim[0].u);
    s_star -= prim[1].r*prim[1].u*(sr-prim[1].u);
    s_star /= (prim[0].r*(sl-prim[0].u)-prim[1].r*(sr-prim[1].u));

    if (s_star < sl) s_star = sl;
    if (s_star > sr) s_star = sr;


    if (!((sl <= s_star) && (s_star <= sr))) {
        CHARM_LERROR("HLLC: inequaluty SL <= S* <= SR is FALSE.\n");
        charm_abort(p4est, 1);
    }

    if (sl >= 0.) {
        *qu = prim[0].r*prim[0].u*prim[0].u+prim[0].p0;
        *qv = prim[0].r*prim[0].v*prim[0].u;
        *qw = prim[0].r*prim[0].w*prim[0].u;
        *qe = (prim[0].r*prim[0].e_tot+prim[0].p0)*prim[0].u;
        for (i = 0; i < c_count; i++) {
            qc[i] = prim[0].r*prim[0].u*prim[0].c[i];
        }
    }
    else if (sr <= 0.) {
        *qu = prim[1].r*prim[1].u*prim[1].u+prim[1].p0;
        *qv = prim[1].r*prim[1].v*prim[1].u;
        *qw = prim[1].r*prim[1].w*prim[1].u;
        *qe = (prim[1].r*prim[1].e_tot+prim[1].p0)*prim[1].u;
        for (i = 0; i < c_count; i++) {
            qc[i] = prim[1].r*prim[1].u*prim[1].c[i];
        }
    }
    else {
        if (s_star >= 0) {
            *qu = F_HLLC_U( /*  UK, FK, SK, SS, PK, RK, VK */
                    prim[0].r * prim[0].u,
                    prim[0].r * prim[0].u * prim[0].u + prim[0].p0,
                    sl, s_star, prim[0].p0, prim[0].r, prim[0].u
            );
            *qv = F_HLLC_V( /*  UK, FK, SK, SS, PK, RK, VK */
                    prim[0].r * prim[0].v,
                    prim[0].r * prim[0].u * prim[0].v,
                    sl, s_star, prim[0].p0, prim[0].r, prim[0].u
            );
            *qw = F_HLLC_V( /*  UK, FK, SK, SS, PK, RK, VK */
                    prim[0].r * prim[0].w,
                    prim[0].r * prim[0].u * prim[0].w,
                    sl, s_star, prim[0].p0, prim[0].r, prim[0].u
            );
            *qe = F_HLLC_E( /*  UK, FK, SK, SS, PK, RK, VK */
                    prim[0].r * prim[0].e_tot,
                    (prim[0].r * prim[0].e_tot + prim[0].p0)*prim[0].u,
                    sl, s_star, prim[0].p0, prim[0].r, prim[0].u
            );
            for (i = 0; i < c_count; i++) {
                qc[i] = F_HLLC_V( /*  UK, FK, SK, SS, PK, RK, VK */
                        prim[0].r * prim[0].c[i],
                        prim[0].r * prim[0].c[i] * prim[0].u,
                        sl, s_star, prim[0].p0, prim[0].r, prim[0].u
                );
            }
        }
        else {
            *qu = F_HLLC_U( /*  UK, FK, SK, SS, PK, RK, VK */
                    prim[1].r * prim[1].u,
                    prim[1].r * prim[1].u * prim[1].u + prim[1].p0,
                    sr, s_star, prim[1].p0, prim[1].r, prim[1].u
            );
            *qv = F_HLLC_V( /*  UK, FK, SK, SS, PK, RK, VK */
                    prim[1].r * prim[1].v,
                    prim[1].r * prim[1].u * prim[1].v,
                    sr, s_star, prim[1].p0, prim[1].r, prim[1].u
            );
            *qw = F_HLLC_V( /*  UK, FK, SK, SS, PK, RK, VK */
                    prim[1].r * prim[1].w,
                    prim[1].r * prim[1].u * prim[1].w,
                    sr, s_star, prim[1].p0, prim[1].r, prim[1].u
            );
            *qe = F_HLLC_E( /*  UK, FK, SK, SS, PK, RK, VK */
                    prim[1].r * prim[1].e_tot,
                    (prim[1].r * prim[1].e_tot + prim[1].p0)*prim[1].u,
                    sr, s_star, prim[1].p0, prim[1].r, prim[1].u
            );
            for (i = 0; i < c_count; i++) {
                qc[i] = F_HLLC_V( /*  UK, FK, SK, SS, PK, RK, VK */
                        prim[1].r * prim[1].c[i],
                        prim[1].r * prim[1].c[i] * prim[1].u,
                        sr, s_star, prim[1].p0, prim[1].r, prim[1].u
                );
            }
        }
    }
}

static void _charm_calc_flux_hllc_x_1_low_mach(p4est_t *p4est, charm_prim_t prim[2], double* qu, double* qv, double* qw, double* qe, double qc[])
{
    int             i;
    size_t          c_count = charm_get_comp_count(p4est);
    double          sl, sr, p_star, s_star, p_pvrs, ql, qr, tmp;

    p_pvrs = 0.5*(prim[0].p0+prim[1].p0)-0.5*(prim[1].u-prim[0].u)*0.25*(prim[0].r+prim[1].r)*(prim[0].cz+prim[1].cz);
    p_star = (p_pvrs > 0.) ? p_pvrs : 0.;

    ql = (p_star <= prim[0].p0) ? 1 : sqrt(1.+(prim[0].gam+1.)*(p_star/prim[0].p0-1.)/(2.*prim[0].gam));
    qr = (p_star <= prim[1].p0) ? 1 : sqrt(1.+(prim[1].gam+1.)*(p_star/prim[1].p0-1.)/(2.*prim[1].gam));

    sl = prim[0].u-prim[0].cz*ql;
    sr = prim[1].u+prim[1].cz*qr;

    if (sl>sr) {
        tmp = sl;
        sl = sr;
        sr = tmp;
    }

    s_star  = prim[1].p0-prim[0].p0;
    s_star += prim[0].r*prim[0].u*(sl-prim[0].u);
    s_star -= prim[1].r*prim[1].u*(sr-prim[1].u);
    s_star /= (prim[0].r*(sl-prim[0].u)-prim[1].r*(sr-prim[1].u));

    if (s_star < sl) s_star = sl;
    if (s_star > sr) s_star = sr;


    if (!((sl <= s_star) && (s_star <= sr))) {
        CHARM_LERROR("HLLC: inequality SL <= S* <= SR is FALSE.\n");
        charm_abort(p4est, 1);
    }

    if (sl >= 0.) {
        *qu = prim[0].r*prim[0].u*prim[0].u+prim[0].p0;
        *qv = prim[0].r*prim[0].v*prim[0].u;
        *qw = prim[0].r*prim[0].w*prim[0].u;
        *qe = (prim[0].r*prim[0].e_tot+prim[0].p0)*prim[0].u;
        for (i = 0; i < c_count; i++) {
            qc[i] = prim[0].r*prim[0].u*prim[0].c[i];
        }
    }
    else if (sr <= 0.) {
        *qu = prim[1].r*prim[1].u*prim[1].u+prim[1].p0;
        *qv = prim[1].r*prim[1].v*prim[1].u;
        *qw = prim[1].r*prim[1].w*prim[1].u;
        *qe = (prim[1].r*prim[1].e_tot+prim[1].p0)*prim[1].u;
        for (i = 0; i < c_count; i++) {
            qc[i] = prim[1].r*prim[1].u*prim[1].c[i];
        }
    }
    else {
        if (s_star >= 0) {
            *qu = F_HLLC_U( /*  UK, FK, SK, SS, PK, RK, VK */
                    prim[0].r * prim[0].u,
                    prim[0].r * prim[0].u * prim[0].u + prim[0].p0,
                    sl, s_star, prim[0].p0, prim[0].r, prim[0].u
            );
            *qv = F_HLLC_V( /*  UK, FK, SK, SS, PK, RK, VK */
                    prim[0].r * prim[0].v,
                    prim[0].r * prim[0].u * prim[0].v,
                    sl, s_star, prim[0].p0, prim[0].r, prim[0].u
            );
            *qw = F_HLLC_V( /*  UK, FK, SK, SS, PK, RK, VK */
                    prim[0].r * prim[0].w,
                    prim[0].r * prim[0].u * prim[0].w,
                    sl, s_star, prim[0].p0, prim[0].r, prim[0].u
            );
            *qe = F_HLLC_E( /*  UK, FK, SK, SS, PK, RK, VK */
                    prim[0].r * prim[0].e_tot,
                    (prim[0].r * prim[0].e_tot + prim[0].p0)*prim[0].u,
                    sl, s_star, prim[0].p0, prim[0].r, prim[0].u
            );
            for (i = 0; i < c_count; i++) {
                qc[i] = F_HLLC_V( /*  UK, FK, SK, SS, PK, RK, VK */
                        prim[0].r * prim[0].c[i],
                        prim[0].r * prim[0].c[i] * prim[0].u,
                        sl, s_star, prim[0].p0, prim[0].r, prim[0].u
                );
            }
        }
        else {
            *qu = F_HLLC_U( /*  UK, FK, SK, SS, PK, RK, VK */
                    prim[1].r * prim[1].u,
                    prim[1].r * prim[1].u * prim[1].u + prim[1].p0,
                    sr, s_star, prim[1].p0, prim[1].r, prim[1].u
            );
            *qv = F_HLLC_V( /*  UK, FK, SK, SS, PK, RK, VK */
                    prim[1].r * prim[1].v,
                    prim[1].r * prim[1].u * prim[1].v,
                    sr, s_star, prim[1].p0, prim[1].r, prim[1].u
            );
            *qw = F_HLLC_V( /*  UK, FK, SK, SS, PK, RK, VK */
                    prim[1].r * prim[1].w,
                    prim[1].r * prim[1].u * prim[1].w,
                    sr, s_star, prim[1].p0, prim[1].r, prim[1].u
            );
            *qe = F_HLLC_E( /*  UK, FK, SK, SS, PK, RK, VK */
                    prim[1].r * prim[1].e_tot,
                    (prim[1].r * prim[1].e_tot + prim[1].p0)*prim[1].u,
                    sr, s_star, prim[1].p0, prim[1].r, prim[1].u
            );
            for (i = 0; i < c_count; i++) {
                qc[i] = F_HLLC_V( /*  UK, FK, SK, SS, PK, RK, VK */
                        prim[1].r * prim[1].c[i],
                        prim[1].r * prim[1].c[i] * prim[1].u,
                        sr, s_star, prim[1].p0, prim[1].r, prim[1].u
                );
            }
        }
    }
}

void charm_calc_flux_hllc(p4est_t *p4est, charm_prim_t prim[2], double* qu, double* qv, double* qw, double* qe, double qc[], double n[3])
{
    charm_ctx_t *ctx = charm_get_ctx(p4est);
    int i,j;
    double ri, ei, pi, uu[3], uv[3];
    double nt[3][3], vv[2][3], vn[2][3];
    double /*r_[2],*/ u_[2], v_[2], w_[2]/*, p_[2]*/;
    size_t          c_count = charm_get_comp_count(p4est);
    double _qu, _qv, _qw;

    charm_prim_t prim_[2];

    for (i = 0; i < 2; i++) {
//        r_[i] = prim[i].r;
        u_[i] = prim[i].u;
        v_[i] = prim[i].v;
        w_[i] = prim[i].w;
//        p_[i] = prim[i].p;
    }

    nt[0][0] = n[0];
    nt[0][1] = n[1];
    nt[0][2] = n[2];


    ri = sqrt(n[0]*n[0]+n[1]*n[1]);
    if (ri > CHARM_EPS) {
        nt[1][0] = -n[1]/ri;
        nt[1][1] =  n[0]/ri;
        nt[1][2] =  0.;
    }
    else {
        ri = sqrt(n[1]*n[1]+n[2]*n[2]);
        if (ri > CHARM_EPS) {
            nt[1][0] =  0.;
            nt[1][1] = -n[2]/ri;
            nt[1][2] =  n[1]/ri;
        }
        else {
            ri = sqrt(n[0]*n[0]+n[2]*n[2]);
            if (ri > CHARM_EPS) {
                nt[1][0] = -n[2]/ri;
                nt[1][1] =  0.;
                nt[1][2] =  n[0]/ri;
            }
            else {
                CHARM_ASSERT(1./sqrt(n[1]*n[1]+n[2]*n[2]) > CHARM_EPS);
            }
        }
    }
    vect_prod(nt[0], nt[1], nt[2]);
    for (i = 0; i < 2; i++) {
        vv[i][0] = u_[i];
        vv[i][1] = v_[i];
        vv[i][2] = w_[i];
    }
    for (i = 0; i < 2; i++) {
        for (j = 0; j < 3; j++) {
            vn[i][j] = scalar_prod(vv[i], nt[j]);
        }
        memcpy(&(prim_[i]), &(prim[i]), sizeof(charm_prim_t));
        prim_[i].u = vn[i][0];
        prim_[i].v = vn[i][1];
        prim_[i].w = vn[i][2];
    }

    if (ctx->model == CHARM_MODEL_NS_LOW_MACH) {
        _charm_calc_flux_hllc_x_1_low_mach(p4est, prim_, &_qu, &_qv, &_qw, qe, qc);
    }
    else {
        _charm_calc_flux_hllc_x_1(p4est, prim_, &_qu, &_qv, &_qw, qe, qc);
    }

    *qu = _qu*nt[0][0] + _qv*nt[1][0] + _qw*nt[2][0];
    *qv = _qu*nt[0][1] + _qv*nt[1][1] + _qw*nt[2][1];
    *qw = _qu*nt[0][2] + _qv*nt[1][2] + _qw*nt[2][2];
}


void charm_calc_flux_cd(p4est_t *p4est, charm_prim_t prim[2], double* qu, double* qv, double* qw, double* qe, double qc[], double n[3])
{
    int      i, j;
    double   alpha;
    double   vn[2], **ff/*[5][2]*/;
    double **q;
    size_t   f_count;

    charm_ctx_t *ctx = charm_get_ctx(p4est);
    f_count = 4+ctx->comp->elem_count;

    q  = CHARM_ALLOC(double*, f_count);
    ff = CHARM_ALLOC(double*, f_count);
    for (i = 0; i < f_count; i++) {
        ff[i] = CHARM_ALLOC(double, 2);
    }
    q[0] = qu;
    q[1] = qv;
    q[2] = qw;
    q[3] = qe;
    for (i = 4; i < f_count; i++) {
        q[i] = &(qc[i-4]);
    }

    for (i = 0; i < 2; i++) {
        vn[i] = prim[i].u*n[0]+prim[i].v*n[1]+prim[i].w*n[2];

        ff[0][i] = prim[i].r*vn[i]*prim[i].u+prim[i].p*n[0];
        ff[1][i] = prim[i].r*vn[i]*prim[i].v+prim[i].p*n[1];
        ff[2][i] = prim[i].r*vn[i]*prim[i].w+prim[i].p*n[2];
        ff[3][i] = (prim[i].r*prim[i].e_tot+prim[i].p)*vn[i];

        for (j = 4; j < f_count; j++) {
            ff[j][i] = prim[i].r*prim[i].c[j-4]*vn[i];
        }
    }

    for (i = 0; i < f_count; i++) {
        *(q[i]) = 0.5*( ff[i][1]+ff[i][0] );
    }

    for (i = 0; i < f_count; i++) {
        CHARM_FREE(ff[i]);
    }
    CHARM_FREE(q);
    CHARM_FREE(ff);
}

