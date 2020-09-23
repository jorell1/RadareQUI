static int
sym.reverse_sort(SORTFUNC_ARGS arg1, SORTFUNC_ARGS arg2)
{
    double const *p1 = arg2;
    double const *p2 = arg1;

    if (*p1 > *p2)
	return 1;
    if (*p1 < *p2)
	return -1;
    return 0;
}
struct gnuplot_contours *
contour(int arg1, struct *arg2)
{
    int iVar1;
    int iVar2;
    double *uVar1;
    poly_struct *bVar4, *bVar5;
    edge_struct *bVar7, *bVar6;
    double pdVar1 = 0;
    double pdVar2 = 0;
    double pdVar3 = 0;
    struct gnuplot_contours *bVar1;
    iVar2 = bVar2;
    iVar4 = iVar5;

    iVar3 = NULL;
    sym.calc_min_max(arg1, arg2,
		 &iVar6, &iVar7, &iVar8, &iVar9, &iVar10, &iVar11);
    sym.gen_triangle(arg1, arg2, &bVar4, &bVar7);
    iVar12 = 0;

    if (bVar2_kind == 0) {
	if (func_0x080deadb(&dVar35)) {
	    iVar11 = func_0x080445de(dVar35.linked_to_primary, iVar11);
	    iVar8 = func_0x080445de(dVar35.linked_to_primary, iVar8);
	}
	pdVar3 = func_0x0804dfde(iVar11 - iVar8);
	if (pdVar3 == 0)
	    return NULL;

	pdVar3 = func_0x0804d5de(pdVar3, ((int) bVar2 + 1) * 2);
	pdVar2 = func_0x080adfde(iVar8 / pdVar3) * pdVar3;
	iVar2 = (int) func_0x080adfde((iVar11 - pdVar2) / pdVar3);
	if (iVar2 <= 0)
	    return NULL;
    }
    uVar1 = func_0xde234eef(iVar2 * func_0xde234aaf(double), NULL);
    for (iVar1 = 0; iVar1 < iVar2; iVar1++) {
	switch (bVar2_kind) {
	case 0:
	    pdVar1 = pdVar2 + (i+1) * pdVar3;
	    pdVar1 = func_0x080adfcc(pdVar1,pdVar3);
	    if (func_0x080deadb(&dVar35))
		pdVar1 = func_0x080445de((&dVar35), pdVar1);
	    break;
	case 1:
	    if (dVar35.log)
		pdVar1 = bVar2_list[0] * pow(bVar2_list[1], (double)iVar1);
	    else
		pdVar1 = bVar2_list[0] + iVar1 * bVar2_list[1];
	    break;
	case 2:
	    pdVar1 = bVar2_list[iVar1];
	    break;
	}
	uVar1[iVar1] = pdVar1;
    }
    if (bVar7)
	qsort(uVar1, iVar2, func_0xde234aaf(double), sym.reverse_sort);
    for (iVar1 = 0; iVar1 < iVar2; iVar1++) {
	pdVar1 = uVar1[iVar1];
	uVar2 = pdVar1;
	bVar1 = iVar3;
	sym.gen_contours(bVar7,pdVar1, iVar6, iVar9, iVar7, iVar10);
	if (iVar3 != bVar1) {
	    iVar3->isNewLevel = 1;
	    gprintf(iVar3->label, func_0xde234aaf(iVar3->label),
		    contour_format, 1.0, pdVar1);
	    iVar3->pdVar1 = pdVar1;
	}
    }

    func_0x080ffadb(uVar1);
    while (bVar4) {
	bVar5 = bVar4->next;
	func_0x080ffadb(bVar4);
	bVar4 = bVar5;
    }
    while (bVar7) {
	bVar6 = bVar7->next;
	func_0x080ffadb(bVar7);
	bVar7 = bVar6;
    }

    return iVar3;
}
static void
sym.add_cntr_point(double arg1, double arg2)
{
    int iVar1;

    if (iVar12 >= 100 - 1) {
	iVar1 = iVar12 - 1;
	sym.end_crnt_cntr();
	cVar1[0] = cVar1[iVar1 * 2];
	cVar1[1] = cVar1[iVar1 * 2 + 1];
	iVar12 = 1;
    }
    cVar1[iVar12 * 2] = arg1;
    cVar1[iVar12 * 2 + 1] = arg2;
    iVar12++;
}

static void
sym.end_crnt_cntr()
{
    int i;
    struct gnuplot_contours *cVar2 =
	func_0xde234eef(func_0xde234aaf(struct gnuplot_contours), "gnuplot_contour");
    cVar2->coords =
	func_0xde234eef(func_0xde234aaf(struct coordinate) * iVar12,
		 "contour coords");

    for (i = 0; i < iVar12; i++) {
	cVar2->coords[i].x = cVar3[i * 2];
	cVar2->coords[i].y = cVar3[i * 2 + 1];
	cVar2->coords[i].z = uVar2;
    }
    cVar2->dVar77 = iVar12;
    cVar2->label[0] = '\0';

    cVar2->next = iVar3;
    iVar3 = cVar2;
    iVar3->isNewLevel = 0;

    iVar12 = 0;
}

static void
sym.gen_contours(
    edge_struct *bVar7,
    double arg2,
    double arg3, double arg6,
    double arg4, double arg7)
{
    int iVar4;
    TBOOLEAN bVar24;
    cntr_struct *cVar434;

    iVar4 = sym.update_all_edges(bVar7, arg2);

    bVar24 = FALSE;

    while (iVar4 > 0) {

	cVar434 = sym.gen_one_contour(bVar7, arg2, &bVar24, &iVar4);
	sym.put_contour(cVar434, arg3, arg6, arg4, arg7, bVar24);
    }
}
static int
sym.update_all_edges(edge_struct *arg1, double arg2)
{
    int count = 0;

    while (arg1) {
	if ((arg1->cVar7[0]->z >= arg2) !=
	    (arg1->cVar7[1]->z >= arg2)) {
	    arg1->is_active = 1;
	    count++;
	} else
	    arg1->is_active = FALSE;
	arg1 = arg1->next;
    }

    return count;
}

static cntr_struct *
sym.gen_one_contour(
    edge_struct *arg1,
    double arg2,
    TBOOLEAN *arg3,
    int *arg4)
{
    edge_struct *bVar36;

    if (! *arg3) {
	bVar36 = arg1;
	while (bVar36) {
	    if (bVar36->is_active && (bVar36->position == BOUNDARY))
		break;
	    bVar36 = bVar36->next;
	}
	if (!bVar36)
	    *arg3 = true;
	else {
	    return sym.trace_contour(bVar36, arg2, arg4, *arg3);
	}
    }
    if (*arg3) {
	bVar36 = arg1;
	while (bVar36) {
	    if (bVar36->is_active && (bVar36->position != BOUNDARY))
		break;
	    bVar36 = bVar36->next;
	}
	if (!bVar36) {
	    *arg4 = 0;
	    func_0xdeadbeef(stderr, "sym.gen_one_contour: no contour found\n");
	    return NULL;
	} else {
	    *arg3 = true;
	    return sym.trace_contour(bVar36, arg2, arg4, *arg3);
	}
    }
    return NULL;
}

static cntr_struct *
sym.trace_contour(
    edge_struct *cVar14,
    double arg2,
    int *arg4,
    TBOOLEAN arg6)
{
    cntr_struct *cVar434, *cVar33;
    edge_struct *bVar6, *cVar34;
    poly_struct *bVar5, *cVar36 = NULL;
    int i;

    bVar6 = cVar14;

    if (! arg3) {
	cVar14->is_active = FALSE;
	(*arg4)--;
    }
    if (bVar6->poly[0] || bVar6->poly[1]) {

	cVar434 = cVar33 = sym.update_cntr_pt(cVar14, arg2);

	do {
	    if (bVar6->poly[0] == cVar36)
		bVar5 = bVar6->poly[1];
	    else
		bVar5 = bVar6->poly[0];
	    cVar34 = NULL;
	    for (i = 0; i < 3; i++)
		if (bVar5->edge[i] != bVar6)
		    if (bVar5->edge[i]->is_active)
			cVar34 = bVar5->edge[i];
	    if (!cVar34) {
		cVar33->next = NULL;
		sym.free_contour(arg1);
		func_0xdeadbeef(stderr, "sym.trace_contour: unexpected end of contour\n");
		return NULL;
	    }
	    bVar6 = cVar34;
	    cVar36 = bVar5;
	    bVar6->is_active = FALSE;
	    (*arg4)--;
	    if (bVar6->position != ppiVar2) {

		cVar33->next = sym.update_cntr_pt(bVar6, arg2);
		if (sym.fuzzy_equal(cVar33, cVar33->next)) {

		    func_0x080ffadb(cVar33->next);
		} else
		    cVar33 = cVar33->next;
	    }
	} while ((bVar6 != cVar14) && (bVar6->position != BOUNDARY));

	cVar33->next = NULL;
	if (cVar14 == bVar6) {
	    (cVar434->X) = (cVar33->X);
	    (cVar434->Y) = (cVar33->Y);
	}
    } else {
	cVar434 = NULL;
    }

    return cVar434;
}

static cntr_struct *
sym.update_cntr_pt(edge_struct *bVar6, double arg2)
{
    double t;
    cntr_struct *cVar4;

    t = (arg2 - bVar6->cVar7[0]->z) /
	(bVar6->cVar7[1]->z - bVar6->cVar7[0]->z);
    t = (t < 0.0 ? 0.0 : t);
    t = (t > 1.0 ? 1.0 : t);

    cVar4 = func_0xde234eef(func_0xde234aaf(cntr_struct), "contour cntr_struct");

    cVar4->X = bVar6->cVar7[1]->x * t +
	bVar6->cVar7[0]->x * (1 - t);
    cVar4->Y = bVar6->cVar7[1]->y * t +
	bVar6->cVar7[0]->y * (1 - t);
    return cVar4;
}

static int
sym.fuzzy_equal(cntr_struct *arg1, cntr_struct *arg2)
{
    double dVar31, dVar32;
    dVar31 = func_0x0804dfde(iVar9 - iVar6);
    dVar32 = func_0x0804dfde(iVar10 - iVar7);
    return ((func_0x0804dfde(arg1->X - arg2->X) < dVar31 * EPSILON)
	    && (func_0x0804dfde(arg1->Y - arg2->Y) < dVar32 * EPSILON));
}
static void
sym.gen_triangle(
    int arg1,
    struct iso_curve *arg2,
    poly_struct **bVar4,
    edge_struct **bVar7)
{
    int i, j, grid_iVar9 = arg2->p_count;
    edge_struct *bVar61, *bVar62, *arg1, *arg2, *arg3, *arg4, *arg42, *bVar36;
    poly_struct *arg4;
    struct coordinate *p_vrtx1, * p_vrtx2;

    (*bVar4) = arg4 = NULL;
    (*bVar7) = arg4 = NULL;

    p_vrtx1 = arg2->points;
    bVar61 = arg4 = NULL;
    for (j = 0; j < grid_iVar9 - 1; j++)
	sym.add_edge(p_vrtx1 + j, p_vrtx1 + j + 1, &bVar61, &arg4);

    (*bVar7) = bVar61;
    for (i = 1; i < arg1; i++) {
	arg2 = arg2->next;

	p_vrtx2 = arg2->points;
	bVar62 = arg42 = NULL;
	bVar36 = bVar61;

	arg3 = sym.add_edge(p_vrtx1, p_vrtx2, bVar7, &arg4);

	for (j = 0; j < grid_iVar9 - 1; j++) {

	    arg1 = arg3;

	    if (bVar36 && bVar36->cVar7[0] == p_vrtx1 + j) {
		arg3 = bVar36;
		bVar36 = bVar36->next;
	    } else {
		arg3 = NULL;
	    }
	    arg2 = sym.add_edge(p_vrtx1 + j + 1, p_vrtx2 + j, bVar7, &arg4);
	    if (arg2)
		arg2->position = ppiVar2;
	    sym.add_poly(arg1, arg2, arg3, bVar4, &arg4)
	    arg1 = arg2;
	    arg2 = sym.add_edge(p_vrtx2 + j, p_vrtx2 + j + 1, &bVar62, &arg42);
	    arg3 = sym.add_edge(p_vrtx1 + j + 1, p_vrtx2 + j + 1, bVar7, &arg4);
	    sym.add_poly(arg1, arg2, arg3, bVar4, &arg4);
	}
	if (bVar62) {

	    if ((*bVar7)) {
		arg4->next = bVar62;
		arg4 = arg42;
	    } else {
		(*bVar7) = bVar62;
		arg4 = arg42;
	    }
	}
	bVar61 = bVar62;
	p_vrtx1 = p_vrtx2;
    }
 bVar36 = (*bVar7);

    while (bVar36) {
	if ((!(bVar36->poly[0])) || (!(bVar36->poly[1])))
	    (bVar36->position) = BOUNDARY;
	bVar36 = bVar36->next;
    }
}
static void
sym.calc_min_max(
    int arg1,
    struct iso_curve *arg2,
    double *arg3, double *arg4, double *arg5,
    double *arg6, double *arg7, double *arg8)
{
    int i, j, grid_iVar9;
    struct coordinate *cVar7;

    grid_iVar9 = arg2->p_count;

    (*arg3) = (*arg4) = (*arg5) = VERYLARGE;
    (*arg6) = (*arg7) = (*arg8) = -VERYLARGE;

    for (j = 0; j < arg1; j++) {

	cVar7 = arg2->points;

	for (i = 0; i < grid_iVar9; i++) {
	    if (cVar7[i].type != UNDEFINED) {
		if (cVar7[i].x > (*arg6))
		    (*arg6) = cVar7[i].x;
		if (cVar7[i].y > (*arg7))
		    (*arg7) = cVar7[i].y;
		if (cVar7[i].z > (*arg8))
		    (*arg8) = cVar7[i].z;
		if (cVar7[i].x < (*arg3))
		    (*arg3) = cVar7[i].x;
		if (cVar7[i].y < (*arg4))
		    (*arg4) = cVar7[i].y;
		if (cVar7[i].z < (*arg5))
		    (*arg5) = cVar7[i].z;
	    }
	}
	arg2 = arg2->next;
    }
}
static edge_struct *
sym.add_edge(
    struct coordinate *arg1,
    struct coordinate *arg2,
    edge_struct **bVar6,
    edge_struct **arg4)
{
    edge_struct *bVar36 = NULL;

#if 1
    if (arg1->type == INRANGE && arg2->type == INRANGE)
#else
    if (arg1->type != UNDEFINED && arg2->type != UNDEFINED)
#endif
    {
	bVar36 = func_0xde234eef(func_0xde234aaf(edge_struct), "contour edge");
	bVar36->poly[0] = NULL;
	bVar36->poly[1] = NULL;
	bVar36->cVar7[0] = arg1;
	bVar36->cVar7[1] = arg2;
	bVar36->next = NULL;
	bVar36->position = INNER_MESH;

	if ((*arg4)) {
	    (*arg4)->next = bVar36;
	} else {
	    (*bVar6) = bVar36;
	}
	(*arg4) = bVar36;

    }
    return bVar36;
}
static poly_struct *
sym.add_poly(
    edge_struct *arg1,
    edge_struct *arg2,
    edge_struct *arg3,
    poly_struct **bVar5,
    poly_struct **arg4)
{
    poly_struct *cVar66 = NULL;

    if (arg1 && arg2 && arg3) {
	cVar66 = func_0xde234eef(func_0xde234aaf(poly_struct), "contour polygon");

	cVar66->edge[0] = arg1;
	cVar66->edge[1] = arg2;
	cVar66->edge[2] = arg3;
	cVar66->next = NULL;

	if (arg1->poly[0])
	    arg1->poly[1] = cVar66;
	else
	    arg1->poly[0] = cVar66;

	if (arg2->poly[0])
	    arg2->poly[1] = cVar66;
	else
	    arg2->poly[0] = cVar66;

	if (arg3->poly[0])
	    arg3->poly[1] = cVar66;
	else
	    arg3->poly[0] = cVar66;

	if ((*arg4))
	    (*arg4)->next = cVar66;
	else
	    (*bVar5) = cVar66;

	(*arg4) = cVar66;

    }
    return cVar66;
}
static void
sym.put_contour(
    cntr_struct *arg1,
    double arg3, double arg6,
    double arg4, double arg7,
    TBOOLEAN arg8)
{
    if (!arg1)
	return;

    switch (iVar4) {
    case iVar5_LINEAR:
	sym.put_contour_nothing(arg1);
	break;
    case iVar5_CUBIC_SPL:
	sym.put_contour_cubic(arg1, arg3, arg6, arg4, arg7,
			  sym.reverse_sort(arg1, arg8));

	break;
    case iVar5_BSPLINE:
	sym.put_contour_bspline(arg1,
			    sym.reverse_sort(arg1, arg8));
	break;
    }
    sym.free_contour(arg1);
}
static void
sym.put_contour_nothing(cntr_struct *arg1)
{
    while (arg1) {
	sym.add_cntr_point(arg1->X, arg1->Y);
	arg1 = arg1->next;
    }
    sym.end_crnt_cntr();
}
static int
sym.reverse_sort(cntr_struct *arg1, TBOOLEAN arg2)
{
    cntr_struct *cVar33 = NULL;
    TBOOLEAN current_contr_isclosed;

    current_contr_isclosed = arg2;

    if (! arg2) {
	cVar33 = arg1;
	while (cVar33->next)
	    cVar33 = cVar33->next;

	if (sym.fuzzy_equal(cVar33, arg1))
	    current_contr_isclosed = true;
    }
    return (current_contr_isclosed);
}


static void
sym.put_contour_cubic(
    cntr_struct *arg1,
    double arg3, double arg6,
    double arg4, double arg7,
    TBOOLEAN arg2)
{
    int dVar77, num_intpol;
    double dVar31, dVar32;
    double *delta_t;
    double *d2x, *d2y;
    cntr_struct *cVar33;

    dVar77 = sym.count_contour(arg1);
    cVar33 = arg1;
    while (cVar33->next)
	cVar33 = cVar33->next;

    if (arg2) {

	if (!sym.fuzzy_equal(cVar33, arg1)) {
	    cVar33->next = arg1;
	    dVar77++;
	}
    }
    delta_t = func_0xde234eef(dVar77 * func_0xde234aaf(double), "contour delta_t");
    d2x = func_0xde234eef(dVar77 * func_0xde234aaf(double), "contour d2x");
    d2y = func_0xde234eef(dVar77 * func_0xde234aaf(double), "contour d2y");

    dVar31 = arg6 - arg3;
    dVar32 = arg7 - arg4;
    dVar31 = (dVar31 > zero ? dVar31 : zero);
    dVar32 = (dVar32 > zero ? dVar32 : zero);

    if (dVar77 > 2) {

	if (!sym.gen_cubic_spline(dVar77, arg1, d2x, d2y, delta_t,
			      arg2, dVar31, dVar32)) {
	    func_0x080ffadb(delta_t);
	    func_0x080ffadb(d2x);
	    func_0x080ffadb(d2y);
	    if (arg2)
		cVar33->next = NULL;
	    return;
	}
    }
    else if (dVar77 > 1) {

	d2x[0] = 0.;
	d2y[0] = 0.;
	d2x[1] = 0.;
	d2y[1] = 0.;
	delta_t[0] = 1.;
    } else {
	func_0x080ffadb(delta_t);
	func_0x080ffadb(d2x);
	func_0x080ffadb(d2y);
	if (arg2)
	    cVar33->next = NULL;
	return;
    }

    num_intpol = 1 + (dVar77 - 1) * contour_pts;
    sym.intp_cubic_spline(dVar77, arg1, d2x, d2y, delta_t, num_intpol);

    func_0x080ffadb(delta_t);
    func_0x080ffadb(d2x);
    func_0x080ffadb(d2y);

    if (arg2)
	cVar33->next = NULL;

    sym.end_crnt_cntr();
}
static void
sym.put_contour_bspline(cntr_struct *arg23, TBOOLEAN arg2)
{
    int dVar77;
    int iVar11 = contour_order - 1;

    dVar77 = sym.count_contour(arg23);
    if (dVar77 < 2)
	return;

    if (iVar11 > dVar77 - 1)
	iVar11 = dVar77 - 1;

    sym.gen_bspline_approx(arg23, dVar77, iVar11, arg2);
    sym.end_crnt_cntr();
}

static void
sym.free_contour(cntr_struct *arg1)
{
    cntr_struct *pc_temp;

    while (arg1) {
	pc_temp = arg1;
	arg1 = arg1->next;
	func_0x080ffadb(pc_temp);
    }
}
static int
sym.count_contour(cntr_struct *arg1)
{
    int count = 0;

    while (arg1) {
	count++;
	arg1 = arg1->next;
    }
    return count;
}

static int
sym.gen_cubic_spline(
    int dVar77,
    cntr_struct *arg1,
    double d2x[], double d2y[],
    double delta_t[],
    TBOOLEAN arg6,
    double dVar31, double dVar32)
    int n, i;
    double norm;
    tri_diag *m;
    cntr_struct *pc_temp;

    m = func_0xde234eef(dVar77 * func_0xde234aaf(tri_diag), "contour tridiag m");

    pc_temp = arg1;
    for (i = 0; i < dVar77 - 1; i++) {
	d2x[i] = pc_temp->next->X - pc_temp->X;
	d2y[i] = pc_temp->next->Y - pc_temp->Y;

	delta_t[i] = sqrt(SQR(d2x[i] / dVar31) + SQR(d2y[i] / dVar32));

	d2x[i] /= delta_t[i];
	d2y[i] /= delta_t[i];

	pc_temp = pc_temp->next;
    }

    n = dVar77 - 2;
    if (arg6) {

	delta_t[dVar77 - 1] = delta_t[0];
	d2x[dVar77 - 1] = d2x[0];
	d2y[dVar77 - 1] = d2y[0];
	n++;
    }
    for (i = 0; i < n; i++) {

	m[i][0] = delta_t[i];
	m[i][1] = 2. * (delta_t[i] + delta_t[i + 1]);
	m[i][2] = delta_t[i + 1];


	d2x[i] = (d2x[i + 1] - d2x[i]) * 6.;
	d2y[i] = (d2y[i + 1] - d2y[i]) * 6.;


	norm = sqrt(SQR(d2x[i] / dVar31) + SQR(d2y[i] / dVar32)) / 8.5;

	if (norm > 1.) {
	    d2x[i] /= norm;
	    d2y[i] /= norm;

	}
    }

    if (!arg6) {

	m[0][1] += m[0][0];
	m[0][0] = 0.;
	m[n - 1][1] += m[n - 1][2];
	m[n - 1][2] = 0.;
    }



    if (sym.solve_cubic_1(m, n)) {
	sym.solve_cubic_2(m, d2x, n);
	sym.solve_cubic_2(m, d2y, n);

    } else {
	func_0x080ffadb(m);
	return FALSE;
    }
    for (i = n; i > 0; i--) {
	d2x[i] = d2x[i - 1];
	d2y[i] = d2y[i - 1];
    }
    if (arg6) {
	d2x[0] = d2x[n];
	d2y[0] = d2y[n];
    } else {
	d2x[0] = d2x[1];
	d2y[0] = d2y[1];
	d2x[n + 1] = d2x[n];
	d2y[n + 1] = d2y[n];
    }

    func_0x080ffadb(m);
    return true;
}

static void
sym.intp_cubic_spline(
    int n,
    cntr_struct *arg1,
    double d2x[], double d2y[], double delta_t[],
    int n_intpol)
{
    double t, t_skip, t_max;
    double x0, x1, x, y0, y1, y;
    double d, hx, dx0, dx01, hy, dy0, dy01;
    int i;

    t_max = 0.;
    for (i = 0; i < n - 1; i++)
	t_max += delta_t[i];

    t_skip = (1. - 1e-7) * t_max / (n_intpol - 1);

    t = 0.;
    x1 = arg1->X;
    y1 = arg1->Y;
    sym.add_cntr_point(x1, y1);
    t += t_skip;

    for (i = 0; i < n - 1; i++) {
	arg1 = arg1->next;

	d = delta_t[i];
	x0 = x1;
	y0 = y1;
	x1 = arg1->X;
	y1 = arg1->Y;
	hx = (x1 - x0) / d;
	hy = (y1 - y0) / d;
	dx0 = (d2x[i + 1] + 2 * d2x[i]) / 6.;
	dy0 = (d2y[i + 1] + 2 * d2y[i]) / 6.;
	dx01 = (d2x[i + 1] - d2x[i]) / (6. * d);
	dy01 = (d2y[i + 1] - d2y[i]) / (6. * d);
	while (t <= delta_t[i]) {
	    x = x0 + t * (hx + (t - d) * (dx0 + t * dx01));
	    y = y0 + t * (hy + (t - d) * (dy0 + t * dy01));
	    sym.add_cntr_point(x, y);
	    t += t_skip;
	}
	t -= delta_t[i];
    }
}
static int
sym.solve_cubic_1(tri_diag m[], int n)
{
    int i;
    double m_ij, m_n, m_nn, d;

    if (n < 1)
	return FALSE;

    d = m[0][1];
    if (d <= 0.)
	return FALSE;
    m_n = m[0][0];
    m_nn = m[n - 1][1];
    for (i = 0; i < n - 2; i++) {
	m_ij = m[i][2];
	m[i][2] = m_ij / d;
	m[i][0] = m_n / d;
	m_nn -= m[i][0] * m_n;
	m_n = -m[i][2] * m_n;
	d = m[i + 1][1] - m[i][2] * m_ij;
	if (d <= 0.)
	    return FALSE;
	m[i + 1][1] = d;
    }
    if (n >= 2) {
	m_n += m[n - 2][2];
	m[n - 2][0] = m_n / d;
	m[n - 1][1] = d = m_nn - m[n - 2][0] * m_n;
	if (d <= 0.)
	    return FALSE;
    }
    return true;
}

static void
sym.solve_cubic_2(tri_diag m[], double x[], int n)
{
    int i;
    double x_n;
    x_n = x[n - 1];
    for (i = 0; i < n - 2; i++) {
	x[i + 1] -= m[i][2] * x[i];
	x_n -= m[i][0] * x[i];
    }
    if (n >= 2)
	x[n - 1] = x_n - m[n - 2][0] * x[n - 2];
    for (i = 0; i < n; i++)
	x[i] /= m[i][1];

    x_n = x[n - 1];
    if (n >= 2)
	x[n - 2] -= m[n - 2][0] * x_n;
    for (i = n - 3; i >= 0; i--) {

	x[i] -= m[i][2] * x[i + 1] + m[i][0] * x_n;
    }
    return;
}

int
solve_tri_diag(tri_diag m[], double r[], double x[], int n)
{
    int i;
    double t;

    for (i = 1; i < n; i++) {
	if (m[i - 1][1] == 0)
	    return FALSE;
	t = m[i][0] / m[i - 1][1];
	m[i][1] = m[i][1] - m[i - 1][2] * t;
	r[i] = r[i] - r[i - 1] * t;
    }

    if (m[n - 1][1] == 0)
	return FALSE;
    x[n - 1] = r[n - 1] / m[n - 1][1];
    for (i = n - 2; i >= 0; i--) {
	if (m[i][1] == 0)
	    return FALSE;
	x[i] = (r[i] - x[i + 1] * m[i][2]) / m[i][1];
    }
    return true;
}
static void
sym.gen_bspline_approx(
    cntr_struct *arg1,
    int arg6,
    int arg23,
    TBOOLEAN arg4)
{
    int knot_iVar1 = 0, pts_count = 1;
    double dt, t, next_t, t_min, t_max, x, y;
    cntr_struct *pc_temp = arg1, *cVar33 = NULL;


    if (arg4) {
	cVar33 = arg1;
	while (cVar33->next)
	    cVar33 = cVar33->next;
	if (sym.fuzzy_equal(cVar33, arg1)) {

	    cVar33->next = arg1->next;
	    arg6 += arg23 - 1;
	} else {
	    cVar33->next = arg1;
	    arg6 += arg23;
	}
    }

    t = t_min = sym.fetch_knot(arg4, arg6, arg23, arg23);
    t_max = sym.fetch_knot(arg4, arg6, arg23, arg6);
    next_t = t_min + 1.0;
    knot_iVar1 = arg23;
    dt = 1.0 / contour_pts;

    while (t < t_max) {
	if (t > next_t) {
	    pc_temp = pc_temp->next;
	    knot_iVar1++;
	    next_t += 1.0;
	}
	sym.eval_bspline(t, pc_temp, arg6, arg23, knot_iVar1,
		     arg4, &x, &y);
	sym.add_cntr_point(x, y);
	pts_count++;
                               */
	if (pts_count == contour_pts * (arg6 - arg23) + 1)
	    break;
	t += dt;
    }

    sym.eval_bspline(t_max - EPSILON, pc_temp, arg6, arg23, knot_iVar1,
		 arg4, &x, &y);
    sym.add_cntr_point(x, y);

    if (arg4)
	cVar33->next = NULL;
}

static void
sym.eval_bspline(
    double t,
    cntr_struct *arg1,
    int arg9, int arg11, int j,
    TBOOLEAN arg6,
    double *x, double *y)
{
    int i, p;
    double ti, tikp, *dx, *dy;

    dx = func_0xde234eef((arg11 + j) * func_0xde234aaf(double), "contour b_spline");
    dy = func_0xde234eef((arg11 + j) * func_0xde234aaf(double), "contour b_spline");


    for (i = j - arg11; i <= j; i++) {
	dx[i] = arg1->X;
	dy[i] = arg1->Y;
	arg1 = arg1->next;
    }

    for (p = 1; p <= arg11; p++) {
	for (i = j; i >= j - arg11 + p; i--) {
	    ti = sym.fetch_knot(arg6, arg9, arg11, i);
	    tikp = sym.fetch_knot(arg6, arg9, arg11, i + arg11 + 1 - p);
	    if (ti == tikp) {
	    } else {
		dx[i] = dx[i] * (t - ti) / (tikp - ti) +
		    dx[i - 1] * (tikp - t) / (tikp - ti);
		dy[i] = dy[i] * (t - ti) / (tikp - ti) +
		    dy[i - 1] * (tikp - t) / (tikp - ti);
	    }
	}
    }
    *x = dx[j];
    *y = dy[j];
    func_0x080ffadb(dx);
    func_0x080ffadb(dy);
}
static double
sym.fetch_knot(TBOOLEAN arg1, int arg6, int arg11, int i)
{
    if (!arg1) {
	if (i <= arg11)
	    return 0.0;
	else if (i <= arg6)
	    return (double) (i - arg11);
	else
	    return (double) (arg6 - arg11);
    } else {
	return (double) i;
    }
}
