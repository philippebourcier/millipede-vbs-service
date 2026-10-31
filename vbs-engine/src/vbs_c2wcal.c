#include <math.h>
#include <stddef.h>
#include <string.h>

#include "vbs/vbs_c2wcal.h"

/* See the header for the measurements behind every choice here. */

int vbs_c2wcal_is_l2_py_code(uint8_t code) {
	return code == CODE_L2W || code == CODE_L2P || code == CODE_L2Y;
}

int vbs_c2wcal_is_l2c_code(uint8_t code) {
	return code == CODE_L2S || code == CODE_L2L || code == CODE_L2X;
}

/* The slots of a GPS record's P(Y) and L2C codes, -1 when absent or empty; 0 for a non-GPS record. */
static int l2_slots(const obsd_t *o, int *jpy, int *jc) {
	int j;

	*jpy = *jc = -1;
	if (satsys(o->sat, NULL) != SYS_GPS)
		return 0;
	for (j = 0; j < NFREQ + NEXOBS; j++) {
		if (o->P[j] == 0.0)
			continue;
		if (*jpy < 0 && vbs_c2wcal_is_l2_py_code(o->code[j]))
			*jpy = j;
		else if (*jc < 0 && vbs_c2wcal_is_l2c_code(o->code[j]))
			*jc = j;
	}
	return 1;
}

/* 2W-2X for one record: 1 and *d when both codes are present and the difference is no blunder. */
static int l2_diff(const obsd_t *o, double *d) {
	int jpy, jc;

	if (!l2_slots(o, &jpy, &jc) || jpy < 0 || jc < 0)
		return 0;
	*d = o->P[jpy] - o->P[jc];
	return fabs(*d) < VBS_C2WCAL_BLUNDER_M;
}

void vbs_c2wcal_probe_reset(struct vbs_c2wcal_probe *p) {
	p->sum_m = 0.0;
	p->n_both = p->n_ident = p->n_py_only = 0;
}

void vbs_c2wcal_probe_add_epoch(struct vbs_c2wcal_probe *p, const obsd_t *obs, int nobs) {
	int i, jpy, jc;
	double d;

	for (i = 0; i < nobs; i++) {
		if (!l2_slots(&obs[i], &jpy, &jc) || jpy < 0)
			continue;
		if (jc < 0) {
			p->n_py_only++;
			continue;
		}
		if (!l2_diff(&obs[i], &d))
			continue;
		p->sum_m += d;
		p->n_both++;
		if (fabs(d) < VBS_C2WCAL_IDENT_M)
			p->n_ident++;
	}
}

int vbs_c2wcal_probe_verdict(const struct vbs_c2wcal_probe *p, double *mean_m) {
	if (mean_m != NULL)
		*mean_m = p->n_both > 0 ? p->sum_m / (double)p->n_both : NAN;
	if (p->n_both < VBS_C2WCAL_MIN_SAMPLES)
		return -1;
	return (double)p->n_ident < VBS_C2WCAL_RELABEL_FRAC * (double)p->n_both;
}

void vbs_c2wcal_groups_reset(struct vbs_c2wcal_groups *g) {
	g->n = 0;
}

struct vbs_c2wcal_probe *vbs_c2wcal_group(struct vbs_c2wcal_groups *g, long build) {
	int k;

	if (build < 0)
		return NULL;
	for (k = 0; k < g->n; k++)
		if (g->build[k] == build)
			return &g->probe[k];
	if (g->n >= VBS_C2WCAL_MAXGROUP)
		return NULL;
	g->build[g->n] = build;
	vbs_c2wcal_probe_reset(&g->probe[g->n]);
	return &g->probe[g->n++];
}

void vbs_c2wcal_est_init(struct vbs_c2wcal_est *e, double window_s) {
	int k;

	e->nbucket = 0;
	if (window_s > 0.0) {
		e->nbucket = (int)ceil(window_s / VBS_C2WCAL_BUCKET_S);
		if (e->nbucket > VBS_C2WCAL_MAXBUCKET)
			e->nbucket = VBS_C2WCAL_MAXBUCKET;
		if (e->nbucket < 1)
			e->nbucket = 1;
	}
	for (k = 0; k < VBS_C2WCAL_MAXBUCKET; k++) {
		e->id[k] = -1;
		e->sum[k] = 0.0;
		e->n[k] = 0;
	}
	e->cum_sum = 0.0;
	e->cum_n = 0;
}

static long minute_of(gtime_t t) {
	return (long)floor(((double)t.time + t.sec) / VBS_C2WCAL_BUCKET_S);
}

void vbs_c2wcal_est_add_epoch(struct vbs_c2wcal_est *e, gtime_t t, const obsd_t *obs, int nobs) {
	long m = minute_of(t);
	int slot = (int)(m % VBS_C2WCAL_MAXBUCKET), i;
	double d;

	if (e->id[slot] != m) {
		e->id[slot] = m;
		e->sum[slot] = 0.0;
		e->n[slot] = 0;
	}
	for (i = 0; i < nobs; i++) {
		if (!l2_diff(&obs[i], &d))
			continue;
		e->sum[slot] += d;
		e->n[slot]++;
		e->cum_sum += d;
		e->cum_n++;
	}
}

int vbs_c2wcal_est_value(const struct vbs_c2wcal_est *e, gtime_t t, double *mean_m, long *n) {
	double s = 0.0;
	long cnt = 0, m = minute_of(t), k;

	if (e->nbucket == 0) {
		s = e->cum_sum;
		cnt = e->cum_n;
	} else {
		for (k = m - e->nbucket + 1; k <= m; k++) {
			int slot = (int)(((k % VBS_C2WCAL_MAXBUCKET) + VBS_C2WCAL_MAXBUCKET) %
			    VBS_C2WCAL_MAXBUCKET);

			if (e->id[slot] == k) {
				s += e->sum[slot];
				cnt += e->n[slot];
			}
		}
	}
	if (n != NULL)
		*n = cnt;
	if (cnt < VBS_C2WCAL_EST_MIN)
		return 0;
	*mean_m = s / (double)cnt;
	return 1;
}

int vbs_c2wcal_apply(obsd_t *obs, int nobs, double bias_m) {
	int i, j, nadj = 0;

	if (bias_m == 0.0)
		return 0;
	for (i = 0; i < nobs; i++) {
		if (satsys(obs[i].sat, NULL) != SYS_GPS)
			continue;
		for (j = 0; j < NFREQ + NEXOBS; j++) {
			if (obs[i].P[j] == 0.0 || !vbs_c2wcal_is_l2_py_code(obs[i].code[j]))
				continue;
			obs[i].P[j] -= bias_m;
			nadj++;
		}
	}
	return nadj;
}

int vbs_c2wcal_apply_l2c(obsd_t *obs, int nobs, double bias_m, int *n_l2c, int *n_py) {
	int i, jpy, jc, a = 0, b = 0;

	for (i = 0; i < nobs; i++) {
		if (!l2_slots(&obs[i], &jpy, &jc) || jpy < 0)
			continue;
		if (jc >= 0) {
			obs[i].P[jpy] = obs[i].P[jc];
			a++;
		} else {
			obs[i].P[jpy] -= bias_m;
			b++;
		}
	}
	if (n_l2c != NULL)
		*n_l2c = a;
	if (n_py != NULL)
		*n_py = b;
	return a + b;
}

void vbs_c2wcal_base_init(struct vbs_c2wcal_base *b, int mode, double window_s) {
	b->mode = mode;
	vbs_c2wcal_probe_reset(&b->probe);
	vbs_c2wcal_est_init(&b->est, window_s);
}

int vbs_c2wcal_base_epoch(struct vbs_c2wcal_base *b, struct vbs_c2wcal_groups *g, long build,
    const obsd_t *obs, int n, obsd_t *out, int nmax, int *verdict, double *offset_m) {
	struct vbs_c2wcal_probe *gp = g != NULL ? vbs_c2wcal_group(g, build) : NULL;
	double off = NAN;
	int v;

	if (offset_m != NULL)
		*offset_m = NAN;
	if (n <= 0) {
		if (verdict != NULL)
			*verdict = -1;
		return 0;
	}
	vbs_c2wcal_probe_add_epoch(&b->probe, obs, n);
	if (gp != NULL)
		vbs_c2wcal_probe_add_epoch(gp, obs, n);
	vbs_c2wcal_est_add_epoch(&b->est, obs[0].time, obs, n);
	v = vbs_c2wcal_probe_verdict(gp != NULL ? gp : &b->probe, NULL);
	if (verdict != NULL)
		*verdict = v;
	if (v != 1 || !vbs_c2wcal_est_value(&b->est, obs[0].time, &off, NULL))
		return 0;
	if (n > nmax)
		n = nmax;
	memcpy(out, obs, (size_t)n * sizeof *obs);
	if (b->mode == VBS_C2WCAL_L2C)
		vbs_c2wcal_apply_l2c(out, n, off, NULL, NULL);
	else
		vbs_c2wcal_apply(out, n, off);
	if (offset_m != NULL)
		*offset_m = off;
	return n;
}
