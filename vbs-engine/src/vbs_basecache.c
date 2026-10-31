#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rtklib.h"
#include "vbs/vbs_basecache.h"

/* ---- a minimal JSON reader ------------------------------------------------------------- */

enum jtype { J_NULL, J_BOOL, J_NUM, J_STR, J_ARR, J_OBJ };

struct jnode {
	enum jtype type;
	double num;
	int boolean;
	char *str;			/* J_STR; also the key of an object member */
	struct jnode *child;		/* first member / element */
	struct jnode *next;		/* next sibling */
	char *key;			/* member name when inside an object */
};

/* far deeper than any cache or coords document; bounds the recursion of jvalue() */
#define JMAXDEPTH	64

struct jparser {
	const char *p;
	int err;
	int depth;
	char msg[96];
};

static struct jnode *jnew(enum jtype t) {
	struct jnode *n = (struct jnode *)calloc(1, sizeof *n);

	if (n)
		n->type = t;
	return n;
}

static void jfree(struct jnode *n) {
	while (n) {
		struct jnode *nx = n->next;

		jfree(n->child);
		free(n->str);
		free(n->key);
		free(n);
		n = nx;
	}
}

static void jskip(struct jparser *ps) {
	while (*ps->p && isspace((unsigned char)*ps->p))
		ps->p++;
}

static void jfail(struct jparser *ps, const char *what) {
	if (!ps->err) {
		ps->err = 1;
		snprintf(ps->msg, sizeof ps->msg, "%s near \"%.20s\"", what, ps->p);
	}
}

static char *jstring(struct jparser *ps) {
	char *out, *w;
	const char *s;
	size_t cap;

	if (*ps->p != '"') {
		jfail(ps, "expected string");
		return NULL;
	}
	ps->p++;
	s = ps->p;
	cap = 16;
	out = (char *)malloc(cap);
	w = out;
	while (*ps->p && *ps->p != '"') {
		char c = *ps->p++;

		if (c == '\\') {
			char e;

			/* a backslash right before the end: reading the escape would step past the NUL */
			if (*ps->p == '\0')
				break;
			e = *ps->p++;

			switch (e) {
			case 'n': c = '\n'; break;
			case 't': c = '\t'; break;
			case 'r': c = '\r'; break;
			case 'b': c = '\b'; break;
			case 'f': c = '\f'; break;
			case 'u':
				/* keep it simple: \u00XX only, else '?' */
				if (isxdigit((unsigned char)ps->p[0]) && isxdigit((unsigned char)ps->p[1]) &&
				    isxdigit((unsigned char)ps->p[2]) && isxdigit((unsigned char)ps->p[3])) {
					char hex[5] = {ps->p[0], ps->p[1], ps->p[2], ps->p[3], 0};
					long v = strtol(hex, NULL, 16);

					c = v < 128 ? (char)v : '?';
					ps->p += 4;
				} else
					c = '?';
				break;
			default: c = e; break;
			}
		}
		if ((size_t)(w - out) + 2 > cap) {
			size_t off = w - out;

			cap *= 2;
			out = (char *)realloc(out, cap);
			w = out + off;
		}
		*w++ = c;
	}
	if (*ps->p != '"') {
		jfail(ps, "unterminated string");
		free(out);
		return NULL;
	}
	ps->p++;
	*w = '\0';
	(void)s;
	return out;
}

static struct jnode *jvalue(struct jparser *ps);

static struct jnode *jarray(struct jparser *ps) {
	struct jnode *a = jnew(J_ARR), **tail = &a->child;

	ps->p++;			/* [ */
	jskip(ps);
	if (*ps->p == ']') {
		ps->p++;
		return a;
	}
	for (;;) {
		struct jnode *v = jvalue(ps);

		if (v == NULL) {
			jfree(a);
			return NULL;
		}
		*tail = v;
		tail = &v->next;
		jskip(ps);
		if (*ps->p == ',') {
			ps->p++;
			jskip(ps);
			continue;
		}
		if (*ps->p == ']') {
			ps->p++;
			return a;
		}
		jfail(ps, "expected , or ]");
		jfree(a);
		return NULL;
	}
}

static struct jnode *jobject(struct jparser *ps) {
	struct jnode *o = jnew(J_OBJ), **tail = &o->child;

	ps->p++;			/* { */
	jskip(ps);
	if (*ps->p == '}') {
		ps->p++;
		return o;
	}
	for (;;) {
		char *key;
		struct jnode *v;

		jskip(ps);
		if ((key = jstring(ps)) == NULL) {
			jfree(o);
			return NULL;
		}
		jskip(ps);
		if (*ps->p != ':') {
			jfail(ps, "expected :");
			free(key);
			jfree(o);
			return NULL;
		}
		ps->p++;
		if ((v = jvalue(ps)) == NULL) {
			free(key);
			jfree(o);
			return NULL;
		}
		v->key = key;
		*tail = v;
		tail = &v->next;
		jskip(ps);
		if (*ps->p == ',') {
			ps->p++;
			continue;
		}
		if (*ps->p == '}') {
			ps->p++;
			return o;
		}
		jfail(ps, "expected , or }");
		jfree(o);
		return NULL;
	}
}

static struct jnode *jvalue(struct jparser *ps) {
	struct jnode *n;

	jskip(ps);
	switch (*ps->p) {
	case '{':
	case '[':
		if (++ps->depth > JMAXDEPTH) {
			jfail(ps, "nesting too deep");
			n = NULL;
		} else
			n = *ps->p == '{' ? jobject(ps) : jarray(ps);
		ps->depth--;
		return n;
	case '"':
		n = jnew(J_STR);
		if ((n->str = jstring(ps)) == NULL) {
			free(n);
			return NULL;
		}
		return n;
	case 't':
		if (!strncmp(ps->p, "true", 4)) {
			ps->p += 4;
			n = jnew(J_BOOL);
			n->boolean = 1;
			return n;
		}
		break;
	case 'f':
		if (!strncmp(ps->p, "false", 5)) {
			ps->p += 5;
			n = jnew(J_BOOL);
			return n;
		}
		break;
	case 'n':
		if (!strncmp(ps->p, "null", 4)) {
			ps->p += 4;
			return jnew(J_NULL);
		}
		break;
	default:
		if (*ps->p == '-' || isdigit((unsigned char)*ps->p)) {
			char *end;

			n = jnew(J_NUM);
			n->num = strtod(ps->p, &end);
			if (end == ps->p) {
				jfail(ps, "bad number");
				free(n);
				return NULL;
			}
			ps->p = end;
			return n;
		}
	}
	jfail(ps, "unexpected token");
	return NULL;
}

static const struct jnode *jget(const struct jnode *obj, const char *key) {
	const struct jnode *c;

	if (obj == NULL || obj->type != J_OBJ)
		return NULL;
	for (c = obj->child; c; c = c->next)
		if (c->key && !strcmp(c->key, key))
			return c;
	return NULL;
}

static void jcopy_str(const struct jnode *obj, const char *key, char *dst, size_t len) {
	const struct jnode *v = jget(obj, key);

	dst[0] = '\0';
	if (v && v->type == J_STR)
		snprintf(dst, len, "%s", v->str);
}

/* ---- the cache ------------------------------------------------------------------------- */

static int bands_from_str(const char *s) {
	int seen[10] = {0}, n = 0;

	for (; s && *s; s++)
		if ((*s == '1' || *s == '2' || *s == '5') && !seen[*s - '0']) {
			seen[*s - '0'] = 1;
			n++;
		}
	return n;
}

/* Read a whole file into a NUL-terminated buffer, or NULL with a message. */
static char *slurp(const char *path, const char *who) {
	FILE *fp;
	char *buf;
	long len;

	if ((fp = fopen(path, "rb")) == NULL) {
		fprintf(stderr, "%s: cannot open %s\n", who, path);
		return NULL;
	}
	fseek(fp, 0, SEEK_END);
	len = ftell(fp);
	fseek(fp, 0, SEEK_SET);
	if (len <= 0 || len > 64L * 1024 * 1024) {
		fclose(fp);
		fprintf(stderr, "%s: %s has an implausible size (%ld)\n", who, path, len);
		return NULL;
	}
	if ((buf = (char *)malloc((size_t)len + 1)) == NULL) {
		fclose(fp);
		return NULL;
	}
	if (fread(buf, 1, (size_t)len, fp) != (size_t)len) {
		fclose(fp);
		free(buf);
		fprintf(stderr, "%s: short read on %s\n", who, path);
		return NULL;
	}
	fclose(fp);
	buf[len] = '\0';
	return buf;
}

static double xyz_dist(const double *a, const double *b) {
	double d0 = a[0] - b[0], d1 = a[1] - b[1], d2 = a[2] - b[2];

	return sqrt(d0 * d0 + d1 * d1 + d2 * d2);
}

static int jnum(const struct jnode *obj, const char *key, double *out) {
	const struct jnode *v = jget(obj, key);

	if (v == NULL || v->type != J_NUM)
		return 0;
	*out = v->num;
	return 1;
}

int vbs_basecache_apply_coords(struct vbs_basecache *bc, const char *path, double max_shift_m,
	int *n_refused, double *worst_applied_m) {
	struct jparser ps;
	struct jnode *root;
	const struct jnode *src, *m;
	char *buf;
	int applied = 0, refused_far = 0, refused_quality = 0, absent = 0;
	double worst = 0.0;

	if (bc == NULL || path == NULL || max_shift_m <= 0.0)
		return -1;
	if ((buf = slurp(path, "coords")) == NULL)
		return -1;
	memset(&ps, 0, sizeof ps);
	ps.p = buf;
	root = jvalue(&ps);
	free(buf);			/* every jnode owns its own strings; see jstring() */
	if (root == NULL || ps.err) {
		fprintf(stderr, "coords: %s: %s\n", path, ps.err ? ps.msg : "empty");
		jfree(root);
		return -1;
	}
	if ((src = jget(root, "sources")) == NULL || src->type != J_OBJ) {
		fprintf(stderr, "coords: %s has no \"sources\" object\n", path);
		jfree(root);
		return -1;
	}
	for (m = src->child; m; m = m->next) {
		struct vbs_basecache_entry *e;
		char quality[VBS_BC_MAXSTR];
		double nxyz[3], shift;
		int idx;

		if (m->key == NULL || m->type != J_OBJ)
			continue;
		if ((idx = vbs_basecache_find(bc, m->key)) < 0) {
			absent++;
			continue;
		}
		if (!jnum(m, "x_m", &nxyz[0]) || !jnum(m, "y_m", &nxyz[1]) ||
		    !jnum(m, "z_m", &nxyz[2]))
			continue;
		jcopy_str(m, "calculation_quality", quality, sizeof quality);
		if (!strcmp(quality, "EXPLORATORY_GLOBAL_BASELINE")) {
			refused_quality++;
			continue;
		}
		e = &bc->e[idx];
		shift = xyz_dist(nxyz, e->xyz);
		/*
		 * The bound is HORIZONTAL. It exists to catch a mountpoint mismatch -- a document
		 * entry keyed to the wrong base -- and a wrong base is somewhere else on the map. A
		 * vertical-only difference is the registry's height being wrong, and the document is
		 * the truth (coordinate rule, 2026-10-02). On the 18 September document, 16 of the 17
		 * entries a 3D bound refused were the same station to 0.15-1.5 m horizontally (one at
		 * 8 m) and differed only in height, by the local geoid: Switzerland +47 to +50 m,
		 * Hungary +44 m -- orthometric heights stored where ellipsoidal ones belong. Only
		 * LAKI, 161 km away, was the case the bound is for.
		 */
		{
			double rel[3] = {nxyz[0] - e->xyz[0], nxyz[1] - e->xyz[1], nxyz[2] - e->xyz[2]};
			double enu[3], horiz;

			ecef2enu(e->llh, rel, enu);
			horiz = sqrt(enu[0] * enu[0] + enu[1] * enu[1]);
			if (horiz > max_shift_m) {
				fprintf(stderr, "coords: %s REFUSED, %.1f m from the registry horizontally "
				    "(bound %.1f m) -- a different place, not a better coordinate\n",
				    e->mp, horiz, max_shift_m);
				refused_far++;
				continue;
			}
			if (fabs(enu[2]) > max_shift_m)
				fprintf(stderr, "coords: %s height corrected by %+.2f m (horizontal %.2f m): "
				    "the registry height is wrong\n", e->mp, enu[2], horiz);
		}
		e->xyz[0] = nxyz[0];
		e->xyz[1] = nxyz[1];
		e->xyz[2] = nxyz[2];
		ecef2pos(e->xyz, e->llh);
		jcopy_str(m, "reference_frame", e->frame, sizeof e->frame);
		/*
		 * Refusals were loud and applications were silent, which is the wrong way round for
		 * the interesting ones. A centimetre correction is the normal case and not worth a
		 * line; a metre means the registry was simply wrong about where that base is, and an
		 * operator should see it go by rather than discover it later.
		 */
		if (shift > VBS_BC_COORD_LOUD_M)
			fprintf(stderr, "coords: %s moved %.3f m -- the registry was wrong about this "
			    "base, not merely stale\n", e->mp, shift);
		e->coord_from_ref = 1;
		if (shift > worst)
			worst = shift;
		applied++;
	}
	jfree(root);
	fprintf(stderr, "coords: %d of %d bases repositioned from %s (worst %.3f m); "
	    "%d refused too far, %d refused on quality, %d not in the cache\n",
	    applied, bc->n, path, worst, refused_far, refused_quality, absent);
	if (n_refused)
		*n_refused = refused_far + refused_quality;
	if (worst_applied_m)
		*worst_applied_m = worst;
	return applied;
}

int vbs_basecache_load(const char *path, struct vbs_basecache *out) {
	FILE *fp;
	char *buf;
	long len;
	struct jparser ps;
	struct jnode *root;
	const struct jnode *arr, *b;
	int n = 0, i;

	memset(out, 0, sizeof *out);
	if ((fp = fopen(path, "rb")) == NULL) {
		fprintf(stderr, "basecache: cannot open %s\n", path);
		return -1;
	}
	fseek(fp, 0, SEEK_END);
	len = ftell(fp);
	fseek(fp, 0, SEEK_SET);
	if (len <= 0 || len > 64L * 1024 * 1024) {
		fclose(fp);
		fprintf(stderr, "basecache: %s has an implausible size (%ld)\n", path, len);
		return -1;
	}
	buf = (char *)malloc((size_t)len + 1);
	if (fread(buf, 1, (size_t)len, fp) != (size_t)len) {
		fclose(fp);
		free(buf);
		fprintf(stderr, "basecache: short read on %s\n", path);
		return -1;
	}
	fclose(fp);
	buf[len] = '\0';

	memset(&ps, 0, sizeof ps);
	ps.p = buf;
	root = jvalue(&ps);
	if (root == NULL || ps.err) {
		fprintf(stderr, "basecache: %s: %s\n", path, ps.err ? ps.msg : "empty");
		jfree(root);
		free(buf);
		return -1;
	}
	arr = jget(root, "bases");
	if (arr == NULL || arr->type != J_ARR) {
		fprintf(stderr, "basecache: %s: no \"bases\" array\n", path);
		jfree(root);
		free(buf);
		return -1;
	}
	jcopy_str(root, "generated", out->generated, sizeof out->generated);
	for (b = arr->child; b; b = b->next)
		n++;
	out->e = (struct vbs_basecache_entry *)calloc(n > 0 ? n : 1, sizeof *out->e);
	for (b = arr->child, i = 0; b; b = b->next) {
		struct vbs_basecache_entry *e = &out->e[i];
		const struct jnode *xyz = jget(b, "xyz"), *v;
		int k;

		if (b->type != J_OBJ)
			continue;
		jcopy_str(b, "mp", e->mp, sizeof e->mp);
		jcopy_str(b, "marker", e->marker, sizeof e->marker);
		jcopy_str(b, "rec", e->rec, sizeof e->rec);
		jcopy_str(b, "ant", e->ant, sizeof e->ant);
		jcopy_str(b, "frame", e->frame, sizeof e->frame);
		/* "fw" is what the caster writes; "firmware" accepted so a hand-made cache works. */
		jcopy_str(b, "fw", e->fw, sizeof e->fw);
		if (e->fw[0] == '\0')
			jcopy_str(b, "firmware", e->fw, sizeof e->fw);
		if (e->mp[0] == '\0' || xyz == NULL || xyz->type != J_ARR)
			continue;
		for (k = 0, v = xyz->child; k < 3 && v; k++, v = v->next)
			e->xyz[k] = v->type == J_NUM ? v->num : 0.0;
		if (k < 3 || (fabs(e->xyz[0]) < 1e5 && fabs(e->xyz[1]) < 1e5))
			continue;
		ecef2pos(e->xyz, e->llh);
		if ((v = jget(b, "gps_bands")) != NULL && v->type == J_STR)
			e->gps_bands = bands_from_str(v->str);
		e->family = vbs_policy_family(e->rec);
		e->live = 1;
		if ((v = jget(b, "live")) != NULL && v->type == J_BOOL)
			e->live = v->boolean;
		e->stream_ok = 0;
		i++;
	}
	out->n = i;
	if (out->n == 0) {
		/* callers treat <= 0 as failure and free nothing */
		free(out->e);
		out->e = NULL;
	}
	jfree(root);
	free(buf);
	return out->n;
}

void vbs_basecache_free(struct vbs_basecache *bc) {
	free(bc->e);
	bc->e = NULL;
	bc->n = 0;
}

int vbs_basecache_find(const struct vbs_basecache *bc, const char *mp) {
	int i;

	for (i = 0; i < bc->n; i++)
		if (!strcmp(bc->e[i].mp, mp))
			return i;
	return -1;
}

int vbs_basecache_merge(struct vbs_basecache *bc, const struct vbs_basecache *nb,
	int *n_updated, int *n_added, int *n_gone) {
	int i, j, nnew = 0, updated = 0, added = 0, gone = 0;
	struct vbs_basecache_entry *grown;

	for (j = 0; j < nb->n; j++)
		if (vbs_basecache_find(bc, nb->e[j].mp) < 0)
			nnew++;
	if (nnew > 0) {
		grown = (struct vbs_basecache_entry *)realloc(bc->e, (bc->n + nnew) * sizeof *bc->e);
		if (grown == NULL)
			return -1;
		bc->e = grown;
	}
	for (i = 0; i < bc->n; i++)
		bc->e[i].live = 0;
	for (j = 0; j < nb->n; j++) {
		const struct vbs_basecache_entry *s = &nb->e[j];
		int k = vbs_basecache_find(bc, s->mp);

		if (k >= 0) {
			struct vbs_basecache_entry *e = &bc->e[k];

			if (memcmp(e->xyz, s->xyz, sizeof e->xyz) != 0 || strcmp(e->rec, s->rec) != 0 ||
			    strcmp(e->frame, s->frame) != 0 || strcmp(e->fw, s->fw) != 0)
				updated++;
			memcpy(e->xyz, s->xyz, sizeof e->xyz);
			memcpy(e->llh, s->llh, sizeof e->llh);
			memcpy(e->rec, s->rec, sizeof e->rec);
			memcpy(e->ant, s->ant, sizeof e->ant);
			memcpy(e->frame, s->frame, sizeof e->frame);
			memcpy(e->marker, s->marker, sizeof e->marker);
			memcpy(e->fw, s->fw, sizeof e->fw);
			e->gps_bands = s->gps_bands;
			e->family = s->family;
			e->live = s->live;
			e->coord_from_ref = s->coord_from_ref;
			e->untrusted = s->untrusted;
		} else {
			struct vbs_basecache_entry *ne = &bc->e[bc->n++];

			*ne = *s;
			ne->stream_ok = 0;
			ne->last_epoch_unix = 0.0;
			added++;
		}
	}
	for (i = 0; i < bc->n; i++)
		if (!bc->e[i].live && vbs_basecache_find(nb, bc->e[i].mp) < 0)
			gone++;
	memcpy(bc->generated, nb->generated, sizeof bc->generated);
	if (n_updated)
		*n_updated = updated;
	if (n_added)
		*n_added = added;
	if (n_gone)
		*n_gone = gone;
	return 0;
}

int vbs_basecache_require_ref(struct vbs_basecache *bc) {
	int i, n = 0;

	for (i = 0; i < bc->n; i++)
		if (!bc->e[i].coord_from_ref) {
			bc->e[i].untrusted = 1;
			bc->e[i].live = 0;
			n++;
		}
	return n;
}

int vbs_basecache_policy_view(const struct vbs_basecache *bc, int use_stream_state,
    struct vbs_policy_base *pb) {
	int i;

	for (i = 0; i < bc->n; i++) {
		pb[i].mp = bc->e[i].mp;
		pb[i].family = bc->e[i].family;
		memcpy(pb[i].xyz, bc->e[i].xyz, sizeof pb[i].xyz);
		pb[i].gps_bands = bc->e[i].gps_bands;
		/* The live string, not a verdict: an upgraded station clears on the next cache read. */
		pb[i].fw = bc->e[i].fw;
		pb[i].live = !bc->e[i].untrusted &&
		    (use_stream_state ? bc->e[i].stream_ok : bc->e[i].live);
	}
	return bc->n;
}

int vbs_basecache_height_at(const struct vbs_basecache *bc, double lat_deg, double lon_deg,
    double radius_km, double *h_out) {
	double pos[3] = {lat_deg * D2R, lon_deg * D2R, 0.0}, xyz[3], wsum = 0.0, hsum = 0.0;
	int i, n = 0;

	pos2ecef(pos, xyz);
	for (i = 0; i < bc->n; i++) {
		double d = sqrt((xyz[0] - bc->e[i].xyz[0]) * (xyz[0] - bc->e[i].xyz[0]) +
		    (xyz[1] - bc->e[i].xyz[1]) * (xyz[1] - bc->e[i].xyz[1]) +
		    (xyz[2] - bc->e[i].xyz[2]) * (xyz[2] - bc->e[i].xyz[2])) / 1000.0;
		double w;

		if (d > radius_km || bc->e[i].untrusted)
			continue;
		/* a base almost under the point should dominate; cap the weight rather than divide by ~0 */
		w = 1.0 / ((d < 0.5 ? 0.5 : d) * (d < 0.5 ? 0.5 : d));
		wsum += w;
		hsum += w * bc->e[i].llh[2];
		n++;
	}
	if (n == 0)
		return -1;
	*h_out = hsum / wsum;
	return 0;
}
