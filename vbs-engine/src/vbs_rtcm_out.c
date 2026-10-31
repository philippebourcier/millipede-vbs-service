#include <stdlib.h>
#include <string.h>

#include "vbs/vbs_rtcm_out.h"

/*
 * RTCM output for a VBS point. See the header for what this owns and what MRTKLIB owns.
 */

struct vbs_rtcm_out {
	rtcm_t rtcm;		/* MRTKLIB encoder state; large, hence always heap-allocated */
	obsd_t obs[MAXOBS];	/* staging area handed to the encoder, one epoch */
	int glo_fcn_n;		/* channel numbers installed by vbs_rtcm_out_set_glo_fcn() */
	long glo_dropped;	/* GLONASS satellites skipped for want of a channel number */
};

/* Constellations we emit, in message order. */
static const struct {
	int sys;
	int type;
} msm_types[] = {
	{ SYS_GPS, VBS_MSM7_GPS },
	{ SYS_GLO, VBS_MSM7_GLO },
	{ SYS_GAL, VBS_MSM7_GAL },
	{ SYS_CMP, VBS_MSM7_BDS },
};
#define NMSM_TYPES	((int)(sizeof msm_types / sizeof msm_types[0]))

struct vbs_rtcm_out *vbs_rtcm_out_new(int staid, const double *xyz) {
	struct vbs_rtcm_out *this = (struct vbs_rtcm_out *)calloc(1, sizeof *this);

	if (this == NULL)
		return NULL;
	if (!init_rtcm(&this->rtcm)) {
		free(this);
		return NULL;
	}
	this->rtcm.staid = staid;
	memcpy(this->rtcm.sta.pos, xyz, 3 * sizeof(double));

	/*
	 * The encoder reads observations from rtcm.obs. Point it at our own staging array: one
	 * epoch is all we ever hold, and this keeps the ownership obvious at free() time.
	 * init_rtcm() has already allocated a buffer there, which would otherwise leak.
	 */
	free(this->rtcm.obs.data);
	this->rtcm.obs.data = this->obs;
	this->rtcm.obs.nmax = MAXOBS;
	this->rtcm.obs.n = 0;
	return this;
}

void vbs_rtcm_out_free(struct vbs_rtcm_out *this) {
	if (this == NULL)
		return;
	/* obs.data is ours, not the encoder's: detach before free_rtcm() so it cannot free it. */
	this->rtcm.obs.data = NULL;
	this->rtcm.obs.nmax = this->rtcm.obs.n = 0;
	free_rtcm(&this->rtcm);
	free(this);
}

/* Append the encoder's current message to buf. Returns bytes appended, or -1. */
int vbs_rtcm_out_set_glo_fcn(struct vbs_rtcm_out *this, const nav_t *nav) {
	int i, prn, n = 0;

	if (this == NULL || nav == NULL)
		return 0;
	for (i = 0; i < nav->ng; i++) {
		if (satsys(nav->geph[i].sat, &prn) != SYS_GLO)
			continue;
		if (prn < 1 || prn > (int)(sizeof this->rtcm.nav.glo_fcn /
		    sizeof this->rtcm.nav.glo_fcn[0]))
			continue;
		/* The encoder's convention is fcn + 8, with 0 meaning "no data". */
		this->rtcm.nav.glo_fcn[prn - 1] = nav->geph[i].frq + 8;
		n++;
	}
	this->glo_fcn_n = n;
	return n;
}

long vbs_rtcm_out_glo_dropped(const struct vbs_rtcm_out *this) {
	return this == NULL ? 0 : this->glo_dropped;
}

static int append_msg(struct vbs_rtcm_out *this, uint8_t *buf, int buflen, int used) {
	if (this->rtcm.nbyte <= 0)
		return -1;
	if (used + this->rtcm.nbyte > buflen)
		return -1;
	memcpy(buf + used, this->rtcm.buff, this->rtcm.nbyte);
	return this->rtcm.nbyte;
}

int vbs_rtcm_out_1005(struct vbs_rtcm_out *this, uint8_t *buf, int buflen) {
	int len;

	if (this == NULL || buf == NULL)
		return -1;
	if (!gen_rtcm3(&this->rtcm, 1005, 0, 0))
		return -1;
	if ((len = append_msg(this, buf, buflen, 0)) < 0)
		return -1;
	return len;
}

void vbs_rtcm_out_set_station(struct vbs_rtcm_out *this, const sta_t *sta) {
	double pos[3];

	if (this == NULL || sta == NULL)
		return;
	/* keep the VBS's own position: the descriptor is the master's, the ARP is the tile's */
	memcpy(pos, this->rtcm.sta.pos, sizeof pos);
	this->rtcm.sta = *sta;
	memcpy(this->rtcm.sta.pos, pos, sizeof pos);
	/* the ARP delta is the master's antenna's business, not the VBS's: it has none */
	this->rtcm.sta.del[0] = this->rtcm.sta.del[1] = this->rtcm.sta.del[2] = 0.0;
	this->rtcm.sta.hgt = 0.0;
	this->rtcm.sta.deltype = 0;
}

int vbs_rtcm_out_station(struct vbs_rtcm_out *this, int glonass, uint8_t *buf, int buflen) {
	int used = 0, len;

	if (this == NULL || buf == NULL)
		return -1;
	if (!gen_rtcm3(&this->rtcm, 1033, 0, 0))
		return -1;
	if ((len = append_msg(this, buf, buflen, used)) < 0)
		return -1;
	used += len;
	if (glonass) {
		if (!gen_rtcm3(&this->rtcm, 1230, 0, 0))
			return -1;
		if ((len = append_msg(this, buf, buflen, used)) < 0)
			return -1;
		used += len;
	}
	return used;
}

int vbs_rtcm_out_msm(struct vbs_rtcm_out *this, gtime_t time,
	const obsd_t *obs, int n, uint8_t *buf, int buflen) {
	int i, k, used = 0, len;
	int present[NMSM_TYPES], last = -1;

	if (this == NULL || obs == NULL || buf == NULL)
		return -1;

	/* Which constellations does this epoch actually carry? */
	for (k = 0; k < NMSM_TYPES; k++) {
		present[k] = 0;
		for (i = 0; i < n; i++)
			if (satsys(obs[i].sat, NULL) == msm_types[k].sys) {
				present[k] = 1;
				break;
			}
		if (present[k])
			last = k;
	}
	if (last < 0)
		return 0;		/* nothing to send this epoch */

	/*
	 * Guard the defect this whole file exists to avoid (D90, and the fifth instance of
	 * "implemented, documented, silently not connected" in this engine). If GLONASS is in
	 * the epoch and NO channel number has ever been installed, the caller forgot
	 * vbs_rtcm_out_set_glo_fcn(): fcn_glo() would return -1 for every satellite and the
	 * encoder would drop the whole constellation without a word. That is a programming
	 * error, not a runtime condition, so it fails here rather than on the wire.
	 *
	 * A *single* satellite missing its channel number is different -- an ephemeris can be
	 * briefly absent in production -- so that is counted (vbs_rtcm_out_glo_dropped) and the
	 * epoch still goes out. Killing the stream for one satellite would be worse than the
	 * silence this guard is here to prevent.
	 */
	for (k = 0; k < NMSM_TYPES; k++) {
		if (msm_types[k].sys != SYS_GLO || !present[k])
			continue;
		if (this->glo_fcn_n == 0)
			return -1;
		for (i = 0; i < n; i++) {
			int prn;

			if (satsys(obs[i].sat, &prn) != SYS_GLO)
				continue;
			if (prn >= 1 && prn <= (int)(sizeof this->rtcm.nav.glo_fcn /
			    sizeof this->rtcm.nav.glo_fcn[0]) &&
			    this->rtcm.nav.glo_fcn[prn - 1] == 0)
				this->glo_dropped++;
		}
	}

	this->rtcm.time = time;

	for (k = 0; k < NMSM_TYPES; k++) {
		if (!present[k])
			continue;
		/*
		 * Stage only this constellation's satellites. The encoder filters by system
		 * itself, but staging per system keeps MAXOBS headroom predictable and makes the
		 * per-message satellite count obvious when debugging a stream.
		 */
		this->rtcm.obs.n = 0;
		for (i = 0; i < n && this->rtcm.obs.n < MAXOBS; i++) {
			if (satsys(obs[i].sat, NULL) != msm_types[k].sys)
				continue;
			this->obs[this->rtcm.obs.n] = obs[i];
			this->obs[this->rtcm.obs.n].time = time;
			this->rtcm.obs.n++;
		}
		if (this->rtcm.obs.n == 0)
			continue;

		/* sync = 1 on every message but the last of the epoch. */
		if (!gen_rtcm3(&this->rtcm, msm_types[k].type, 0, k != last)) {
			this->rtcm.obs.n = 0;
			return -1;
		}
		if ((len = append_msg(this, buf, buflen, used)) < 0) {
			this->rtcm.obs.n = 0;
			return -1;
		}
		used += len;
	}
	this->rtcm.obs.n = 0;
	return used;
}
