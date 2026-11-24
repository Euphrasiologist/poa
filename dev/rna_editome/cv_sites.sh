#!/bin/bash
# cv_sites.sh KIND SEED TAG SPECIES... : cv.py for angiosperms, genus held out, with orfedit/transsplice that know editing
# sites (final3/bin). KIND ds: dataset species, both organelles; gb: RefSeq test genomes, mito only. Settings from the
# environment: EP, PP (edit cost mito, plastid), KN, PKN (cost at a column edited in the references), NV, PNV (at one not edited), SG (genera a site needs)
set -uo pipefail
R=/lustre/scratch122/tol/teams/blaxter/users/mb39/plant_organelle_annotator; H=$R/build/rnaedit
K=$1; SEED=$2; TAG=$3; shift 3
X="CV_NONE=1"; [ "$K" == "ds" ] && X="$X USE_DATASET=1" || X="$X CV_NO_PLASTID=1"
[ -n "${EP:-}" ] && X="$X CV_EDIT_PENALTY=$EP"; [ -n "${PP:-}" ] && X="$X CV_PLTD_EDIT_PENALTY=$PP"
[ -n "${KN:-}" ] && X="$X CV_KNOWN=$KN"; [ -n "${PKN:-}" ] && X="$X CV_PLTD_KNOWN=$PKN"; [ -n "${NV:-}" ] && X="$X CV_NOVEL=$NV"; [ -n "${PNV:-}" ] && X="$X CV_PLTD_NOVEL=$PNV"; [ -n "${SG:-}" ] && X="$X CV_SITE_GROUPS=$SG"
env -i HOME=$HOME LSB_JOBID=$LSB_JOBID PYTHONNOUSERSITE=1 TMPDIR=/tmp PYTHONPATH=$R/src PATH=$H/final3/bin:$R/build/pkgenv/bin:/usr/bin:/bin $X python $H/cv.py angiosperm $SEED $TAG genus "$@" > $H/pilot/cv.$TAG.$LSB_JOBID.log 2>&1 || { echo "cv failed"; tail -5 $H/pilot/cv.$TAG.$LSB_JOBID.log; exit 1; }
grep -c "held out" $H/pilot/cv.$TAG.$LSB_JOBID.log
