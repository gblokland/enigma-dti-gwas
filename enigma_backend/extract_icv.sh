#!/usr/bin/env bash
# =====================================================================
#  Universal ICV extraction script for T1 MRI data
#  Supports FreeSurfer, FSL, ANTs, and Python/nibabel fallback.
#  Author: ChatGPT (GPT-5)
# =====================================================================

set -euo pipefail

SUBJLIST=${1:-subjlist.txt}
ROOTDIR=${2:-$PWD}
OUTCSV="icv_summary.csv"

# --- detect available tools ---------------------------------------------------
has_freesurfer=$(command -v freesurfer || command -v recon-all >/dev/null && echo 1 || echo 0)
has_fsl=$(command -v fslstats >/dev/null && echo 1 || echo 0)
has_ants=$(command -v antsBrainExtraction.sh >/dev/null && echo 1 || echo 0)
has_python=$(command -v python >/dev/null && echo 1 || echo 0)

echo "Detected tools:"
echo "  FreeSurfer: $has_freesurfer"
echo "  FSL:         $has_fsl"
echo "  ANTs:        $has_ants"
echo "  Python:      $has_python"
echo "==================================================================="

echo "Subject,Method,ICV_mm3,ICV_mL" > "$OUTCSV"

while read -r subj; do
  [ -z "$subj" ] && continue
  echo "Processing $subj ..."

  subjdir="${ROOTDIR}/${subj}"
  icv_mm3=""
  method=""

  # --- 1) FreeSurfer eTIV -----------------------------------------------------
  if [ "$has_freesurfer" -eq 1 ] && [ -f "${subjdir}/stats/aseg.stats" ]; then
    aseg="${subjdir}/stats/aseg.stats"
    icv_mm3=$(grep -i "EstimatedTotalIntraCranialVol" "$aseg" | awk '{print $NF}' || true)
    [ -z "$icv_mm3" ] && icv_mm3=$(grep -i "IntraCranialVol" "$aseg" | awk '{print $NF}' || true)
    if [ -n "$icv_mm3" ]; then
      method="FreeSurfer_eTIV"
    fi
  fi

  # --- 2) FSL fallback --------------------------------------------------------
  if [ -z "$icv_mm3" ] && [ "$has_fsl" -eq 1 ]; then
    t1=$(find "$subjdir" -maxdepth 2 -type f -name "*T1*.nii*" | head -n1)
    if [ -n "$t1" ]; then
      outbase="${t1%.nii*}_brain"
      mask="${outbase}_mask.nii.gz"
      echo "  Using FSL BET for ${subj}..."
      bet "$t1" "$outbase" -R -f 0.5 -g 0 || true
      fslmaths "${outbase}.nii.gz" -bin "$mask"
      read nvox icv_mm3 <<< $(fslstats "$mask" -V)
      method="FSL_BET"
    fi
  fi

  # --- 3) ANTs fallback -------------------------------------------------------
  if [ -z "$icv_mm3" ] && [ "$has_ants" -eq 1 ]; then
    t1=$(find "$subjdir" -maxdepth 2 -type f -name "*T1*.nii*" | head -n1)
    if [ -n "$t1" ]; then
      echo "  Using ANTs BrainExtraction for ${subj}..."
      antsBrainExtraction.sh -d 3 -a "$t1" \
        -e $ANTSPATH/tpl/OASIS_template.nii.gz \
        -m $ANTSPATH/tpl/OASIS_brainmask.nii.gz \
        -o "${subjdir}/ants_"
      mask="${subjdir}/ants_BrainExtractionBrainMask.nii.gz"
      read nvox icv_mm3 <<< $(fslstats "$mask" -V)
      method="ANTs_BrainExtraction"
    fi
  fi

  # --- 4) Python nibabel fallback --------------------------------------------
  if [ -z "$icv_mm3" ] && [ "$has_python" -eq 1 ]; then
    mask=$(find "$subjdir" -maxdepth 2 -type f -name "*mask*.nii*" | head -n1)
    if [ -n "$mask" ]; then
      echo "  Using Python/nibabel for ${subj}..."
      icv_mm3=$(python - <<'PYCODE' "$mask"
import nibabel as nib, numpy as np, sys
img = nib.load(sys.argv[1])
data = img.get_fdata()
voxvol = np.prod(img.header.get_zooms())
print(np.sum(data > 0.5) * voxvol)
PYCODE
)
      method="Python_mask_volume"
    fi
  fi

  # --- results ---------------------------------------------------------------
  if [ -n "$icv_mm3" ]; then
    icv_ml=$(awk -v v="$icv_mm3" 'BEGIN{printf "%.2f", v/1000}')
    echo "$subj,$method,$icv_mm3,$icv_ml" >> "$OUTCSV"
    echo "  -> $method: $icv_ml mL"
  else
    echo "$subj,NA,NA,NA" >> "$OUTCSV"
    echo "  !! No ICV could be determined for $subj"
  fi
done < "$SUBJLIST"

echo "==================================================================="
echo "✅ Done. Summary written to: $OUTCSV"
echo "==================================================================="
