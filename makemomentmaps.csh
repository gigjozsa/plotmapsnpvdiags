#! /usr/bin/tcsh

set prefix = eso149
set incube_original = eso149_in.fits
set incube_model = final_eso149_out.fits
# Calculations will be done on this sub-region
set region = "box(19,14,65,60)"
set conv_kernel = 60
set threshold_original = 0.003
set threshold_model = 0.0009
set mom1thresh = 0.1
set optical_1 = dss2_red_gj.fits
set ellipsepos_x = 358.01
set ellipsepos_y = -52.5550
set ellipse_maj = 0.0375
set ellipse_min = 0.0275
set ellipse_pa = 30
set distance = 7.01

#Script below:
#${prefix}_original.fits: original data cube reduced to size as given in region
#${prefix}_model.fits: model data cube reduced to size as given in region
# ${prefix}_normal_mom0.fits: Convolve measured cube to conv_kernel arcsec beam and clip at threshold_original Jy/beam, then blink, create mom0
# ${prefix}_normal_mom0_sol.fits: Convolve measured cube to conv_kernel arcsec beam and clip at threshold_original Jy/beam, then blink, create mom0, solar masses per square parsec
# ${prefix}_normal_mom1.fits: Convolve measured cube to conv_kernel arcsec beam and clip at threshold_original Jy/beam, then blink, create mom1 and clip that at mom1thresh using the moment 0
# ${prefix}_cdiff_mom0.fits: Subtract model cube from measured cube, then convolve to conv_kernel arcsec beam, then clip at 0.003 Jy/beam, blink, create mom0
# ${prefix}_cdiff_mom0_sol.fits op=xyout: Subtract model cube from measured cube, then convolve to conv_kernel arcsec beam, then clip at 0.003 Jy/beam, blink, create mom0, solar masses per square parsec
# ${prefix}_model_normal_mom0.fits: Convolve measured cube to conv_kernel arcsec beam and clip at threshold_model Jy/beam, then blink, create mom0
# ${prefix}_model_normal_mom0_sol.fits: Convolve measured cube to conv_kernel arcsec beam and clip at threshold_model Jy/beam, then blink, create mom0, solar masses per square parsec
# ${prefix}_model_normal_mom1.fits: Convolve measured cube to conv_kernel arcsec beam and clip at threshold_model Jy/beam, then blink, create mom0
# ${prefix}_mdiff_normal_mom0.fits: Subtract ${prefix}_model_normal_mom0.fits from ${prefix}_normal_mom0.fits
# ${prefix}_mdiff_normal_mom0_sol.fits: Subtract ${prefix}_model_normal_mom0.fits from ${prefix}_normal_mom0.fits in solar masses per square parsec
# ${prefix}_not_mask_mom0.fits: Exclude the mask from the model created for ${prefix}_mdiff_normal_mom0.fits from the mask for ${prefix}_normal_mom0.fits and calculate mom0 from that.
# ${prefix}_not_mask_mom0_sol.fits: ${prefix}_not_mask_mom0.fits in solar masses per square parsec
# ${prefix}_or_mask_mom0.fits: combine masks from the model created for ${prefix}_mdiff_normal_mom0.fits and for ${prefix}_normal_mom0.fits and calculate mom0 from that mask applied on the measured data
# ${prefix}_or_mask_mom0_sol.fits: Same in solar masses per square parsec
# ${prefix}_or_mask_mom1.fits: combine masks from the model created for ${prefix}_mdiff_normal_mom0.fits and for ${prefix}_normal_mom0.fits and calculate mom1 from that mask applied on the measured data, then clip to mom1thresh on the moment0
# ${prefix}_model_or_mask_mom0.fits: Same for model
# ${prefix}_model_or_mask_mom0_sol.fits: Same for model in solar masses per square parsec
# ${prefix}_model_or_mask_mom1.fits: Same thing for model moment 1

rm -r ${prefix}_raw
rm -r ${prefix}_model_raw
rm -r ${prefix}
rm -r ${prefix}_model
fits in=${incube_original} out=${prefix}_raw op=xyin
fits in=${incube_model} out=${prefix}_model_raw op=xyin
imsub in=${prefix}_raw region=${region} out=${prefix}
imsub in=${prefix}_model_raw region=${region} out=${prefix}_model
rm -r ${prefix}_original.fits
rm -r ${prefix}_model.fits
fits in=${prefix} out=${prefix}_original.fits op=xyout
fits in=${prefix}_model out=${prefix}_model.fits op=xyout

set bmaj = `itemize in=eso149 | grep bmaj | awk '{print $3}'`
set bmin = `itemize in=eso149 | grep bmin | awk '{print $3}'`

# Create normal moment0
rm -r ${prefix}_con
rm -r ${prefix}_masked
rm -r ${prefix}_masked_mom0
rm -r ${prefix}_masked_mom1_raw
rm -r ${prefix}_masked_mom1
rm -r ${prefix}_normal_mom1.fits
rm -r ${prefix}_masked_mom0_blr
rm -r ${prefix}_normal_mom0.fits
rm -r ${prefix}_masked_mom0_blr_sol
rm -r ${prefix}_normal_mom0_sol.fits
convol map=${prefix} fwhm=${conv_kernel} options=final out=${prefix}_con
maths exp="<${prefix}>" mask="<${prefix}_con>.gt.${threshold_original}" out=${prefix}_masked
mafia in=${prefix}_masked options=mask,3dim,blink
moment in=${prefix}_masked out=${prefix}_masked_mom0
moment in=${prefix}_masked out=${prefix}_masked_mom1_raw mom=1
maths exp="<${prefix}_masked_mom1_raw>" mask="<${prefix}_masked_mom0>.gt.${mom1thresh}" out=${prefix}_masked_mom1
fits in=${prefix}_masked_mom1 out=${prefix}_normal_mom1.fits op=xyout
imblr in=${prefix}_masked_mom0 out=${prefix}_masked_mom0_blr
fits in=${prefix}_masked_mom0_blr out=${prefix}_normal_mom0.fits op=xyout
maths exp="<${prefix}_masked_mom0_blr>*1.248683E24*8.01325e-21/(1.1330900354567984*(180*3600/(4*ATAN(1.d0)))**2*${bmaj}*${bmin})" out=${prefix}_masked_mom0_blr_sol
puthd in=${prefix}_masked_mom0_blr_sol/bunit value='MSOL/PC**2'
fits in=${prefix}_masked_mom0_blr_sol out=${prefix}_normal_mom0_sol.fits op=xyout

# Create mom0 from difference cube
rm -r ${prefix}_cdiff
rm -r ${prefix}_cdiff_con
rm -r ${prefix}_cdiff_masked
rm -r ${prefix}_cdiff_masked_mom0
rm -r ${prefix}_cdiff_masked_mom0_blr
rm -r ${prefix}_cdiff_mom0.fits
rm -r ${prefix}_cdiff_masked_mom0_blr_sol
rm -r ${prefix}_cdiff_mom0_sol.fits
# Create difference cube
maths exp="<${prefix}>-<${prefix}_model>" out=${prefix}_cdiff
convol map=${prefix}_cdiff fwhm=${conv_kernel} options=final out=${prefix}_cdiff_con
maths exp="<${prefix}_cdiff>" mask="<${prefix}_cdiff_con>.gt.${threshold_original}" out=${prefix}_cdiff_masked
mafia in=${prefix}_cdiff_masked options=mask,3dim,blink
moment in=${prefix}_cdiff_masked out=${prefix}_cdiff_masked_mom0
imblr in=${prefix}_cdiff_masked_mom0 out=${prefix}_cdiff_masked_mom0_blr
fits in=${prefix}_cdiff_masked_mom0_blr out=${prefix}_cdiff_mom0.fits op=xyout
maths exp="<${prefix}_cdiff_masked_mom0_blr>*1.248683E24*8.01325e-21/(1.1330900354567984*(180*3600/(4*ATAN(1.d0)))**2*${bmaj}*${bmin})" out=${prefix}_cdiff_masked_mom0_blr_sol
puthd in=${prefix}_cdiff_masked_mom0_blr_sol/bunit value='MSOL/PC**2'
fits in=${prefix}_cdiff_masked_mom0_blr_sol out=${prefix}_cdiff_mom0_sol.fits op=xyout

# Create difference mom0 with model mom0 expanded until in most areas normal mom0 is masked out
rm -r ${prefix}_model_con
rm -r ${prefix}_model_masked
rm -r ${prefix}_model_masked_mom0
rm -r ${prefix}_model_masked_mom1_raw
rm -r ${prefix}_model_masked_mom1
rm -r ${prefix}_model_masked_mom0_blr
rm -r ${prefix}_model_normal_mom0.fits
rm -r ${prefix}_model_masked_mom0_blr_sol
rm -r ${prefix}_model_normal_mom0_sol.fits
rm -r ${prefix}_mdiff_normal_mom0
rm -r ${prefix}_mdiff_normal_mom0.fits
rm -r ${prefix}_mdiff_normal_mom0_sol
rm -r ${prefix}_mdiff_normal_mom0_sol.fits
convol map=${prefix}_model fwhm=${conv_kernel} options=final out=${prefix}_model_con
maths exp="<${prefix}_model>" mask="<${prefix}_model_con>.gt.${threshold_model}" out=${prefix}_model_masked
#maths exp="<${prefix}_model>" mask="<${prefix}_con>.gt.${threshold_original}" out=${prefix}_model_masked
mafia in=${prefix}_model_masked options=mask,3dim,blink
moment in=${prefix}_model_masked out=${prefix}_model_masked_mom0
moment in=${prefix}_model_masked out=${prefix}_model_masked_mom1_raw mom=1
maths exp="<${prefix}_model_masked_mom1_raw>" mask="<${prefix}_model_masked_mom0>.gt.${mom1thresh}" out=${prefix}_model_masked_mom1
fits in=${prefix}_model_masked_mom1 out=${prefix}_model_normal_mom1.fits op=xyout
imblr in=${prefix}_model_masked_mom0 out=${prefix}_model_masked_mom0_blr
fits in=${prefix}_model_masked_mom0_blr out=${prefix}_model_normal_mom0.fits op=xyout
maths exp="<${prefix}_model_masked_mom0_blr>*1.248683E24*8.01325e-21/(1.1330900354567984*(180*3600/(4*ATAN(1.d0)))**2*${bmaj}*${bmin})" out=${prefix}_model_masked_mom0_blr_sol
puthd in=${prefix}_model_masked_mom0_blr_sol/bunit value='MSOL/PC**2'
fits in=${prefix}_model_masked_mom0_blr_sol out=${prefix}_model_normal_mom0_sol.fits op=xyout
maths exp="<${prefix}_masked_mom0_blr>-<${prefix}_model_masked_mom0_blr>" out=${prefix}_mdiff_normal_mom0
fits in=${prefix}_mdiff_normal_mom0 out=${prefix}_mdiff_normal_mom0.fits op=xyout
maths exp="<${prefix}_mdiff_normal_mom0>*1.248683E24*8.01325e-21/(1.1330900354567984*(180*3600/(4*ATAN(1.d0)))**2*${bmaj}*${bmin})" out=${prefix}_mdiff_normal_mom0_sol
puthd in=${prefix}_mdiff_normal_mom0_sol/bunit value='MSOL/PC**2'
fits in=${prefix}_mdiff_normal_mom0_sol out=${prefix}_mdiff_normal_mom0_sol.fits op=xyout

# Mask out further regions in normal mask using model mask
rm -r ${prefix}_model_masked_blr
rm -r ${prefix}_not_mask
rm -r ${prefix}_not_mask_mom0
rm -r ${prefix}_not_mask_mom0_blr
rm -r ${prefix}_not_mask_mom0.fits
rm -r ${prefix}_not_mask_mom0_blr_sol
rm -r ${prefix}_not_mask_mom0_sol.fits
imblr in=${prefix}_model_masked out=${prefix}_model_masked_blr
maths exp="<${prefix}_masked>" mask="<${prefix}_model_masked_blr>.eq.0" out=${prefix}_not_mask
moment in=${prefix}_not_mask out=${prefix}_not_mask_mom0
imblr in=${prefix}_not_mask_mom0 out=${prefix}_not_mask_mom0_blr
fits in=${prefix}_not_mask_mom0_blr out=${prefix}_not_mask_mom0.fits op=xyout
maths exp="<${prefix}_not_mask_mom0_blr>*1.248683E24*8.01325e-21/(1.1330900354567984*(180*3600/(4*ATAN(1.d0)))**2*${bmaj}*${bmin})" out=${prefix}_not_mask_mom0_blr_sol
puthd in=${prefix}_not_mask_mom0_blr_sol/bunit value='MSOL/PC**2'
fits in=${prefix}_not_mask_mom0_blr_sol out=${prefix}_not_mask_mom0_sol.fits op=xyout

##
# Create combined mask for model and original data
rm -r blam1
rm -r blam1b
rm -r blam2
rm -r blam2b
rm -r ${prefix}_or_mask
rm -r ${prefix}_model_or_mask
maths exp="<${prefix}>*0+1" mask="<${prefix}_con>.gt.${threshold_original}" out=blam1
mafia in=blam1 options=mask,3dim,blink
imblr in=blam1 out=blam1b
maths exp="<${prefix}>*0+1" mask="<${prefix}_model_con>.gt.${threshold_model}" out=blam2
mafia in=blam2 options=mask,3dim,blink
imblr in=blam2 out=blam2b
maths exp="<${prefix}>" mask="(<blam1b>.gt.0).or.(<blam2b>.gt.0)" out=${prefix}_or_mask
maths exp="<${prefix}_model>" mask="(<blam1b>.gt.0).or.(<blam2b>.gt.0)" out=${prefix}_model_or_mask

# Create moment0 on model and original data
# Data
rm -r ${prefix}_or_mask_mom1_raw
rm -r ${prefix}_or_mask_mom1
rm -r ${prefix}_or_mask_mom1.fits
rm -r ${prefix}_or_mask_mom0
rm -r ${prefix}_or_mask_mom0_blr
rm -r ${prefix}_or_mask_mom0.fits
rm -r ${prefix}_or_mask_mom0_blr_sol
rm -r ${prefix}_or_mask_mom0_sol.fits
moment in=${prefix}_or_mask out=${prefix}_or_mask_mom0
moment in=${prefix}_or_mask out=${prefix}_or_mask_mom1_raw mom=1
maths exp="<${prefix}_or_mask_mom1_raw>" mask="<${prefix}_or_mask_mom0>.gt.${mom1thresh}" out=${prefix}_or_mask_mom1
fits in=${prefix}_or_mask_mom1 out=${prefix}_or_mask_mom1.fits op=xyout
imblr in=${prefix}_or_mask_mom0 out=${prefix}_or_mask_mom0_blr
fits in=${prefix}_or_mask_mom0_blr out=${prefix}_or_mask_mom0.fits op=xyout
maths exp="<${prefix}_or_mask_mom0_blr>*1.248683E24*8.01325e-21/(1.1330900354567984*(180*3600/(4*ATAN(1.d0)))**2*${bmaj}*${bmin})" out=${prefix}_or_mask_mom0_blr_sol
puthd in=${prefix}_or_mask_mom0_blr_sol/bunit value='MSOL/PC**2'
fits in=${prefix}_or_mask_mom0_blr_sol out=${prefix}_or_mask_mom0_sol.fits op=xyout

# Models
rm -r ${prefix}_model_or_mask_mom0
rm -r ${prefix}_model_or_mask_mom1_raw
rm -r ${prefix}_model_or_mask_mom1
rm -r ${prefix}_model_or_mask_mom1.fits
rm -r ${prefix}_model_or_mask_mom0_blr
rm -r ${prefix}_model_or_mask_mom0.fits
rm -r ${prefix}_model_or_mask_mom0_blr_sol
rm -r ${prefix}_model_or_mask_mom0_sol.fits
moment in=${prefix}_model_or_mask out=${prefix}_model_or_mask_mom0
moment in=${prefix}_model_or_mask out=${prefix}_model_or_mask_mom1_raw mom=1
maths exp="<${prefix}_model_or_mask_mom1_raw>" mask="<${prefix}_model_or_mask_mom0>.gt.${mom1thresh}" out=${prefix}_model_or_mask_mom1
fits in=${prefix}_model_or_mask_mom1 out=${prefix}_model_or_mask_mom1.fits op=xyout
imblr in=${prefix}_model_or_mask_mom0 out=${prefix}_model_or_mask_mom0_blr
fits in=${prefix}_model_or_mask_mom0_blr out=${prefix}_model_or_mask_mom0.fits op=xyout
maths exp="<${prefix}_model_or_mask_mom0_blr>*1.248683E24*8.01325e-21/(1.1330900354567984*(180*3600/(4*ATAN(1.d0)))**2*${bmaj}*${bmin})" out=${prefix}_model_or_mask_mom0_blr_sol
puthd in=${prefix}_model_or_mask_mom0_blr_sol/bunit value='MSOL/PC**2'
fits in=${prefix}_model_or_mask_mom0_blr_sol out=${prefix}_model_or_mask_mom0_sol.fits op=xyout

## Subtract model-moment-0 from original on combined mask
rm -r ${prefix}_or_mdiff_mom0
rm -r ${prefix}_or_mdiff_mom0.fits
rm -r ${prefix}_or_mdiff_mom0_sol
rm -r ${prefix}_or_mdiff_mom0_sol.fits
maths exp="<${prefix}_or_mask_mom0_blr>-<${prefix}_model_or_mask_mom0_blr>" out=${prefix}_or_mdiff_mom0
fits in=${prefix}_or_mdiff_mom0 out=${prefix}_or_mdiff_mom0.fits op=xyout
maths exp="<${prefix}_or_mdiff_mom0>*1.248683E24*8.01325e-21/(1.1330900354567984*(180*3600/(4*ATAN(1.d0)))**2*${bmaj}*${bmin})" out=${prefix}_or_mdiff_mom0_sol
puthd in=${prefix}_or_mdiff_mom0_sol/bunit value='MSOL/PC**2'
fits in=${prefix}_or_mdiff_mom0_sol out=${prefix}_or_mdiff_mom0_sol.fits op=xyout

# Some trick to get an optical image of the size of the radio image
rm -r ${prefix}_optical_1
rm -r ${prefix}_model_masked_mom0_blr_sol_reg
rm -r ${prefix}_optical_1_masked
rm -r ${prefix}_preoptical
rm -r ${prefix}_preoptical_blr
rm -r ${prefix}_optical
rm -r ${prefix}_optical_dss_r.fits
fits in=${optical_1} out=${prefix}_optical_1 op=xyin
regrid in=${prefix}_model_masked_mom0_blr_sol tin=${prefix}_optical_1 axes=1,2 out=${prefix}_model_masked_mom0_blr_sol_reg
maths exp="<${prefix}_optical_1>" mask="<${prefix}_model_masked_mom0_blr_sol_reg>.gt.-100000000" out=${prefix}_optical_1_masked
imsub in=${prefix}_optical_1 region=mask"(${prefix}_optical_1_masked)" out=${prefix}_preoptical
imblr in=${prefix}_preoptical out=${prefix}_preoptical_blr
regrid in=${prefix}_optical_1 tin=${prefix}_preoptical axes=1,2 out=${prefix}_optical
fits in=${prefix}_optical out=${prefix}_optical_dss_r.fits op=xyout
convfits.py eso149_optical_dss_r.fits 3 3

# Measurements
rm -r ${prefix}_or_mask_mom0_blr_ellipse
rm -r ${prefix}_not_mask_mom0_blr_emasked
rm -r ${prefix}_or_mdiff_mom0_emasked
rm -r ${prefix}_cdiff_masked_mom0_blr_emasked
rm -r ${prefix}_not_mask_mom0_blr_emasked.fits
rm -r ${prefix}_or_mdiff_mom0_emasked.fits
rm -r ${prefix}_cdiff_masked_mom0_blr_emasked.fits

set ellipse_maj_as = `calc ${ellipse_maj}"*"3600`
set ellipse_min_as = `calc ${ellipse_min}"*"3600`

set xoffset = `impos in=${prefix}_or_mask_mom0_blr coord=${ellipsepos_x},${ellipsepos_y} type=absdeg,absdeg | grep -A 2 "Offset world coordinates" | awk '{if (match($2,"1:")) {printf("%.6f",$5)}}'`
set yoffset = `impos in=${prefix}_or_mask_mom0_blr coord=${ellipsepos_x},${ellipsepos_y} type=absdeg,absdeg | grep -A 2 "Offset world coordinates" | awk '{if (match($2,"2:")) {printf("%.6f",$5)}}'`

imgen in=${prefix}_or_mask_mom0_blr out=${prefix}_or_mask_mom0_blr_ellipse factor=0 object=disk spar=1,${xoffset},${yoffset},${ellipse_maj_as},${ellipse_min_as},${ellipse_pa}

maths exp="<${prefix}_not_mask_mom0_blr>*<${prefix}_or_mask_mom0_blr_ellipse>" out=${prefix}_not_mask_mom0_blr_emasked
maths exp="<${prefix}_or_mdiff_mom0>*<${prefix}_or_mask_mom0_blr_ellipse>" out=${prefix}_or_mdiff_mom0_emasked
maths exp="<${prefix}_cdiff_masked_mom0_blr>*<${prefix}_or_mask_mom0_blr_ellipse>" out=${prefix}_cdiff_masked_mom0_blr_emasked

fits in=${prefix}_not_mask_mom0_blr_emasked      op=xyout out=${prefix}_not_mask_mom0_blr_emasked.fits     
fits in=${prefix}_or_mdiff_mom0_emasked         op=xyout out=${prefix}_or_mdiff_mom0_emasked.fits        
fits in=${prefix}_cdiff_masked_mom0_blr_emasked op=xyout out=${prefix}_cdiff_masked_mom0_blr_emasked.fits

getfluxnmass.py ${prefix}_or_mask_mom0.fits 7.01 ${distance}
getfluxnmass.py ${prefix}_not_mask_mom0_blr_emasked.fits ${distance}
getfluxnmass.py ${prefix}_or_mdiff_mom0_emasked.fits ${distance}
getfluxnmass.py ${prefix}_cdiff_masked_mom0_blr_emasked.fits ${distance}

rm -r ${prefix}_raw
rm -r ${prefix}_model_raw
rm -r ${prefix}
rm -r ${prefix}_model
rm -r ${prefix}_con
rm -r ${prefix}_masked
rm -r ${prefix}_masked_mom0
rm -r ${prefix}_masked_mom1_raw
rm -r ${prefix}_masked_mom1
#rm -r ${prefix}_normal_mom1.fits
rm -r ${prefix}_masked_mom0_blr
#rm -r ${prefix}_normal_mom0.fits
rm -r ${prefix}_masked_mom0_blr_sol
#rm -r ${prefix}_normal_mom0_sol.fits
rm -r ${prefix}_cdiff
rm -r ${prefix}_cdiff_con
rm -r ${prefix}_cdiff_masked
rm -r ${prefix}_cdiff_masked_mom0
rm -r ${prefix}_cdiff_masked_mom0_blr
#rm -r ${prefix}_cdiff_mom0.fits
rm -r ${prefix}_cdiff_masked_mom0_blr_sol
#rm -r ${prefix}_cdiff_mom0_sol.fits
rm -r ${prefix}_model_con
rm -r ${prefix}_model_masked
rm -r ${prefix}_model_masked_mom0
rm -r ${prefix}_model_masked_mom1_raw
rm -r ${prefix}_model_masked_mom1
rm -r ${prefix}_model_masked_mom0_blr
#rm -r ${prefix}_model_normal_mom0.fits
rm -r ${prefix}_model_masked_mom0_blr_sol
#rm -r ${prefix}_model_normal_mom0_sol.fits
rm -r ${prefix}_mdiff_normal_mom0
#rm -r ${prefix}_mdiff_normal_mom0.fits
rm -r ${prefix}_mdiff_normal_mom0_sol
#rm -r ${prefix}_mdiff_normal_mom0_sol.fits
rm -r ${prefix}_model_masked_blr
rm -r ${prefix}_not_mask
rm -r ${prefix}_not_mask_mom0
rm -r ${prefix}_not_mask_mom0_blr
#rm -r ${prefix}_not_mask_mom0.fits
rm -r ${prefix}_not_mask_mom0_blr_sol
#rm -r ${prefix}_not_mask_mom0_sol.fits
rm -r blam1
rm -r blam1b
rm -r blam2
rm -r blam2b
rm -r ${prefix}_or_mask
rm -r ${prefix}_model_or_mask
rm -r ${prefix}_or_mask_mom1_raw
rm -r ${prefix}_or_mask_mom1
#rm -r ${prefix}_or_mask_mom1.fits
rm -r ${prefix}_or_mask_mom0
rm -r ${prefix}_or_mask_mom0_blr
#rm -r ${prefix}_or_mask_mom0.fits
rm -r ${prefix}_or_mask_mom0_blr_sol
#rm -r ${prefix}_or_mask_mom0_sol.fits
rm -r ${prefix}_model_or_mask_mom0
rm -r ${prefix}_model_or_mask_mom1_raw
rm -r ${prefix}_model_or_mask_mom1
#rm -r ${prefix}_model_or_mask_mom1.fits
rm -r ${prefix}_model_or_mask_mom0_blr
#rm -r ${prefix}_model_or_mask_mom0.fits
rm -r ${prefix}_model_or_mask_mom0_blr_sol
#rm -r ${prefix}_model_or_mask_mom0_sol.fits
rm -r ${prefix}_or_mdiff_mom0
#rm -r ${prefix}_or_mdiff_mom0.fits
rm -r ${prefix}_or_mdiff_mom0_sol
#rm -r ${prefix}_or_mdiff_mom0_sol.fits
rm -r ${prefix}_optical_1
rm -r ${prefix}_model_masked_mom0_blr_sol_reg
rm -r ${prefix}_optical_1_masked
rm -r ${prefix}_optical
#rm -r ${prefix}_optical.fits
rm -r ${prefix}_or_mask_mom0_blr_ellipse
rm -r ${prefix}_not_mask_mom0_blr_emasked
rm -r ${prefix}_or_mdiff_mom0_emasked
rm -r ${prefix}_cdiff_masked_mom0_blr_emasked
#rm -r ${prefix}_not_mask_mom0_blr_emasked.fits
#rm -r ${prefix}_or_mdiff_mom0_emasked.fits
#rm -r ${prefix}_cdiff_masked_mom0_blr_emasked.fits
