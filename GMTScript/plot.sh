#!/bin/bash

#===============================================================================
# Author: Tche L., USTC, seistche@gmail.com
# Created at: Thu 14 Apr 2022 04:07:31 PM CST
#-------------------------------------------------------------------------------

ter=ter.nc
pgv=PGV.nc
ann=ann.txt
cnt=cnt.lev
fau=faultLocationOnSurf.txt
R=101.7200/103.9023/34.9473/36.7096
R1=101.7200/103.9023/34.9473/36.7096
J=M10c

J1=m5c

illu=illu.nc
#vext=`gmt grdinfo $pgv -T`
#vext=${vext##-T}
vext=-6/-0.6

B="a0.5f0.25"

grd=/data0/home/wangwq/earth_relief_15s.grd

gmt begin pdf
  gmt set PS_MEDIA A0
  gmt set FONT_ANNOT_PRIMARY=10p,4
  gmt set MAP_FRAME_WIDTH=0.1c

  gmt grdgradient $ter -A20 -Nt0.99 -G$illu

  gmt basemap -J$J -R$R -BWSen -B$B
  gmt grd2cpt $ter -CgrayC -Z
  gmt grdimage -J$J -R$R $ter -I$illu

  gmt makecpt -Cjet -D -Z -T$vext/0.01  #-T-9/1.5/0.1 #
  gmt grdimage -J$J -R$R $pgv -I$illu -t20


  gmt grdcontour -J$J -R$R $pgv -W0.5p,black -C$cnt

  gmt colorbar -DjBC+w8c/0.2c+o0c/-1.2c+m+h+e -Bpxc$ann -By+l"  m/s" \
    --MAP_TICK_LENGTH_PRIMARY=2.5p -G-6/-0.6
gmt end #show
