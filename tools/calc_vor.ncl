begin
;************************************************
; variable and file handling
;************************************************
  in  = addfile("mz.u_pole.dt60.360x180.h0.nc","r")                             ; open netcdf file
  u   = in->ua                                            ; pull u off file
  v   = in->va                                            ; pull v off file
  lat = in->lat
  lon = in->lon
;************************************************
; calculate vorticity on a Gaussian Grid
; scale for visual convenience
;************************************************
  scale = 1.e05
  vrt   = u                                      ; retain coordinates
  vrt   = uv2vr_cfd(u,v,lat,lon,3) * scale
  vrt@long_name = "vorticity"
  vrt@units     = "10^-5"
;************************************************
; create plot
;************************************************
  wks  = gsn_open_wks("pdf","vort")              ; send graphics to PNG file
  gsn_define_colormap(wks, "ncl_default") 
  
  res                   = True                   ; plot mods desired
;  res@mpProjection      = "Orthographic"
  res@gsnPolar          = "NH"
  res@mpLimitMode       = "LatLon"
  res@mpCenterLonF      = 0
  res@mpCenterLatF      = 90
  res@mpMinLatF         = 70
  res@mpGridAndLimbOn   = True
  res@mpOutlineOn       = False
  res@cnFillOn          = True                   ; turn on color
;  res@cnFillPalette     = "BlWhRe"               ; set color map
  res@cnLinesOn         = True                  ; no contour lines

  res@cnLevelSelectionMode = "ManualLevels"       ; set manual contour levels
  res@cnMinLevelValF       = -2.7                   ; set min contour level
  res@cnMaxLevelValF       =  3.4                   ; set max contour level
  res@cnLevelSpacingF      =  0.2                 ; set contour spacing


  plot = gsn_csm_contour_map_polar(wks,vrt(24,:,:),res) ; create plot  
;************************************************
end
