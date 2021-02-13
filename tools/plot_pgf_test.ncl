;*******************************
; plot_pgf_test.ncl
;*******************************

begin

  fi0 = "../build/pgf.360x180/pgf.dt60.360x181.l30_vor_damp_lin97.h0.nc"
  fi1 = "../build/pgf.360x180/pgf.dt60.360x181.l30_vor_damp_lin97.h1.nc"
  
  out_name = "pgf_test_case_lev30_reduce_vor_damp"
  wks = gsn_open_wks("eps", out_name)
  gsn_define_colormap(wks, "gui_default")  
  
  ;*******************************
  ; Panel plot
  ;*******************************
  plot = new(2, graphic)
  
  df0  = addfile(fi0, "r")

  u    = df0->u(24,:,90,:)
  zs   = df0->zs

  lon  = df0->lon
  ilon = df0->ilon
  ilev = df0->ilev
  lev  = df0->lev

  res_u                   = True
  res_u@gsnDraw           = False
  res_u@gsnFrame          = False
  res_u@gsnAddCyclic      = True
  res_u@cnFillOn          = True
  res_u@cnLinesOn         = False
  res_u@lbLabelAutoStride = True
  
  res_u@cnLevelSelectionMode = "ManualLevels"
  res_u@cnLevelSpacingF      = 0.02
  res_u@cnMinLevelValF       = -0.1
  res_u@cnMaxLevelValF       = 0.1
  res_u@trYReverse           = True  
  res_u@vpWidthF             = 0.6
  res_u@vpHeightF            = 0.3
  res_u@gsnCenterStringFontHeightF = 0.015
  res_u@gsnCenterString      = "a) U"
  res_u@gsnRightStringFontHeightF  = 0.015
  res_u@gsnRightString       = "m/s"
  res_u@gsnLeftString        = ""
  res_u@tiYAxisString        = "Height (m)" 
  res_u@tiXAxisString        = "Longitude"

  res_u@tmYLMode             = "Explicit"
  res_u@tmYLValues           = ilev(5:25:5) 
  res_u@tmYLLabels           = ispan(10000,2000,2000)

  res_u@tmYLMinorOn          = True
  res_u@tmYLMinorValues      = ilev

  plot(0) = gsn_csm_contour(wks, u, res_u)
  
;  res_zs  = True
;  res_zs@trYMinF   = 0
;  res_zs@trYMaxF   = 12000
;  res_zs@gsnYRefLine = 0.0
;  res_zs@gsnAboveYRefLineColor = "White"
;  plot = gsn_csm_xy(wks, lon, zs(0,90,:), res_zs)
;  overlay(plot(0), plot)
;  draw(plot(0))
  ;------------
  df1   = addfile(fi1,"r")
  omega = df1->wp(24,:,90,:)

  res_omega                   = True
  res_omega@gsnDraw           = False
  res_omega@gsnFrame          = False
  res_omega@gsnAddCyclic      = True
  res_omega@cnFillOn          = True
  res_omega@cnLinesOn         = False
  res_omega@lbLabelAutoStride = True

  res_omega@cnLevelSelectionMode = "ManualLevels"
  res_omega@cnLevelSpacingF      =  0.0001
  res_omega@cnMinLevelValF       = -0.0003
  res_omega@cnMaxLevelValF       =  0.0003
  res_omega@trYReverse           = True
  res_omega@vpWidthF             = 0.6
  res_omega@vpHeightF            = 0.3
  res_omega@gsnCenterStringFontHeightF = 0.015
  res_omega@gsnCenterString      = "b) OMEGA"
  res_omega@gsnRightStringFontHeightF  = 0.015
  res_omega@gsnRightString       = "Pa/s"
  res_omega@gsnLeftString        = ""
  res_omega@tiYAxisString        = "Height (m)"
  res_omega@tiXAxisString        = "Longitude"
  
  res_omega@tmYLMode             = "Explicit"
  res_omega@tmYLValues           = ilev(5:25:5)
  res_omega@tmYLLabels           = ispan(10000,2000,2000)
  res_omega@tmYLMinorOn          = True
  res_omega@tmYLMinorValues      = ilev

  plot(1) = gsn_csm_contour(wks, omega, res_omega)
  ;------------
  res                            = True
  res@gsnFrame                   = False
;  res@gsnPanelYWhiteSpacePercent = 5
;  res@gsnPanelBottom             = 0.2
  res@gsnPanelMainFontHeightF     = 0.02 
  res@gsnPanelMainString         = "GMCORE, Test 200, t = 6 days, 30 level"
  gsn_panel(wks, plot, (/2,1/), res)
  
  frame(wks)
  system("convert -trim -density 300 " + out_name + ".eps" + " " + out_name + ".png")
end
