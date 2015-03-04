;+
; NAME:
;   xis_plot
;
; PURPOSE:
;   To plot xi(s)
;
; CALLING SEQUENCE:
;    IDL> xis_plot 
;
; INPUTS:
;
; REVISION HISTORY:
;   10-Sep-2014  v0.0.1     NPR
;-

;; eg. Fig. 6, Ross et al., 2007, MNRAS, 381, 573
readcol, 'k_output_jack_perl_full_newcor_v2_ed.dat', $
         log_s, s, xi, error, dd, dr, xi_ls, xi_ham, rr
log_xi_ls = alog10(xi_ls)

readcol, 'k_output_jack_perl_full_temp.dat', $
         log_s_temp, s_temp, xi_temp, error_temp, dd_temp, dr_temp, xi_ls_temp, xi_ham_temp, rr_temp
log_xi_ls_temp = alog10(xi_ls_temp)

;;
;; Colour Table
;; http://ham.space.umn.edu/johnd/ct/ct-names.html
clr_table =13
loadct, clr_table

;; Colours for clr_table =13
black      =   0
purple     =  32
deep_blue  =  48
blue       =  64
light_blue =  80
turquiose  = 128
green      = 150
yellow     = 210
orange     = 232
red        = 254

charsize=2.6
charthick=4.8
thick=3.8
xthick=4
ythick=4
XTICKLEN  = 0.03
YTICKLEN  = 0.03

;; positions...
xpos_min = 0.20
xpos_max = 0.98
ypos_min = 0.20
ypos_max = 0.98

;; x-ranges
xmin = -0.6  ;; Doing everything in log_space...
xmax =  2.3  ;; Doing everything in log_space...

;; y-ranges
ymin = -2.0  ;; Doing everything in log_space...
ymax =  2.3  ;; Doing everything in log_space...


;; 0 - circle, 3 - 5 pointed star, 4 - triangle, 8 - square
;
plotsym, 8, 2.4, /fill
plot, log_s, log_xi_ls, $
      psym=8, $
      xrange=[xmin, xmax], yrange=[ymin, ymax], xstyle=1, ystyle=1, $
      thick=thick/1.2, charsize=charsize*1.2, $
      xtitle='!6log!D10!N(s / Mpc h!E-1!N!6', $
      ytitle='!6lgo!D10!N(!7n(s))!6'

oplot, log_s, log_xi_ls, thick=thick, linestyle=0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; FULL N(z)'s
;;
set_plot, 'ps'
device, filename='xis_plot_temp.eps', $
        xsize=8.0, ysize=8.0, $
        xoffset=0.2, yoffset=0.2, $
        /inches, /color, /encapsulated

plotsym, 8, 2.4, /fill
plot, log_s, log_xi_ls, $
      psym=8, $
      position=[xpos_min, ypos_min, xpos_max, ypos_max], $
      xrange=[xmin, xmax], $
      yrange=[ymin, ymax], $
      xstyle=1, ystyle=1, $
      xthick=xthick, ythick=ythick, $
      XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
      charsize=charsize, $
      charthick=charthick, $
      thick=thick,$ 
      /nodata, $
      xtitle='!6log!D10!N(s / Mpc h!E-1!N!6)', $
      ytitle='!6log!D10!N(!7n!6(s))!6', $
      color=black

oplot, log_s, log_xi_ls, psym=8, color=red
oplot, log_s, log_xi_ls, thick=thick, linestyle=0, color=red

;; The temp file you've got going on...
oplot, log_s_temp, log_xi_ls_temp, thick=thick, linestyle=0, color=black


xyouts, xmax*0.40, ymax*0.80, '2SLAQ LRGs',  charsize=charsize, charthick=charthick, color=color
;xyouts, xmax*0.80, ymax*0.80, No_of_words,  charsize=charsize, charthick=charthick, color=color

legend, pos=[-22.2, -7.5], ' ',   box=0, thick=14, linestyle = 0, charsize=1.2
xyouts,      -25.0, -7.8,  'PLE', charsize=2.2, charthick=8.


device, /close
set_plot, 'X'     



end
