# Research
This will contain the code for Detecting Binaries Via Cross Correlation Function Subtraction 

Below is the IDL formatted code thus far:

PRO test_cc_subtraction, spectra1file, spectra2file
cd, '/Volumes/coveydata/APOGEE_Spectra'
;read in spectra1, assuming I'm starting from the APOGEE_Spectra directory.
apload, 'APOGEE_DR12_ApStar/'+spectra1file, spectra1
HELP, spectra1.rv
;spectra1 = MRDFITS(spectra1file, 2)
;make my life a **lot** easier by stripping off the combined CCFs
all_ccfs = N_ELEMENTS(spectra1.rv.ccf [0,*])
n_visits = all_ccfs-2
print, n_visits

;set guess for APOGEE_resolution
resolution = 3.
cspeed = 2.99792458d5
;set factor for converting lag space to velocities. ccfdw = 6e-6
ccfdw = spectra1.rv.ccfdw

---------------------------------------------------------------------------------

Code for plotting the ccf subtraction

make a loop that runs through the first n-2 visits anddifferences the CCFs 
from the last (which is n-1 since we start with 0, not 1).

start the loop, which goes from 0 to n-2
FOR i=0,n_visits-2 DO BEGIN
This next FOR loop will need to take each component and subtract it from the last so that it counts down.
 i.e : 0=15, 1=14, 2=13,......,9=6, ..., 13=1
  
  ;Generalize the saving of files 
  file_name = string('2M03430679+3148204,')
  
  This will store the iterations from a given i into folders in a directory Reyna
  FILE_MKDIR, '/Volumes/coveydata/APOGEE_Spectra/Reyna/' + file_name + strtrim(string(i),1)
  
    FOR j =i+1, n_visits-1 DO BEGIN

  difference the n-th CCF from the last (which is n_visits-1)

  this_diff = spectra1.rv.ccf[*,i] - spectra1.rv.ccf[*,j]

  plot the difference CCF along with the two that I started with to check if it looks the way I expect.
   PLOT, spectra1.rv.ccf[*,i], YRANGE = [-0.5 * MAX(spectra1.rv.ccf[*,i]), MAX(spectra1.rv.ccf[*,i])], /YSTY

   /buffer command allows the plots to not be shown

  plot1=plot(spectra1.rv.ccf[*,i], YRANGE = [-0.5 * MAX(spectra1.rv.ccf[*,i]), MAX(spectra1.rv.ccf[*,i])], /buffer, NAME='    CCF 1') (solid line)
          plot2=plot(spectra1.rv.ccf[*,n_visits-1], LINESTYLE=2 ,/overplot, /buffer, NAME='--- CCF 2') ; plot last CCF as dashed line

  Label the visit and run number (number of times it goes through the loop)
          plot3=plot(this_diff, /overplot, /buffer, 'r2', NAME='CCF Diff',TITLE = ' FITS NAME: '+ file_name + ' '+ 'Visit: '+ '  ' + strtrim(string(i),1) + '   '+ 'Iteration: ' + strtrim(string(j),1)) 


          
          ;Creating a legend

          leg=LEGEND(TARGET=[plot1, plot2, plot3], POSITION = [0.9,0.8], sample_width = 0, /NORMAL, /AUTO_TEXT_COLOR) 

_______________________________________________________________________

Detecting the change. Need to integrate over the residuals and divide by the "window of interest" 

Finding the RMS.
These are the elements that are equivalent to N
n = 75.
m = 80.


This array will take the sum of the first section of 'wiggles'
x = spectra1.rv.ccf [0.0:75.0]
This array will take the sum of the second section of 'wiggles' 
y = spectra1.rv.ccf[375.0:400.0]
;Testing if input of math effects the reading of code
z = (this_diff)^2.

Sigma_1 = SQRT((1/n)*TOTAL(x ^ 2.))
Sig = SQRT((1/n)*TOTAL(y ^ 2.))
sums = TSUM( this_diff^2.0 / Sig)
sums2 = TSUM( this_diff / Sig)

print, 'Outer Integration: ' , sums
print, ' Outer Integ w/o square: ', sums2


Sigma_2 = SQRT((1/m)*TOTAL(x^2.))

 
 Now we need to find the summation of the residuals divided by the RMS: Try squaring and not squaring the residuals.
 
  integrated_res2 = TSUM( z / Sigma_1)
  integrated_res1 = TSUM( this_diff / Sigma_1)
 
 
 Now, let's check the other possibilities. What if we don't divide by Sigma. Try squaring and not squaring the residuals.
  integrated_resid1 = TSUM( this_diff)
  integrated_resid2 = TSUM( this_diff^2)
 
 Using Sigma_2

  integ_res2 = TSUM( z / Sigma_2)
  integ_res1 = TSUM( this_diff / Sigma_2)
 
  print, 'IRNS&D2 by RMS: ' , integ_res2 
  print, 'IRS&D2 by RMS: ' , integ_res1 
 

  print, 'IRNS&D by RMS:', integrated_res1
  print,  'IRS&D by RMS:', integrated_res2


  print, 'IRNS:',integrated_resid1         
  print, 'IRS:',integrated_resid2

      print,'Visit number:' , i
      ;Date of the visits
      
      print, 'Visit i:' , spectra1.rv.jd[i]-2.45D6 , '   Visit j:' ,  spectra1.rv.jd[j]-2.45D6

Call the area under the residuals and print the value of them:
      integral = TSUM(this_diff)
      print,'Integral: ', integral

Call the squared of the area under the residual to see the differences summed
      integral_squared = TSUM(this_diff ^ 2.0)
      print, 'Squared:', integral_squared


Call the max value/ possible index of the result
      maximum = MAX(this_diff)
      print, 'Maximum: ' ,maximum
 __________________________________________________________________________________________________________________________
   XYOUTS for plotting important values

Integ = TEXT( 25.,0.4, 'Integration : '+ STRMID(STRTRIM(STRING(integral),2),0,10), /DATA, FONT_SIZE = 10, FONT_NAME = 'Helvetica')
IntegS = TEXT( 25,0.37, '$( Integrat )^2$: '+ STRMID(STRTRIM(STRING(integral_squared),2),0,10), /DATA, FONT_SIZE = 10, FONT_NAME = 'Helvetica')
Maxs = TEXT( 25.0,0.34, 'Max: '+ STRMID(STRTRIM(STRING(MAX(this_diff)),2),0,6), /DATA, FONT_SIZE = 10, FONT_NAME = 'Helvetica')
RootMS = TEXT( 25.0,0.31, 'Rms: '+ STRMID(STRTRIM(STRING(Sigma_1),2),0,6), /DATA, FONT_SIZE = 10, FONT_NAME = 'Helvetica')
Integ_resd2 = TEXT(25.0,0.28, 'diff/$\sigma$:'+ STRMID(STRTRIM(STRING(integ_res1),2),0,6), /DATA, FONT_SIZE = 10, FONT_NAME = 'Helvetica')
Integ_resd1 = TEXT(25.0,0.25, '$diff^2$/$\sigma$:'+ STRMID(STRTRIM(STRING(integ_res2),2),0,6), /DATA, FONT_SIZE = 10, FONT_NAME = 'Helvetica')

________________________________________________________________________

;Saving individual plots to Reyna Folder as titles of data
plot1.save, '/Volumes/coveydata/APOGEE_Spectra/Reyna/'+file_name + strtrim(string(i),1) + '/CCF_Diff' + strtrim(string(j),1)+'.png'
----------------------------------------------------------------------
     Code for making a CSV file
    
  CLOSE, 2
  OPENW, 2, 'CCF_subtraction_values.tex'

  columnalign = 'cccccccccccccc'
  columnnames = ['2MASS ID', 'Integration', 'Integration Squared', $
    'Maximum', 'RMS', $
      'Integ. not squared/RMS ', 'Integ. squared/RMS', $
      'Confirmed Binary?', 'Sigma 1', $
      'Sigma 2']

    tabletitle = 'CCF Subtraction Values'
    tablelabel = 'tab:Subtraction Values'
   , columnnames, columnalign, tabletitle, tablelabel, '\tiny', '0', 2, /LANDSCAPE

  FOR i=0,n_visits-1 DO BEGIN
 PRINTF, 2, FORMAT = '(A20, 2x,A1, 2x, D13.8,2x,A1,2x,D13.8,2x,A1,2x,F6.2,2x,A1,2x,A10,2x,A1,2x,F6.2,2x,A1,2x,F8.4,2x,A1,2x,A10,2x,A1,2x,F6.2,2x,A1,2x,F6.2,2x,A1,2x,I1,2x,A1,2x,F5.2,2x,A1,2x,F6.2,2x,A1,2x,F7.2,2x,A1,2x,I3,2x,A1,2x,F7.2,2x,A2)', spectra1.rv.ccf[*,i], ',',$
       spectra1.rv.ccf[*,i], ',', this_diff[*,i], ',', $
       
   ENDFOR 
close the loop (so let the loop go back to the top.)
ENDFOR
END  
