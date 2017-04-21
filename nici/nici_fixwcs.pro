;+
; NAME:
;	NICI_FIXWCS
;
; PURPOSE:
;	This procedure corrects the WCS in NICI FITS files
;       for IAA errors, and sets the 1.1 deg rotation between
;       holmes and watson channels.  The required corrections are
;       stored in a lookup table and selected automatically
;       based on the observation date.
;
; CATEGORY:
;	NICI
;
; CALLING SEQUENCE:
;	NICI_FIXWCS, Prefix, N1
;
; INPUTS:
;       Prefix     : String containing DHS prefix
;       N1         : Integer representing first file number
;
; OPTIONAL INPUTS: 
;       N2         : Integer representing last file number
;
; KEYWORD PARAMETERS:
;       SUFFIX     : Suffix to add to name of output file.
;                    Default is no suffix, meaning that input 
;                    file will be overwritten.
;       NOWRITEFITS: Set to inhibit writing of output files.
;       SIMDATE    : String with date; used for testing.
;       HELP       : Print summary of command arguments.
;
; OUTPUTS:
;	The output MEF files have the WCS records in the
;       holmes and watson extensions corrected.
;
; RESTRICTIONS:
;       Test version, only partially tested at this time.
;
;       The supplied LUT's are subject to change at any time
;       as astrometry is updated and NICI is changed from
;       port to port!  Make sure to use the most recent version
;       of NICI_FIXWCS!
;
; ROUTINES USED:
;       From the IDL astrometry library: 
;       fits_open, fits_read, fits_write, sxpar, sxaddpar
;       extast to extract FITS astrometry records
;       The rotation matrix algorithm is taken from 
;       the hrot procedure.
;
; METHOD:
;       - Reads the IAA from the primary header, and finds
;         the "true" IAA in the LUT in the code.
;
;       - For some ranges of dates, the holmes WCS was mistakenly
;         set so it does not correspond to the header IAA.  In
;         these cases, the IAA represented by the WCS is supplied
;         in a separate LUT.
;
;       - Corrects the holmes extension WCS based on the
;         IAA difference. 
;
;       - If necessary, recomputes the watson WCS to be
;         rotated 1.1 degree from holmes.
;
;       - All corrections are made solely to the WCS.
;         Image data are not rotated or modified in any way.
;
; EXAMPLES:
;       nici_fixwcs, 'S20101031S', 1
;       - Processes the file S20101031S0001.fits.
;       - Input file is overwritten with updated headers.
;
;	nici_fixwcs, 'S20101031S', 1, 10
;       - Processes files S20101031S0001 - S0010.fits
;
;	nici_fixwcs, 'S20101031S', 1, 10, suffix='_1'
;       - Instead of overwriting, new files are created with
;         suffix '_1.fits'
;
;       nici_fixwcs, 'S20101031S', 1, 10, /nowritefits
;       - Diagnostic messages are printed to show the WCS
;         corrections, but the output files are not written.
;
; MODIFICATION HISTORY
;       2010 Nov  8: Original NICI_FIXIAA, T. Hayward
;       2010 Nov 18: test version
;       2010 Dec 03: removed forcwcs, added iaawcsh param.
;                    corrected watson rotation sign error for bottom port
;       2010 Dec 21: Corrected IAA parsing error
;                    Updated watson relative rotation to 1.100 deg.
;
;       2011 Aug 20: NICI_FIXIAA2 created from NICI_FIXIAA to add 
;                    date LUT. (TLH)
;       2012 Jan 23: Bug fixes, initial release (TLH).
;
;       2012 Jan 26: NICI_FIXWCS adapted from NICI_FIXIAA2.
;                    - Add loop over multiple files.  
;                    - Fixed bug in drot computation.
;                    - Improved messages.
;
;       2012 Jan 27: Fixed bug with selection of last date in LUT.
;
;       2012 Jan 28: Initial CD matrix now saved as CDX_Y_O keywords.
;                    Added NOWRITEFITS keyword to inhibit writing of
;                    output files.
;                    Write output files only if WCS has been updated.
;
;       2013 Jan 21: Significant rewrite to compute WCS from first
;                    principles instead of rotating the existing
;                    matrix, which allows incorporation of mean
;                    holmes & watson pixel scales.  The change in
;                    scale requires that the WCS _always_ be
;                    regenerated even if the rotation angle is
;                    identical.  Updated LUT's through 
;                    mounting 10, 2012-07-13. 
;
;       2013 Nov 12: Updated LUT's for mountings 11 and 12 through
;                    2013 June.
;-

function r180, theta
   while theta GT  180.D0 do theta -= 360.D0
   while theta LT -180.D0 do theta += 360.D0
   return, theta
end

PRO NICI_FIXWCS, prefix, n1, n2, suffix=suffix, simdate=simdate, $
   gzip=gzip, nowritefits=nowritefits, debug=debug, help=help

   on_error, 2

   if (n_params() LT 2 OR keyword_set(Help)) then begin
      print, 'Syntax: NICI_FIXWCS, prefix, N1 [, N2, SUFFIX=, SIMDATE=, '
      print, '          /GZIP, /NOWRITEFITS, /DEBUG, /HELP]'
      return
   endif

   if n_elements(n2) eq 0 then n2 = n1
   if n_elements(suffix) eq 0 then suffix = ''
   nowritefits = keyword_set(nowritefits)
   wcsmodified = 0
   cdr   = !DPI/180.0D

   ;; Initialize NICI remounting lookup tables.
   ;;
   ;; dates    : Dates of NICI mountings or remountings on the telescope.
   ;;
   ;; iaas_init: IAA value initially measured for a particular mounting.
   ;;            The header IAA should equal this value for dates after 
   ;;            the IAA was measured.  The header IAA may be different 
   ;;            for dates immediately after the mounting date, before
   ;;            IAA was measured and updated.
   ;;
   ;; iaas_true: The final "true" IAA value determined for the mounting
   ;;            after accounting for all known errors.  Unless
   ;;            iaas_wcs>0, the holmes WCS will be rotated by the 
   ;;            difference between the header IAA and iaa_true.
   ;;
   ;; iaas_wcs : If < 0, the holmes WCS is assumed to correspond to
   ;;            the header IAA.
   ;;            If > 0, equals the IAA that the holmes WCS represents,
   ;;            instead of the header IAA.  In this case the holmes 
   ;;            WCS will be rotated by the difference between
   ;;            iaas_wcs and iaas_true.  Used to correct one error on
   ;;            2010-11-02 when the IAA was set to 112.20, but the
   ;;            holmes WCS was not updated to correspond to the IAA.

   mounting  = [         1,          2,          3,          4,          4,          4,          4,          5,          6,          7,          8,          9,         10,         11,         12]
   dates     = ['20080727', '20091019', '20100331', '20101020', '20101102', '20101110', '20101214', '20110114', '20110311', '20110427', '20110608', '20120215', '20120713', '20130507', '20130619']
   epochs    = [ 2008.5695,  2009.7967,  2010.2437,  2010.7995,  2010.8350,  2010.8569,  2010.9509,  2011.0356,  2011.1889,  2011.3176,  2011.4326,  2012.1232,  2012.5311,  2013.3450,  2013.4627]
   ports     = [         5,          1,          5,          5,          5,          5,          5,          1,          1,          1,          5,          5,          5,          5,         1 ]
   iaas_init = [    112.30,     247.50,     111.70,     111.32,     112.20,     112.20,     112.61,     247.50,     247.50,     247.50,     112.82,     112.64,     112.42,     112.25,    247.50 ]
   iaas_true = [    112.30,     247.50,     112.20,     112.61,     112.61,     112.61,     112.61,     247.50,     247.50,     247.50,     112.82,     112.64,     112.42,     112.25,    247.50 ]
   iaas_wcs  = [     -1.00,      -1.00,      -1.00,      -1.00,     111.32,      -1.00,      -1.00,      -1.00,      -1.00,      -1.00,      -1.00,      -1.00,      -1.00,      -1.00,     -1.00 ]
   lastvaliddate = '20130901'

   ; Holmes and Watson mean scales, arcsec/pixel
   cdelt_h = 0.017958D
   cdelt_w = cdelt_h * 1.002D

   ;; Print Version Message
   message, 'Version 2013 Jan 21', /info

   ;; Loop through files with suffix n1 to n2
   for n = n1, n2 do begin

      ;; Construct the filenames
      infilename  = string(format='(%"%s%04d.fits")'  , prefix, n)
      if keyword_set(gzip) then infilename = infilename + '.gz'
      outfilename = string(format='(%"%s%04d%s.fits")', prefix, n, suffix)
      message, infilename, /inf

      ;; Open the MEF file
      fits_open, infilename, fcb, update=0, message=message
      if (message NE '') then return

      ;; If the Open was successful, read the primary header
      fits_read, fcb, im, hdr0, message=message, exten_no=0, /header_only
      if (message NE '') then begin
         message, message, /cont
         fits_close, fcb
         return
      endif

      ;; Get the HDR0 parameters
      instrume = strupcase(strtrim(sxpar(hdr0, 'INSTRUME'), 2)) ; remove lead & trail blanks
      inport   = sxpar(hdr0, 'INPORT')
      iaa_hdr  = sxpar(hdr0, 'IAA')
      iaa_wcsh = sxpar(hdr0, 'IAAWCS', count=counth)
      if n_elements(simdate) GT 0 then dateobs = simdate $
      else dateobs = strjoin(strsplit(sxpar(hdr0, 'DATE-OBS'), '-', /extract), '')

      ;; Proceed only if a NICI file
      if instrume NE 'NICI' then begin
         message, 'Error: Not a NICI file', /cont
         fits_close, fcb
         return
      endif

      ;; Check ISS port
      if      (inport EQ 1) then dirsign = -1. $
      else if (inport EQ 5) then dirsign =  1. $
      else begin
         message, 'Error: unknown port', /inf
         fits_close, fcb
         return
      endelse

      ;; Check that Date of Observation is within valid range
      if dateobs GT lastvaliddate then begin
         message, 'Error: File date is after last valid date of ' + lastvaliddate, /inf
         fits_close, fcb
         return
      endif

      ;; Check for IAAWCS keyword
      if counth GT 0 then begin
         message, string(format='(%"Found IAAWCS = %7.2f from previous run of FIXWCS")', $ 
                      iaa_wcsh), /inf
      endif

      ;; Read hdr and image data from each extension
      fits_read, fcb, im_h, hdr_h, message=message, exten_no=1, /no_pdu
      fits_read, fcb, im_w, hdr_w, message=message, exten_no=2, /no_pdu
      fits_close, fcb
      if (message NE '') then begin
         message, infilename + ': ' + message, /cont
         return
      endif

      ;; Set new IAA depending on date
      nlast = n_elements(dates)-1
      for i = 0, nlast do begin
         if i EQ nlast then break
         if dates[i+1] GT dateobs then break
      endfor
      ; print, 'Selected: ', i, ' ', dates[i]

      if counth gt 0 then begin
         iaa1 = iaa_wcsh    ; Use value from previous run of fixwcs
      endif else begin
         if iaas_wcs[i] GT 0 then iaa1=iaas_wcs[i] else iaa1=iaa_hdr
      endelse
      iaa2 = iaas_true[i]

      message, string(format='(%"Date = %s  IAAs: HDR = %7.2f, WCS_H = %7.2f, True = %7.2f")',$ 
                   dateobs, iaa_hdr, iaa1, iaa2), /inf

      ;; Compute angle difference in radians
      dtheta_h = dirsign * (iaa2 - iaa1)
      msg = string(format='(%"IAA = %7.2f, should be %7.2f, error = %7.2f deg")',$
                      iaa1, iaa2, dtheta_h)
      message, msg, /inf

      ;; Extract existing WCS's
      extast, hdr_h, astrom_h1
      extast, hdr_w, astrom_w1
      getrot, astrom_h1, theta_h1, cdelt_h1, /silent
      getrot, astrom_w1, theta_w1, cdelt_w1, /silent

      ;; Construct new holmes and watson WCS's with new scale and rotation.
      ;; Initially both coordinate systems will be left-handed.
      theta_h2 = theta_h1 + dtheta_h
      th = theta_h2 * cdr
      cd_h2 = (cdelt_h / 3600.)*[ [-cos(th),-sin(th)],[-sin(th), cos(th)] ]

      theta_w2 = theta_h2 + dirsign * 1.1
      th = theta_w2 * cdr
      cd_w2 = (cdelt_w / 3600.)*[ [-cos(th),-sin(th)],[-sin(th), cos(th)] ]

      ;; Reflect one WCS L-R depending on port.
      ;; Taken from hrotate.pro for dir=5
      rot_mat_trans = transpose([ [-1,0], [0,-1] ])
      if inport EQ 5 then begin
         cd_h2 = cd_h2 # rot_mat_trans
         cd_h2[*,1] = -cd_h2[*,1]
      endif else if inport EQ 1 then begin
         cd_w2 = cd_w2 # rot_mat_trans
         cd_w2[*,1] = -cd_w2[*,1]
      endif

      ;; Update headers
      sxaddpar, hdr0,  'IAAWCS' , iaa2

      sxaddpar, hdr_h, 'CD1_1_O', astrom_h1.cd[0,0]
      sxaddpar, hdr_h, 'CD1_2_O', astrom_h1.cd[0,1]
      sxaddpar, hdr_h, 'CD2_1_O', astrom_h1.cd[1,0]
      sxaddpar, hdr_h, 'CD2_2_O', astrom_h1.cd[1,1]
      sxaddpar, hdr_h, 'CD1_1'  , cd_h2[0,0]
      sxaddpar, hdr_h, 'CD1_2'  , cd_h2[0,1]
      sxaddpar, hdr_h, 'CD2_1'  , cd_h2[1,0]
      sxaddpar, hdr_h, 'CD2_2'  , cd_h2[1,1]
      sxaddpar, hdr_h, 'HISTORY', 'CD records updated by NICI_FIXWCS'
      sxaddpar, hdr_h, 'HISTORY', 'Original CD saved as CDX_Y_O'

      sxaddpar, hdr_w, 'CD1_1_O', astrom_w1.cd[0,0]
      sxaddpar, hdr_w, 'CD1_2_O', astrom_w1.cd[0,1]
      sxaddpar, hdr_w, 'CD2_1_O', astrom_w1.cd[1,0]
      sxaddpar, hdr_w, 'CD2_2_O', astrom_w1.cd[1,1]
      sxaddpar, hdr_w, 'CD1_1'  , cd_w2[0,0]
      sxaddpar, hdr_w, 'CD1_2'  , cd_w2[0,1]
      sxaddpar, hdr_w, 'CD2_1'  , cd_w2[1,0]
      sxaddpar, hdr_w, 'CD2_2'  , cd_w2[1,1]
      sxaddpar, hdr_w, 'HISTORY', 'CD records updated by NICI_FIXWCS'
      sxaddpar, hdr_w, 'HISTORY', 'Original CD saved as CDX_Y_O'

      wcsmodified = 1

      getrot, hdr_h, theta_h2, cdelt_h2, /silent
      getrot, hdr_w, theta_w2, cdelt_w2, /silent

      theta_h1 = r180(theta_h1)
      theta_w1 = r180(theta_w1)
      theta_h2 = r180(theta_h2)
      theta_w2 = r180(theta_w2)
      drot1 = -theta_w1 - theta_h1
      drot2 = -theta_w2 - theta_h2

      ;; Print the before and after CD matrices and CDELT values
      ;; if the DEBUG keyword is set.
      if (keyword_set(debug)) then begin
         handed = ['left-handed', 'right-handed']

         print, ''
         print, '---HOLMES---'
         print, 'CD1: ', astrom_h1.cd[*]
         print, 'CD2: ', cd_h2[*]
         print, 'CDELT1: ', cdelt_h1[0]*3600., cdelt_h1[1]*3600., ': ', handed[cdelt_h1[0] GT 0]
         print, 'CDELT2: ', cdelt_h2[0]*3600., cdelt_h2[1]*3600., ': ', handed[cdelt_h2[0] GT 0]
         print, ''

         print, '---WATSON---'
         print, 'CD1: ', astrom_w1.cd[*]
         print, 'CD2: ', cd_w2[*]
         print, 'CDELT1; ', cdelt_w1[0]*3600., cdelt_w1[1]*3600., ': ', handed[cdelt_w1[0] GT 0]
         print, 'CDELT2: ', cdelt_w2[0]*3600., cdelt_w2[1]*3600., ': ', handed[cdelt_w2[0] GT 0]
         print, ''

      endif

      message, string(format='(%"holmes WCS Rot = %7.2f -> %7.2f, rotated by %7.2f deg")',$
                       theta_h1, theta_h2, theta_h2-theta_h1), /info
      message, string(format='(%"watson WCS Rot = %7.2f -> %7.2f, rotated by %7.2f deg")',$
                       theta_w1, theta_w2, theta_w2-theta_w1), /info
      message, string(format='(%"h-w   WCS dRot = %7.2f -> %7.2f deg")', $
                       drot1, drot2), /info

      ;; Write output file
      if wcsmodified and not nowritefits then begin
         fits_open, outfilename, outfcb, /write
         fits_write, outfcb, 0, hdr0
         fits_write, outfcb, im_h, hdr_h
         fits_write, outfcb, im_w, hdr_w
         fits_close, outfcb
         message, string(format='(%"Wrote to %s")', outfilename), /inf
      endif else begin
         message, 'No output file', /inf
      endelse
      print, ''

   endfor
   
end
