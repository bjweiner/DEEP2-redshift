
c  use a set of previously computed eigentemplates
c  to find redshifts

      program pcredshift

c  read the eigentemplates as FITS files

c  read the spectrum as FITS file

c knock out regions around sky lines?

c rebin into log space if necessary

c do weights as fn of wavelength in template-wl space?

c compute a cross-correlation of the spectrum
c with each eigentemplate

c compute chi-squared

c  find the chi-squared minimum

c  output redshift and PCA-compressed spectrum
c  (coeffs of lin comb of eigentemplates)

c  lather, rinse, repeat


      end
