program optimal_phot

!=============================================================================== 
! VERSION FOR RTPhoS - No User INPUT, all input controlled from run_RTPhos.py
!===============================================================================
! 
! Program to calculate the flux of stars in CCD images using Tim Naylor's optimal 
! extraction routines. For a thorough explanation of the algorithms read the 
! following two papers:
!
! "An Optimal Extraction Algorithm for Imaging Photometry"
! Naylor Tim, MNRAS, 296, 399, 1998
!
! "Optimal Photometry for Colour-Magnitude Diagrams and its Application on NGC 2547"
! Naylor, Tim et al. MNRAS, 335, 291, 2002
!
! Requirements:
! opt_extr.f90 - All the required routines for optimal and aperature photometry
! marq.f90     - Computes Marquardt routines needed for profile fitting
! subs.f90     - Various array and file sorting routines
! cfitsio      - NASA's FITS file manipulation library needs to be installed.
!
! Compiles with:
! edit optimal.csh and source it. (Edit needed to provide path for cfitsio library)
!
! Alternatively (I use g95 but any Fortran compiler will do):
! g95 -c opt_extr.f90 marq.f90
! g95 -o optimal optimal.f90 marq.o opt_extr.o -L/<Local Path>/cfitsio -lcfitsio
!
! Required improvements:
! 1. Include flagging of bad pixels, bad stars, bad sky etc. This feature was removed from
!    The original code by Tim in order to make the code independent of the ARK
!    software package. 

  use opt_extr

  implicit none

!----------------------------------------------------------------------------

! Variable declarations
! ------------------------------------------------------------------------------

! Various counters & test variables
  integer :: iunit1, iostat, eof1, i, j, rows, cols, blines, fnum, num, ounit
  integer :: progcnt, stepcnt

  real :: searchrad

  logical :: there, verbose, filelist

  character (len=1)  :: verbosein, filelistin
  character (len=3)  :: suf  
  character (len=80) :: wformat
  character (len=80) :: dummy
  character (len=80) :: input(13)

! Variables for opening FITS files and reading the image data
  integer :: status, unit, readwrite, blocksize, naxis, group, nfound
  integer :: naxes(2),fpix(2), lpix(2), inc(2)
  integer :: low(2), high(2)

  real :: datamin, datamax, nullval
  real, allocatable :: data(:,:) ! **** This is the image data array

  logical :: anynul

  character (len=50) :: filename
  character, allocatable, dimension(:,:) :: pix_flg ! **** This is the image flag array

! Variables regarding PSF stars
  integer :: npsf, ipsf, nfit

  real, dimension(3) :: shape_par !The 3 params defining the shape of the PSF 
  real, allocatable, dimension(:) :: dpsf, xpos0, ypos0  
  real, allocatable, dimension(:,:,:) :: newpsfpos

  double precision, allocatable, dimension(:,:) :: psfstars
  
  character (len=50) :: psfpos
  character (len=15), allocatable, dimension(:) :: psfnames

! Variables used for clipping the PSF star
  real :: clip_fwhm, cliprad, fwhm
  
  character(len=10) :: clipanswer, distanswer 

! Variables regarding Target stars
  integer :: istar, nstar, nframes, ntimes 
     
  logical :: comp, iftimes
  real :: xcomp, ycomp
 
  real :: xpos, ypos, badpix, posfix
  real :: optflux, opterror, optnrm, peak, xfit, yfit, xerr, yerr
  real, allocatable, dimension(:) :: dpos, seeing
  real, allocatable, dimension(:,:) :: offsets
  real, allocatable, dimension(:,:,:) :: optres, apres, newxypos

  double precision, allocatable, dimension(:,:) :: stars, times
    
  character (len=50) :: starpos, offsetfile, timesfile

  character (len=15), allocatable, dimension(:) :: starnames
  character (len=50), allocatable, dimension(:) :: datafiles, timefiles

! Varialbes used for sky estimation
  real :: bad_sky_chi, bad_sky_skw, skynos, skycnt, skyerr
  
  character :: cflag
  character(len=10) :: bskyanswer 

  logical :: sky_stuffed

! Variables for aperature photometry
  integer :: ibox

  real :: aprad, adu
  real :: apx, apy
  real :: apflux, aperror

  character(len=10) :: apanswer

  logical :: aperature

! Variables for optimal photometry
  integer :: iopt  ! The ID# of the star whose values are to be optimized

  character(len=10) :: optanswer 

  logical :: optimal

! ***************************************************************
! Actual program starts here.
! ***************************************************************

! Read input parameters
  read*, filename, psfpos, starpos, verbosein
  read*, bad_sky_skw, bad_sky_chi, fwhm, clip_fwhm, aprad, iopt, searchrad, adu

! J1753.fits psf.cat stars.cat N
! -1 -1 2 5 10 1 5 6.2

! Initialize variables, set global constants

  rows = 0    ! Number of rows in a text file    (used in subs.f90)
  cols = 0    ! Number of columns in a text file (used in subs.f90)
  blines = 3  ! Number of top lines to ignore in a text file (used in subs.f90)
  npsf = 0    ! Number of PSF stars
  status = 0  ! Required by CFITSIO

  optimal=.true.
  aperature=.true.
!  verbosein = "N" 

! Set the program output type
  verbose = .False.
  if (verbosein=='y'.or.verbosein=='Y') verbose=.True.

! ***************************************************************
! Read the co-ordinates for each PSF star.
! ***************************************************************

! Open the file with the position(s) of the psf star(s) and
! read the positions of the psf star(s). 
! The format of the input PSF star file should have 3 comment or
! blank lines at the top and then 3 columns:
! ID#, X-coordinate, Y-Coordinate

! Check to see that the file exists
  inquire(file=psfpos, exist=there)
  if (.not. there) then
     print*, '* PSF Positions File Not Found.'
     print*, '* Proceeding with Aperature Photometry only!'
     print*
     optimal=.False.
     aperature=.True.
     clip_fwhm = -1.0
  end if 
  if (optimal) then
     call ftgiou(iunit1, status)            
     open (unit=iunit1, file=psfpos, status='old', access='sequential')
     call xsize(iunit1, cols, rows, blines)
     if (verbose) print*, "* Found ", rows, "PSF star location(s): "
     npsf = rows
     allocate(psfstars(cols,rows))
     allocate(psfnames(rows))
     rewind(iunit1)
     call readit(iunit1, cols, rows, blines, psfstars, psfnames)
     close(iunit1)
     call ftfiou(iunit1, status)
     if (verbose) then
        do i=1, rows
           print*, "* PSF Star ",psfnames(i)," at position (X,Y):", &
                   real(psfstars(2,i)), real(psfstars(3,i))    
        end do
     end if 
     if (verbose) print*
  end if

  allocate(xpos0(npsf))
  allocate(ypos0(npsf))
  allocate(dpsf(npsf))

! *******************************************************************
! Read the co-ordinates for each target star.
! *******************************************************************

! Open the file with the position(s) of the target star(s) and
! read the positions of the target star(s). 
! The format of the input file should have 3 comment or
! blank lines at the top and then 3 columns:
! ID#, X-coordinate, Y-Coordinate  

! Initialize the file opening variables
  rows = 0
  cols = 0
  blines = 3 

! Check to see that the file exists.
  inquire(file=starpos, exist=there)
  if (.not. there) then
     print*, '* WARNING Target Star Positions Not Found. Nothing to do. Exiting....'
     stop
  end if

  call ftgiou(iunit1, status)           
  open (unit=iunit1, file=starpos, status='old', access='sequential')
  call xsize(iunit1, cols, rows, blines)
  if (verbose) print*, "* Found ", rows, "star location(s): "
  nstar = rows
  allocate(stars(cols,rows))
  allocate(starnames(rows))
  allocate(dpos(nstar))
  rewind(iunit1)
  call readit(iunit1, cols, rows, blines, stars, starnames)
  close(iunit1)
  call ftfiou(iunit1, status)
  if (verbose) then
     do i=1, rows
        print*, "* Object ",starnames(i),"at position (X,Y):", &
                real(stars(2,i)), real(stars(3,i))    
     end do
  end if 

! Allocate all the arrays.
! Since this is the pipeline version it runs every time for a single file.
! Therefore ntimes must equal 1.
  ntimes = 1
  allocate(optres(nstar,2,ntimes))
  allocate(apres(nstar,2,ntimes))
  allocate(newxypos(nstar,2,ntimes))
  allocate(newpsfpos(npsf,2,ntimes))
  allocate(seeing(ntimes))

! ***************************************************************
! Setup all the optimal and aperature photometry options
! ***************************************************************

  ! If not set, set the clipping radius for the PSF.
  ! Most of the signal-to-noise is obtained by setting the clipping
  ! radius to be 2*fwhm.  But there is little gain in speed by
  ! making it less than 3 pixels.
  if (clip_fwhm < 0.0) then 
     optimal = .False.
     clip_fwhm=-1.0
  else
     if (clip_fwhm == 0.0 .or. clip_fwhm < fwhm) clip_fwhm=max(3.0/fwhm, 2.0)
  end if

  ! If not set, set the default aperture radius to the optimum (as shown in the paper),
  ! but don't let it be less than 2 pixels.
  if (aprad<0.0) then
     aperature = .False.
     aprad=-1.0
  else
     if (aprad == 0.0 .or. aprad < fwhm) aprad=max(2.0, 2.0*fwhm/3.0)
  end if

  ! Set the stage for the photometry to be optimized for a particular star or
  ! for a sky limited case.
  if (optimal) then
     optstar: do
              if (iopt > 0) then
                 do i=1, nstar
                    if (int(stars(1,i)) == iopt) then
                       iopt=i
                       exit optstar
                    end if
      	         end do
	             print*, "* WARNING: The star to optimize for is not in the target list."
                 print*, " (Optimal) Will proceed with sky limited optimization."
                 print*
                 iopt=-1
                 exit optstar
              else
                 iopt=-1
                 exit optstar
	      end if   
     end do optstar
  end if

! Set the search radius. If the radius is negative then the positions will be 
! considered as fixed. If not supplied or if the supplied one is less than 2*FWHM 
! then search radius is set to 2*FWHM
  dpsf=searchrad
  dpos=searchrad
  if (searchrad < 2.0*fwhm) then
     dpsf=2.0*fwhm
     dpos=2.0*fwhm
  end if
  if (searchrad <=0.0) then
     print*, "* WARNING: Search radius is zero or negative. No centroiding will be performed"
     print*, "*(Optimal) and the positions will be fixed to the input values."
     print*
  end if

! Summary of all the setup options.
  if (verbose) then
     print*, 'Photometry will be performed with the following options:'
     print*, '========================================================'
     print*, 'Input PSF file                      :',trim(psfpos)
     print*, 'Stars listed in the PSF file        :',npsf
     print*, 'Input target star file              :',trim(starpos)
     print*, 'Stars listed in the target star file:',nstar
     print*, 'Sky skew flag limit                 :',bad_sky_skw
     print*, 'Sky Chi^2 flat limit                :',bad_sky_chi
     print*, 'FWHM of image                       :',fwhm
     print*, 'PSF clipping radius                 :',clip_fwhm
     print*, 'Aperature photometry radius         :',aprad
     print*, 'Star to be optimized (-ve if all)   :',iopt
     print*, 'Centroid search radius              :',dpsf(1)
     print*, 'Detector gain (e-/ADU)              :',adu
     print*, '======================================================='
     print*
  end if

! If both the optimal and photometry options are negative exit the program.
  if (.not.optimal) print*, '* Optical Photometry will not be performed'
  if (.not.aperature) print*, '* Aperature Photometry will not be performed'
  if ((.not.optimal).and.(.not.aperature)) then
     print*, '* Nothing to do! Exiting...'
     stop
  end if

! ***************************************************************
! Begin the Photometry Processes
! ***************************************************************

  frame: do fnum=1, ntimes

         ! First read in the image data into a data array.
         ! Initialize CFITSIO required variables. 
         if (verbose) print*, "Filename: ", filename         
         nfound=0
         naxis=2
         group=1
         inc=1
         nullval=-999
         status=0
         readwrite=0

         call ftgiou(unit, status)
         call ftopen(unit, filename, readwrite, blocksize, status)
         call ftgknj(unit, 'NAXIS', 1, 2, naxes, nfound, status)
         if (nfound /= 2)then
            if (verbose) print *,'Could not read image. Moving to the next frame...'
            cycle frame
         end if
         ! Allocate the data array and determine the first and last pixels of the image.
         fpix=1
         lpix=naxes
	     low=fpix
         high=lpix
         allocate(data(1:naxes(1),1:naxes(2)))
         ! Read the data into the data array
         call ftgsve(unit, group, naxis, naxes, fpix, lpix, inc, nullval, data, anynul, status )
         ! Print somethings from the FITS file to confirm that it read properly.
!         print*, naxis, naxes(1), naxes(2)
!         print*, fpix(1), fpix(2)
!         print*, lpix(1), lpix(2)
!         print*, "array value at x=100, y=100: ", data(100,100)
!         print*, "maximum value position: ", maxloc(data)
! 	      print*, "*********************"
         call ftclos(unit, status)
         call ftfiou(unit, status)
         ! Check for any error, and if so print out error messages.
         ! The PRINTERROR subroutine is listed at the end of this file.
         if (status .gt. 0 .and. verbose) call printerror(status)

!        Setup a dummy flag array - It must be properly set later  ********* 
	     allocate(pix_flg(1:naxes(1),1:naxes(2)))
         pix_flg="O"

!        Estimate the position of the PSF stars
         if (clip_fwhm > 0.0) then
            xpos0 = real(psfstars(2,:))
            ypos0 = real(psfstars(3,:))

            do i=1, npsf              
               if (verbose) print*, 'Estimated position of PSF star ',i, 'is ', &
                            xpos0(i), ypos0(i)
            end do

            if (verbose) print*, "@ Starting FWHM and Clipping radius is: ", fwhm, clip_fwhm
!           Fit the PSF stars
            call psf_calc(data, pix_flg, npsf, xpos0, ypos0, dpsf, &
                          adu, high, low, fwhm, shape_par, ipsf, nfit, verbose)

            if (ipsf == -1) then
               if (verbose) print*, 'No good PSF star could be found within frame.'
               deallocate(data)
	       deallocate(pix_flg)
               cycle frame
            end if
	    
            ! We now have a better estimate of fwhm, so use this instead.
            fwhm=1.665*sqrt(shape_par(1)*shape_par(2))
!            cliprad = clip_fwhm*1.665*sqrt(shape_par(1)*shape_par(2))  
            cliprad = 3.0*fwhm            

            ! Save the estimate of the seeing to an array.
            seeing(fnum)=sqrt(1.665*shape_par(1)*1.665*shape_par(2))
	    if (verbose) print*, "@ New FWHM and Clipping radius is: ", fwhm, cliprad
	    if (verbose) print*, "@ Seeing is: ", seeing(fnum)

	    if (verbose) then
               print*, 'Fitted PSF star ', int(psfstars(1,ipsf)), &
                       ' which gave FWsHM ', 1.665*shape_par(1), 1.665*shape_par(2)
               print*, 'Rotated at an angle of ', 57.29*shape_par(3), &
                       ' degrees from the vertical.'
               print*, 'Chosen using ', nfit, ' stars.'
               print*, 'Will use a clipping radius of ', cliprad, ' pixels.'
	    end if

            ! Now get a handle on the flux in the star to be extracted optimally.
            if (iopt > 0) then
               xpos = real(stars(2,iopt))
               ypos = real(stars(3,iopt))
               ! Call the extraction, with the normalisation explicitly set to one, and
               ! return the peak flux in optnrm.
               call extr(data, pix_flg, xpos, ypos, dpos(1), adu, high, low, fwhm, &
                         cliprad, shape_par, 0.0, .false., xcomp, ycomp, optflux, opterror, &
                         xfit, yfit, xerr, yerr, optnrm, cflag, skynos, verbose)
               if (cflag /= 'O') then
                  if (verbose) print*, 'Star to be optimised has flag ', cflag, &
                                       'going to next frame.'
                  deallocate(data)
	          deallocate(pix_flg)
                  cycle frame
               end if
               if (verbose) print*, 'Extractions will be optimised for stars with ',&
                                     optnrm,' peak counts.'
               optnrm=optnrm/(skynos*skynos)
            else                            
               optnrm=0.0
            end if

            if (verbose) then
               print*
               print*
            end if

         else
            ! Don't do optimal photometry, but we will need a FWHM for fitting
            ! the star to find its position.
            shape_par(1)=fwhm/1.665
            shape_par(2)=fwhm/1.665
            shape_par(3)=0.0
            ! And a cliprad.
            cliprad=2.0*fwhm
            ! And (believe it or not) a star flux to be optimised.
            optnrm = 0.0
         end if

         ! ***********************************************************
         ! Go around this loop for each star.
         ! ***********************************************************

         each_star: do istar=1, nstar

	            comp=.False.

                    ! Calculate approximate position of star.
                    xpos = real(stars(2,istar))
                    ypos = real(stars(3,istar))

                    if (optimal) then

                   ! Call the optimal extraction routine.
                     call extr(data, pix_flg, xpos, ypos, dpos(1), adu, & 
                              high, low, fwhm, cliprad, shape_par, optnrm, &
                              comp, xcomp, ycomp, optflux, opterror, &
                              xfit, yfit, xerr, yerr, peak, cflag, skynos,verbose)
		    if (verbose) then
		       print*, 'New position of star',starnames(istar),':', xfit, yfit
		       print*, '1st pass flux is: ', optflux, opterror
		    end if
                    ! Put the flux and its error in the output array.	    
       	            optres(istar,1,fnum)=optflux
	                optres(istar,2,fnum)=opterror

                    ! If the star was not in the frame go to the next one
		            ! otherwise perform a second pass this time fixing the
                    ! position of the star.	    
		    if (opterror<0.0) then 
		       newxypos(istar,1,fnum)=-1.0
		       newxypos(istar,2,fnum)=-1.0
	           apres(istar,1,fnum)=0.0 
		       apres(istar,2,fnum)=-1.0
                       cycle each_star
		    else
                  ! *** Make the code keep a track of the star position in all frames.		    		       
		       newxypos(istar,1,fnum)=xfit
		       newxypos(istar,2,fnum)=yfit
       	       optres(istar,1,fnum)=optflux
	           optres(istar,2,fnum)=opterror
!              Do a second pass but this time with the star positions fixed.
               xpos = xfit
               ypos = yfit
               posfix=-1.0    ! Do not centroid

		       if (verbose) print*, 'Perfoming 2nd pass for star',starnames(istar),&
                                            'with its position fixed to', xfit, yfit
                       call extr(data, pix_flg, xpos, ypos, posfix, adu, & 
                                 high, low, fwhm, cliprad, shape_par, optnrm, &
                                 comp, xcomp, ycomp, optflux, opterror, &
                                 xfit, yfit, xerr, yerr, peak, cflag, skynos,verbose)
                       ! Update the flux and its error in the output array.	    
          	       optres(istar,1,fnum)=optflux
	               optres(istar,2,fnum)=opterror
		    end if

		    if (verbose) print*, 'Opphot flux for object', &
                                 starnames(istar),'is:',optflux,'+/-',opterror
		    end if

		    ! In addition to Optimal photometry do aperature photometry as well
		    ! using the centroided positions that we got from optimal photometry
		    ! xfit, yfit.
            if (aperature) then

			   if (.not.optimal) then
                  xfit = xpos
                  yfit = ypos
               end if

                       sky_stuffed = .false.
                       ! The box size is determined in the way described in the paper.
                       ibox=int(sqrt(628.4*aprad*aprad+ 4.0*fwhm*fwhm))
                       call skyfit(data, xfit, yfit, fwhm, ibox, low, high, &
                                   pix_flg, skycnt, skyerr, skynos, cflag)
!                       if (cflag == 'I') sky_stuffed = .true.
!                       if (cflag == 'B') then
!                          ap%col(1)%data = 0.0
!                          ap%col(1)%err = 0.0
!                       else
                          apflux=sum_aper(data, aprad, skycnt, skyerr, skynos, &
                                          xfit, yfit, adu, low, high, pix_flg, &
                                          aperror, cflag)
!                          if (sky_stuffed) cflag = 'I'
!                          if (ap%col(1)%flg/='ON') ap%col(1)%flg='O'//cflag
!                          if (clip_fwhm < 0.0) then   
!                             if (ap%col(1)%flg == 'OO') then
                        if (verbose) print*, 'Apphot flux for star', stars(1,istar), 'is', &
                                         apflux, 'counts +/- ',aperror
!                             else
                                ! print*, 'Target star ', ap%id,  &
                                !         ' was flagged with ', ap%col(1)%flg
!                             end if
!                          else  
                             ! if (ap%col(1)%flg /= 'OO') &
!                             ! print*, 'Aperture photometry was flagged with ', ap%col(1)%flg
!                          end if
!                       end if
      	            apres(istar,1,fnum)=apflux
	                apres(istar,2,fnum)=aperror
                    end if

		    if (verbose) print*

         end do each_star
 

         deallocate(data)
         deallocate(pix_flg)

  end do frame

  if (verbose) then
     print*, "*****************************"
     print*
  end if

! ***************************************************************
! Write the output files.
! ***************************************************************
!  wformat='(I5," ",F13.4," ",F10.4," ",F7.4)'
! Write the optimal photometry output files.
  if (optimal) then
     do i=1, nstar
        do j=1, ntimes
           print*, starnames(i), optres(i,1,j), optres(i,2,j), seeing(j), newxypos(i,1,j), newxypos(i,2,j) 
        end do   
     end do
  end if

! Write the aperature photometry output files.
  if (aperature) then
     do i=1, nstar
        do j=1, ntimes
           print*, starnames(i), apres(i,1,j), apres(i,2,j), seeing(j), newxypos(i,1,j), newxypos(i,2,j) 
        end do   
     end do
  end if

end program optimal_phot

!******************************************************************************
! printerror - Prints an error if the input FITS file is not properly read.
!******************************************************************************

subroutine printerror(status)

! This subroutine prints out the descriptive text corresponding to the
! error status value and prints out the contents of the internal
! error message stack generated by FITSIO whenever an error occurs.

! Variable declarations
! ---------------------
  integer :: status

  character (len=30) :: errtext
  character (len=80) :: errmessage*80
! ---------------------

! Check if status is OK (no error); if so, simply return
  if (status <= 0)return

! The FTGERR subroutine returns a descriptive 30-character text string that
! corresponds to the integer error status number.  A complete list of all
! the error numbers can be found in the back of the FITSIO User's Guide.
  call ftgerr(status, errtext)
  print *,'* FITSIO Error Status =',status,': ',errtext

! FITSIO usually generates an internal stack of error messages whenever
! an error occurs.  These messages provide much more information on the
! cause of the problem than can be provided by the single integer error
! status value.  The FTGMSG subroutine retrieves the oldest message from
! the stack and shifts any remaining messages on the stack down one
! position.  FTGMSG is called repeatedly until a blank message is
! returned, which indicates that the stack is empty.  Each error message
! may be up to 80 characters in length.  Another subroutine, called
! FTCMSG, is available to simply clear the whole error message stack in
! cases where one is not interested in the contents.
  call ftgmsg(errmessage)
  do while (errmessage /= ' ')
     print *,errmessage
     call ftgmsg(errmessage)
  end do

end subroutine printerror

!******************************************************************************
! xsize - gets the number of columns and rows in an ascii file.
!******************************************************************************

subroutine xsize(unit, cols, rows, blines)

!       This subroutine finds the number of columns and rows in an ASCII file.
!       The program can handle up to 3 blank lines at the top of the file and
!       the file line length must be 80 characters.

!       Variable definitions.
!       *********************************************************************

        character(len=80) :: line

        integer :: i, j, nline, unit, eof1, cols, rows, flag, blines 

        double precision :: dummy

!       *********************************************************************

        cols = 0
        rows = 0
        nline = 4

        if (blines/=0) then
           do j=1, blines
              read(unit,*)
           end do
           nline = 1
        end if

        blanks: do j=1, nline

                   read(unit,'(a)')line

                   columns: do i=1, 80
 
                            if (line(i:i)==' ') then
                               flag = 0
                            else
                               flag = 1
                            end if

                            if (i==1.and.flag==1) then
                               cols = cols + 1
                               cycle columns 
                            end if
                            if (flag==1.and.line(i-1:i-1)==' ') then
                               cols = cols + 1
                            end if

                   end do columns

                   if (cols==0) then
                      blines = blines + 1
                      cycle blanks
                   else
                      exit blanks
                   end if

        end do blanks

        if (cols==0) then
           print*
           print*, '*** WARNING: Failed to find any columns, please check file!'
           print*
           stop
        end if

!       Now measure the number of rows in the file.
        rewind(unit)
        eof1 = 0

        if (blines/=0) then
           do i=1, blines
              read(unit,*)
           end do
        end if

        do
          read (unit,*,iostat=eof1)dummy
          if (eof1/=0) exit
          rows = rows + 1
        end do

end subroutine xsize

!******************************************************************************
! readit - reads in the coordinates and names of the PSF and target stars.
!******************************************************************************
subroutine readit(unit, cols, rows, blines, data, starname)

!       This subroutine reads in files into a data array.
!       The program can handle arrays with up to 5 columns.

!       Variable definitions.
!       *********************************************************************

        integer :: i, unit, eof1, cols, rows, blines

        double precision :: data(cols,rows)   

        character (len=15):: starname(rows)

!       *********************************************************************

        eof1 = 0

!       Go past any blank lines that might exist at the top of the file.
        if (blines/=0) then
           do i=1, blines
              read(unit,*)
           end do
        end if

        do i=1, rows
           if (eof1/=0) exit
           if (cols==4) then
              read(unit,*,iostat=eof1)data(1,i), data(2,i), data(3,i), starname(i)
           else
              print*
              print*, '*** WARNING: Something wrong with the stars or PSF files,'
              print*, '             Check their format and try again!! Exiting...'
              print*
              stop
           end if
        end do

end subroutine readit
