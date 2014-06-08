***********************************************************************
*** nulike_amap returns the index of the analysis identified by the
*** passed string.
***
*** input:  analysis_name  name of the analysis
*** output:                index of the analysis
***       
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: Jun 6, 2014
***********************************************************************

      integer function nulike_amap(analysis_name)

      implicit none
      include 'nulike.h'
      
      character(len=*) analysis_name
      integer i

      do i = 1, nAnalyses
        if (trim(analysis_name_array(i)) .eq. trim(analysis_name)) then
          nulike_amap = i
          return
        endif
      enddo  
      nulike_amap = 0

      end function nulike_amap     
