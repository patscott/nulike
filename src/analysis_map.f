***********************************************************************
*** nulike_amap returns the index of the analysis identified by the
*** passed string.
***
*** input:  analysis_name  name of the analysis
*** output:                index of the analysis
***       
*** Author: Pat Scott (p.scott@imperial.ac.uk)
*** Date: Jun 6, 2014
***********************************************************************

      integer function nulike_amap(analysis_name)

      use iso_c_binding, only: c_ptr

      implicit none
      include 'nulike_internal.h'
      
      character(len=nulike_clen) analysis_name
      integer i

      do i = 1, max_analyses
        if (trim(analysis_name_array(i)) .eq. trim(analysis_name)) then
          nulike_amap = i
          return
        endif
      enddo  
      nulike_amap = 0

      end function nulike_amap     
