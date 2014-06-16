C**********************************************************************
C Initialization of common block elements must appear in block data
C routines
C**********************************************************************
      block data nulike_init_flagblock
      logical nulike_init_called    
      common /nulike_init_flag/ nulike_init_called
      data nulike_init_called/.false./
      end

      block data nulike_credit_flagblock
      logical credits_rolled    
      common /nulike_credit_flag/ credits_rolled
      data credits_rolled/.false./
      end

