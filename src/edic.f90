Program test
  USE environment,   ONLY : environment_start, environment_end

  Implicit None

  

    


  call initialization()


  call calcmdefect_noncolin()
  call environment_end('BANDS')
End Program test






