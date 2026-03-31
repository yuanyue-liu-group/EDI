MODULE openfilepw_mod
  IMPLICIT NONE
  SAVE
  INTEGER :: iuwfc_epw = 20
  INTEGER :: lrwfc_epw = 0
END MODULE openfilepw_mod

SUBROUTINE openfilepw()
  USE io_files, ONLY : prefix, diropn
  USE wvfct, ONLY : nbnd, npwx
  USE noncollin_module, ONLY : npol
  USE openfilepw_mod, ONLY : iuwfc_epw, lrwfc_epw
  IMPLICIT NONE
  LOGICAL :: exst

  IF (LEN_TRIM(prefix) == 0) CALL errore('openfilepw', 'wrong prefix', 1)
  iuwfc_epw = 20
  lrwfc_epw = 2 * nbnd * npwx * npol
  CALL diropn(iuwfc_epw, 'wfc', lrwfc_epw, exst)
  IF (.NOT. exst) CALL errore('openfilepw', &
       'file ' // TRIM(prefix) // '.wfc not found', 1)
END SUBROUTINE openfilepw
