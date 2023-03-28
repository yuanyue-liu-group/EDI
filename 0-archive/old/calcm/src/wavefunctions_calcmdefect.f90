Module wavefunctions_calcmdefect
    Use kinds,   Only : dp

    Implicit None
    save

    complex(dp), allocatable, target :: &
        evc1(: , :), &
        evc2(: , :), &
        evc3(: , :), &
        evc4(: , :) 

    complex(dp), allocatable, target :: &
        psic1(:), &
        psic2(:), &
        psic3(:), &
        psic4(:) 

End Module wavefunctions_calcmdefect