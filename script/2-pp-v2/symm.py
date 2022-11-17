

    mob = 0.0
    for ik in range(nktot):
        ikirr = ibz2bz[bz2ibz[ik]]
        for ibnd in range(nbnd):
            if scat_kirr[ibnd, ikirr] <= 0.0: #1e-6 a nice try
                continue
            if (not (ecut is None)) and _bande[ik,ibnd] > abs(ecut):
                continue

            # if _bande[ik,ibnd]-np.min(_bande) < 0.005: # FIXME
            #     continue

            mob += 1./scat_kirr[ibnd, ikirr] * mob_w[ik,ibnd]

    mob = mob/ev2inv_s

